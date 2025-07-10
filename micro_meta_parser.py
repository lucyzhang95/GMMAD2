import asyncio
import csv
import json
import os
import pickle
import tarfile
import time
import uuid
import zipfile
from collections import Counter
from collections.abc import Iterator
from pathlib import Path
from typing import Dict, List

import aiohttp
import biothings_client as bt
import pandas as pd
import requests
from dotenv import load_dotenv
from ete3 import NCBITaxa
from tqdm.asyncio import tqdm_asyncio

"""
column names with index:
{0: 'id', -> file index
 1: 'g_micro',
 2: 'organism',
 3: 'g_meta',
 4: 'metabolic',
 5: 'pubchem_compound',
 6: 'pubchem_id',
 7: 'formula',
 8: 'kegg_id',
 9: 'tax_id',  -> current taxid
 10: 'phylum',
 11: 'class',
 12: 'order',
 13: 'family',
 14: 'genus',
 15: 'species',
 16: 'species_id', -> species taxid
 17: 'source',
 18: 'smiles_sequence',
 19: 'HMDBID',
 20: 'Origin'}
"""

load_dotenv()
CACHE_DIR = os.path.join(os.getcwd(), "cache")
os.makedirs(CACHE_DIR, exist_ok=True)


def save_pickle(obj, f_name):
    """
    :param obj: data to be saved as a pickle file
    :param f_name: files should only be existing in the cache directory
    :return:
    """
    with open(os.path.join(CACHE_DIR, f_name), "wb") as out_f:
        pickle.dump(obj, out_f)


def load_pickle(f_name):
    path = os.path.join(CACHE_DIR, f_name)
    return (
        pickle.load(open(path, "rb")) if os.path.exists(path) else print("The file does not exist.")
    )


def save_json(obj, f_name):
    with open(os.path.join(CACHE_DIR, f_name), "w") as out_f:
        json.dump(obj, out_f, indent=4)


def line_generator(in_file: str | os.PathLike, delimiter=",") -> Iterator[list]:
    """generates lines from a CSV file, yielding each line as a list of strings
    This function opens the specified CSV file, skips the header row, and yields each subsequent line as a list of strings.

    :param in_file: The path to the CSV file.
    :param delimiter: The character used to separate values in the CSV file (default is comma).
    :return: An iterator that yields each line of the CSV file as a list of strings.
    """
    with open(in_file, "r") as in_f:
        reader = csv.reader(in_f, delimiter=delimiter)
        next(reader)
        for line in reader:
            yield line


def get_primary_id(line, bigg_map):
    pubchem_id = line[6]
    kegg_id = line[8]
    smiles = line[18]
    hmdb_id = line[19]
    bigg_id = bigg_map.get(line[5].lower())

    # (line, prefix) pairs for ID hierarchy
    id_hierarchy = [
        (pubchem_id, "PUBCHEM.COMPOUND"),
        (kegg_id, None),
        (smiles, ""),
        (hmdb_id, "HMDB"),
        (bigg_id, "BIGG.METABOLITE"),
    ]

    def classify_kegg(val):
        if val.startswith("C"):
            return "KEGG.COMPOUND"
        elif val.startswith("G"):
            return "KEGG.GLYCAN"
        elif val.startswith("D"):
            return "KEGG.DRUG"
        return "KEGG"

    xrefs = {}
    primary_id = None

    for val, prefix in id_hierarchy:
        if not val or val.strip().lower() == "not available":
            continue
        if prefix is None and val == kegg_id:
            prefix = classify_kegg(val)

        curie = f"{prefix}:{val}" if prefix else val
        if prefix == "PUBCHEM.COMPOUND":
            key = "pubchem_cid"
        elif prefix == "HMDB":
            key = "hmdb"
        elif prefix is not None and prefix.startswith("KEGG"):
            key = "kegg"
        elif prefix == "":
            key = "smiles"
        elif prefix == "BIGG.METABOLITE":
            key = "bigg"
        else:
            continue

        if primary_id is None:
            primary_id = curie
        xrefs.setdefault(key, curie)

    return primary_id, xrefs


def get_taxon_info(taxids: list) -> list:
    """retrieves taxonomic information for a given list of taxon IDs from disease_species.csv

    This function reads taxon IDs, removes duplicates, and queries taxonomic info from biothings_client
    to retrieve detailed taxonomic information including scientific name, parent taxid, lineage, and rank.

    :param taxids:
    :return: A list of dictionaries containing taxonomic information.
    """
    taxids = sorted(set(taxids))
    t = bt.get_client("taxon")
    taxon_info = t.gettaxa(taxids, fields=["scientific_name", "parent_taxid", "lineage", "rank"])
    return taxon_info


async def pug_query_pubchem_description(
    cid: int,
    session: aiohttp.ClientSession,
    sem: asyncio.Semaphore,
    max_retries: int = 3,
    delay: float = 1.0,
):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
    for attempt in range(max_retries):
        async with sem:
            try:
                async with session.get(url, timeout=30) as resp:
                    resp.raise_for_status()
                    data = await resp.json()
                break
            except aiohttp.ClientError as e:
                if attempt < max_retries - 1:
                    await asyncio.sleep(delay * (2**attempt))
                    continue
                else:
                    raise e

    description = None
    synonyms = []

    for sec in data.get("Record", {}).get("Section", []):
        heading = sec.get("TOCHeading", "")
        if heading == "Names and Identifiers" and description is None:
            for sub in sec.get("Section", []):
                if sub.get("TOCHeading") == "Record Description":
                    for info in sub.get("Information", []):
                        for mark in info.get("Value", {}).get("StringWithMarkup", []):
                            text = mark.get("String", "").strip()
                            if text:
                                description = text
                                break
                        if description:
                            break
                elif sub.get("TOCHeading") == "Synonyms":
                    for sec in sub.get("Section", []):
                        for info in sec.get("Information", []):
                            for mark in info.get("Value", {}).get("StringWithMarkup", []):
                                text = mark.get("String", "").strip()
                                if text:
                                    synonyms.append(text.lower().strip())
            break

    seen = set()
    synonyms = [s for s in synonyms if not (s in seen or seen.add(s))]

    return cid, {
        "id": f"PUBCHEM.COMPOUND:{cid}",
        "description": f"{description}[PUBCHEM]" if description else "",
        "synonyms": synonyms,
    }


async def get_batch_pubchem_descriptions_async(
    cids: List[int],
    workers: int = 5,
) -> Dict[int, Dict[str, str]]:
    cids = list(set(cids))

    workers = min(workers, 5)
    sem = asyncio.Semaphore(workers)
    connector = aiohttp.TCPConnector(limit_per_host=workers)
    cids = sorted(list(set(cids)))

    results = {}
    async with aiohttp.ClientSession(connector=connector) as session:
        tasks = []
        last_launch = 0.0

        for cid in cids:
            elapsed = time.perf_counter() - last_launch
            if elapsed < 1.0 / 5:
                await asyncio.sleep(1.0 / 5 - elapsed)
            last_launch = time.perf_counter()

            task = asyncio.create_task(pug_query_pubchem_description(cid, session, sem))
            tasks.append(task)

        for cid, payload in await tqdm_asyncio.gather(*tasks, total=len(tasks)):
            if payload:
                results[cid] = payload

    return results


def get_pubchem_descriptions(cids: List[int], workers: int = 5):
    return asyncio.run(get_batch_pubchem_descriptions_async(cids, workers))


def get_bigg_metabolite_mapping(in_f):
    bigg_map = {line[2].lower(): line[1] for line in line_generator(in_f, delimiter="\t")}
    return bigg_map


def get_organism_type(node) -> str:
    """
    Inspect node['lineage'] for known taxids.
    Return the matching biolink CURIE, or Other if no match.
    Types include: 3 domains of life (Bacteria, Archaea, Eukaryota) and Virus.
    """
    taxon_map = {
        2: "biolink:Bacterium",
        2157: "Archaeon",
        2759: "Eukaryote",
        10239: "biolink:Virus",
    }

    lineage = node.get("lineage")
    if not isinstance(lineage, list):
        return "Other"
    for taxid, organism_type in taxon_map.items():
        if taxid in lineage:
            return organism_type
    return "Other"


def bt_get_mw_logp(cids: list):
    cids = sorted(set(cids))
    t = bt.get_client("chem")
    get_chem = t.querymany(
        cids,
        scopes="pubchem.cid",
        fields=["pubchem.molecular_weight", "pubchem.monoisotopic_weight", "pubchem.xlogp"],
    )

    q_out = {}
    for q in get_chem:
        if "notfound" in q:
            continue

        pub = q.get("pubchem")
        if isinstance(pub, list):
            info = pub[0]
        elif isinstance(pub, dict):
            info = pub
        else:
            continue

        mw = info.get("molecular_weight")
        mono = info.get("monoisotopic_weight")
        logp = info.get("xlogp")
        q_out[q["query"]] = {
            "molecular_weight": {
                "average_molecular_weight": mw,
                "monoisotopic_molecular_weight": mono,
            },
            "xlogp": logp,
        }

    return q_out


def load_merged_from_tar(tar_gz_path, f_name="merged.dmp"):
    """Parse 'merged.dmp' of taxdump.tar.gz downloaded by ete3.
    Returns a dict {old_taxid: new_taxid}.
    """
    if not os.path.exists(tar_gz_path):
        ncbi = NCBITaxa()
        ncbi.update_taxonomy_database()

    taxid_mapping = {}
    with tarfile.open(tar_gz_path, "r:gz") as tar:
        f = tar.getmember(f_name)
        with tar.extractfile(f) as fp:
            for line in fp:
                parts = line.decode("utf-8").split("\t")
                old, new = parts[0], parts[2]
                taxid_mapping[old] = new
    return taxid_mapping


def get_current_taxid(old_taxids: list, merged_mapping: dict) -> dict[str, str]:
    taxid_mapping = {}
    for old_taxid in old_taxids:
        taxid_mapping[old_taxid] = merged_mapping.get(old_taxid, None)
    return taxid_mapping


def get_taxon_names(taxon_info: dict) -> list[str]:
    """Extracts biothings names from the taxon_info dictionary."""
    taxon_names = set()
    for _, taxon in taxon_info.items():
        if "scientific_name" in taxon:
            taxon_names.add(taxon["scientific_name"].lower())
    return list(taxon_names)


def dump_ncit_source(f_name, out_path):
    url = f"https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/{f_name}"
    try:
        resp = requests.get(url, stream=True)
        resp.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"{f_name} downloaded successfully.")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred while downloading the file: {e}")


def build_ncit_organism_mapping(f_name, path):
    if not Path(path).exists():
        dump_ncit_source(f_name, path)

    with zipfile.ZipFile(path) as zf:
        with zf.open("Thesaurus.txt") as f:
            df = pd.read_csv(f, sep="\t", dtype=str).fillna("")

    df.columns = [
        "code",
        "concept IRI",
        "parents",
        "synonyms",
        "definition",
        "display name",
        "concept status",
        "semantic type",
        "concept in subset",
    ]

    relevant_types = {
        "Bacterium",
        "Virus",
        "Fungus",
        "Eukaryote",
        "Organism",
        "Animal",
        "Group|Organism",
    }

    df_organism = df[df["semantic type"].isin(relevant_types)].copy()

    mapping_result = {}
    for _, row in df_organism.iterrows():
        name = str(row["display name"]).strip().lower()
        description = str(row["definition"]).strip()
        ncit_code = str(row["code"]).strip()
        if not name or pd.isna(description):
            continue
        mapping_result[name] = {"description": description, "xrefs": {"ncit": ncit_code}}
    return mapping_result


async def fetch_ncit_description(
    session: aiohttp.ClientSession, name: str, sem: asyncio.Semaphore
) -> tuple[str, dict]:
    NCIT_API_KEY = os.getenv("NCIT_API_KEY")
    SEARCH_URL = "https://data.bioontology.org/search"
    params = {
        "q": name,
        "ontologies": "NCIT",
        "apikey": NCIT_API_KEY,
    }
    async with sem:
        async with session.get(SEARCH_URL, params=params) as resp:
            resp.raise_for_status()
            data = await resp.json()
    for result in data.get("collection", []):
        if not result:
            continue
        pref_label = result.get("prefLabel", "").lower()
        if pref_label != name:
            continue
        definition = result.get("definition", [])
        ncit_id = result.get("@id", "").split("#")[-1]
        return name, {
            "description": f"{definition[0]} [NCIT]" if definition else "",
            "xrefs": {"ncit": ncit_id},
        }


async def get_ncit_taxon_description_async(taxon_names, max_concurrent=5):
    unique_names = {n.lower() for n in taxon_names}
    sem = asyncio.Semaphore(max_concurrent)
    connector = aiohttp.TCPConnector(limit_per_host=max_concurrent)
    async with aiohttp.ClientSession(connector=connector) as session:
        tasks = [fetch_ncit_description(session, name, sem) for name in unique_names]
        results = await asyncio.gather(*tasks)
    return {name: data for item in results if item for name, data in (item,)}


def get_ncit_taxon_description(taxon_names):
    """

    :param taxon_names:
    :return:
    {
       "550":{
          "query":"550",
          "_id":"550",
          "_version":1,
          "lineage":[
             550,
             354276,
             547,
             543,
             91347,
             1236,
             1224,
             3379134,
             2,
             131567,
             1
          ],
          "parent_taxid":354276,
          "rank":"species",
          "scientific_name":"enterobacter cloacae",
          "description":"A species of facultatively anaerobic, Gram negative, rod shaped bacterium in the phylum Proteobacteria. This species is motile by peritrichous flagella, oxidase, urease and indole negative, catalase positive, reduces nitrate, does not degrade pectate and produces acid from sorbitol. E. cloacae is associated with hospital-acquired urinary and respiratory tract infections and is used in industry for explosives biodegradation. [NCIT]",
          "xrefs": {
             "ncit":"C86360"
          }
       }
    }
    """
    return asyncio.run(get_ncit_taxon_description_async(taxon_names))


def add_description2taxon_info(taxon_info: dict, descriptions: dict) -> dict:
    for info in taxon_info.values():
        name = info.get("scientific_name").lower()
        descr_info = descriptions.get(name, {})
        info.update(descr_info)
    return taxon_info


def remove_empty_none_values(obj):
    if isinstance(obj, dict):
        cleaned = {}
        for k, v in obj.items():
            v_clean = remove_empty_none_values(v)
            if v_clean not in (None, {}, []):
                cleaned[k] = v_clean
        return cleaned

    if isinstance(obj, list):
        cleaned_list = []
        for v in obj:
            v_clean = remove_empty_none_values(v)
            if v_clean not in (None, {}, []):
                cleaned_list.append(v_clean)
        return cleaned_list
    return obj


def cache_data(
    f_path,
    gzip_path=None,
    bigg_path=None,
):
    if gzip_path is None:
        gzip_path = os.path.join("downloads", "taxdump.tar.gz")
    if bigg_path is None:
        bigg_path = os.path.join("downloads", "bigg_models_metabolites.txt")

    # cache metabolite descriptions
    pubchem_cids = [
        line[6] for line in line_generator(f_path) if line[6] and line[6] != "not available"
    ]
    print(f"Total unique pubchem_cids: {len(pubchem_cids)}")
    pubchem_descr = get_pubchem_descriptions(pubchem_cids)
    print(f"Total pubchem_cid with descriptions: {len(pubchem_descr)}")
    save_pickle(pubchem_descr, "gmmad2_micro_meta_description.pkl")

    # cache molecular weight and xlogp
    pubchem_mw_logp = bt_get_mw_logp(pubchem_cids)
    print(f"Total pubchem_cid with molecular weight and xlogp: {len(pubchem_mw_logp)}")
    save_pickle(pubchem_mw_logp, "gmmad2_micro_meta_pubchem_mw.pkl")

    # cache bigg metabolite mapping
    bigg_mapping = get_bigg_metabolite_mapping(bigg_path)
    save_pickle(bigg_mapping, "gmmad2_micro_meta_bigg_mapping.pkl")

    # cache taxon info
    taxids = [
        line[9]
        if (line[9] and line[9] != "not available")
        else line[16]
        if (line[16] and line[16] != "not available")
        else None
        for line in line_generator(f_path)
    ]  #
    taxids = [t for t in taxids if t]
    print(f"Total unique taxids: {len(set(taxids))}")
    taxon_info_q = get_taxon_info(taxids)
    taxon_info = {t["query"]: t for t in taxon_info_q if "notfound" not in t.keys()}
    print(f"Biothings mapped {len(taxon_info)} taxids.")
    notfound = sorted(set([t["query"] for t in taxon_info_q if "notfound" in t.keys()]))
    print(f"Not found taxids: {len(notfound)}")
    taxid_mapping = load_merged_from_tar(gzip_path)
    new_taxid_map = get_current_taxid(notfound, taxid_mapping)
    print(f"NCBI mapped taxids: {len(new_taxid_map)}")
    new_taxids = [new for old, new in new_taxid_map.items()]
    print(f"taxids need to be queried: {len(new_taxids)}")
    new_taxon_q = {t["query"]: t for t in get_taxon_info(new_taxids) if "notfound" not in t}
    new_taxon_info = {
        old: new_taxon_q[new] for old, new in new_taxid_map.items() if new in new_taxon_q
    }
    print(f"Biothings mapped new taxids: {len(new_taxon_info)}")
    taxon_info.update(new_taxon_info)
    print(f"Merged taxon info: {len(taxon_info)}")
    save_pickle(taxon_info, "gmmad2_micro_meta_taxon_info.pkl")

    # cache taxon descriptions from NCIT
    taxon_names = get_taxon_names(taxon_info)
    print(f"Total unique taxon names: {len(taxon_names)}")
    ncit_descriptions = get_ncit_taxon_description(taxon_names)
    print(f"NCIT descriptions found for {len(ncit_descriptions)} taxon names.")
    taxon_info_w_descr = add_description2taxon_info(taxon_info, ncit_descriptions)
    save_pickle(taxon_info_w_descr, "gmmad2_micro_meta_taxon_info_w_descr.pkl")


def get_suffix(identifier: str) -> str:
    return identifier.split(":", 1)[1].strip() if ":" in identifier else identifier.strip()


def get_node_info(file_path: str | os.PathLike) -> Iterator[dict]:
    """generates node information from micro_metabolic.csv.
    This function reads data from micro_metabolic.csv, processes taxonomic information,
    and generates nodes representing metabolites and microbes, including their relationships.

    :param file_path: Path to the micro_metabolic.csv file
    :return: An iterator of dictionaries containing node information.
    """
    if not os.path.exists(Path(os.path.join("cache", "gmmad2_micro_meta_taxon_info_w_descr.pkl"))):
        print("Taxon info not found in cache, running cache_data()...")
        cache_data(file_path)

    # load cached data
    taxon_info = load_pickle("gmmad2_micro_meta_taxon_info_w_descr.pkl")
    pubchem_descr = load_pickle("gmmad2_micro_meta_description.pkl")
    pubchem_mw = load_pickle("gmmad2_micro_meta_pubchem_mw.pkl")
    bigg_mapping = load_pickle("gmmad2_micro_meta_bigg_mapping.pkl")

    for line in line_generator(file_path):
        # create object node (metabolites)
        object_node = {
            "id": None,
            "name": line[5].lower() if line[5] else line[4].lower(),
            "synonym": pubchem_descr.get(line[6], {}).get("synonyms", []),
            "description": pubchem_descr.get(line[6], {}).get("description", ""),
            "chemical_formula": line[7] if line[7] else None,
            "molecular_weight": pubchem_mw.get(line[6], {}).get("molecular_weight", {}),
            "xlogp": pubchem_mw.get(line[6], {}).get("xlogp", None),
            "type": "biolink:SmallMolecule",
            "xrefs": {},
        }

        id_xrefs = get_primary_id(line, bigg_mapping)
        if id_xrefs:
            primary_id, xrefs = id_xrefs
        else:
            primary_id, xrefs = None, {}
        object_node["id"] = primary_id if primary_id else str(uuid.uuid4())
        object_node["xrefs"] = xrefs
        object_node = remove_empty_none_values(object_node)

        # create subject node (microbes)
        taxid = (
            line[9]
            if line[9] and line[9] != "not available"
            else line[16]
            if line[16] and line[16] != "not available"
            else None
        )
        if taxid and taxid in taxon_info:
            subject_node = {
                "id": f"NCBITaxon:{taxon_info[taxid].get('_id')}",
                "taxid": int(taxon_info[taxid].get("_id")),
                "name": taxon_info[taxid].get("scientific_name").lower(),
                "original_name": line[2].lower(),
                "description": taxon_info[taxid].get("description"),
                "parent_taxid": taxon_info[taxid].get("parent_taxid"),
                "lineage": taxon_info[taxid].get("lineage", []),
                "rank": taxon_info[taxid].get("rank"),
                "type": "biolink:OrganismTaxon",
                "organism_type": None,
                "xrefs": taxon_info[taxid].get("xrefs", {}),
            }
            subject_node["organism_type"] = get_organism_type(subject_node)
            subject_node = remove_empty_none_values(subject_node)
        else:
            subject_node = {
                "id": str(uuid.uuid4()),
                "original_name": line[2].lower(),
                "type": "biolink:OrganismTaxon",
                "organism_type": "Other",
            }
            subject_node = remove_empty_none_values(subject_node)

        # association node has the habitat and source of the interaction
        evidence_map = {
            "infores:wom": "ECO:0001230",  # mass spec + manual
            "infores:vmh": "ECO:0000218",  # manual assertion
            "infores:gutMGene": "ECO:0000218",  # manual assertion
            "infores:Metabolomics": "ECO:0001230",  # mass spec + manual
        }

        src = f"infores:{line[17].strip()}"
        habitat = (
            [h.strip().lower() for h in line[20].split(";")]
            if line[20] and line[20] != "Unknown"
            else None
        )

        association_node = {
            "predicate": "biolink:OrganismTaxonToChemicalEntityAssociation",
            "type": "has_metabolic_interaction_with",
            "association_habitat": habitat,
            "primary_knowledge_source": src,
            "aggregator_knowledge_source": "infores:GMMAD2",
            "evidence_type": evidence_map.get(src, "ECO:0000000"),
        }

        association_node = remove_empty_none_values(association_node)

        output_dict = {
            "_id": None,
            "association": association_node,
            "object": object_node,
            "subject": subject_node,
        }
        subject_suffix = get_suffix(subject_node["id"])
        object_suffix = get_suffix(object_node["id"])
        output_dict["_id"] = f"{subject_suffix}_has_metabolic_interaction_with_{object_suffix}"

        yield output_dict


def load_micro_meta_data(f_path) -> Iterator[dict]:
    """loads and yields unique microbe-metabolite data records from micro_metabolic.csv file
    This function constructs the file path to the micro_metabolic.csv file,
    retrieves node information using the `get_node_info` function, and yields unique records based on `_id`.

    :return: An iterator of dictionaries containing microbe-metabolite data.
    """
    assert os.path.exists(Path(f_path))

    dup_ids = set()
    recs = get_node_info(f_path)
    for rec in recs:
        if rec["_id"] not in dup_ids:
            dup_ids.add(rec["_id"])
            yield rec


if __name__ == "__main__":
    file_path = os.path.join("downloads", "micro_metabolic.csv")
    micro_meta_data = [obj for obj in load_micro_meta_data(file_path)]
    type_list = [obj["subject"].get("organism_type") for obj in micro_meta_data]
    type_counts = Counter(type_list)
    print(type_counts)

    _ids = []
    for obj in micro_meta_data:
        # print(obj)
        _ids.append(obj["_id"])
    print(f"total records: {len(_ids)}")
    print(f"total records without duplicates: {len(set(_ids))}")
