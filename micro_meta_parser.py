import asyncio
import csv
import json
import os
import pickle
import tarfile
import time
import urllib.parse
import uuid
from collections import Counter
from collections.abc import Iterator
from pathlib import Path
from typing import Dict, List, Tuple

import aiohttp
import biothings_client as bt
from dotenv import load_dotenv
from ete3 import NCBITaxa
from tqdm.asyncio import tqdm_asyncio

"""
column names with index:
{0: 'id',
 1: 'g_micro',
 2: 'organism',
 3: 'g_meta',
 4: 'metabolic',
 5: 'pubchem_compound',
 6: 'pubchem_id',
 7: 'formula',
 8: 'kegg_id',
 9: 'tax_id',
 10: 'phylum',
 11: 'class',
 12: 'order',
 13: 'family',
 14: 'genus',
 15: 'species',
 16: 'species_id',
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


def assign_col_val_if_available(node: dict, key: str, val: str | int, transform=None):
    """assigns a value to a specified key in a dictionary if the value is available and not equal to "not available"
    This function updates the given dictionary with the provided key and value.
    It can also transform the value using a specified function before assigning it.

    :param node: The dictionary to be updated.
    :param key: The key to be assigned the value in the dictionary.
    :param val: The value to be assigned to the key in the dictionary.
    :param transform: An optional function to transform the value before assignment. If None, the value is assigned as is.
    :return: None
    """
    if val and val != "not available":
        node[key] = transform(val) if transform else val


def assign_to_xrefs_if_available(node: dict, key: str, val: str | int, transform=None):
    """assigns a value to the 'xrefs' sub-dictionary of a given dictionary
    This function checks if the 'xrefs' key exists in the given dictionary.
    If not, it initializes 'xrefs' as an empty dictionary.
    Then, it assigns the provided value to the specified key within the 'xrefs' dictionary
    It can also transform the value using a specified function.

    :param node: The dictionary to be updated.
    :param key: The key to be assigned the value in the 'xrefs' sub-dictionary.
    :param val: The value to be assigned to the key in the 'xrefs' sub-dictionary.
    :param transform: An optional function to transform the value before assignment. If None, the value is assigned as is.
    :return: None
    """
    if val and val != "not available":
        if "xrefs" not in node:
            node["xrefs"] = {}

        node["xrefs"][key] = transform(val) if transform else val


def get_taxon_info(taxids: list) -> list:
    """retrieves taxonomic information for a given list of taxon IDs from disease_species.csv

    This function reads taxon IDs, removes duplicates, and queries taxonomic info from biothings_client
    to retrieve detailed taxonomic information including scientific name, parent taxid, lineage, and rank.

    :param taxids:
    :return: A list of dictionaries containing taxonomic information.
    """
    taxids = set(taxids)
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

    sections = data.get("Record", {}).get("Section", [])
    for sec in sections:
        if sec.get("TOCHeading") == "Names and Identifiers":
            for sub in sec.get("Section", []):
                for info in sub.get("Information", []):
                    markup = info.get("Value", {}).get("StringWithMarkup", [])
                    if markup:
                        descr = markup[0].get("String")
                        if descr:
                            description = descr.strip()
        elif sec.get("TOCHeading") == "Synonyms":
            for sub in sec.get("Section", []):
                for info in sub.get("Information", []):
                    markup = info.get("Value", {}).get("StringWithMarkup", [])
                    if markup:
                        synonym = markup[0].get("String")
                        if synonym:
                            synonyms.append(synonym)

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

    for taxid, organism_type in taxon_map.items():
        if taxid in node.get("lineage", []):
            return organism_type

    return "Other"


def bt_get_mw_logp(cids: list):
    cids = set(cids)
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
        taxid_mapping[old_taxid] = merged_mapping[old_taxid]
    return taxid_mapping


def get_taxon_names(taxon_info: dict) -> list[str]:
    """Extracts biothings names from the taxon_info dictionary."""
    taxon_names = set()
    for _, taxon in taxon_info.items():
        if "scientific_name" in taxon:
            taxon_names.add(taxon["scientific_name"].lower())
    return list(taxon_names)


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


def get_node_info(file_path: str | os.PathLike) -> Iterator[dict]:
    """generates node information from micro_metabolic.csv.
    This function reads data from micro_metabolic.csv, processes taxonomic information,
    and generates nodes representing metabolites and microbes, including their relationships.

    :param file_path: Path to the micro_metabolic.csv file
    :return: An iterator of dictionaries containing node information.
    """
    taxon_info = {
        int(taxon["query"]): taxon
        for taxon in get_taxon_info(file_path)
        if "notfound" not in taxon.keys()
    }

    for line in line_generator(file_path):
        # create object node (metabolites)
        object_node = {
            "id": None,
            "name": line[5].lower(),
            "type": "biolink:SmallMolecule",
        }

        assign_col_val_if_available(object_node, "pubchem_cid", line[6], int)
        assign_col_val_if_available(object_node, "chemical_formula", line[7])
        assign_col_val_if_available(object_node, "smiles", line[18])

        if "pubchem_cid" in object_node:
            assign_to_xrefs_if_available(object_node, "kegg_compound", line[8])
        else:
            assign_col_val_if_available(object_node, "kegg_compound", line[8])
        if "pubchem_cid" not in object_node and "kegg_compound" not in object_node:
            assign_col_val_if_available(object_node, "hmdb", line[19])
        else:
            assign_to_xrefs_if_available(object_node, "hmdb", line[19])

        if "pubchem_cid" in object_node:
            object_node["id"] = f"PUBCHEM.COMPOUND:{object_node['pubchem_cid']}"
        elif "kegg_compound" in object_node:
            object_node["id"] = f"KEGG.COMPOUND:{object_node['kegg_compound']}"
        elif "hmdb" in object_node:
            object_node["id"] = f"HMDB:{object_node['hmdb']}"
        else:
            object_node["id"] = str(uuid.uuid4())

        # create subject node (microbes)
        subject_node = {"id": None, "name": line[2].lower(), "type": "biolink:OrganismalEntity"}

        assign_col_val_if_available(subject_node, "taxid", line[9], int)
        if "taxid" in subject_node:
            subject_node["id"] = f"taxid:{subject_node['taxid']}"
        else:
            subject_node["id"] = str(uuid.uuid4())

        if subject_node.get("taxid") in taxon_info:
            subject_node["scientific_name"] = taxon_info[subject_node["taxid"]]["scientific_name"]
            subject_node["parent_taxid"] = taxon_info[subject_node["taxid"]]["parent_taxid"]
            subject_node["lineage"] = taxon_info[subject_node["taxid"]]["lineage"]
            subject_node["rank"] = taxon_info[subject_node["taxid"]]["rank"]

        # association node has the reference and source of metabolites
        association_node = {
            "predicate": "biolink:associated_with",
            "infores": line[17],
        }
        if line[20] and line[20] != "Unknown":
            association_node["sources"] = [src.strip().lower() for src in line[20].split(";")]

        output_dict = {
            "_id": None,
            "association": association_node,
            "object": object_node,
            "subject": subject_node,
        }
        if ":" in object_node["id"] and ":" in subject_node["id"]:
            output_dict[
                "_id"
            ] = f"{subject_node['id'].split(':')[1].strip()}_associated_with_{object_node['id'].split(':')[1].strip()}"
        elif ":" not in object_node["id"] and ":" in subject_node["id"]:
            output_dict[
                "_id"
            ] = f"{subject_node['id'].split(':')[1].strip()}_associated_with_{object_node['id']}"
        elif ":" in object_node["id"] and ":" not in subject_node["id"]:
            output_dict[
                "_id"
            ] = f"{subject_node['id']}_associated_with_{object_node['id'].split(':')[1].strip()}"
        else:
            output_dict["_id"] = f"{subject_node['id']}_associated_with_{object_node['id']}"

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


# if __name__ == "__main__":
#     path = os.getcwd()
#     file_path = os.path.join(path, "data", "micro_metabolic.csv")
# micro_meta_data = load_micro_meta_data()
# type_list = [obj["subject"]["type"] for obj in micro_meta_data]
# type_counts = Counter(type_list)
# for value, count in type_counts.items():
#     print(f"{value}: {count}")

# _ids = []
# for obj in micro_meta_data:
#     print(obj)
#     _ids.append(obj["_id"])
# print(f"total records: {len(_ids)}")
# print(f"total records without duplicates: {len(set(_ids))}")
