import asyncio
import json
import os
import pickle
import tarfile
from collections import Counter
from collections.abc import Iterator
from pathlib import Path

import aiohttp
import biothings_client
from dotenv import load_dotenv
from ete3 import NCBITaxa

"""updated on 07/07/2025
column names with index:
{0: 'id',
 1: 'disease_id',
 2: 'disease',
 3: 'organism',
 4: 'level',
 5: 'species_id',
 6: 'disease_samples',
 7: 'disease_mean',
 8: 'disease_median',
 9: 'disease_sd',
 10: 'health_id',
 11: 'health',
 12: 'health_samples',
 13: 'health_mean',
 14: 'health_median',
 15: 'health_sd',
 16: 'change',
 17: 'alteration',
 18: 'disease_info',
 19: 'phylum',
 20: 'class',
 21: 'order',
 22: 'family',
 23: 'genus'}
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


def line_generator_4_midi(in_file: str | os.PathLike) -> Iterator[list[str]]:
    """Yield each CSV line as a list of exactly 24 fields,
    rejoining misaligned columns with commas for the 'disease' and 'disease_info' columns.
    """
    EXPECTED_COUNT = 24
    HEAD_COUNT = 2
    TAIL_COUNT = 5
    FIXED_COUNT = 15

    with open(in_file, "r", encoding="utf-8") as f:
        next(f)
        for line in f:
            parts = [part.strip() for part in line.split(",")]
            if len(parts) == EXPECTED_COUNT:
                yield parts
            else:
                alteration_idx = next(
                    i for i, p in enumerate(parts) if p in ("Increase", "Decrease")
                )
                head = parts[:HEAD_COUNT]
                fixed_start = alteration_idx - (FIXED_COUNT - 1)
                disease = ",".join(parts[2:fixed_start])
                fixed = parts[fixed_start : alteration_idx + 1]
                disease_info = ",".join(parts[alteration_idx + 1 : len(parts) - TAIL_COUNT])
                tail = parts[-TAIL_COUNT:]
                new_line = head + [disease] + fixed + [disease_info] + tail
                assert (
                    len(new_line) == EXPECTED_COUNT
                ), f"Expected {EXPECTED_COUNT} cols, got {len(new_line)}"
                yield new_line


def get_taxon_info(taxids: list) -> list:
    """retrieves taxonomic information for a given list of taxon IDs from disease_species.csv

    This function reads taxon IDs, removes duplicates, and queries taxonomic info from biothings_client
    to retrieve detailed taxonomic information including scientific name, parent taxid, lineage, and rank.

    :param taxids:
    :return: A list of dictionaries containing taxonomic information.
    """
    taxids = set(taxids)
    t = biothings_client.get_client("taxon")
    taxon_info = t.gettaxa(taxids, fields=["scientific_name", "parent_taxid", "lineage", "rank"])
    return taxon_info


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


def cache_data(f_path, gzip_path="taxdump.tar.gz"):
    # cache all taxon info including lineage, rank, parent_taxid from disease_species.csv
    taxids = sorted(set([line[5] for line in line_generator_4_midi(f_path)]))
    print(f"Total unique taxids: {len(taxids)}")
    taxon_info_q = get_taxon_info(taxids)
    taxon_info = {t["query"]: t for t in taxon_info_q if "notfound" not in t.keys()}
    print(f"Biothings mapped {len(taxon_info)} taxids.")
    notfound = sorted(set([t["query"] for t in taxon_info_q if "notfound" in t.keys()]))
    print(f"Not found taxids: {len(notfound)}")
    taxid_mapping = load_merged_from_tar(gzip_path)
    new_taxid_map = get_current_taxid(notfound, taxid_mapping)
    print(f"NCBI mapped taxids: {len(new_taxid_map)}")
    new_taxids = sorted(set([new for old, new in new_taxid_map.items()]))
    print(f"taxids need to be queried: {len(new_taxids)}")
    new_taxon_q = {t["query"]: t for t in get_taxon_info(new_taxids) if "notfound" not in t}
    new_taxon_info = {
        old: new_taxon_q[new] for old, new in new_taxid_map.items() if new in new_taxon_q
    }
    print(f"Biothings mapped new taxids: {len(new_taxon_info)}")
    taxon_info.update(new_taxon_info)
    print(f"Merged taxon info: {len(taxon_info)}")
    save_pickle(taxon_info, "gmmad2_microbe_disease_taxon_info.pkl")

    # cache taxon descriptions from NCIT
    taxon_names = get_taxon_names(taxon_info)
    print(f"Total unique taxon names: {len(taxon_names)}")
    ncit_descriptions = get_ncit_taxon_description(taxon_names)
    print(f"NCIT descriptions found for {len(ncit_descriptions)} taxon names.")
    taxon_info_w_descr = add_description2taxon_info(taxon_info, ncit_descriptions)
    save_pickle(taxon_info_w_descr, "gmmad2_microbe_disease_taxon_info_w_descr.pkl")


def get_node_info(f_path: str | os.PathLike) -> Iterator[dict]:
    """generates node dictionaries and parse through the disease_species.csv file
    This function reads data, processes taxonomic information,
    and generates subject, object and association nodes,
    representing diseases and microbes, as well as their relationships.

    :param f_path: Path to the disease_species.csv file
    :return: An iterator of dictionaries containing node information.
    """
    if not os.path.exists(
        Path(os.path.join("cache", "gmmad2_microbe_disease_taxon_info_w_descr.pkl"))
    ):
        print("Taxon info not found in cache, running cache_data()...")
        cache_data(f_path)

    taxon_info = load_pickle("gmmad2_microbe_disease_taxon_info_w_descr.pkl")
    for line in line_generator_4_midi(f_path):
        # create object node (diseases)
        object_node = {
            "id": f"MESH:{line[1]}",
            "name": line[2].lower(),
            "type": "biolink:Disease",
            "description": line[18],
            "xrefs": {"mesh": f"MESH:{line[1]}"},
        }
        object_node = remove_empty_none_values(object_node)

        # create subject node (microbes)
        taxid = line[5]
        if taxid in taxon_info:
            subject_node = {
                "id": f"NCBITaxon:{taxon_info[taxid].get('_id')}",
                "taxid": taxon_info[taxid].get("_id"),
                "name": taxon_info[taxid].get("scientific_name").lower(),
                "original_name": line[3].lower(),
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

        # association node
        # includes disease and health sample sizes, microbial abundance mean, median, sd, qualifier
        association_node = {
            "predicate": "OrganismalEntityAsAModelOfDiseaseAssociation",
            "control_name": "healthy control",
            "qualifier": line[17].lower(),
            "qualifier_ratio": line[16],
            "disease_sample_size": line[6],
            "disease_abundance_mean": line[7],
            "disease_abundance_median": line[8],
            "disease_abundance_sd": line[9],
            "healthy_sample_size": line[12],
            "healthy_abundance_mean": line[13],
            "healthy_abundance_median": line[14],
            "healthy_abundance_sd": line[15],
            "infores": "GMMAD2-GMrepo",  # original knowledge source
        }
        association_node = remove_empty_none_values(association_node)

        output_dict = {
            "_id": f"{subject_node['id'].split(':')[1]}_OrganismalEntityAsAModelOfDiseaseAssociation_{object_node['id'].split(':')[1]}",
            "association": association_node,
            "object": object_node,
            "subject": subject_node,
        }

        yield output_dict


def load_micro_disease_data(f_path) -> Iterator[dict]:
    """loads and yields microbe-disease data records from disease_species.csv file
    This function constructs the file path to the disease_species.csv file,
    retrieves node information using the `get_node_info` function, and yields each record.

    :return: An iterator of dictionaries containing microbe-disease data.
    """
    assert os.path.exists(f_path), f"The file {file_path} does not exist."

    recs = get_node_info(file_path)
    for rec in recs:
        # exclude records with sample size == 0
        if (
            int(rec["association"].get("disease_sample_size")) != 0
            and int(rec["association"].get("healthy_sample_size")) != 0
        ):
            yield rec


if __name__ == "__main__":
    file_path = os.path.join("downloads", "disease_species.csv")

    data = load_micro_disease_data(file_path)
    recs = [rec for rec in data]
    save_pickle(recs, "gmmad2_microbe_disease.pkl")
    save_json(recs, "gmmad2_microbe_disease.json")

    _ids = [rec["_id"] for rec in recs]
    disease_ids = []

    for obj in recs:
        if "id" in obj["object"]:
            disease_ids.append(obj["object"]["id"])

    unique_ids = set(_ids)
    unique_disease_ids = set(disease_ids)

    print(f"total records: {len(_ids)}")
    print(f"total records without duplicates: {len(unique_ids)}")
    print(f"Number of unique disease IDs: {len(unique_disease_ids)}")

    type_list = [obj["subject"]["organism_type"] for obj in recs]
    type_counts = Counter(type_list)  # count organism types
    print(type_counts)
    rank_list = [obj["subject"]["rank"] for obj in recs if "rank" in obj["subject"]]
    rank_counts = Counter(rank_list)
    print(rank_counts)
