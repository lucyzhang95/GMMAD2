import json
import os
import pickle
import tarfile
from collections.abc import Iterator

import biothings_client
import requests
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


def get_taxon_info(f_path) -> list:
    """retrieves taxonomic information for a given list of taxon IDs from disease_species.csv

    This function reads taxon IDs, removes duplicates, and queries taxonomic info from biothings_client
    to retrieve detailed taxonomic information including scientific name, parent taxid, lineage, and rank.

    :param f_path: Path to disease_species.csv containing the taxids.
    :return: A list of dictionaries containing taxonomic information.
    """
    taxids = [line[5] for line in line_generator_4_midi(f_path)]
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


def get_current_taxid(old_taxids, merged_mapping):
    taxid_mapping = {}
    for old_taxid in old_taxids:
        taxid_mapping[old_taxid] = merged_mapping[old_taxid]
    return taxid_mapping


def get_taxon_names(taxon_info: dict):
    """Extracts biothings names from the taxon_info dictionary."""
    taxon_names = set()


def get_ncit_taxon_description(taxon_names):
    """

    :param taxon_names:
    :return:
    {'serratia': {'description':
    'A genus of small motile peritrichous bacteria in the Enterobacteriacaea family
    consisting of Gram-negative rods.
    [NCIT]',
    'xrefs': {'ncit': 'C86010', }} ...}
    """
    NCIT_API_KEY = os.getenv("NCIT_API_KEY")
    search_url = "https://data.bioontology.org/search"
    taxon_names = set(taxon_names)
    mapping_result = {}
    for name in taxon_names:
        params = {
            "q": name,
            "ontologies": "NCIT",
            "apikey": NCIT_API_KEY,
        }
        response = requests.get(search_url, params=params)
        data = response.json()
        for result in data.get("collection", []):
            if result:
                ncit_output = {
                    "name": result.get("prefLabel").lower(),
                    "description": f"{result.get('definition')[0]} [NCIT]"
                    if "definition" in result
                    else "",
                    "xrefs": {"ncit": result.get("@id").split("#")[1]},
                }
                if ncit_output["name"] == name:
                    mapping_result[name] = ncit_output
                    del mapping_result[name]["name"]
    return mapping_result


def add_description2taxon_info(taxon_info: dict, descriptions: dict) -> dict:
    for _, info in taxon_info.items():
        name = info.get("name")
        if name in descriptions:
            info.update(descriptions[name])
        else:
            pass

    return taxon_info


def get_full_taxon_info(mapped_taxon_names: dict, taxon_info: dict) -> dict:
    """A dictionary of taxon metadata keyed by taxon name,
    using taxid lookup from `taxon_info`.

    :param mapped_taxon_names: Mapping of taxon name to a dict with taxon taxid, lineage, rank...
    {'clostridia propionicum': {'taxid': 28446, 'mapping_tool': 'manual'} ...}
    :param taxon_info: Mapping of taxid (as str) to a dictionary with mapped taxon names and taxid.
    '28450': {'id': 'taxid:28450', 'taxid': 28450, 'name': 'burkholderia pseudomallei', 'parent_taxid': 111527, 'lineage': [28450, 111527, 32008, 119060, 80840, 28216, 1224, 3379134, 2, 131567, 1], 'rank': 'species', 'description': 'A species of aerobic, Gram-negative, rod shaped bacteria assigned to the phylum Proteobacteria.
    This species is motile, non-spore forming, oxidase and catalase positive and indole negative.
    B. pseudomallei is found in contaminated soil, water,
    and produce and causes melioidosis in humans with the highest rate of disease occurring in southeast Asia.
    [NCIT]', 'xrefs': {'ncit': 'C86010' }} ...}

    :return:

    """
    full_taxon = {
        name: taxon_info[str(taxon_entry["taxid"])]
        for name, taxon_entry in mapped_taxon_names.items()
        if str(taxon_entry["taxid"]) in taxon_info
    }
    return full_taxon


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


def get_node_info(f_path: str | os.PathLike) -> Iterator[dict]:
    """generates node dictionaries and parse through the disease_species.csv file
    This function reads data, processes taxonomic information,
    and generates subject, object and association nodes,
    representing diseases and microbes, as well as their relationships.

    :param f_path: Path to the disease_species.csv file
    :return: An iterator of dictionaries containing node information.
    """
    taxon_info = {
        int(taxon["query"]): taxon
        for taxon in get_taxon_info(f_path)
        if "notfound" not in taxon.keys()
    }
    for line in line_generator_4_midi(f_path):
        # create object node (diseases)
        object_node = {
            "id": f"MESH:{line[0]}",
            "name": line[1].lower(),
            "mesh": line[0],
            "type": "biolink:Disease",
            "description": line[17],
        }

        # create subject node (microbes)
        taxid = int(line[4])
        subject_node = {
            "id": f"taxid:{taxid}",
            "taxid": taxid,
            "name": None,
            "original_name": line[2].lower(),
            "type": "biolink:OrganismTaxon",
        }
        if subject_node["taxid"] in taxon_info:
            subject_node["name"] = taxon_info[subject_node["taxid"]]["scientific_name"]
            subject_node["parent_taxid"] = taxon_info[subject_node["taxid"]]["parent_taxid"]
            subject_node["lineage"] = taxon_info[subject_node["taxid"]]["lineage"]
            subject_node["rank"] = taxon_info[subject_node["taxid"]]["rank"]
            subject_node["organism_type"] = get_organism_type(subject_node)

        # association node
        # includes disease and health sample sizes, microbial abundance mean, median, sd, qualifier
        association_node = {
            "predicate": "OrganismalEntityAsAModelOfDiseaseAssociation",
            "control_name": "healthy control",
            "qualifier": line[16].lower(),
            "qualifier_ratio": line[15],
            "disease_sample_size": line[5],
            "disease_abundance_mean": line[6],
            "disease_abundance_median": line[7],
            "disease_abundance_sd": line[8],
            "healthy_sample_size": line[11],
            "healthy_abundance_mean": line[12],
            "healthy_abundance_median": line[13],
            "healthy_abundance_sd": line[14],
            "infores": "GMMAD2-GMrepo",  # knowledge source
        }

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
    _ids = [obj["_id"] for obj in data]
    disease_ids = []
    disease_names = []

    for obj in data:
        print(obj)
        if "id" in obj["object"]:
            disease_ids.append(obj["object"]["id"])
        else:
            disease_names.append(obj["object"]["name"])

    unique_ids = set(_ids)
    unique_disease_ids = set(disease_ids)
    unique_disease_names = set(disease_names)

    print(f"total records: {len(_ids)}")
    print(f"total records without duplicates: {len(unique_ids)}")
    print(f"Number of unique disease IDs: {len(unique_disease_ids)}")
    print(f"Number of disease names: {len(unique_disease_names)}")

    from collections import Counter

    type_list = [obj["subject"]["type"] for obj in data]
    type_counts = Counter(type_list)  # count microorganism type

    for value, count in type_counts.items():
        print(f"{value}: {count}")

    rank_list = [obj["subject"]["rank"] for obj in data if "rank" in obj["subject"]]
    rank_counts = Counter(rank_list)
    for value, count in rank_counts.items():
        print(f"{value}: {count}")
