import csv
import os
from collections.abc import Iterator

import biothings_client

"""
column names with index: 
{
    0: 'disease_id', 
    1: 'disease', 
    2: 'organism', 
    3: 'level', 
    4: 'species_id', 
    5: 'disease_samples', 
    6: 'disease_mean', 
    7: 'disease_median', 
    8: 'disease_sd', 
    9: 'health_id', 
    10: 'health', 
    11: 'health_samples', 
    12: 'health_mean', 
    13: 'health_median', 
    14: 'health_sd', 
    15: 'change', 
    16: 'alteration', 
    17: 'disease_info', 
    18: 'phylum', 
    19: 'class', 
    20: 'order', 
    21: 'family', 
    22: 'genus'
}
"""


def line_generator(in_file: str | os.PathLike) -> Iterator[list]:
    """generates lines from a CSV file, yielding each line as a list of strings
    This function opens the specified CSV file, skips the header row, and yields each subsequent line as a list of strings.

    :param in_file: The path to the CSV file.
    :return: An iterator that yields each line of the CSV file as a list of strings.
    """
    with open(in_file, "r") as in_f:
        reader = csv.reader(in_f)
        next(reader)
        for line in reader:
            yield line


def get_taxon_info(file_path) -> list:
    """retrieves taxonomic information for a given list of taxon IDs from disease_species.csv

    This function reads taxon IDs, removes duplicates, and queries taxonomic info from biothings_client
    to retrieve detailed taxonomic information including scientific name, parent taxid, lineage, and rank.

    :param file_path: Path to disease_species.csv containing the taxids.
    :return: A list of dictionaries containing taxonomic information.
    """
    taxids = [line[4] for line in line_generator(file_path)]
    taxids = set(taxids)
    t = biothings_client.get_client("taxon")
    taxon_info = t.gettaxa(taxids, fields=["scientific_name", "parent_taxid", "lineage", "rank"])
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

    for taxid, biolink_type in taxon_map.items():
        if taxid in node.get("lineage", []):
            return biolink_type

    return "Other"


def get_node_info(file_path: str | os.PathLike) -> Iterator[dict]:
    """generates node dictionaries and parse through the disease_species.csv file
    This function reads data, processes taxonomic information,
    and generates subject, object and association nodes,
    representing diseases and microbes, as well as their relationships.

    :param file_path: Path to the disease_species.csv file
    :return: An iterator of dictionaries containing node information.
    """
    taxon_info = {
        int(taxon["query"]): taxon
        for taxon in get_taxon_info(file_path)
        if "notfound" not in taxon.keys()
    }
    for line in line_generator(file_path):
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
