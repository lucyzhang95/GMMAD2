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


def line_generator(in_file):
    with open(in_file, "r") as in_f:
        reader = csv.reader(in_f)
        next(reader)
        for line in reader:
            yield line


def get_taxon_info(file_path) -> list:
    taxids = [line[4] for line in line_generator(file_path)]
    taxids = set(taxids)
    t = biothings_client.get_client("taxon")
    taxon_info = t.gettaxa(taxids, fields=["scientific_name", "parent_taxid", "lineage", "rank"])
    return taxon_info


def get_node_info(file_path):
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
            "name": line[2].lower(),
            "type": "biolink:OrganismalEntity",
        }
        if subject_node["taxid"] in taxon_info:
            subject_node["scientific_name"] = taxon_info[subject_node["taxid"]]["scientific_name"]
            subject_node["parent_taxid"] = taxon_info[subject_node["taxid"]]["parent_taxid"]
            subject_node["lineage"] = taxon_info[subject_node["taxid"]]["lineage"]
            subject_node["rank"] = taxon_info[subject_node["taxid"]]["rank"]

        # categorize subject microbial super kingdom type
        if "lineage" in subject_node:
            if 2 in subject_node["lineage"]:
                subject_node["type"] = "biolink:Bacterium"
            elif 10239 in subject_node["lineage"]:
                subject_node["type"] = "biolink:Virus"
            elif 4751 in subject_node["lineage"]:
                subject_node["type"] = "biolink:Fungus"
            else:
                subject_node["type"] = "biolink:OrganismalEntity"

        # create association node
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
        }

        output_dict = {
            "_id": f"{subject_node['id'].split(':')[1]}_OrganismalEntityAsAModelOfDiseaseAssociation_{object_node['id'].split(':')[1]}",
            "association": association_node,
            "object": object_node,
            "subject": subject_node,
        }

        yield output_dict


def load_micro_disease_data():
    path = os.getcwd()
    file_path = os.path.join(path, "data", "disease_species.csv")
    assert os.path.exists(file_path), f"The file {file_path} does not exist."

    recs = get_node_info(file_path)
    for rec in recs:
        yield rec


# if __name__ == "__main__":
#     from collections import Counter
#
#     data = load_micro_disease_data()
#
#     type_list = [obj["subject"]["type"] for obj in data]
#     type_counts = Counter(type_list)
#
#     for value, count in type_counts.items():
#         print(f"{value}: {count}")
#
#     rank_list = [obj["subject"]["rank"] for obj in data if "rank" in obj["subject"]]
#     rank_counts = Counter(rank_list)
#     for value, count in rank_counts.items():
#         print(f"{value}: {count}")
#
#     _ids = []
#     for obj in data:
#         # print(obj)
#         _ids.append(obj["_id"])
#     print(f"total records: {len(_ids)}")
#     print(f"total records without duplicates: {len(set(_ids))}")
