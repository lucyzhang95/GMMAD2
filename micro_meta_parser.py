import csv
import os
import uuid
from collections.abc import Iterator

import biothings_client

"""
column names with index: 
{
    0: 'id', 
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
    20: 'Origin'
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
        # skip the header
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


def get_taxon_info(file_path: str | os.PathLike) -> list:
    """retrieves taxonomic information for a given list of taxon IDs from micro_metabolic.csv
    This function reads taxon IDs, removes duplicates, and queries taxonomic info from biothings_client
    to retrieve detailed taxonomic information including scientific name, parent taxid, lineage, and rank.

    :param file_path: Path to micro_metabolic.csv containing the taxids.
    :return: A list of dictionaries containing taxonomic information.
    """
    taxids = [line[9] for line in line_generator(file_path)]
    taxids = set(taxids)
    t = biothings_client.get_client("taxon")
    taxon_info = t.gettaxa(taxids, fields=["scientific_name", "parent_taxid", "lineage", "rank"])
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
            assign_to_xrefs_if_available(object_node, "kegg", line[8])
        else:
            assign_col_val_if_available(object_node, "kegg", line[8])
        if "pubchem_cid" not in object_node and "kegg" not in object_node:
            assign_col_val_if_available(object_node, "hmdb", line[19])
        else:
            assign_to_xrefs_if_available(object_node, "hmdb", line[19])

        if "pubchem_cid" in object_node:
            object_node["id"] = f"PUBCHEM.COMPOUND:{object_node['pubchem_cid']}"
        elif "kegg" in object_node:
            object_node["id"] = f"KEGG.COMPOUND:{object_node['kegg']}"
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
            output_dict["_id"] = (
                f"{subject_node['id'].split(':')[1].strip()}_associated_with_{object_node['id'].split(':')[1].strip()}"
            )
        elif ":" not in object_node["id"] and ":" in subject_node["id"]:
            output_dict["_id"] = (
                f"{subject_node['id'].split(':')[1].strip()}_associated_with_{object_node['id']}"
            )
        elif ":" in object_node["id"] and ":" not in subject_node["id"]:
            output_dict["_id"] = (
                f"{subject_node['id']}_associated_with_{object_node['id'].split(':')[1].strip()}"
            )
        else:
            output_dict["_id"] = f"{subject_node['id']}_associated_with_{object_node['id']}"

        yield output_dict


def load_micro_meta_data() -> Iterator[dict]:
    """loads and yields unique microbe-metabolite data records from micro_metabolic.csv file
    This function constructs the file path to the micro_metabolic.csv file,
    retrieves node information using the `get_node_info` function, and yields unique records based on `_id`.

    :return: An iterator of dictionaries containing microbe-metabolite data.
    """
    path = os.getcwd()
    file_path = os.path.join(path, "data", "micro_metabolic.csv")
    assert os.path.exists(file_path), f"The file {file_path} does not exist."

    dup_ids = set()
    recs = get_node_info(file_path)
    for rec in recs:
        if rec["_id"] not in dup_ids:
            dup_ids.add(rec["_id"])
            yield rec


# if __name__ == "__main__":
# from collections import Counter
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
