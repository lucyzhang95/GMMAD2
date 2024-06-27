import csv
import os
import re
import uuid
from collections.abc import Iterator

import biothings_client

"""
column names with index:
{
    0: 'id', internal id
    1: 'g_meta', internal id
    2: 'compound', name
    3: 'pubchem_id', not available exists
    4: 'formula', chemical formula
    5: 'kegg_id', 
    6: 'HMDBID', not available exists
    7: 'drug_id', not available exists
    8: 'drug_name', not available exists
    9: 'Origin', list of items
    10: 'smiles_sequence', 
    11: 'gene_id', internal id
    12: 'gene', symbol
    13: 'ensembl_id', 
    14: 'NCBI',
    15: 'HGNC',
    16: 'UniProt',
    17: 'protein_size', 
    18: 'annonation', description
    19: 'score', ?
    20: 'alteration', qualifier, Unknown exists
    21: 'PMID', not available exists
    22: 'source', infores
}
"""


def line_generator(in_file: str | os.PathLike) -> Iterator[list]:
    """generates lines from a CSV file, yielding each line as a list of strings
    This function opens the specified CSV file, skips the header row, and yields each subsequent line as a list of strings.

    :param in_file: The path to the CSV file.
    :return: An iterator that yields each line of the CSV file as a list of strings.
    """
    with open(in_file) as in_f:
        reader = csv.reader(in_f)
        next(reader)
        for line in reader:
            yield line


def assign_col_val_if_available(node: dict, key: str, val: str | int | float, transform=None):
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


def assign_to_xrefs_if_available(node: dict, key: str, val: str | int | float, transform=None):
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


def get_gene_name(gene_ids: list) -> list:
    """retrieves gene names for a given list of gene IDs using biothings_client
        This function takes a list of gene IDs, removes any duplicates by converting the list to a set
        Queries a gene database client for the gene names associated with these IDs.
        The IDs are searched across multiple scopes: "entrezgene", "ensembl.gene", and "uniprot".

        :param gene_ids: A list of gene IDs to be queried.
        :type gene_ids: list
        :return: A list of dictionaries containing the gene names and associated information.
        """
    gene_ids = set(gene_ids)
    t = biothings_client.get_client("gene")
    gene_names = t.querymany(
        gene_ids, scopes=["entrezgene", "ensembl.gene", "uniprot"], fields=["name"]
    )
    return gene_names


def get_node_info(file_path: str | os.PathLike) -> Iterator[dict]:
    """generates node dictionaries from meta_gene_net.csv file
    This function reads gene and metabolite data and processes it.
    It generates object, subject, and association nodes.

    :param file_path: path to meta_gene_net.csv file
    :return: An iterator of dictionaries containing node information.
    """

    # gather gene ids from the file
    entrezgene_ids = [
        line[14] for line in line_generator(file_path) if "not available" not in line[14]
    ]
    ensembl_ids = [
        line[13]
        for line in line_generator(file_path)
        if "not available" in line[14] and "not available" not in line[13]
    ]
    uniprot_ids = [
        line[16]
        for line in line_generator(file_path)
        if "not available" in line[14] and "not available" in line[13]
    ]
    gene_ids = entrezgene_ids + ensembl_ids + uniprot_ids

    # get gene name using get_gene_name() function
    gene_name = {
        gene_id["query"]: gene_id
        for gene_id in get_gene_name(gene_ids)
        if "notfound" not in gene_id.keys() and "name" in gene_id.keys()
    }

    # parse the data
    for line in line_generator(file_path):
        # create object node (genes)
        object_node = {"id": None, "symbol": line[12], "type": "biolink:Gene"}

        assign_col_val_if_available(object_node, "entrezgene", line[14])
        assign_col_val_if_available(object_node, "protein_size", line[17], int)

        # add gene id via a hierarchical order: 1.entrezgene, 2.ensembl, 3.hgnc, and 4.uniportkb
        if "entrezgene" in object_node:
            assign_to_xrefs_if_available(object_node, "ensembl", line[13])
        else:
            assign_col_val_if_available(object_node, "ensembl", line[13])
        if "entrezgene" not in object_node and "ensembl" not in object_node:
            assign_col_val_if_available(object_node, "hgnc", line[15], int)
        else:
            assign_to_xrefs_if_available(object_node, "hgnc", line[15], int)
        if (
            "entrezgene" not in object_node
            and "ensembl" not in object_node
            and "hgnc" not in object_node
        ):
            assign_col_val_if_available(object_node, "uniprotkb", line[16])
        else:
            assign_to_xrefs_if_available(object_node, "uniprotkb", line[16])

        # assign ids via a hierarchical order: 1.entrezgene, 2.ensembl, 3.hgnc, and 4.uniprotkb
        if "entrezgene" in object_node:
            object_node["id"] = f"NCBIGene:{object_node['entrezgene']}"
        elif "ensembl" in object_node:
            object_node["id"] = f"ENSEMBL:{object_node['ensembl']}"
        elif "hgnc" in object_node:
            object_node["id"] = f"HGNC:{object_node['hgnc']}"
        else:
            object_node["id"] = f"UniProtKG:{object_node['uniprotkb']}"

        # assign gene names by using biothings_client
        if "entrezgene" in object_node and object_node["entrezgene"] in gene_name:
            object_node["name"] = gene_name[object_node["entrezgene"]]["name"]
        elif "ensembl" in object_node and object_node["ensembl"] in gene_name:
            object_node["name"] = gene_name[object_node["ensembl"]].get("name")
        elif "uniprotkb" in object_node and object_node["uniprotkb"] in gene_name:
            object_node["name"] = gene_name[object_node["uniprotkb"]]["name"]

        # divide annotation `line[18]` to description and reference
        # some entries have both, and some entries only have description
        descr_match = re.search(r"(.+?)\s*\[(.+)\]$", line[18])
        if "[" in line[18]:
            if descr_match:
                object_node["description"] = descr_match.group(1).strip()
                object_node["ref"] = descr_match.group(2).strip()
        else:
            object_node["description"] = line[18].strip()

        # convert entrezgene to integers
        if "entrezgene" in object_node:
            object_node["entrezgene"] = int(object_node["entrezgene"])

        # create subject node (metabolites)
        subject_node = {
            "id": None,
            "name": line[2].lower(),
            "type": "biolink:SmallMolecule",
        }

        assign_col_val_if_available(subject_node, "pubchem_cid", line[3], int)
        assign_col_val_if_available(subject_node, "drugbank", line[7])
        assign_col_val_if_available(subject_node, "drug_name", line[8].lower())
        assign_col_val_if_available(subject_node, "chemical_formula", line[4])
        assign_col_val_if_available(subject_node, "smiles", line[10])

        # add chemicals via a hierarchical order: 1.pubchem_cid, 2.kegg, 3.hmdb, and 4.drugbank
        if "pubchem_cid" in subject_node:
            assign_to_xrefs_if_available(subject_node, "kegg", line[5])
        else:
            assign_col_val_if_available(subject_node, "kegg", line[5])
        if "pubchem_cid" not in subject_node and "kegg" not in subject_node:
            assign_col_val_if_available(subject_node, "hmdb", line[6])
        else:
            assign_to_xrefs_if_available(subject_node, "hmdb", line[6])

        # assign chemical id via a hierarchical order: 1.pubchem_cid, and 2.kegg
        if "pubchem_cid" in subject_node:
            subject_node["id"] = f"PUBCHEM.COMPOUND:{subject_node['pubchem_cid']}"
        else:
            subject_node["id"] = str(uuid.uuid4())

        # association node has the qualifier, reference and source of metabolites
        association_node = {"predicate": "biolink:associated_with"}

        assign_col_val_if_available(association_node, "score", line[19], float)
        assign_col_val_if_available(association_node, "pmid", line[21], int)

        if line[9] and line[9] != "Unknown":
            association_node["sources"] = [src.strip().lower() for src in line[9].split(";")]
        if line[22] and line[22] != "Unknown":
            association_node["infores"] = [src.strip().lower() for src in line[22].split(",")]
        if line[20] and line[20] != "Unknown":
            association_node["qualifier"] = line[20].lower()
        if "elevated" in association_node.get("qualifier", ""):
            association_node["qualifier"] = association_node["qualifier"].replace(
                "elevated", "increase"
            )
        if "reduced" in association_node.get("qualifier", ""):
            association_node["qualifier"] = association_node["qualifier"].replace(
                "reduced", "decrease"
            )

        # combine all the nodes together
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
        else:
            output_dict["_id"] = (
                f"{subject_node['id']}_associated_with_{object_node['id'].split(':')[1].strip()}"
            )
        yield output_dict


def load_meta_gene_data() -> Iterator[dict]:
    """loads and yields unique meta gene data records from meta_gene_net.csv file.

    This function constructs the file path to meta_gene_net.csv file,
    retrieves node information using the `get_node_info` function,
    and yields unique records based on `_id`.

    :return: An iterator of unique dictionaries containing meta gene data.
    """
    path = os.getcwd()
    file_path = os.path.join(path, "data", "meta_gene_net.csv")
    assert os.path.exists(file_path), f"The file {file_path} does not exist."

    dup_ids = set()
    recs = get_node_info(file_path)
    for rec in recs:
        if rec["_id"] not in dup_ids:
            dup_ids.add(rec["_id"])
            yield rec


# if __name__ == "__main__":
#     _ids = []
#     meta_gene_data = load_meta_gene_data()
#     for obj in meta_gene_data:
#         print(obj)
#         _ids.append(obj["_id"])
#     print(f"total records: {len(_ids)}")
#     print(f"total records without duplicates: {len(set(_ids))}")
