import csv
import os
import uuid
from collections.abc import Iterator

import biothings_client

"""
{
 0: 'id',   # '1001'
 1: 'g_meta',   # 'meta195'
 2: 'compound', # '1,4-Dihydronicotinamide adenine dinucleotide'
 3: 'pubchem_id',   # '439153'
 4: 'formula',  # 'C21H29N7O14P2'
 5: 'kegg_id',  # 'C00004' or 'not available'
 6: 'HMDBID',   # 'HMDB0001487' or 'not available'
 7: 'drug_id',  # 'DB00157' or 'Not available'
 8: 'drug_name',    # 'NADH' or 'Not available'
 9: 'Origin',   # 'Microbiota; Food related; Drug related' or 'Unknown'
 10: 'smiles_sequence', # 'C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O'
 11: 'gene_id', # 'g5070'
 12: 'gene',    # 'NSDHL' or 'Unknow'
 13: 'ensembl_id',  # 'ENSG00000147383' or 'Unknow'
 14: 'NCBI',    # '50814' or 'Unknow'
 15: 'HGNC',    # '13398' or 'Unknow'
 16: 'UniProt', # 'Q15738' or 'Unknow'
 17: 'protein_size',    # '373' or 'Unknow'
 18: 'annonation',  # "The protein encoded by this gene is localized in the endoplasmic reticulum ...' or 'Unknow'
 19: 'score',   # '0.95' or 'Unknow'
 20: 'alteration',  # ['Unknown', 'elevated', 'reduced', 'target', 'Inhibitor', 'Activator']
 21: 'PMID',   # 'Not available' or 'Unknow'
 22: 'source'   # 'stitch, drugbank' or 'Unknow'
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

        # add gene id via a hierarchical order: 1.entrezgene, 2.ensemblgene, 3.hgnc, and 4.uniportkb
        if "entrezgene" in object_node:
            assign_to_xrefs_if_available(object_node, "ensemblgene", line[13])
        else:
            assign_col_val_if_available(object_node, "ensemblgene", line[13])
        if "entrezgene" not in object_node and "ensemblgene" not in object_node:
            assign_col_val_if_available(object_node, "hgnc", line[15], int)
        else:
            assign_to_xrefs_if_available(object_node, "hgnc", line[15], int)
        if (
            "entrezgene" not in object_node
            and "ensemblgene" not in object_node
            and "hgnc" not in object_node
        ):
            assign_col_val_if_available(object_node, "uniprotkb", line[16])
        else:
            assign_to_xrefs_if_available(object_node, "uniprotkb", line[16])

        # assign ids via a hierarchical order: 1.entrezgene, 2.ensemblgene, 3.hgnc, and 4.uniprotkb
        if "entrezgene" in object_node:
            object_node["id"] = f"NCBIGene:{object_node['entrezgene']}"
        elif "ensemblgene" in object_node:
            object_node["id"] = f"ENSEMBL:{object_node['ensemblgene']}"
        elif "hgnc" in object_node:
            object_node["id"] = f"HGNC:{object_node['hgnc']}"
        else:
            object_node["id"] = f"UniProtKG:{object_node['uniprotkb']}"

        # assign gene names by using biothings_client
        for key in ("entrezgene", "ensemblgene", "uniprotkb"):
            if key in object_node and object_node[key] in gene_name:
                object_node["name"] = gene_name[object_node[key]].get("name")
                break

        # add gene summary to the object_node
        if line[18]:
            object_node["summary"] = line[18]

        # convert entrezgene to integers
        if "entrezgene" in object_node:
            object_node["entrezgene"] = int(object_node["entrezgene"])

        # change object_node type to biolink:Protein if there is only uniprot exists
        if "uniportkb" in object_node:
            object_node["type"] = "biolink:Protein"

        # create subject node (metabolites)
        subject_node = {
            "id": None,
            "name": line[2].lower(),
            "type": "biolink:SmallMolecule",
        }

        assign_col_val_if_available(subject_node, "pubchem_cid", line[3], int)
        assign_col_val_if_available(subject_node, "drug_name", line[8].lower())
        assign_col_val_if_available(subject_node, "chemical_formula", line[4])
        assign_col_val_if_available(subject_node, "smiles", line[10])

        # add chemicals via a hierarchical order: 1.pubchem_cid, 2.kegg_compound, 3.hmdb, and 4.drugbank
        if "pubchem_cid" in subject_node:
            assign_to_xrefs_if_available(subject_node, "kegg_compound", line[5])
        else:
            assign_col_val_if_available(subject_node, "kegg_compound", line[5])
        if "pubchem_cid" not in subject_node and "kegg_compound" not in subject_node:
            assign_col_val_if_available(subject_node, "hmdb", line[6])
        else:
            assign_to_xrefs_if_available(subject_node, "hmdb", line[6])
        if (
            "pubchem_cid" not in subject_node
            and "kegg_compound" not in subject_node
            and "drugbank" not in subject_node
        ):
            assign_col_val_if_available(subject_node, "drugbank", line[7])
        else:
            assign_to_xrefs_if_available(subject_node, "drugbank", line[7])

        # assign chemical id via a hierarchical order: 1.pubchem_cid, and 2.kegg_compound
        if "pubchem_cid" in subject_node:
            subject_node["id"] = f"PUBCHEM.COMPOUND:{subject_node['pubchem_cid']}"
        elif "kegg_compound" in subject_node:
            subject_node["id"] = f"KEGG.COMPOUND:{subject_node['kegg_compound']}"
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
                "elevated", "increased"
            )
            association_node["category"] = "biolink:ChemicalAffectsGeneAssociation"
        if "reduced" in association_node.get("qualifier", ""):
            association_node["qualifier"] = association_node["qualifier"].replace(
                "reduced", "decreased"
            )
            association_node["category"] = "biolink:ChemicalAffectsGeneAssociation"

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
