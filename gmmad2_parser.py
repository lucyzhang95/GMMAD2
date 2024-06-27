import csv
import os
import uuid
from collections.abc import Iterator

import biothings_client


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


def get_taxon_info(taxids: list) -> list:
    """retrieves taxonomic information for a given list of taxon IDs from micro_metabolic.csv
    This function reads taxon IDs, removes duplicates, and queries taxonomic info from biothings_client
    to retrieve detailed taxonomic information including scientific name, parent taxid, lineage, and rank.

    :param taxids: a list containing the taxids
    :return: A list of dictionaries containing taxonomic information.
    """
    taxids = set(taxids)
    t = biothings_client.get_client("taxon")
    taxon_info = t.gettaxa(taxids, fields=["scientific_name", "parent_taxid", "lineage", "rank"])
    return taxon_info


# microbe-metabolite parser
def get_micro_meta_node_info(file_path: str | os.PathLike) -> Iterator[dict]:
    """generates node information from micro_metabolic.csv.
    This function reads data from micro_metabolic.csv, processes taxonomic information,
    and generates nodes representing metabolites and microbes, including their relationships.

    :param file_path: Path to the micro_metabolic.csv file
    :return: An iterator of dictionaries containing node information.
    """
    taxids = [line[9] for line in line_generator(file_path)]

    taxon_info = {
        int(taxon["query"]): taxon
        for taxon in get_taxon_info(taxids)
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

        # categorize subject microbial super kingdom type
        if "lineage" in subject_node:
            if 2 in subject_node["lineage"]:
                subject_node["type"] = "biolink:Bacterium"
            elif 10239 in subject_node["lineage"]:
                subject_node["type"] = "biolink:Virus"
            elif 4751 in subject_node["lineage"]:
                subject_node["type"] = "biolink:Fungus"
            elif 2157 in subject_node["lineage"]:
                subject_node["type"] = "biolink:Archaea"
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
    recs = get_micro_meta_node_info(file_path)
    for rec in recs:
        if rec["_id"] not in dup_ids:
            dup_ids.add(rec["_id"])
            yield rec


# microbe-disease parser
def get_micro_disease_node_info(file_path: str | os.PathLike) -> Iterator[dict]:
    """generates node dictionaries and parse through the disease_species.csv file
    This function reads data, processes taxonomic information,
    and generates subject, object and association nodes,
    representing diseases and microbes, as well as their relationships.

    :param file_path: Path to the disease_species.csv file
    :return: An iterator of dictionaries containing node information.
    """
    taxids = [line[4] for line in line_generator(file_path)]

    taxon_info = {
        int(taxon["query"]): taxon
        for taxon in get_taxon_info(taxids)
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
            elif 2157 in subject_node["lineage"]:
                subject_node["type"] = "biolink:Archaea"
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


def load_micro_disease_data() -> Iterator[dict]:
    """loads and yields microbe-disease data records from disease_species.csv file
    This function constructs the file path to the disease_species.csv file,
    retrieves node information using the `get_node_info` function, and yields each record.

    :return: An iterator of dictionaries containing microbe-disease data.
    """
    path = os.getcwd()
    file_path = os.path.join(path, "data", "disease_species.csv")
    assert os.path.exists(file_path), f"The file {file_path} does not exist."

    recs = get_micro_disease_node_info(file_path)
    for rec in recs:
        yield rec


# metabolite-gene parser
def get_meta_gene_node_info(file_path: str | os.PathLike) -> Iterator[dict]:
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
    recs = get_meta_gene_node_info(file_path)
    for rec in recs:
        if rec["_id"] not in dup_ids:
            dup_ids.add(rec["_id"])
            yield rec


def load_data() -> Iterator[dict]:
    """loads and combines all 3 relationships among microbes, diseases, and metabolites in GMMAD2 database
    sorts them by their '_id' field, and yields each document in sorted order.

    This function utilizes the following data loading functions:
    - `load_micro_meta_data()`
    - `load_micro_disease_data()`
    - `load_meta_gene_data()`

    The documents from these sources are combined into a single list,
    sorted by the '_id' field, and then yielded one by one.

    :return: An iterator that yields each sorted document as a dictionary.
    """
    from itertools import chain
    from operator import itemgetter

    doc_list = list(
        chain(
            load_micro_meta_data(),
            load_micro_disease_data(),
            load_meta_gene_data(),
        )
    )
    sorted_docs = sorted(doc_list, key=itemgetter("_id"))

    for doc in sorted_docs:
        yield doc


# if __name__ == "__main__":
#     data_list = []
#     data = load_data()
#     for obj in data:
#         data_list.append(obj["_id"])
#
#     print("total records", len(data_list))
#     print("total records without duplications", len(set(data_list)))
