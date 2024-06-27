import csv
import os
from collections.abc import Iterator

import biothings_client
import re

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


def line_generator(in_file):
    with open(in_file) as in_f:
        reader = csv.reader(in_f)
        next(reader)
        for line in reader:
            yield line


def assign_col_val_if_available(node, key, val, transform=None):
    if val and val != "not available":
        node[key] = transform(val) if transform else val


def assign_to_xrefs_if_available(node, key, val, transform=None):
    if val and val != "not available":
        if "xrefs" not in node:
            node["xrefs"] = {}

        node["xrefs"][key] = transform(val) if transform else val


def get_taxon_info(file_path) -> Iterator[dict]:
    taxids = [line[9] for line in line_generator(file_path)]
    taxids = set(taxids)
    t = biothings_client.get_client("taxon")
    taxon_info = t.gettaxa(taxids, fields=["scientific_name", "parent_taxid", "lineage", "rank"])
    return taxon_info


def get_gene_name(gene_ids):
    gene_ids = set(gene_ids)
    t = biothings_client.get_client("gene")
    gene_names = t.querymany(gene_ids, scopes=["entrezgene", "ensembl.gene", "uniprot"], fields=["name"])
    return gene_names


path = os.getcwd()
file_path = os.path.join(path, "data", "meta_gene_net.csv")
assert os.path.exists(file_path), f"The file {file_path} does not exist."

entrezgene_ids = [line[14] for line in line_generator(file_path) if "not available" not in line[14]]
ensembl_ids = [line[13] for line in line_generator(file_path) if "not available" in line[14] and "not available" not in line[13]]
uniprot_ids = [line[16] for line in line_generator(file_path) if "not available" in line[14] and "not available" in line[13]]
gene_ids = entrezgene_ids + ensembl_ids + uniprot_ids

# get gene name using get_gene_name() function
gene_name = {
    gene_id["query"]: gene_id
    for gene_id in get_gene_name(gene_ids)
    if "notfound" not in gene_id.keys() and "name" in gene_id.keys()
}

# taxon info using get_taxon_info() function
taxon_info = {
        taxon["query"]: taxon
        for taxon in get_taxon_info(file_path)
        if "notfound" not in taxon.keys()
}

# parse the data
for line in line_generator(file_path):
    # create object node (genes)
    object_node = {
        "id": None,
        "symbol": line[12],
        "type": "biolink:Gene"
    }

    assign_col_val_if_available(object_node, "entrezgene", line[14])
    assign_col_val_if_available(object_node, "protein_size", line[17], int)

    # add gene id via a hierarchical order: 1.entrezgene, 2.ensembl, 3.hgnc, and 4.uniportkb
    if "entrezgene" in object_node:
        assign_to_xrefs_if_available(object_node, "ensembl", line[13])
    else:
        assign_col_val_if_available(object_node, "ensembl", line[13])
    if "entrezgene" not in object_node and "ensembl" not in object_node:
        assign_col_val_if_available(object_node, "hgnc", line[15])
    else:
        assign_to_xrefs_if_available(object_node, "hgnc", line[15])
    if "entrezgene" not in object_node and "ensembl" not in object_node and "hgnc" not in object_node:
        assign_col_val_if_available(object_node, "uniprotkb", line[16])
    else:
        assign_to_xrefs_if_available(object_node, "uniprotkb", line[15])

    # assign ids via a hierarchical order: 1.entrezgene, 2.ensembl, 3.hgnc, and 4.uniprotkb
    if "entrezgene" in object_node:
        object_node["id"] = f"NCBIGene: {object_node['entrezgene']}"
    elif "ensembl" in object_node:
        object_node["id"] = f"ENSEMBL: {object_node['ensembl']}"
    elif "hgnc" in object_node:
        object_node["id"] = f"HGNC: {object_node['hgnc']}"
    else:
        object_node["id"] = f"UniProtKG: {object_node['uniprotkb']}"

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

    print(object_node)







