import os
import uuid

import csv
import biothings_client
from collections.abc import Iterator

"""
column names with index: {0: 'id', 1: 'g_micro', 2: 'organism', 3: 'g_meta', 4: 'metabolic', 5: 'pubchem_compound', 6: 'pubchem_id', 7: 'formula', 8: 'kegg_id', 9: 'tax_id', 10: 'phylum', 11: 'class', 12: 'order', 13: 'family', 14: 'genus', 15: 'species', 16: 'species_id', 17: 'source', 18: 'smiles_sequence', 19: 'HMDBID', 20: 'Origin'}
"""

path = os.getcwd()
file_path = os.path.join(path, "data", "micro_metabolic.csv")
assert os.path.exists(file_path), f"The file {file_path} does not exist."


def line_generator(file_path):
    with open(file_path, "r") as in_f:
        reader = csv.reader(in_f)
        # skip the header
        next(reader)
        for row in reader:
            yield row


def assign_col_val_if_available(node, key, val, transform=None):
    if val and val != "not available":
        node[key] = transform(val) if transform else val


def create_object_node():
    for ls in line_generator(file_path):
        object_node = {
            "id": None,
            "name": ls[5].lower()
        }

        assign_col_val_if_available(object_node, "pubchem_cid", ls[6], int)
        assign_col_val_if_available(object_node, "kegg", ls[8])
        assign_col_val_if_available(object_node, "hmdb", ls[19])
        assign_col_val_if_available(object_node, "chemical_formula", ls[7])
        assign_col_val_if_available(object_node, "smiles", ls[18])

        if "pubchem_cid" in object_node:
            object_node["id"] = f"PUBCHEM.COMPOUND:{object_node['pubchem_cid']}"
        elif "kegg" in object_node:
            object_node["id"] = f"KEGG.COMPOUND:{object_node['kegg']}"
        elif "hmdb" in object_node:
            object_node["id"] = f"HMDB:{object_node['hmdb']}"
        else:
            object_node["id"] = str(uuid.uuid4())

        yield object_node







