import os
import time
import uuid

import pandas as pd
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


for ls in line_generator(file_path):
    object_node = {
        "id": None,
        "name": ls[5]
    }
    if ls[6] and ls[6] == "not available":
        if ls[8] and ls[8] != "not available":
            object_node["id"] = f"KEGG.COMPOUND:{ls[8]}"
            object_node["kegg"] = ls[8]
        elif ls[19] and ls[19] != "not available":
            object_node["id"] = f"HMDB:{ls[19]}"
            object_node["hmdb"] = ls[19]
        else:
            object_node["id"] = int(uuid.uuid4())
    else:
        object_node["id"] = f"PUBCHEM.COMPOUND:{int(ls[6])}"
        object_node["pubchem_cid"] = int(ls[6])

    if ls[7] and ls[7] != "not available":
        object_node["chemical_formula"] = ls[7]
    if ls[18] and ls[18] != "not available":
        object_node["smiles"] = ls[18]

    print(object_node)





