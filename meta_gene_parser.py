import csv
import os
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
    11: 'gene_id', symbol
    12: 'gene', name
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




