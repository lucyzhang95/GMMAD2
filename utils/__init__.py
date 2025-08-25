#%%
import os
import pandas as pd
from ete3 import NCBITaxa
import biothings_client as bt
import tarfile
import csv
import gzip

from micro_disease_parser import line_generator_4_midi, get_taxon_info, load_merged_from_tar, get_current_taxid
from micro_meta_parser import get_bigg_metabolite_mapping
from data_utils import line_generator, check_line_fields, print_misalignment_report, check_missing_data
#%%
os.getcwd()
#%% md
## Microbe-Disease
#%%
micro_disease_path = os.path.join("downloads", "disease_species.csv")
micro_disease = [line for line in line_generator(micro_disease_path)]
#%%
header = micro_disease[0]
print(header)
print(len(header))
#%%
fields = check_line_fields(micro_disease)
print_misalignment_report(fields, len(micro_disease))
#%%
micro_disease = [line for line in line_generator_4_midi(micro_disease_path, skip_header=False)]
for i, line in enumerate(micro_disease):
    if i == 68162:
        print(len(line))
        print(line)
        break
#%%
header = micro_disease[0]
data_rows = micro_disease[1:]
micro_disease_df = pd.DataFrame(data_rows, columns=header)
#%%
print(micro_disease_df.columns)
#%%
check_missing_data(micro_disease_df)
#%%
taxids = [line[5] for line in line_generator_4_midi(micro_disease_path)]
print(len(set(taxids)))
#%%
notfound = [
    taxon["query"]
    for taxon in get_taxon_info(micro_disease_path)
    if "notfound" in taxon.keys()
]
#%%
print(len(notfound))
print(len(set(notfound)))
#%%
mapping = load_merged_from_tar("taxdump.tar.gz")
mapped_taxid = get_current_taxid(notfound, mapping)
#%%
if "194866" in mapping:
    print("194866 is in the mapping")
    print(f"current taxid for 194866: {mapping['194866']}")
else:
    print("194866 is not in the mapping")
#%%
len(mapped_taxid)
#%%
new_taxids = [new for old, new in mapped_taxid.items()]
print(len(new_taxids))
print(len(set(new_taxids)))
#%%
print(set(new_taxids))
#%%
still_notfound = [taxid for taxid in notfound if taxid not in mapped_taxid.keys()]
print(len(still_notfound))
#%%
new_taxid_mapped = get_taxon_info(new_taxids)
new_taxid_mapped
#%%
_ids = [taxon["query"] for taxon in new_taxid_mapped]
still_notfound = [taxid for taxid in new_taxids if taxid not in _ids]
#%%
len(set(_ids))
#%%
taxids = sorted([line[5] for line in line_generator_4_midi(micro_disease_path)])
print(len(taxids))
print(len(set(taxids)))
#%%

#%% md
## Microbe-Metabolite
#%%
micro_meta_df = pd.read_csv(os.path.join("downloads", "micro_metabolic.csv"), low_memory=False)
#%%
micro_meta_col_map = dict(enumerate(micro_meta_df.columns))
micro_meta_col_map
#%%
micro_meta = [line for line in line_generator(os.path.join("downloads", "micro_metabolic.csv"), skip_header=False)]
#%%
fields = check_line_fields(micro_meta)
print_misalignment_report(fields, len(micro_meta))
#%%
check_missing_data(micro_meta_df)
#%%
no_chem_id = []
for line in micro_meta:
    if "not available" in line[6] and "not available" in line[19]:
        no_chem_id.append(line[4])
print(len(no_chem_id))
#%%
print(len(set(no_chem_id)))
#%%
pubchem_cids = [line[6] for line in micro_meta if "not available" not in line[6]]
print(len(set(pubchem_cids)))
#%%
hmdb_ids = [line[19] for line in micro_meta if "not available" not in line[19] and "not available" in line[6]]
print(len(set(hmdb_ids)))
#%% md

#%%
total_metabolites = [line[4] for line in micro_meta]
print(len(set(total_metabolites)))
#%%

#%%
bigg_mapped = get_bigg_metabolite_mapping(os.path.join("downloads", "bigg_models_metabolites.txt"))
#%%
len(bigg_mapped)
#%%
bigg_mapped_ids = [name.lower() for name in set(no_chem_id) if name.lower() in bigg_mapped]
#%%
len(bigg_mapped_ids)
#%%
set(no_chem_id)
#%%

#%% md
## Metabolite - Gene
#%%
meta_gene_df = pd.read_csv(os.path.join("downloads", "meta_gene_net.csv"), low_memory=False)
#%%
dict(enumerate(meta_gene_df.columns))
#%%
meta_gene_data = [line for line in line_generator(os.path.join("downloads", "meta_gene_net.csv"), skip_header=False)]
#%%
fields = check_line_fields(meta_gene_data)
print_misalignment_report(fields, len(meta_gene_data))
#%%
for i, line in enumerate(meta_gene_data):
    if i == 2:
        print(len(line))
        print(line)
        break
#%%

#%%
print(len(check_missing_data(meta_gene_df)))
check_missing_data(meta_gene_df)
#%%
meta_gene_df["pubchem_id"].unique()
#%%
meta_gene_df["Origin"].unique()
#%%
meta_gene_df["alteration"].unique()
#%%
meta_gene_df["source"].unique()
#%%
meta_gene_df["PMID"].unique()
#%%
