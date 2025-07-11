import asyncio
import csv
import json
import os
import pickle
import time
import uuid
from collections import Counter
from collections.abc import Iterator
from pathlib import Path
from typing import Dict, List

import aiohttp
import biothings_client as bt
import pandas as pd
import requests
from dotenv import load_dotenv
from tqdm.asyncio import tqdm_asyncio

"""
{
 0: 'id',   # '1001'
 1: 'g_meta',   # 'meta195'
 2: 'compound', # '1,4-Dihydronicotinamide adenine dinucleotide'
 3: 'pubchem_id',   # 'Not available' or '439153'
 4: 'formula',  # 'not available' or 'C21H29N7O14P2'
 5: 'kegg_id',  # 'not available' or 'C00004'
 6: 'HMDBID',   # 'not available' or 'HMDB0001487'
 7: 'drug_id',  # 'Not available' or 'DB00157'
 8: 'drug_name',    # 'Not available' or 'NADH'
 9: 'Origin',   # 'Unknown' or 'Microbiota; Food related; Drug related'
 10: 'smiles_sequence', # 'not available' or 'C1C=CN(C=C1C(=O)N)C2C(C(C(O2)COP(=O)(O)OP(=O)(O)OCC3C(C(C(O3)N4C=NC5=C(N=CN=C54)N)O)O)O)O'
 11: 'gene_id', # 'g5070'
 12: 'gene',    # 'NSDHL'
 13: 'ensembl_id',  # 'Not available' or 'ENSG00000147383'
 14: 'NCBI',    # 'Not available' or '50814'
 15: 'HGNC',    # 'Not available' or '13398'
 16: 'UniProt', # 'Not available' or 'Q15738'
 17: 'protein_size',    # '373'
 18: 'annonation',  # 'Not available' or 'The protein encoded by this gene is localized in the endoplasmic reticulum ...'
 19: 'score',   # 'Not available' or '0.95'
 20: 'alteration',  # 'Unknown' or ['elevated', 'reduced', 'target', 'Inhibitor', 'Activator']
 21: 'PMID',   # 'Not available' or '31142855'
 22: 'source'   # 'stitch, drugbank'
}
"""

load_dotenv()
CACHE_DIR = os.path.join(os.getcwd(), "cache")
os.makedirs(CACHE_DIR, exist_ok=True)


def save_pickle(obj, f_name):
    """
    :param obj: data to be saved as a pickle file
    :param f_name: files should only be existing in the cache directory
    :return:
    """
    with open(os.path.join(CACHE_DIR, f_name), "wb") as out_f:
        pickle.dump(obj, out_f)


def load_pickle(f_name):
    path = os.path.join(CACHE_DIR, f_name)
    return (
        pickle.load(open(path, "rb")) if os.path.exists(path) else print("The file does not exist.")
    )


def save_json(obj, f_name):
    with open(os.path.join(CACHE_DIR, f_name), "w") as out_f:
        json.dump(obj, out_f, indent=4)


def line_generator(in_file: str | os.PathLike, delimiter=",", skip_header=True) -> Iterator[list]:
    """generates lines from a CSV file, yielding each line as a list of strings
    This function opens the specified CSV file, skips the header row, and yields each subsequent line as a list of strings.

    :param skip_header:
    :param in_file: The path to the CSV file.
    :param delimiter: The character used to separate values in the CSV file (default is comma).
    :return: An iterator that yields each line of the CSV file as a list of strings.
    """
    with open(in_file, "r") as in_f:
        reader = csv.reader(in_f, delimiter=delimiter)
        if skip_header:
            next(reader)
        else:
            pass

        for line in reader:
            yield line


def get_gene_name(gene_ids: list) -> list:
    """
    Retrieves gene names for a given list of gene IDs using biothings_client
    The IDs are searched across multiple scopes: "ensembl.gene" and "uniprot".

    :param gene_ids: A list of gene IDs to be queried.
    :type gene_ids: List
    :return: A list of dictionaries containing the gene names and associated information.
    """
    gene_ids = set(gene_ids)
    t = bt.get_client("gene")
    gene_names = t.querymany(
        gene_ids, scopes=["uniprot", "ensembl.gene"], fields=["name"]
    )
    return gene_names


def get_bigg_metabolite_mapping(in_f):
    bigg_map = {line[2].lower(): line[1] for line in line_generator(in_f, delimiter="\t")}
    return bigg_map


def get_primary_chem_id(line, bigg_map):
    pubchem_id = line[3]
    drug_id = line[7]
    kegg_id = line[5]
    smiles = line[10]
    hmdb_id = line[6]
    bigg_id = bigg_map.get(line[2].lower())

    # (line, prefix) pairs for ID hierarchy
    id_hierarchy = [
        (pubchem_id, "PUBCHEM.COMPOUND"),
        (drug_id, "DRUGBANK"),
        (kegg_id, None),
        (smiles, ""),
        (hmdb_id, "HMDB"),
        (bigg_id, "BIGG.METABOLITE"),
    ]

    def classify_kegg(val):
        if val.startswith("C"):
            return "KEGG.COMPOUND"
        elif val.startswith("G"):
            return "KEGG.GLYCAN"
        elif val.startswith("D"):
            return "KEGG.DRUG"
        return "KEGG"

    xrefs = {}
    primary_id = None

    for val, prefix in id_hierarchy:
        if not val or val.strip().lower() == "not available":
            continue
        if prefix is None and val == kegg_id:
            prefix = classify_kegg(val)

        curie = f"{prefix}:{val}" if prefix else val
        if prefix == "PUBCHEM.COMPOUND":
            key = "pubchem_cid"
        elif prefix == "DRUGBANK":
            key = "drugbank"
        elif prefix == "HMDB":
            key = "hmdb"
        elif prefix is not None and prefix.startswith("KEGG"):
            key = "kegg"
        elif prefix == "":
            key = "smiles"
        elif prefix == "BIGG.METABOLITE":
            key = "bigg"
        else:
            continue

        if primary_id is None:
            primary_id = curie
        xrefs.setdefault(key, curie)

    return primary_id, xrefs


def get_primary_gene_id(line):
    ensembl_id = line[13]
    ncbi_id = line[14]
    hgnc_id = line[15]
    uniprotkb = line[16]

    # (line, prefix) pairs for ID hierarchy
    id_hierarchy = [
        (uniprotkb, "UniProtKB"),
        (ensembl_id, "ENSEMBL"),
        (ncbi_id, "NCBIGene"),
        (hgnc_id, "HGNC"),
    ]

    xrefs = {}
    primary_id = None

    for val, prefix in id_hierarchy:
        if not val or val.strip().lower() == "not available":
            continue

        curie = f"{prefix}:{val}" if prefix else val
        if prefix == "NCBIGene":
            key = "entrezgene"
        else:
            key = prefix.lower()

        if primary_id is None:
            primary_id = curie
        xrefs.setdefault(key, curie)

    return primary_id, xrefs


async def pug_query_pubchem_description(
    cid: int,
    session: aiohttp.ClientSession,
    sem: asyncio.Semaphore,
    max_retries: int = 3,
    delay: float = 1.0,
):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
    for attempt in range(max_retries):
        async with sem:
            try:
                async with session.get(url, timeout=30) as resp:
                    resp.raise_for_status()
                    data = await resp.json()
                break
            except aiohttp.ClientError as e:
                if attempt < max_retries - 1:
                    await asyncio.sleep(delay * (2**attempt))
                    continue
                else:
                    raise e

    description = None
    synonyms = []

    for sec in data.get("Record", {}).get("Section", []):
        heading = sec.get("TOCHeading", "")
        if heading == "Names and Identifiers" and description is None:
            for sub in sec.get("Section", []):
                if sub.get("TOCHeading") == "Record Description":
                    for info in sub.get("Information", []):
                        for mark in info.get("Value", {}).get("StringWithMarkup", []):
                            text = mark.get("String", "").strip()
                            if text:
                                description = text
                                break
                        if description:
                            break
                elif sub.get("TOCHeading") == "Synonyms":
                    for sec in sub.get("Section", []):
                        for info in sec.get("Information", []):
                            for mark in info.get("Value", {}).get("StringWithMarkup", []):
                                text = mark.get("String", "").strip()
                                if text:
                                    synonyms.append(text.lower().strip())
            break

    seen = set()
    synonyms = [s for s in synonyms if not (s in seen or seen.add(s))]

    return cid, {
        "id": f"PUBCHEM.COMPOUND:{cid}",
        "description": f"{description}[PUBCHEM]" if description else "",
        "synonyms": synonyms,
    }


async def get_batch_pubchem_descriptions_async(
    cids: List[int],
    workers: int = 5,
) -> Dict[int, Dict[str, str]]:
    cids = list(set(cids))

    workers = min(workers, 5)
    sem = asyncio.Semaphore(workers)
    connector = aiohttp.TCPConnector(limit_per_host=workers)
    cids = sorted(list(set(cids)))

    results = {}
    async with aiohttp.ClientSession(connector=connector) as session:
        tasks = []
        last_launch = 0.0

        for cid in cids:
            elapsed = time.perf_counter() - last_launch
            if elapsed < 1.0 / 5:
                await asyncio.sleep(1.0 / 5 - elapsed)
            last_launch = time.perf_counter()

            task = asyncio.create_task(pug_query_pubchem_description(cid, session, sem))
            tasks.append(task)

        for cid, payload in await tqdm_asyncio.gather(*tasks, total=len(tasks)):
            if payload:
                results[cid] = payload

    return results


def get_pubchem_descriptions(cids: List[int], workers: int = 5):
    return asyncio.run(get_batch_pubchem_descriptions_async(cids, workers))


def get_organism_type(node) -> str:
    """
    Inspect node['lineage'] for known taxids.
    Return the matching biolink CURIE, or Other if no match.
    Types include: 3 domains of life (Bacteria, Archaea, Eukaryota) and Virus.
    """
    taxon_map = {
        2: "biolink:Bacterium",
        2157: "Archaeon",
        2759: "Eukaryote",
        10239: "biolink:Virus",
    }

    lineage = node.get("lineage")
    if not isinstance(lineage, list):
        return "Other"
    for taxid, organism_type in taxon_map.items():
        if taxid in lineage:
            return organism_type
    return "Other"


def bt_get_mw_logp(cids: list):
    cids = sorted(set(cids))
    t = bt.get_client("chem")
    get_chem = t.querymany(
        cids,
        scopes="pubchem.cid",
        fields=["pubchem.molecular_weight", "pubchem.monoisotopic_weight", "pubchem.xlogp"],
    )

    q_out = {}
    for q in get_chem:
        if "notfound" in q:
            continue

        pub = q.get("pubchem")
        if isinstance(pub, list):
            info = pub[0]
        elif isinstance(pub, dict):
            info = pub
        else:
            continue

        mw = info.get("molecular_weight")
        mono = info.get("monoisotopic_weight")
        logp = info.get("xlogp")
        q_out[q["query"]] = {
            "molecular_weight": {
                "average_molecular_weight": mw,
                "monoisotopic_molecular_weight": mono,
            },
            "xlogp": logp,
        }

    return q_out


def cache_data(in_f):
    # cache pubchem descriptions
    pubchem_cids = [
        line[3] for line in line_generator(in_f) if line[3] and line[3] != "Not available"
    ]
    print(f"Total unique pubchem_cids: {len(set(pubchem_cids))}")
    pubchem_descr = get_pubchem_descriptions(pubchem_cids, workers=5)
    print(f"Total pubchem_cid with descriptions: {len(pubchem_descr)}")
    save_pickle(pubchem_descr, "gmmad2_meta_gene_description.pkl")

    # cache molecular weight and xlogp
    pubchem_mw_logp = bt_get_mw_logp(pubchem_cids)
    print(f"Total pubchem_cid with molecular weight and xlogp: {len(pubchem_mw_logp)}")
    save_pickle(pubchem_mw_logp, "gmmad2_meta_gene_pubchem_mw.pkl")


def get_node_info(file_path: str | os.PathLike) -> Iterator[dict]:
    """generates node dictionaries from meta_gene_net.csv file
    This function reads gene and metabolite data and processes it.

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

        # association node
        association_node = {"predicate": "biolink:associated_with"}

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

        # combine all the nodes
        output_dict = {
            "_id": None,
            "association": association_node,
            "object": object_node,
            "subject": subject_node,
        }

        if ":" in object_node["id"] and ":" in subject_node["id"]:
            output_dict[
                "_id"
            ] = f"{subject_node['id'].split(':')[1].strip()}_associated_with_{object_node['id'].split(':')[1].strip()}"
        else:
            output_dict[
                "_id"
            ] = f"{subject_node['id']}_associated_with_{object_node['id'].split(':')[1].strip()}"
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
