import asyncio
import csv
import json
import os
import pickle
import time
import uuid
from io import BufferedWriter
from typing import Dict, Iterator, List

import aiohttp
import biothings_client as bt
from Bio import Entrez
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
 12: 'gene',    # 'NSDHL' -> as original_name: from GENECARD
 13: 'ensembl_id',  # 'Not available' or 'ENSG00000147383'
 14: 'NCBI',    # 'Not available' or '50814'
 15: 'HGNC',    # 'Not available' or '13398'
 16: 'UniProt', # 'Not available' or 'Q15738'
 17: 'protein_size',    # '373'
 18: 'annonation',  # 'Not available' or 'The protein encoded by this gene is localized in the endoplasmic reticulum ...'
 19: 'score',   # 'Not available' or '0.95'
 20: 'alteration',  # 'Unknown' or ['elevated', 'reduced', 'target', 'Inhibitor', 'Activator']
 21: 'PMID',   # 'Not available' or '31142855'
 22: 'source'   # ['stitch', 'gutMGene', 'stitch, gutMGene', 'stitch, drugbank', 'drugbank']
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
        pickle.load(open(path, "rb"))
        if os.path.exists(path)
        else print(f"The pickle file {path} does not exist.")
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


def get_gene_name(gene_ids: list) -> dict:
    """
    Retrieves gene names for a given list of gene IDs using biothings_client
    The IDs are searched across multiple scopes: "ensembl.gene".

    :param gene_ids: A list of gene IDs to be queried.
    :type gene_ids: List
    :return: A dictionary containing the gene names and associated information.
    """
    gene_ids = set(gene_ids)
    t = bt.get_client("gene")
    gene_q = t.querymany(gene_ids, scopes=["ensembl.gene"], fields=["name", "symbol"])
    gene_names = {
        d["query"]: {"name": d.get("symbol"), "full_name": d.get("name")}
        for d in gene_q
        if "notfound" not in d
    }
    return gene_names


async def uniprot_query_protein_info(
    uniprot_id: str,
    session: aiohttp.ClientSession,
    sem: asyncio.Semaphore,
    max_retries: int = 3,
    delay: float = 1.0,
    timeout: float = 30.0,
) -> (str, Dict[str, str]):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    for attempt in range(max_retries):
        async with sem:
            try:
                async with session.get(url, timeout=timeout) as resp:
                    resp.raise_for_status()
                    data = await resp.json()
                break
            except aiohttp.ClientError:
                if attempt < max_retries - 1:
                    await asyncio.sleep(delay * (2**attempt))
                    continue
                raise

    rec_name = data.get("proteinDescription", {}).get("recommendedName", {})
    raw_full = rec_name.get("fullName", {}).get("value")
    full_name = raw_full.lower() if isinstance(raw_full, str) else None

    raw_short = next(
        (sn.get("value") for sn in rec_name.get("shortNames", []) if sn.get("value")), None
    )
    if raw_short is None:
        alter_list = data.get("proteinDescription", {}).get("alternativeNames", [])
        raw_short = next(
            (
                sn.get("value")
                for alt in alter_list
                for sn in alt.get("shortNames", [])
                if sn.get("value")
            ),
            None,
        )
    name = raw_short if isinstance(raw_short, str) else None

    desc = None
    for comment in data.get("comments", []):
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts", [])
            if texts and texts[0].get("value"):
                desc = texts[0]["value"]
                break
    if desc is None:
        desc = "Function not found"

    return uniprot_id, {
        "name": name,
        "full_name": full_name,
        "description": desc,
    }


async def get_batch_protein_info_async(
    uniprot_ids: List[str],
    workers: int = 5,
    rate_limit: float = 5.0,
) -> Dict[str, Dict[str, str]]:
    ids = sorted(set(uniprot_ids))
    sem = asyncio.Semaphore(workers)
    connector = aiohttp.TCPConnector(limit_per_host=workers)
    results: Dict[str, Dict[str, str]] = {}

    async with aiohttp.ClientSession(connector=connector) as session:
        tasks = []
        last_launch = 0.0
        for uid in ids:
            elapsed = time.perf_counter() - last_launch
            min_interval = 1.0 / rate_limit
            if elapsed < min_interval:
                await asyncio.sleep(min_interval - elapsed)
            last_launch = time.perf_counter()

            task = asyncio.create_task(uniprot_query_protein_info(uid, session, sem))
            tasks.append(task)

        for uid, payload in await tqdm_asyncio.gather(*tasks, total=len(tasks)):
            results[uid] = payload

    return results


def get_protein_info(uniprot_ids: List[str], workers: int = 5) -> Dict[str, Dict[str, str]]:
    return asyncio.run(get_batch_protein_info_async(uniprot_ids, workers))


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


def get_pubmed_metadata(pmids):
    """Get title, DOI, and abstract for a list of pmids using Entrez.

    :param pmids: a list of pmids obtained from core_table.txt
    :return: a dictionary with pmid as key and a dictionary with title, abstract, and doi as value.
    """
    Entrez.email = os.getenv("EMAIL_ADDRESS")
    pmids = set(pmids)
    handle = Entrez.efetch(db="pubmed", id=",".join(map(str, pmids)), retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    result = {}
    for article in records["PubmedArticle"]:
        try:
            pmid = str(article["MedlineCitation"]["PMID"])
            article_data = article["MedlineCitation"]["Article"]

            title = article_data.get("ArticleTitle", "")
            abstract = ""
            if "Abstract" in article_data:
                abstract_parts = article_data["Abstract"].get("AbstractText", [])
                abstract = " ".join(str(part) for part in abstract_parts)

            doi = ""
            elist = article.get("PubmedData", {}).get("ArticleIdList", [])
            for el in elist:
                if el.attributes.get("IdType") == "doi":
                    doi = str(el)
                    break

            result[pmid] = {
                "name": title,
                "summary": f"{abstract} [abstract]",
                "doi": doi,
            }
        except Exception as e:
            print(f"Failed to parse article: {e}")
    return result


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

    # cache gene/protein information
    gene_ids = [
        line[16]
        if line[16] and line[16] != "Not available"
        else line[13]
        if line[13] and line[13] != "Not available"
        else None
        for line in line_generator(file_path)
    ]
    gene_ids = list(set(gene_ids))
    print(f"Total unique gene/protein ids: {len(gene_ids)}")
    uniprot_ids = [_id for _id in gene_ids if "ENSG" not in _id]
    print(f"Total unique uniprot_ids: {len(set(uniprot_ids))}")
    prot_info = get_protein_info(uniprot_ids)
    ensembl_gene_ids = [_id for _id in gene_ids if "ENSG" in _id]
    print(f"Total unique gene_ids: {len(set(ensembl_gene_ids))}")
    gene_info = get_gene_name(ensembl_gene_ids)
    full_gene_info = prot_info | gene_info
    print(f"Total unique gene/protein info: {len(full_gene_info)}")
    save_pickle(full_gene_info, "gmmad2_meta_gene_uniprot_ensemble_info.pkl")

    # cache pmid metadata
    pmids = [
        line[21] for line in line_generator(file_path) if line[21] and line[21] != "Not available"
    ]
    print(f"Total unique pmids: {len(set(pmids))}")
    pmid_metadata = get_pubmed_metadata(pmids)
    save_pickle(pmid_metadata, "gmmad2_meta_gene_pmid_metadata.pkl")


def remove_empty_none_values(obj):
    if isinstance(obj, dict):
        cleaned = {}
        for k, v in obj.items():
            v_clean = remove_empty_none_values(v)
            if v_clean not in (None, {}, []):
                cleaned[k] = v_clean
        return cleaned

    if isinstance(obj, list):
        cleaned_list = []
        for v in obj:
            v_clean = remove_empty_none_values(v)
            if v_clean not in (None, {}, []):
                cleaned_list.append(v_clean)
        return cleaned_list
    return obj


def get_suffix(identifier: str) -> str:
    return identifier.split(":", 1)[1].strip() if ":" in identifier else identifier.strip()


def get_node_info(file_path: str | os.PathLike) -> Iterator[dict]:
    """generates node dictionaries from meta_gene_net.csv file
    This function reads gene and metabolite data and processes it.

    :param file_path: path to meta_gene_net.csv file
    :return: An iterator of dictionaries containing node information.
    """
    if not os.path.exists(os.path.join(CACHE_DIR, "gmmad2_meta_gene_description.pkl")):
        cache_data(file_path)

    # load cached data
    pubchem_descr = load_pickle("gmmad2_meta_gene_description.pkl")
    pubchem_mw = load_pickle("gmmad2_meta_gene_pubchem_mw.pkl")
    gene_info = load_pickle("gmmad2_meta_gene_uniprot_ensemble_info.pkl")
    bigg_mapping = load_pickle("gmmad2_micro_meta_bigg_mapping.pkl")
    pmid_metadata = load_pickle("gmmad2_meta_gene_pmid_metadata.pkl")

    for line in line_generator(file_path):
        uniprot_id = line[16] if line[16] and line[16] != "Not available" else None
        ensemble = line[13] if line[13] and line[13] != "Not available" else None

        gene_record = (gene_info.get(uniprot_id) if uniprot_id else None) or (
            gene_info.get(ensemble) if ensemble else None
        )

        # object node (genes/proteins)
        object_node = {
            "id": f"UniProtKB:{uniprot_id}" if uniprot_id else f"ENSEMBL:{ensemble}",
            "name": gene_record.get("name") if gene_record else None,
            "full_name": gene_record.get("full_name") if gene_record else None,
            "original_name": line[12],
            "description": (
                line[18]
                if line[18] and line[18] != "Not available"
                else gene_record.get("description")
                if gene_record
                else None
            ),
            "residue_num": int(line[17]) if line[17] else None,
            "type": "biolink:Protein",
            "xrefs": {},
        }
        _, gene_xrefs = get_primary_gene_id(line)
        object_node["xrefs"] = gene_xrefs
        object_node = remove_empty_none_values(object_node)

        # subject node (metabolites)
        subject_node = {
            "id": None,
            "name": line[2].lower(),
            "synonym": pubchem_descr.get(line[3], {}).get("synonyms", []),
            "description": pubchem_descr.get(line[3], {}).get("description", ""),
            "chemical_formula": line[4] if line[4] and line[4] != "not available" else None,
            "molecular_weight": pubchem_mw.get(line[3], {}).get("molecular_weight", {}),
            "xlogp": pubchem_mw.get(line[3], {}).get("xlogp", None),
            "type": "biolink:SmallMolecule",
            "xrefs": {},
        }
        chem_id_xrefs = get_primary_chem_id(line, bigg_mapping)
        if chem_id_xrefs:
            chem_primary_id, chem_xrefs = chem_id_xrefs
        else:
            chem_primary_id, chem_xrefs = None, {}
        subject_node["id"] = chem_primary_id if chem_primary_id else str(uuid.uuid4())
        subject_node["xrefs"] = chem_xrefs
        subject_node = remove_empty_none_values(subject_node)

        # association node
        habitat = (
            [h.strip().lower() for h in line[9].split(";")]
            if line[9] and line[9] != "Unknown"
            else None
        )
        srcs = (
            [f"infores:{src.strip()}" for src in line[22].split(",")]
            if "," in line[22]
            else f"infores:{line[22].lower().strip()}"
        )
        evidence_map = {
            "infores:stitch": "ECO:0007669",  # computational evidence used in automatic assertion
            "infores:drugbank": "ECO:0000269",  # experimental evidence used in manual assertion
            "infores:gutMGene": "ECO:0000218",  # manual assertion
        }
        if isinstance(srcs, list):
            ecos = [evidence_map.get(s, "ECO:0000000") for s in srcs]
            evidence_type = ecos if len(ecos) > 1 else ecos[0]
        elif isinstance(srcs, str):
            evidence_type = evidence_map.get(srcs, "ECO:0000000")
        else:
            evidence_type = "ECO:0000000"

        pmid_key = line[21]
        metadata = pmid_metadata.get(pmid_key) if pmid_key and pmid_key != "Not available" else None

        association_node = {
            "id": "RO:0002434",
            "predicate": "biolink:ChemicalGeneInteractionAssociation",
            "type": "interacts_with",
            "association_habitat": habitat,
            "score": float(line[19]) if line[19] and line[19] != "Not available" else None,
            "qualifier": (
                "decrease"
                if line[20] and "reduced" in line[20].lower()
                else "increase"
                if line[20] and "elevated" in line[20].lower()
                else line[20].lower()
                if line[20]
                else None
            ),
            "primary_knowledge_source": srcs,
            "aggregator_knowledge_source": "infores:gmmad2",
            "evidence_type": evidence_type,
            "publications": {
                "pmid": int(line[21]) if line[21] and line[21] != "Not available" else None,
                "type": "biolink:Publication",
                "summary": metadata.get("summary") if metadata else None,
                "name": metadata.get("name") if metadata else None,
                "doi": metadata.get("doi") if metadata else None,
            },
        }
        association_node = remove_empty_none_values(association_node)

        output_dict = {
            "_id": None,
            "association": association_node,
            "object": object_node,
            "subject": subject_node,
        }

        subject_suffix = get_suffix(subject_node["id"])
        object_suffix = get_suffix(object_node["id"])
        output_dict["_id"] = f"{subject_suffix}_interacts_with_{object_suffix}"

        yield output_dict


def load_meta_gene_data(f_path) -> Iterator[dict]:
    """loads and yields unique meta gene data records from meta_gene_net.csv file.

    This function constructs the file path to meta_gene_net.csv file,
    retrieves node information using the `get_node_info` function,
    and yields unique records based on `_id`.

    :return: An iterator of unique dictionaries containing meta gene data.
    """
    assert os.path.exists(f_path), f"The file {f_path} does not exist."
    recs = get_node_info(f_path)
    for rec in recs:
        yield rec

    # dup_ids = set()
    # for rec in recs:
    #     if rec["_id"] not in dup_ids:
    #         dup_ids.add(rec["_id"])
    #         yield rec


if __name__ == "__main__":
    file_path = os.path.join("downloads", "meta_gene_net.csv")
    cache_data(file_path)

    meta_gene_data = [line for line in load_meta_gene_data(file_path)]
    _ids = []
    for obj in meta_gene_data:
        # print(obj)
        _ids.append(obj["_id"])
    print(f"total records: {len(_ids)}")
    print(f"total records without duplicates: {len(set(_ids))}")
