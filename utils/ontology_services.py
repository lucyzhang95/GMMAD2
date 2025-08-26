import asyncio
import os
import tarfile
import zipfile
from pathlib import Path
from typing import Dict, List

import aiohttp
import biothings_client as bt
import pandas as pd
import requests
from Bio import Entrez
from dotenv import load_dotenv
from ete3 import NCBITaxa
from tqdm.asyncio import tqdm


class NCBITaxonomyService:
    """NCBI Taxonomy lookups and ID remapping."""

    def __init__(self):
        """Initializes the service and the ETE3 NCBI Taxa object."""
        self.ncbi = NCBITaxa()
        print("NCBI Taxonomy Service initialized.")

    def _ensure_taxdump_exists(self, taxdump_path="taxdump.tar.gz") -> Path:
        """
        Ensures the taxdump file exists, downloading it if necessary.
        Returns the valid path to the file.
        """
        taxdump_path = Path(taxdump_path)
        if taxdump_path.exists():
            print(f"{taxdump_path} is found, using existing file.")
            return taxdump_path

        print(f"NCBI taxdump not found at {taxdump_path} – downloading via ETE3…")
        self.ncbi.update_taxonomy_database()
        return Path(taxdump_path)

    def get_taxid_mapping_from_ncbi_merged_dmp(
        self, tar_gz_path: str = "taxdump.tar.gz", f_name: str = "merged.dmp"
    ) -> dict:
        """
        Parses 'merged.dmp' from the provided NCBI taxdump tarball.
        Returns a dictionary mapping old taxids to new taxids.
        """
        taxid_mapping = {}
        with tarfile.open(tar_gz_path, "r:gz") as tar:
            member = tar.getmember(f_name)
            with tar.extractfile(member) as fp:
                for line in fp:
                    old_taxid, _, new_taxid, *_ = line.decode("utf-8").split("\t")
                    taxid_mapping[old_taxid] = new_taxid
        return taxid_mapping

    def get_current_taxid_mapping(self, old_taxids: list, merged_mapping: dict) -> dict[str, str]:
        """

        :param old_taxids: outdated taxids in a database
        :param merged_mapping: output of get_merged_taxid_mapping
        :return:
        """
        taxid_mapping = {}
        for old_taxid in old_taxids:
            taxid_mapping[old_taxid] = merged_mapping.get(old_taxid, None)
        return taxid_mapping

    def query_taxon_info_from_biothings(self, taxids: list) -> list:
        """retrieves taxonomic information for a given list of taxon IDs from disease_species.csv

        This function reads taxon IDs, removes duplicates, and queries taxonomic info from biothings_client
        to retrieve detailed taxonomic information including scientific name, parent taxid, lineage, and rank.

        :param taxids:
        :return: A list of dictionaries containing taxonomic information.
        """
        taxids = sorted(set(taxids))
        t = bt.get_client("taxon")
        taxon_info = t.gettaxa(
            taxids, fields=["scientific_name", "parent_taxid", "lineage", "rank"]
        )
        print(f"(BioThings API Call) Querying for {len(taxids)} taxids: {taxids[:5]}...")
        return taxon_info


class NCItTaxonomyService:
    """NCIT-based organism mappings and descriptions."""

    load_dotenv()

    def __init__(self):
        self.NCIT_API_KEY = os.getenv("NCIT_API_KEY")

    def download_ncit_source(self, f_name, out_path):
        """Downloads the NCI Thesaurus_25.06e.FLAT.zipfile from the NCI FTP server."""
        url = f"https://evs.nci.nih.gov/ftp1/NCI_Thesaurus/{f_name}"
        try:
            resp = requests.get(url, stream=True)
            resp.raise_for_status()
            with open(out_path, "wb") as f:
                for chunk in resp.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"{f_name} downloaded successfully.")

        except requests.exceptions.RequestException as e:
            print(f"An error occurred while downloading the file: {e}")

    def get_ncit_organism_mapping(self, f_name, path):
        if not Path(path).exists():
            self.download_ncit_source(f_name, path)

        with zipfile.ZipFile(path) as zf:
            with zf.open("Thesaurus.txt") as f:
                df = pd.read_csv(f, sep="\t", dtype=str).fillna("")

        df.columns = [
            "code",
            "concept IRI",
            "parents",
            "synonyms",
            "definition",
            "display name",
            "concept status",
            "semantic type",
            "concept in subset",
        ]

        relevant_types = {
            "Bacterium",
            "Virus",
            "Fungus",
            "Eukaryote",
            "Organism",
            "Animal",
            "Group|Organism",
        }

        df_organism = df[df["semantic type"].isin(relevant_types)].copy()

        mapping_result = {}
        for _, row in df_organism.iterrows():
            name = str(row["display name"]).strip().lower()
            description = str(row["definition"]).strip()
            ncit_code = str(row["code"]).strip()
            if not name or pd.isna(description):
                continue
            mapping_result[name] = {"description": description, "xrefs": {"ncit": ncit_code}}
        return mapping_result

    async def async_query_ncit_taxon_description(
        self, session: aiohttp.ClientSession, name: str, sem: asyncio.Semaphore
    ) -> tuple[str, dict] or None:
        SEARCH_URL = "https://data.bioontology.org/search"
        params = {
            "q": name,
            "ontologies": "NCIT",
            "apikey": self.NCIT_API_KEY,
        }
        async with sem:
            async with session.get(SEARCH_URL, params=params) as resp:
                resp.raise_for_status()
                data = await resp.json()
        for result in data.get("collection", []):
            if not result:
                continue
            pref_label = result.get("prefLabel", "").lower()
            if pref_label != name:
                continue

            definition = result.get("definition", [])

            uri = result.get("@id", "")
            ncit_id = uri.rsplit("/", 1)[-1] or uri.split("#")[-1]

            return name, {
                "description": f"{definition[0]} [NCIT]" if definition else "",
                "xrefs": {"ncit": ncit_id},
            }
        return None

    async def async_query_ncit_taxon_descriptions(
        self, taxon_names, max_concurrent=5
    ) -> dict[str, dict]:
        unique_names = sorted(
            list(set(str(n.lower().strip()) for n in taxon_names if n is not None))
        )
        print(
            f"(NCIt BioPortal Call) Querying for {len(unique_names)} taxon names: {unique_names[:5]}..."
        )
        sem = asyncio.Semaphore(max_concurrent)
        connector = aiohttp.TCPConnector(limit_per_host=max_concurrent)
        async with aiohttp.ClientSession(connector=connector) as session:
            tasks = [
                self.async_query_ncit_taxon_description(session, name, sem) for name in unique_names
            ]
            results = await tqdm.gather(*tasks, desc="Querying NCIT Descriptions")
        return {name: data for item in results if item for name, data in (item,)}

    def run_async_query_ncit_taxon_descriptions(self, taxon_names):
        """

        :param taxon_names:
        :return:
        {
           "enterobacter cloacae":{
              "query":"550",
              "_id":"550",
              "_version":1,
              "lineage":[
                 550,
                 354276,
                 547,
                 543,
                 91347,
                 1236,
                 1224,
                 3379134,
                 2,
                 131567,
                 1
              ],
              "parent_taxid":354276,
              "rank":"species",
              "scientific_name":"enterobacter cloacae",
              "description": "A species of facultatively anaerobic, Gram negative, rod shaped bacterium in the phylum Proteobacteria. This species is motile by peritrichous flagella, oxidase, urease and indole negative, catalase positive, reduces nitrate, does not degrade pectate and produces acid from sorbitol. E. cloacae is associated with hospital-acquired urinary and respiratory tract infections and is used in industry for explosives biodegradation. [NCIT]",
              "xrefs": {
                 "ncit":"C86360"
              }
           }
        }
        """
        return asyncio.run(self.async_query_ncit_taxon_descriptions(taxon_names))


class PubChemService:
    """PubChem PUG-REST description fetchers."""

    async def async_query_pug_pubchem_description(
        self,
        cid: int,
        session: aiohttp.ClientSession,
        sem: asyncio.Semaphore,
        max_retries: int = 3,
        delay: float = 1.0,
    ):
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
        data = {}
        for attempt in range(max_retries):
            async with sem:
                try:
                    async with session.get(url, timeout=30) as resp:
                        resp.raise_for_status()
                        data = await resp.json()
                    break
                except aiohttp.ClientError as e:
                    print(
                        f"Warning: Description query for pubchem_cid: {cid} failed on attempt {attempt + 1}. Error: {e}"
                    )
                    if attempt < max_retries - 1:
                        await asyncio.sleep(delay * (2**attempt))
                        continue
                    else:
                        print(f"Error: All retries failed for querying pubchem_cid: {cid}.")
                        return None
                finally:
                    await asyncio.sleep(2)

        if not data:
            return None

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

    async def async_query_pug_pubchem_descriptions(
        self,
        cids: List[int],
        workers: int = 5,
    ) -> Dict[int, Dict[str, str]]:
        """Queries PubChem for descriptions of a list of compound IDs (CIDs)."""

        unique_cids = sorted(list(set(cids)))
        sem = asyncio.Semaphore(workers)
        connector = aiohttp.TCPConnector(limit_per_host=workers)

        async with aiohttp.ClientSession(connector=connector) as session:
            tasks = [
                self.async_query_pug_pubchem_description(cid, session, sem) for cid in unique_cids
            ]

            print(f"Querying PubChem for {len(unique_cids)} CIDs...")
            results = await tqdm.gather(*tasks, desc="Querying PubChem Metabolite Descriptions")
        output = {cid: payload for item in results if item for cid, payload in (item,)}

        return output

    def run_async_query_pug_pubchem_descriptions(self, cids: List[int], workers: int = 5):
        return asyncio.run(self.async_query_pug_pubchem_descriptions(cids, workers))


class UniProtService:
    """UniProt REST API helpers."""

    async def async_query_uniprot_name_and_function(
        self,
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

    async def async_query_uniprot_names_and_functions(
        self,
        uniprot_ids: List[str],
        workers: int = 5,
    ) -> Dict[str, Dict[str, str]]:

        unique_ids = sorted(list(set(uniprot_ids)))

        sem = asyncio.Semaphore(workers)
        connector = aiohttp.TCPConnector(limit_per_host=workers)

        async with aiohttp.ClientSession(connector=connector) as session:
            tasks = [
                self.async_query_uniprot_name_and_function(uid, session, sem) for uid in unique_ids
            ]

            print(f"Querying UniProt for {len(unique_ids)} IDs...")
            results = await tqdm.gather(*tasks, desc="Querying UniProt")

        output = {uid: payload for item in results if item for uid, payload in (item,)}
        return output

    def run_async_query_uniprot_names_and_functions(
        self, uniprot_ids: List[str], workers: int = 5
    ) -> Dict[str, Dict[str, str]]:
        return asyncio.run(self.async_query_uniprot_names_and_functions(uniprot_ids, workers))


class GeneOntologyService:
    def query_gene_name_from_biothings(self, gene_ids: list) -> dict:
        """
        Retrieves gene names for a given list of gene IDs using biothings_client
        The IDs are searched across multiple scopes: "ensembl.gene".

        :param gene_ids: A list of gene IDs to be queried.
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


class ChemPropertyServices:
    def query_mw_and_xlogp_from_biothings(self, cids: list):
        cids = sorted(list(set(cids)))
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


class PubMedService:
    def __init__(self):
        self.EMAIL_ADDRESS = os.getenv("EMAIL_ADDRESS")

    def query_pubmed_metadata(self, pmids):
        """Get title, DOI, and abstract for a list of pmids using Entrez.

        :param pmids: a list of pmids obtained from core_table.txt
        :return: a dictionary with pmid as key and a dictionary with title, abstract, and doi as value.
        """
        Entrez.email = self.EMAIL_ADDRESS
        pmids = set(pmids)
        handle = Entrez.efetch(db="pubmed", id=",".join(map(str, pmids)), retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        result = {}
        for article in tqdm(records["PubmedArticle"], desc="Processing PubMed Articles"):
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
