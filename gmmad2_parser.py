import asyncio
import csv
import itertools
import json
import logging
import os
import pickle
import tarfile
import uuid
import zipfile
from pathlib import Path
from typing import Dict, Iterator, List

import aiohttp
import biothings_client as bt
import pandas as pd
import requests
from Bio import Entrez
from dotenv import load_dotenv
from ete3 import NCBITaxa
from tqdm.asyncio import tqdm


class CSVParser:
    """Lightweight generator for large CSV/TXT tables."""

    def line_generator(
        self, in_file: str | os.PathLike, delimiter=",", skip_header=True
    ) -> Iterator[list]:
        """Generates lines from a CSV file, yielding each line as a list of strings
        This function opens the specified CSV file, skips the header row, and yields each following line as a list of strings.

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

    def line_generator_for_microbe_disease(
        self, in_file: str | os.PathLike, skip_header=True
    ) -> Iterator[list[str]]:
        """Yield each CSV line as a list of exactly 24 fields, specially handling disease_species.csv
        rejoining misaligned columns with commas for the 'disease' and 'disease_info' columns.
        """
        EXPECTED_COUNT = 24
        HEAD_COUNT = 2
        TAIL_COUNT = 5
        FIXED_COUNT = 15

        with open(in_file, "r", encoding="utf-8") as f:
            if skip_header:
                next(f)
            else:
                pass

            for line in f:
                parts = [part.strip() for part in line.split(",")]
                if len(parts) == EXPECTED_COUNT:
                    yield parts
                else:
                    alteration_idx = next(
                        i for i, p in enumerate(parts) if p in ("Increase", "Decrease")
                    )
                    head = parts[:HEAD_COUNT]
                    fixed_start = alteration_idx - (FIXED_COUNT - 1)
                    disease = ",".join(parts[2:fixed_start])
                    fixed = parts[fixed_start : alteration_idx + 1]
                    disease_info = ",".join(parts[alteration_idx + 1 : len(parts) - TAIL_COUNT])
                    tail = parts[-TAIL_COUNT:]
                    new_line = head + [disease] + fixed + [disease_info] + tail
                    assert (
                        len(new_line) == EXPECTED_COUNT
                    ), f"Expected {EXPECTED_COUNT} cols, got {len(new_line)}"
                    yield new_line


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

    def filter_taxon_info(self, taxon_info: list) -> dict and list:
        """
        Filters the taxon_info dictionary to remove not found entries.

        :param taxon_info: Output of query_taxon_info_from_biothings.
        :return:
        """
        taxon_info_filtered = {t["query"]: t for t in taxon_info if "notfound" not in t.keys()}
        taxid_notfound = set(t["query"] for t in taxon_info if "notfound" in t.keys())
        return taxon_info_filtered, list(taxid_notfound)

    def fetch_taxon_names_from_taxon_info(self, taxon_info: dict) -> list[str]:
        """Extracts biothings names from the taxon_info dictionary."""
        taxon_names = set()
        for _, taxon in taxon_info.items():
            if "scientific_name" in taxon:
                taxon_names.add(taxon["scientific_name"].lower())
        return list(taxon_names)


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
            ncit_id = result.get("@id", "").split("#")[-1]
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


class BiGGParser:
    """BiGG metabolite mapping helper."""

    def __init__(self):
        self.csv_parser = CSVParser()

    def get_bigg_metabolite_mapping(
        self, in_f: str | os.PathLike, delimiter: str = "\t", skip_header: bool = True
    ) -> Dict[str, str]:
        bigg_map: Dict[str, str] = {}
        for fields in self.csv_parser.line_generator(
            in_f, delimiter=delimiter, skip_header=skip_header
        ):
            bigg_id = fields[1].strip()
            metab_key = fields[2].strip().lower()
            bigg_map[metab_key] = bigg_id

        return bigg_map


class ChemPropertyUtils:
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


class CacheHelper:
    """Generic on-disk cache and serialization helpers."""

    DEFAULT_CACHE_DIR = os.path.join(os.getcwd(), "cache")

    def __init__(self, cache_dir=None):
        """Initializes the CacheManager and ensures the cache directory exists."""
        self.cache_dir = cache_dir or self.DEFAULT_CACHE_DIR
        os.makedirs(self.cache_dir, exist_ok=True)
        print(f"CacheManager initialized. Cache directory is {self.cache_dir}")

    def _get_path(self, f_name: str) -> str:
        """Helper method to construct the full path for a cache file."""
        return os.path.join(self.cache_dir, f_name)

    def save_pickle(self, obj, f_name: str):
        """Saves an object to a pickle file in the cache directory."""
        path = self._get_path(f_name)
        try:
            with open(path, "wb") as out_f:
                pickle.dump(obj, out_f)
            print(f"Saved pickle to: {path}")
        except (IOError, pickle.PicklingError) as e:
            print(f"Error saving pickle file {path}: {e}")

    def load_pickle(self, f_name: str):
        """Loads an object from a pickle file. Returns None if the file doesn't exist."""
        path = self._get_path(f_name)
        if not os.path.exists(path):
            return None
        try:
            with open(path, "rb") as in_f:
                return pickle.load(in_f)
        except (pickle.UnpicklingError, EOFError) as e:
            print(f"Error loading pickle file {path}: {e}")
            return None

    def save_json(self, obj, f_name: str, indent=4) -> str:
        """Saves an object to a JSON file in the cache directory."""
        path = self._get_path(f_name)
        try:
            with open(path, "w", encoding="utf-8") as out_f:
                json.dump(obj, out_f, indent=indent)
            print(f"Saved JSON to: {path}")
        except (IOError, TypeError) as e:
            print(f"Error saving JSON file {path}: {e}")

    def load_json(self, f_name: str):
        """Loads an object from a JSON file. Returns None if it fails."""
        path = self._get_path(f_name)
        if not os.path.exists(path):
            return None
        try:
            with open(path, "r", encoding="utf-8") as in_f:
                return json.load(in_f)
        except (IOError, json.JSONDecodeError) as e:
            print(f"Error loading JSON file {path}: {e}")
            return None


class CacheManager(CacheHelper):
    """CacheManager for managing on-disk caches and serialization."""

    def __init__(self, cache_dir=None):
        """Initializes the CacheManager and ensures the cache directory exists."""
        super().__init__(cache_dir)
        self.COMBINED_GENE_CACHE_F_NAME = "gmmad2_protein_gene_combined_info.pkl"

    def cache_entity(self, entity_type: str, **kwargs):
        """
        Build and cache data for a single entity type, updating the existing
        cache file instead of overwriting it.
        """
        entity_map = {
            "taxon_info": self._cache_taxon_info,
            "taxon_description": self.update_taxon_info_with_descriptions,
            "pubchem_description": self._cache_pubchem_description,
            "pubchem_mw": self._cache_pubchem_mw,
            "bigg_mapping": self._cache_bigg_mapping,
            "uniprot_info": self._cache_uniprot_info,
            "gene_info": self._cache_gene_info,
            "pubmed_metadata": self._cache_pubmed_metadata,
        }

        handler = entity_map.get(entity_type)
        if not handler:
            raise ValueError(f"Unknown entity type: {entity_type!r}")

        print(f"\n---Caching entity: '{entity_type}' ---")
        f_name = f"gmmad2_{entity_type}.pkl"
        existing_data = self.load_pickle(f_name) or {}
        if not isinstance(existing_data, dict):
            print(f"*Warning: Data in {f_name} is not a dictionary. Overwriting.")
            existing_data = {}

        kwargs["existing_data"] = existing_data
        new_data = handler(**kwargs)
        if new_data:
            print(f"-> Received {len(new_data)} new items to cache.")
            existing_data.update(new_data)
            self.save_pickle(existing_data, f_name)
        else:
            print(f"-> No new data to cache for '{entity_type}'.")

    def _cache_taxon_info(self, **kwargs):
        taxids = kwargs.get("taxids", [])
        if not taxids:
            return None

        f_name = "gmmad2_taxon_info.pkl"
        existing_data = self.load_pickle(f_name) or {}
        taxids_to_query = [tid for tid in taxids if tid not in existing_data]
        if not taxids_to_query:
            print("All requested taxids are already in the cache.")
            return None

        ncbi_service = NCBITaxonomyService()
        taxon_info = ncbi_service.query_taxon_info_from_biothings(taxids=taxids_to_query)
        filtered_taxon_info, taxid_notfound = ncbi_service.filter_taxon_info(taxon_info)

        if taxid_notfound:
            print(f"Found {len(taxid_notfound)} taxids that were not found: {taxid_notfound[:5]}")
            notfound_f_name = "gmmad2_taxon_info_notfound.pkl"
            existing_notfound = self.load_pickle(notfound_f_name) or []
            updated_notfound = sorted(list(set(existing_notfound + taxid_notfound)))
            self.save_pickle(updated_notfound, notfound_f_name)

        return filtered_taxon_info

    def update_taxon_info_with_current_taxids(self):
        """
        Loads not-found taxids, finds their current IDs, queries for them,
        and updates the main taxon_info cache.
        """
        print("\n---Updating cache with re-mapped taxids---")

        notfound_f_name = "gmmad2_taxon_info_notfound.pkl"
        main_cache_f_name = "gmmad2_taxon_info.pkl"
        taxid_notfound = self.load_pickle(notfound_f_name) or []
        main_taxon_info_cache = self.load_pickle(main_cache_f_name) or {}

        old_taxids_to_process = [tid for tid in taxid_notfound if tid not in main_taxon_info_cache]
        if not old_taxids_to_process:
            print("All 'not found' taxids are already in the cache.")
            return

        ncbi_service = NCBITaxonomyService()
        taxid_map = ncbi_service.get_taxid_mapping_from_ncbi_merged_dmp()
        new_taxid_map = ncbi_service.get_current_taxid_mapping(
            old_taxids=old_taxids_to_process, merged_mapping=taxid_map
        )

        # remap old taxids to existing taxon info with the new taxids
        # "93930" is in taxid_notfound, "93930" -> "93929" is already in main_taxon_info_cache
        # reverse map "93930" to the taxon info of "93929"
        info_to_add_from_remap = {}
        new_taxids_to_query = []
        for old_id, new_id in new_taxid_map.items():
            if not new_id:
                continue
            if new_id in main_taxon_info_cache:
                info_to_add_from_remap[old_id] = main_taxon_info_cache[new_id]
            else:
                new_taxids_to_query.append(new_id)

        new_taxon_info_to_add = {}  # new taxids after remapping
        if new_taxids_to_query:
            new_taxon_info = ncbi_service.query_taxon_info_from_biothings(
                taxids=new_taxids_to_query
            )
            valid_new_taxon_info, _ = ncbi_service.filter_taxon_info(new_taxon_info)

            reverse_map = {v: k for k, v in new_taxid_map.items()}
            for new_id, new_info in valid_new_taxon_info.items():
                old_id = reverse_map.get(new_id)
                if old_id:
                    new_taxon_info_to_add[old_id] = new_info
                new_taxon_info_to_add[new_id] = new_info

        if not info_to_add_from_remap and not new_taxon_info_to_add:
            print("No new information to add to the cache.")
            return

        print(
            f"Adding {len(info_to_add_from_remap) + len(new_taxon_info_to_add)} new/remapped entries to the main cache."
        )
        main_taxon_info_cache.update(info_to_add_from_remap)
        main_taxon_info_cache.update(new_taxon_info_to_add)
        self.save_pickle(main_taxon_info_cache, main_cache_f_name)

    def update_taxon_info_with_descriptions(self):
        """
        Loads the main taxon info cache, finds entries missing a description,
        queries for them, and saves the updated data back to the same file.
        """
        print("\n---Updating taxon info with NCIT descriptions---")

        main_cache_f_name = "gmmad2_taxon_info.pkl"
        desc_cache_f_name = "gmmad2_ncit_descriptions.pkl"

        taxon_info_cache = self.load_pickle(main_cache_f_name)
        if not taxon_info_cache:
            print(f"Main taxon info cache not found at '{main_cache_f_name}'. Cannot update.")
            return

        print(f"Loaded '{main_cache_f_name}' with {len(taxon_info_cache)} keys.")

        desc_cache = self.load_pickle(desc_cache_f_name) or {}
        names_need_desc = set()
        for info in taxon_info_cache.values():
            if "description" not in info or not info.get("description"):
                name = info["scientific_name"]
                names_need_desc.add(name)

        if not names_need_desc:
            print("All taxon entries already have descriptions.")
            return

        names_to_query = sorted([name for name in names_need_desc if name not in desc_cache])

        if names_to_query:
            ncit_service = NCItTaxonomyService()
            new_descriptions = ncit_service.run_async_query_ncit_taxon_descriptions(
                taxon_names=names_to_query
            )
            if new_descriptions:
                print(f"Found {len(new_descriptions)} new descriptions from API.")
                desc_cache.update(new_descriptions)
                self.save_pickle(desc_cache, desc_cache_f_name)
        else:
            print("All needed descriptions were already in the local description cache.")

        update_count = 0
        for taxid, info in taxon_info_cache.items():
            if ("description" not in info or not info.get("description")) and info.get(
                "scientific_name"
            ):
                name = info["scientific_name"]
                if name in desc_cache:
                    taxon_info_cache[taxid]["description"] = desc_cache[name]
                    update_count += 1
        if update_count > 0:
            print(f"Updated {update_count} entries in the main taxon info cache with descriptions.")
            print(f"Saving '{main_cache_f_name}' with {len(taxon_info_cache)} keys.")
            self.save_pickle(taxon_info_cache, main_cache_f_name)
        else:
            print("-> No entries in the main taxon info cache needed an update.")

    def _cache_pubchem_description(self, **kwargs):
        cids = kwargs.get("cids", [])
        if not cids:
            return
        f_name = "gmmad2_pubchem_descriptions.pkl"
        existing_data = self.load_pickle(f_name) or {}
        cids_to_query = [cid for cid in cids if cid not in existing_data]
        if not cids_to_query:
            print("All requested pubchem_cids are already in the cache.")
            return

        pubchem_service = PubChemService()
        pubchem_desc = pubchem_service.run_async_query_pug_pubchem_descriptions(cids=cids_to_query)

        if not pubchem_desc:
            print("-> API query did not return any new descriptions.")
            return
        print(f"Received {len(pubchem_desc)} new PubChem descriptions to cache.")
        existing_data.update(pubchem_desc)

        self.save_pickle(existing_data, f_name)

    def _cache_pubchem_mw(self, **kwargs):
        print("\n---Caching PubChem Metabolite Physical Properties---")
        cids = kwargs.get("cids", [])
        if not cids:
            return

        existing_data = self.load_pickle("gmmad2_pubchem_mw.pkl") or {}
        cids_to_query = [cid for cid in cids if cid not in existing_data]
        if not cids_to_query:
            print("All requested pubchem_cids are already in the cache.")
            return

        chem_property_service = ChemPropertyUtils()
        new_chem_properties = chem_property_service.query_mw_and_xlogp_from_biothings(
            cids=cids_to_query
        )
        if not new_chem_properties:
            print("-> API query did not return any new properties.")
            return
        print(f"Received {len(new_chem_properties)} new PubChem properties to cache.")
        existing_data.update(new_chem_properties)
        self.save_pickle(existing_data, "gmmad2_pubchem_mw.pkl")

    def _cache_bigg_mapping(self, **kwargs):
        print("\n---Caching BiGG Metabolite Mapping---")
        in_f = kwargs.get("downloads", "bigg_models_metabolites.txt")
        if not os.path.exists(in_f):
            print(f"BiGG mapping file '{in_f}' does not exist. Skipping caching.")
            return

        bigg_parser = BiGGParser()
        bigg_map = bigg_parser.get_bigg_metabolite_mapping(in_f=in_f)

        if not bigg_map:
            print("No BiGG mapping data found.")
            return

        self.save_pickle(bigg_map, "gmmad2_bigg_metabolite_mapping.pkl")
        print(f"BiGG mapping cached with {len(bigg_map)} entries.")

    def _cache_uniprot_info(self, **kwargs):
        print("\n---Caching UniProt Information---")
        uniprot_ids = kwargs.get("uniprots", [])
        if not uniprot_ids:
            return
        existing_data = self.load_pickle(self.COMBINED_GENE_CACHE_F_NAME) or {}
        uids_to_query = [uid for uid in uniprot_ids if uid not in existing_data]
        if not uids_to_query:
            print("All requested uniprot ids are already in the cache.")
            return

        uniprot_service = UniProtService()
        new_protein_info = uniprot_service.run_async_query_uniprot_names_and_functions(
            uniprot_ids=uids_to_query
        )

        if not new_protein_info:
            print("-> API query did not return any new descriptions.")
            return
        print(f"Received {len(new_protein_info)} new Uniprot entires to cache.")
        existing_data.update(new_protein_info)
        self.save_pickle(existing_data, self.COMBINED_GENE_CACHE_F_NAME)

    def _cache_gene_info(self, **kwargs):
        print("\n---Caching Gene Information---")
        gene_ids = kwargs.get("gene_ids", [])
        if not gene_ids:
            return

        existing_data = self.load_pickle(self.COMBINED_GENE_CACHE_F_NAME) or {}
        gids_to_query = [gid for gid in gene_ids if gid not in existing_data]
        if not gids_to_query:
            print("All requested gene ids are already in the cache.")
            return

        go_service = GeneOntologyService()
        new_gene_info = go_service.query_gene_name_from_biothings(gene_ids=gids_to_query)

        if not new_gene_info:
            print("-> API query did not return any new descriptions.")
            return
        print(f"Received {len(new_gene_info)} new Gene descriptions to cache.")
        existing_data.update(new_gene_info)
        self.save_pickle(existing_data, self.COMBINED_GENE_CACHE_F_NAME)

    def _cache_pubmed_metadata(self, **kwargs):
        print("\n---Caching PubMed Metadata---")
        pmids = kwargs.get("pmids", [])
        if not pmids:
            return

        existing_data = self.load_pickle("gmmad2_pubmed_metadata.pkl") or {}
        pmids_to_query = [pmid for pmid in pmids if pmid not in existing_data]
        if not pmids_to_query:
            print("All requested PMIDs are already in the cache.")
            return

        pubmed_service = PubMedService()
        new_pubmed_metadata = pubmed_service.query_pubmed_metadata(pmids=pmids_to_query)

        if not new_pubmed_metadata:
            print("-> API query did not return any new metadata.")
            return
        print(f"Received {len(new_pubmed_metadata)} new PubMed entries to cache.")
        existing_data.update(new_pubmed_metadata)
        self.save_pickle(existing_data, "gmmad2_pubmed_metadata.pkl")


class DataCachePipeline:
    def __init__(self, cache_dir="cache", downloads_dir="downloads"):
        self.downloads_dir = downloads_dir
        self.cache_dir = cache_dir
        self.midi_path = os.path.join(self.downloads_dir, "disease_species.csv")
        self.mime_path = os.path.join(self.downloads_dir, "micro_metabolic.csv")
        self.mege_path = os.path.join(self.downloads_dir, "meta_gene_net.csv")
        self.csv_parser = CSVParser()
        self.cache_manager = CacheManager(cache_dir=self.cache_dir)
        print("Data Cache Pipeline initialized...")

    def _cache_midi_taxon_info(self):
        taxids = [
            line[5] for line in self.csv_parser.line_generator_for_microbe_disease(self.midi_path)
        ]
        self.cache_manager.cache_entity("taxon_info", taxids=taxids)

    def _cache_mime_taxon_info(self):
        taxids = [
            line[9]
            if (line[9] and line[9] != "not available")
            else line[16]
            if (line[16] and line[16] != "not available")
            else None
            for line in self.csv_parser.line_generator(self.mime_path)
        ]
        taxids = [t for t in taxids if t]
        self.cache_manager.cache_entity("taxon_info", taxids=taxids)

    def _update_taxon_info(self):
        self.cache_manager.update_taxon_info_with_current_taxids()

    def _update_taxon_info_with_ncit_descriptions(self):
        self.cache_manager.update_taxon_info_with_descriptions()

    def _verify_taxon_info_cache(self):
        taxon_info_cache = self.cache_manager.load_pickle("gmmad2_taxon_info.pkl")
        if taxon_info_cache:
            print(f"Final cache contains {len(taxon_info_cache.keys())} unique taxids.")
        else:
            print("Final cache is empty or could not be loaded.")

    def _cache_mime_pubchem_descriptions(self):
        pubchem_cids = [
            line[6]
            for line in self.csv_parser.line_generator(self.mime_path)
            if line[6] and line[6] != "not available"
        ]
        self.cache_manager.cache_entity("pubchem_description", cids=pubchem_cids)

    def _cache_mege_pubchem_descriptions(self):
        pubchem_cids = [
            line[3]
            for line in self.csv_parser.line_generator(self.mege_path)
            if line[3] and line[3] != "Not available"
        ]
        self.cache_manager.cache_entity("pubchem_description", cids=pubchem_cids)

    def _verify_pubchem_cache(self):
        pubchem_desc_cache = self.cache_manager.load_pickle("gmmad2_pubchem_descriptions.pkl")
        if pubchem_desc_cache:
            print(f"PubChem cache contains {len(pubchem_desc_cache.keys())} unique CIDs.")
        else:
            print("PubChem cache is empty or could not be loaded.")

    def _cache_mege_gene_and_protein_info(self):
        gene_ids = [
            line[16]
            if line[16] and line[16] != "Not available"
            else line[13]
            if line[13] and line[13] != "Not available"
            else None
            for line in self.csv_parser.line_generator(self.mege_path)
        ]
        gene_ids = list(set(gene_ids))
        uniprot_ids = [_id for _id in gene_ids if "ENSG" not in _id]
        ensembl_gene_ids = [_id for _id in gene_ids if "ENSG" in _id]
        self.cache_manager.cache_entity("uniprot_info", uniprots=uniprot_ids)
        self.cache_manager.cache_entity("gene_info", gene_ids=ensembl_gene_ids)

    def _verify_gene_and_protein_info_cache(self):
        gene_and_protein_info = self.cache_manager.load_pickle(
            self.cache_manager.COMBINED_GENE_CACHE_F_NAME
        )
        if gene_and_protein_info:
            print(
                f"Gene and protein info cache contains {len(gene_and_protein_info.keys())} unique entries."
            )
        else:
            print("Gene and protein info cache is empty or could not be loaded.")

    def _cache_mege_pubmed_metadata(self):
        pmids = [
            line[21]
            for line in self.csv_parser.line_generator(self.mege_path)
            if line[21] and line[21] != "Not available"
        ]
        self.cache_manager.cache_entity("pubmed_metadata", pmids=pmids)

    def _verify_pubmed_metadata_cache(self):
        pubmed_metadata = self.cache_manager.load_pickle("gmmad2_pubmed_metadata.pkl")
        if pubmed_metadata:
            print(f"PubMed metadata cache contains {len(pubmed_metadata.keys())} unique PMIDs.")
        else:
            print("PubMed metadata cache is empty or could not be loaded.")

    def _cache_mime_chem_properties(self):
        pubchem_cids = [
            line[6]
            for line in self.csv_parser.line_generator(self.mime_path)
            if line[6] and line[6] != "not available"
        ]
        self.cache_manager.cache_entity("pubchem_mw", cids=pubchem_cids)

    def _cache_mege_chem_properties(self):
        pubchem_cids = [
            line[3]
            for line in self.csv_parser.line_generator(self.mege_path)
            if line[3] and line[3] != "Not available"
        ]
        self.cache_manager.cache_entity("pubchem_mw", cids=pubchem_cids)

    def _verify_pubchem_mw_cache(self):
        pubchem_mw_cache = self.cache_manager.load_pickle("gmmad2_pubchem_mw.pkl")
        if pubchem_mw_cache:
            print(f"PubChem MW cache contains {len(pubchem_mw_cache.keys())} unique CIDs.")
        else:
            print("PubChem MW cache is empty or could not be loaded.")

    def _cache_bigg_mapping(self):
        bigg_parser = BiGGParser()
        bigg_map = bigg_parser.get_bigg_metabolite_mapping(
            in_f=os.path.join(self.downloads_dir, "bigg_models_metabolites.txt")
        )
        if not bigg_map:
            print("No BiGG mapping data found.")
            return
        self.cache_manager.save_pickle(bigg_map, "gmmad2_bigg_metabolite_mapping.pkl")
        print(f"BiGG mapping cached with {len(bigg_map)} entries.")

    def run_cache_pipeline(self):
        print("Running data cache pipeline...")
        self._cache_midi_taxon_info()
        self._cache_mime_taxon_info()
        self._cache_mime_pubchem_descriptions()
        self._cache_mege_pubchem_descriptions()
        self._cache_mege_gene_and_protein_info()
        self._cache_mege_pubmed_metadata()
        self._cache_mime_chem_properties()
        self._cache_mege_chem_properties()
        self._cache_bigg_mapping()

        self._update_taxon_info()
        self._update_taxon_info_with_ncit_descriptions()

        self._verify_taxon_info_cache()
        self._verify_pubchem_cache()
        self._verify_gene_and_protein_info_cache()
        self._verify_pubmed_metadata_cache()
        self._verify_pubchem_mw_cache()


class ParserHelper:
    def get_organism_type(self, node) -> str:
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

        for taxid, organism_type in taxon_map.items():
            if taxid in node.get("lineage", []):
                return organism_type

        return "Other"

    def get_suffix(self, identifier: str) -> str:
        return identifier.split(":", 1)[1].strip() if ":" in identifier else identifier.strip()

    def assign_primary_metabolite_id(self, id_hierarchy: List[tuple]) -> tuple[str, dict] or None:
        def classify_kegg(val):
            if val.startswith("C"):
                return "KEGG.COMPOUND"
            if val.startswith("G"):
                return "KEGG.GLYCAN"
            if val.startswith("D"):
                return "KEGG.DRUG"
            return "KEGG"

        xrefs = {}
        primary_id = None
        for val, prefix in id_hierarchy:
            if not val or str(val).strip().lower() == "not available":
                continue
            if prefix == "KEGG" and isinstance(val, str):
                prefix = classify_kegg(val)

            curie = f"{prefix}:{val}" if prefix else str(val)
            key_map = {
                "PUBCHEM.COMPOUND": "pubchem_cid",
                "DRUGBANK": "drugbank",
                "HMDB": "hmdb",
                "SMILES": "smiles",
                "BIGG.METABOLITE": "bigg",
                "KEGG.COMPOUND": "kegg",
                "KEGG.GLYCAN": "kegg",
                "KEGG.DRUG": "kegg",
            }
            key = key_map.get(prefix)
            if not key:
                continue

            if primary_id is None:
                primary_id = curie
            xrefs.setdefault(key, curie)

        return primary_id, xrefs

    def assign_primary_gene_id(self, line):
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

    def remove_empty_none_values(self, obj):
        if isinstance(obj, dict):
            cleaned = {}
            for k, v in obj.items():
                v_clean = self.remove_empty_none_values(v)
                if v_clean not in (None, {}, []):
                    cleaned[k] = v_clean
            return cleaned

        if isinstance(obj, list):
            cleaned_list = []
            for v in obj:
                v_clean = self.remove_empty_none_values(v)
                if v_clean not in (None, {}, []):
                    cleaned_list.append(v_clean)
            return cleaned_list
        return obj


class GMMAD2Parser(CacheHelper):
    """
    Parses and merges microbe-disease, microbe-metabolite,
    metabolite-gene associations.
    """

    def __init__(
        self,
        csv_parser: CSVParser,
        cache_pipeline: DataCachePipeline,  # require CacheManager and CacheHelper
        parser_helpers: ParserHelper,
        cache_dir="cache",
        downloads_dir="downloads",
    ):
        super().__init__(cache_dir)  # inherit from CacheHelper
        self.csv_parser = csv_parser
        self.cache_pipeline = cache_pipeline
        self.parser_helpers = parser_helpers

        self.midi_path = os.path.join(downloads_dir, "disease_species.csv")
        self.mime_path = os.path.join(downloads_dir, "micro_metabolic.csv")
        self.mege_path = os.path.join(downloads_dir, "meta_gene_net.csv")

        self.logger = logging.getLogger(__name__)
        logging.basicConfig(level=logging.INFO)

    def _get_disease_node(self, line: list) -> dict:
        node = {
            "id": f"MESH:{line[1]}",
            "name": line[2].lower(),
            "type": "biolink:Disease",
            "description": line[18],
            "xrefs": {"mesh": f"MESH:{line[1]}"},
        }
        return self.parser_helpers.remove_empty_none_values(node)

    def _get_microbe_node(self, taxid: str, original_name: str, taxon_info: dict) -> dict:
        if taxid in taxon_info:
            info = taxon_info[taxid]
            node = {
                "id": f"NCBITaxon:{info.get('_id')}",
                "taxid": int(info.get("_id")) if info.get("_id") else None,
                "name": info.get("scientific_name", "").lower(),
                "original_name": original_name.lower(),
                "description": info.get("description"),
                "parent_taxid": info.get("parent_taxid"),
                "lineage": info.get("lineage", []),
                "rank": info.get("rank"),
                "type": "biolink:OrganismTaxon",
                "xrefs": info.get("xrefs", {}),
            }
            node["organism_type"] = self.parser_helpers.get_organism_type(node)
        else:
            if taxid:
                self.logger.warning(
                    f"NCBI taxid '{taxid}' not found in cache. Creating a uuid for microbe node."
                )
            else:
                self.logger.warning(
                    f"NCBI taxid missing for '{original_name.lower()}'. Using uuid for node id."
                )

            node = {
                "id": f"uuid:{str(uuid.uuid4())}",
                "original_name": original_name.lower(),
                "type": "biolink:OrganismTaxon",
                "organism_type": "Other",
            }
        return self.parser_helpers.remove_empty_none_values(node)

    def _get_metabolite_node(
        self,
        cid: str,
        name: str,
        formula: str,
        primary_id: str,
        xrefs: dict,
        pubchem_desc: dict,
        pubchem_mw: dict,
    ) -> dict:

        node = {
            "id": primary_id
            if primary_id
            else f"uuid:{str(uuid.uuid4())}",  # use uuid if no primary id but has compound name
            "name": name.lower(),
            "synonym": pubchem_desc.get(cid, {}).get("synonyms", []),
            "description": pubchem_desc.get(cid, {}).get("description", ""),
            "chemical_formula": formula if formula and formula != "not available" else None,
            "molecular_weight": pubchem_mw.get(cid, {}).get("molecular_weight", {}),
            "xlogp": pubchem_mw.get(cid, {}).get("xlogp", None),
            "type": "biolink:SmallMolecule",
            "xrefs": xrefs,
        }
        return self.parser_helpers.remove_empty_none_values(node)

    def _get_gene_node(self, line: list, protein_and_gene_info: dict) -> dict:
        uniprot_id = line[16] if line[16] and line[16] != "Not available" else None
        ensemble = line[13] if line[13] and line[13] != "Not available" else None
        gene_record = (protein_and_gene_info.get(uniprot_id) if uniprot_id else None) or (
            protein_and_gene_info.get(ensemble) if ensemble else None
        )

        node = {
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
        }
        _, gene_xrefs = self.parser_helpers.assign_primary_gene_id(line)
        node["xrefs"] = gene_xrefs
        return self.parser_helpers.remove_empty_none_values(node)

    def _get_midi_association_node(self, line: list) -> dict:
        node = {
            "predicate": "OrganismalEntityAsAModelOfDiseaseAssociation",
            "type": "biolink:associated_with",
            "qualifier": line[17].lower(),
            "qualifier_ratio": float(line[16]) if line[16] else None,
            "disease_sample_size": int(line[6]) if line[6] else None,
            "disease_abundance_mean": float(line[7]) if line[7] else None,
            "disease_abundance_median": float(line[8]) if line[8] else None,
            "disease_abundance_sd": float(line[9]) if line[9] else None,
            "control_id": line[10],
            "control_name": line[11].lower(),
            "healthy_sample_size": int(line[12]) if line[12] else None,
            "healthy_abundance_mean": float(line[13]) if line[13] else None,
            "healthy_abundance_median": float(line[14]) if line[14] else None,
            "healthy_abundance_sd": float(line[15]) if line[15] else None,
            "primary_knowledge_source": "infores:GMrepo",
            "aggregator_knowledge_source": "infores:GMMAD2",
            "evidence_type": "ECO:0000221",  # high throughput nucleotide sequencing assay evidence
        }
        return self.parser_helpers.remove_empty_none_values(node)

    def _get_mime_association_node(self, line: list) -> dict:
        evidence_map = {
            "infores:wom": "ECO:0001230",  # mass spec + manual
            "infores:vmh": "ECO:0000218",  # manual assertion
            "infores:gutMGene": "ECO:0000218",  # manual assertion
            "infores:Metabolomics": "ECO:0001230",  # mass spec + manual
        }
        src = f"infores:{line[17].strip()}"
        habitat = (
            [h.strip().lower() for h in line[20].split(";")]
            if line[20] and line[20] != "Unknown"
            else None
        )
        node = {
            "predicate": "biolink:OrganismTaxonToChemicalEntityAssociation",
            "type": "has_metabolic_interaction_with",
            "association_habitat": habitat,
            "primary_knowledge_source": src,
            "aggregator_knowledge_source": "infores:GMMAD2",
            "evidence_type": evidence_map.get(src, "ECO:0000000"),
        }
        return self.parser_helpers.remove_empty_none_values(node)

    def _get_mege_association_node(self, line: list, pmid_metadata: dict) -> dict:
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
        else:
            evidence_type = evidence_map.get(srcs, "ECO:0000000")

        pmid_key = line[21]
        metadata = pmid_metadata.get(pmid_key) if pmid_key and pmid_key != "Not available" else None

        node = {
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
            "aggregator_knowledge_source": "infores:GMMAD2",
            "evidence_type": evidence_type,
            "publications": {
                "pmid": int(line[21]) if line[21] and line[21] != "Not available" else None,
                "type": "biolink:Publication",
                "summary": metadata.get("summary") if metadata else None,
                "name": metadata.get("name") if metadata else None,
                "doi": metadata.get("doi") if metadata else None,
            },
        }
        return self.parser_helpers.remove_empty_none_values(node)

    def parse_microbe_disease(self) -> Iterator[dict]:
        print("\n--- Parsing Microbe-Disease Data ---")
        required_cache = "gmmad2_taxon_info.pkl"
        if not os.path.exists(self._get_path(required_cache)):
            self.cache_pipeline.run_cache_pipeline()
        taxon_info = self.load_pickle(required_cache)

        rec_ct = 0
        for line in self.csv_parser.line_generator_for_microbe_disease(self.midi_path):
            subject_node = self._get_microbe_node(line[5], line[3], taxon_info)
            object_node = self._get_disease_node(line)
            association_node = self._get_midi_association_node(line)
            yield {
                "_id": f"{self.parser_helpers.get_suffix(subject_node['id'])}_associated_with_{self.parser_helpers.get_suffix(object_node['id'])}",
                "association": association_node,
                "object": object_node,
                "subject": subject_node,
            }
            rec_ct += 1
        print(f"--- Finished Parsing Microbe-Disease. Generated {rec_ct} records. ---\n")

    def parse_microbe_metabolite(self) -> Iterator[dict]:
        print("\n--- Parsing Microbe-Metabolite Data ---")
        taxon_cache = self.load_pickle("gmmad2_taxon_info.pkl")
        pubchem_desc = self.load_pickle("gmmad2_pubchem_descriptions.pkl")
        pubchem_mw = self.load_pickle("gmmad2_pubchem_mw.pkl")
        bigg_mapping = self.load_pickle("gmmad2_bigg_metabolite_mapping.pkl")

        rec_ct = 0
        for line in self.csv_parser.line_generator(self.mime_path):
            taxid = (
                line[9]
                if line[9] and line[9] != "not available"
                else line[16]
                if line[16] and line[16] != "not available"
                else None
            )
            subject_node = self._get_microbe_node(taxid, line[2], taxon_cache)

            mime_id_hierarchy = [
                (line[6], "PUBCHEM.COMPOUND"),
                (line[8], "KEGG"),
                (line[18], "SMILES"),
                (line[19], "HMDB"),
                (bigg_mapping.get(line[5].lower()), "BIGG.METABOLITE"),
            ]
            primary_id, xrefs = self.parser_helpers.assign_primary_metabolite_id(mime_id_hierarchy)
            object_node = self._get_metabolite_node(
                cid=line[6],
                name=(line[5] if line[5] else line[4]),
                formula=line[7],
                primary_id=primary_id,
                xrefs=xrefs,
                pubchem_desc=pubchem_desc,
                pubchem_mw=pubchem_mw,
            )

            association_node = self._get_mime_association_node(line)

            yield {
                "_id": f"{self.parser_helpers.get_suffix(subject_node['id'])}_has_metabolic_interaction_with_{self.parser_helpers.get_suffix(object_node['id'])}",
                "association": association_node,
                "object": object_node,
                "subject": subject_node,
            }
            rec_ct += 1
        print(f"--- Finished Parsing Microbe-Metabolite. Generated {rec_ct} records. ---\n")

    def parse_metabolite_gene(self) -> Iterator[dict]:
        print("\n--- Parsing Metabolite-Gene Data ---")
        pubchem_descr = self.load_pickle("gmmad2_pubchem_descriptions.pkl")
        pubchem_mw = self.load_pickle("gmmad2_pubchem_mw.pkl")
        gene_info = self.load_pickle("gmmad2_protein_gene_combined_info.pkl")
        bigg_mapping = self.load_pickle("gmmad2_bigg_metabolite_mapping.pkl")
        pmid_metadata = self.load_pickle("gmmad2_pubmed_metadata.pkl")

        rec_ct = 0
        for line in self.csv_parser.line_generator(self.mege_path):
            mege_id_hierarchy = [
                (line[3], "PUBCHEM.COMPOUND"),
                (line[7], "DRUGBANK"),
                (line[5], "KEGG"),
                (line[10], "SMILES"),
                (line[6], "HMDB"),
                (bigg_mapping.get(line[2].lower()), "BIGG.METABOLITE"),
            ]
            primary_id, xrefs = self.parser_helpers.assign_primary_metabolite_id(mege_id_hierarchy)

            subject_node = self._get_metabolite_node(
                cid=line[3],
                name=line[2],
                formula=line[4],
                primary_id=primary_id,
                xrefs=xrefs,
                pubchem_desc=pubchem_descr,
                pubchem_mw=pubchem_mw,
            )

            object_node = self._get_gene_node(line, gene_info)

            association_node = self._get_mege_association_node(line, pmid_metadata)

            yield {
                "_id": f"{self.parser_helpers.get_suffix(subject_node['id'])}_interacts_with_{self.parser_helpers.get_suffix(object_node['id'])}",
                "association": association_node,
                "object": object_node,
                "subject": subject_node,
            }
            rec_ct += 1
        print(f"--- Finished Parsing Metabolite-Gene. Generated {rec_ct} records. ---\n")


class RecordCacheManager(CacheHelper):
    """Manages the creation of a combined cache from individual relationship files."""

    COMBINED_FILENAME = "gmmad2_parsed_records.pkl"

    def __init__(self, cache_dir=None):
        super().__init__(cache_dir)
        print("RecordCacheManager initialized.")

    def cache_combined_associations(self, data_loader):
        """
        Uses the DataLoader to parse each relationship and combines them into one file.
        The final format is a dictionary grouped by relationship type.
        """
        print("\nCreating full GMMAD2 association record cache...")

        combined_data = {
            "microbe-disease": list(data_loader.load_microbe_disease_data()),
            "microbe-metabolite": list(data_loader.load_microbe_metabolite_data()),
            "metabolite-gene": list(data_loader.load_metabolite_gene_data()),
        }

        self.save_pickle(combined_data, self.COMBINED_FILENAME)
        total_records = sum(len(v) for v in combined_data.values())
        print(
            f"Full GMMAD2 association record cache with {total_records} records created at {self.COMBINED_FILENAME}"
        )
        return combined_data


class DataLoader:
    def __init__(self, cache_dir="cache", downloads_dir="downloads"):
        csv_parser = CSVParser()
        cache_pipeline = DataCachePipeline()
        parser_helpers = ParserHelper()

        self.parser = GMMAD2Parser(
            csv_parser=csv_parser,
            cache_pipeline=cache_pipeline,
            parser_helpers=parser_helpers,
            cache_dir=cache_dir,
            downloads_dir=downloads_dir,
        )
        print("🚀 DataLoader initialized.")

    def load_microbe_disease_data(self) -> Iterator[dict]:
        """Loads and yields data from the microbe-disease file."""
        return self.parser.parse_microbe_disease()

    def load_microbe_metabolite_data(self) -> Iterator[dict]:
        """Loads and yields data from the microbe-metabolite file."""
        return self.parser.parse_microbe_metabolite()

    def load_metabolite_gene_data(self) -> Iterator[dict]:
        """Loads and yields data from the metabolite-gene file."""
        return self.parser.parse_metabolite_gene()

    def load_entire_gmmad2_data(self) -> Iterator[dict]:
        """Loads and yields all GMMAD2 data sources."""
        print("\n--- Loading all GMMAD2 data sources ---")
        return itertools.chain(
            self.load_microbe_disease_data(),
            self.load_microbe_metabolite_data(),
            self.load_metabolite_gene_data(),
        )


if __name__ == "__main__":
    data_loader = DataLoader()
    all_gmmad2_recs = data_loader.load_entire_gmmad2_data()

    record_cache_manager = RecordCacheManager()
    record_cache_manager.cache_combined_associations(data_loader=DataLoader())
