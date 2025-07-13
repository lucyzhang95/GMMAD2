import asyncio
import csv
import json
import os
import pickle
import tarfile
import time
import uuid
import zipfile
from collections import Counter
from pathlib import Path
from typing import Dict, Iterator, List

import aiohttp
import biothings_client as bt
import pandas as pd
import requests
from Bio import Entrez
from dotenv import load_dotenv
from ete3 import NCBITaxa
from tqdm.asyncio import tqdm_asyncio


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

    def get_merged_taxid_mapping(self, tar_gz_path: str, f_name: str = "merged.dmp") -> dict:
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

    def get_current_taxid(self, old_taxids: list, merged_mapping: dict) -> dict[str, str]:
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
        return taxon_info

    def fetch_taxon_names_from_taxon_info(self, taxon_info: dict) -> list[str]:
        """Extracts biothings names from the taxon_info dictionary."""
        taxon_names = set()
        for _, taxon in taxon_info.items():
            if "scientific_name" in taxon:
                taxon_names.add(taxon["scientific_name"].lower())
        return list(taxon_names)


class NCITTaxonomyService:
    """NCIT-based organism mappings and descriptions."""

    def __init__(self):
        self.NCIT_API_KEY = os.getenv("NCIT_API_KEY")

    def dump_ncit_source(self, name, out_path):
        # TODO: dump NCIT source
        pass

    def build_ncit_organism_mapping(self, f_name, path):
        # TODO: build mapping
        pass

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
        unique_names = {n.lower() for n in taxon_names}
        sem = asyncio.Semaphore(max_concurrent)
        connector = aiohttp.TCPConnector(limit_per_host=max_concurrent)
        async with aiohttp.ClientSession(connector=connector) as session:
            tasks = [
                self.async_query_ncit_taxon_description(session, name, sem) for name in unique_names
            ]
            results = await asyncio.gather(*tasks)
        return {name: data for item in results if item for name, data in (item,)}

    def run_async_query_ncit_taxon_descriptions(self, taxon_names):
        """

        :param taxon_names:
        :return:
        {
           "550":{
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

    def update_taxon_info_with_ncit_description(self, taxon_info: dict, descriptions: dict) -> dict:
        """

        :param taxon_info: output of NCBITaxonomyService.query_taxon_info_from_biothings
        :param descriptions: output of NCITTaxonomyService.run_async_query_ncit_taxon_descriptions
        :return:
        """
        for info in taxon_info.values():
            name = info.get("scientific_name").lower()
            descr_info = descriptions.get(name, {})
            info.update(descr_info)
        return taxon_info


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

    async def async_query_pubchem_descriptions(
        self,
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

                task = asyncio.create_task(
                    self.async_query_pug_pubchem_description(cid, session, sem)
                )
                tasks.append(task)

            for cid, payload in await tqdm_asyncio.gather(*tasks, total=len(tasks)):
                if payload:
                    results[cid] = payload

        return results

    def run_async_query_pubchem_descriptions(self, cids: List[int], workers: int = 5):
        return asyncio.run(self.async_query_pubchem_descriptions(cids, workers))


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

                task = asyncio.create_task(
                    self.async_query_uniprot_name_and_function(uid, session, sem)
                )
                tasks.append(task)

            for uid, payload in await tqdm_asyncio.gather(*tasks, total=len(tasks)):
                results[uid] = payload

        return results

    def run_async_query_uniprot_names_and_functions(
        self, uniprot_ids: List[str], workers: int = 5
    ) -> Dict[str, Dict[str, str]]:
        return asyncio.run(self.async_query_uniprot_names_and_functions(uniprot_ids, workers))


class BiGGParser:
    """BiGG metabolite mapping helper."""

    def __init__(self, csv_parser):
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


class PubMedService:
    def query_pubmed_metadata(self, pmids):
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
            print(f"{path} does not exist.")
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

    CACHE_FILENAME = "gmmad2_all_association_cache.pkl"

    def __init__(self, cache_dir=None):
        """Initializes the CacheManager and ensures the cache directory exists."""
        super().__init__(cache_dir)
        print(f"CacheManager initialized. Cache directory is {self.cache_dir}")

    def cache_entity(self, entity_type: str, **kwargs):
        """
        Build and cache data for a single entity type, updating the existing
        cache file instead of overwriting it.
        """
        entity_map = {
            "taxon_info": self._cache_taxon_info,
            "taxon_description": self._cache_taxon_description,
            "pubchem_description": self._cache_pubchem_description,
            "pubchem_mw": self._cache_pubchem_mw,
            "bigg_mapping": self._cache_bigg_mapping,
            "uniprot": self._cache_uniprot,
            "pubmed_metadata": self._cache_pubmed_metadata,
        }

        handler = entity_map.get(entity_type)
        if not handler:
            raise ValueError(f"Unknown entity type: {entity_type!r}")

        print(f"\nCaching entity: '{entity_type}'...")
        new_data = handler(**kwargs)

        if new_data:
            f_name = f"{entity_type}.pkl"
            existing_data = self.load_pickle(f_name) or {}

            if not isinstance(existing_data, dict):
                print(f"Warning: Data in {f_name} is not a dictionary. Overwriting.")
                existing_data = {}
            print(f"Updating cache for {f_name}...")
            existing_data.update(new_data)
            self.save_pickle(existing_data, f_name)
        else:
            print(f"No new data generated for '{entity_type}'. Cache not updated.")

    # TODO: Need futher detailed implementations for each caching method
    def _cache_taxon_info(self, **kwargs):
        print("Generating taxon info...")
        return NCBITaxonomyService.query_taxon_info_from_biothings(**kwargs)

    def _cache_taxon_description(self, **kwargs):
        print("Generating taxon descriptions...")
        pass

    def _cache_pubchem_description(self, **kwargs):
        print("Generating PubChem descriptions...")
        pass

    def _cache_pubchem_mw(self, **kwargs):
        print("Generating PubChem molecular weights...")
        pass

    def _cache_bigg_mapping(self, **kwargs):
        print("Generating BiGG mappings...")
        pass

    def _cache_uniprot(self, **kwargs):
        print("Generating UniProt data...")
        pass

    def _cache_pubmed_metadata(self, **kwargs):
        print("Generating PMID data...")
        pass


class CombinedCacheManager(CacheHelper):
    """Manages the creation of a combined cache from individual relationship files."""

    COMBINED_FILENAME = "gmmad2_combined_associations.pkl"

    def __init__(self, cache_dir=None):
        super().__init__(cache_dir)
        print("CombinedCacheManager initialized.")

    def create_combined_cache(self):
        """
        Loads individual relationship caches and combines them into one file.
        """
        print("\nCreating combined relationship cache...")

        self._cache_relationship("microbe-disease")
        self._cache_relationship("microbe-metabolite")
        self._cache_relationship("metabolite-gene")

        microbe_disease_data = self.load_pickle("microbe-disease.pkl") or {}
        microbe_metabolite_data = self.load_pickle("microbe-metabolite.pkl") or {}
        metabolite_gene_data = self.load_pickle("metabolite-gene.pkl") or {}

        combined_data = {
            "microbe_disease": microbe_disease_data,
            "microbe_metabolite": microbe_metabolite_data,
            "metabolite_gene": metabolite_gene_data,
        }

        self.save_pickle(combined_data, self.COMBINED_FILENAME)
        print(f"Combined cache created at {self.COMBINED_FILENAME}")
        return combined_data

    def _cache_relationship(self, relationship_name: str):
        """Placeholder to generate and save a single relationship cache."""
        print(f"-> Generating temporary cache for '{relationship_name}'...")
        if relationship_name == "microbe-disease":
            data = {"microbeA": ["diseaseX"], "microbeB": ["diseaseY"]}
        elif relationship_name == "microbe-metabolite":
            data = {"microbeA": ["metabolite1"], "microbeC": ["metabolite2"]}
        elif relationship_name == "metabolite-gene":
            data = {"metabolite1": ["geneZ"], "metabolite2": ["geneW"]}
        else:
            data = {}

        self.save_pickle(data, f"{relationship_name}.pkl")


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

    def assign_metabolite_primary_id(self, line, bigg_map):
        pubchem_id = line[6]
        kegg_id = line[8]
        smiles = line[18]
        hmdb_id = line[19]
        bigg_id = bigg_map.get(line[5].lower())

        # (line, prefix) pairs for ID hierarchy
        id_hierarchy = [
            (pubchem_id, "PUBCHEM.COMPOUND"),
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


class GMMAD2Parser:
    """
    Parses and merges microbe-disease, microbe-metabolite,
    metabolite-gene associations.
    """

    def __init__(
        self,
        csv_parser: CSVParser,
        cache_mgr: CacheManager,
        taxonomy_svc: NCBITaxonomyService,
        ncit_svc: NCITTaxonomyService,
        parser_helpers: ParserHelper,
    ):
        load_dotenv()
        self.csv = csv_parser
        self.cache = cache_mgr
        self.taxonomy = taxonomy_svc
        self.ncit = ncit_svc
        self.parser_helpers = parser_helpers

    def parse_microbe_disease(self, in_file):
        # TODO: implement parsing
        pass

    def parse_microbe_metabolite(self, in_file):
        # TODO: implement parsing
        pass

    def parse_metabolite_gene(self, in_file):
        # TODO: implement parsing
        pass


class DataLoader:
    """Data loaders."""

    def load_micro_disease_data(self, f_path):
        # TODO: parse microbe-disease table
        pass

    def load_micro_meta_data(self, f_path):
        # TODO: parse microbe-metabolite table
        pass

    def load_entire_gmmad2_data(self, f_path):
        # TODO: concatenate all generators into a single iterator
        pass
