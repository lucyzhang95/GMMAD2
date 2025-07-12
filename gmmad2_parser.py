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
from tqdm.asyncio import tqdm, tqdm_asyncio


class CacheManager:
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

    def save_pickle(self, obj, f_name: str) -> str:
        """Saves an object to a pickle file in the cache directory."""
        path = self._get_path(f_name)
        with open(path, "wb") as out_f:
            pickle.dump(obj, out_f)
        print(f"{path} saved to cache.")
        return path

    def load_pickle(self, f_name: str):
        """Loads an object from a pickle file. Returns None if the file doesn't exist."""
        path = self._get_path(f_name)
        if not os.path.exists(path):
            print(f"Cache miss: {path} does not exist.")
            return None

        try:
            with open(path, "rb") as in_f:
                return pickle.load(in_f)
        except (pickle.UnpicklingError, EOFError) as e:
            print(f"Error loading pickle file {path}: {e}")
            return None

    def save_json(self, obj, f_name: str) -> str:
        """Saves an object to a JSON file in the cache directory."""
        path = self._get_path(f_name)
        with open(path, "w", encoding="utf-8") as out_f:
            json.dump(obj, out_f, indent=4)
        print(f"Saved: {path}")
        return path

    def cache_data(self, f_path, gzip_path=None, bigg_path=None):
        if gzip_path is None:
            gzip_path =





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

    def pug_query_pubchem_description(self, cid, session, sem, max_retries=3, delay=1.0):
        # TODO: implement PUG-REST query
        pass

    async def get_batch_pubchem_descriptions_async(self, cids, workers=5):
        # TODO: async batch fetch
        pass

    def get_pubchem_descriptions(self, cids, workers=5):
        # TODO: sync batch fetch
        pass


class UniProtService:
    """UniProt REST API helpers."""

    async def uniprot_query_protein_info(
        self, uniprot_id, session, sem, max_retries=3, delay=1.0, timeout=30.0
    ):
        # TODO: implement UniProt REST call
        pass

    async def get_batch_protein_info_async(self, uniprot_ids, workers=5, rate_limit=5.0):
        # TODO: async batch UniProt queries
        pass

    def get_protein_info(self, uniprot_ids, workers=5):
        # TODO: sync wrapper
        pass


class BiGGService:
    """BiGG metabolite mapping helper."""

    def get_bigg_metabolite_mapping(self, in_f):
        # TODO: parse BiGG mapping file
        pass


class CheminformaticsUtils:
    def bt_get_mw_logp(self, pubchem_cids):
        # TODO: compute molecular weight, logP
        pass

    def get_organism_type(self, node):
        # TODO: determine organism type
        pass

    def get_suffix(self, identifier):
        # TODO: parse identifier suffix
        pass


class ParserHelpers:
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
        parser_helpers: ParserHelpers,
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
