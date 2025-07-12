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
        # TODO: implement generic caching logic
        pass


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

    def ensure_ncbi_taxdump(self, tar_gz_path: str | Path) -> Path:
        tar_gz_path = Path(tar_gz_path)  # TODO: need to specify the path later
        if tar_gz_path.exists():
            print(f"NCBI taxdump found at {tar_gz_path}, using existing file.")
            return tar_gz_path

        print("NCBI taxdump not found – downloading via ETE3 …")
        ncbi = NCBITaxa()
        ncbi.update_taxonomy_database()
        return tar_gz_path

    def get_ncbi_taxdump_path(self, tar_gz_path=None) -> str:
        """Ensures the NCBI taxdump.tar.gz file exists and returns its path.
        Downloads the file if it's not found.
        """
        return ensure_ncbi_taxdump(tar_gz_path)

    def parse_taxid_from_ncbi_taxdump_merged(
        self, tar_gz_path: str, f_name: str = "merged.dmp"
    ) -> dict:
        """
        Parses 'merged.dmp' from the provided taxdump tarball.
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

    def get_current_taxid(self, old_taxids, merged_mapping):
        # TODO: map old IDs to current
        pass

    def get_taxon_info(self, taxids):
        # TODO: fetch taxon info
        pass

    def get_taxon_names(self, taxon_info):
        # TODO: extract names
        pass


class NCITTaxonomyService:
    """NCIT-based organism mappings and descriptions."""

    def dump_ncit_source(self, name, out_path):
        # TODO: dump NCIT source
        pass

    def build_ncit_organism_mapping(self, f_name, path):
        # TODO: build mapping
        pass

    def fetch_ncit_description(self, session, name, sem):
        # TODO: synchronous fetch
        pass

    async def get_ncit_taxon_description_async(self, taxon_names, max_concurrent=5):
        # TODO: async fetch
        pass

    def get_ncit_taxon_description(self, taxon_names):
        # TODO: batch fetch
        pass

    def add_description_to_taxon_info(self, taxon_info, descriptions):
        # TODO: annotate taxon info
        pass


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


class GMMAD2Parser:
    """
    Parses and merges microbe-disease, microbe-metabolite,
    metabolite-gene associations.
    """

    def __init__(
        self,
        csv_parser: CSVParser,
        cache_mgr: CacheManager,
        taxonomy_svc: TaxonomyService,
        ncit_svc: NCITService,
    ):
        load_dotenv()
        self.csv = csv_parser
        self.cache = cache_mgr
        self.taxonomy = taxonomy_svc
        self.ncit = ncit_svc

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
