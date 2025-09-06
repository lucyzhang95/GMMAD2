import json
import os
import pickle
from datetime import datetime
from typing import Any, Dict, List

from .ontology_mapper import BiGGMapper
from .ontology_services import (
    ChemPropertyServices,
    GeneOntologyService,
    NCBITaxonomyService,
    NCItTaxonomyService,
    PubChemService,
    PubMedService,
    UniProtService,
)
from .parser_helper import ParserHelper


class CacheHelper:
    """Generic on-disk cache and serialization helpers."""

    DEFAULT_CACHE_DIR = os.path.join(os.getcwd(), "cache")

    def __init__(self, cache_dir=None):
        """Initializes the CacheManager and ensures the cache directory exists."""
        self.cache_dir = cache_dir or self.DEFAULT_CACHE_DIR
        os.makedirs(self.cache_dir, exist_ok=True)
        print(f"▶️CacheManager initialized. Cache directory is {self.cache_dir}")

    def _get_path(self, f_name: str) -> str:
        """Helper method to construct the full path for a cache file."""
        return os.path.join(self.cache_dir, f_name)

    def save_pickle(self, obj, f_name: str):
        """Saves an object to a pickle file in the cache directory."""
        path = self._get_path(f_name)
        try:
            with open(path, "wb") as out_f:
                pickle.dump(obj, out_f)
            print(f"✅Saved pickle to: {path}")
        except (IOError, pickle.PicklingError) as e:
            print(f"‼️Error saving pickle file {path}: {e}")

    def load_pickle(self, f_name: str):
        """Loads an object from a pickle file. Returns None if the file doesn't exist."""
        path = self._get_path(f_name)
        if not os.path.exists(path):
            return None
        try:
            with open(path, "rb") as in_f:
                return pickle.load(in_f)
        except (pickle.UnpicklingError, EOFError) as e:
            print(f"‼️Error loading pickle file {path}: {e}")
            return None


class RecordHelper:
    """Handles saving and loading of records to/from JSON and JSONL files."""

    RECORD_DIR = os.path.join("records")

    def __init__(self, rec_dir=RECORD_DIR):
        self.rec_dir = os.path.join(os.getcwd(), rec_dir)
        os.makedirs(self.rec_dir, exist_ok=True)

    def save_json(self, obj, f_name):
        """Saves an object to a JSON file."""
        with open(os.path.join(self.rec_dir, f_name), "w") as out_f:
            json.dump(obj, out_f, indent=4)

    def load_json(self, f_name):
        """Loads an object from a JSON file."""
        path = os.path.join(self.rec_dir, f_name)
        if os.path.exists(path):
            with open(path, "r") as in_f:
                return json.load(in_f)
        return None

    def save_jsonl(self, records: List[Dict[str, Any]], f_name: str):
        """Saves records to a JSONL file with standardized formatting."""
        if not records:
            print("!!! Warning: No records provided for JSONL export.")
            return None

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base_name = f_name.replace(".jsonl", "") if f_name.endswith(".jsonl") else f_name
        jsonl_filename = f"{base_name}_{timestamp}.jsonl"
        jsonl_path = os.path.join(self.rec_dir, jsonl_filename)

        exported_count = 0
        with open(jsonl_path, "w", encoding="utf-8") as f:
            for record in records:
                json.dump(record, f, ensure_ascii=False, separators=(",", ":"))
                f.write("\n")
                exported_count += 1

        print(f"-> Exported {exported_count} records to {jsonl_path}")
        return jsonl_path

    def load_jsonl(self, f_name: str) -> List[Dict[str, Any]]:
        """Load records from JSONL file."""
        path = os.path.join(self.rec_dir, f_name)

        if not os.path.exists(path):
            return []

        records = []
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line:
                    try:
                        records.append(json.loads(line))
                    except json.JSONDecodeError:
                        continue

        print(f"-> Loaded {len(records)} records from {f_name}")
        return records


class CacheManager(CacheHelper):
    """CacheManager for managing on-disk caches and serialization."""

    def __init__(self, cache_dir=None):
        """Initializes the CacheManager and ensures the cache directory exists."""
        super().__init__(cache_dir)
        self.COMBINED_GENE_CACHE_F_NAME = "gmmad2_protein_gene_combined_info.pkl"
        self.parser_helper = ParserHelper()
        self.ncbi_service = NCBITaxonomyService()

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

        print(f"\n▶️---Caching entity: '{entity_type}' ---")
        f_name = f"gmmad2_{entity_type}.pkl"
        existing_data = self.load_pickle(f_name) or {}
        if not isinstance(existing_data, dict):
            print(f"‼️Warning: Data in {f_name} is not a dictionary. Overwriting.")
            existing_data = {}

        kwargs["existing_data"] = existing_data
        new_data = handler(**kwargs)
        if new_data:
            print(f"▶️Received {len(new_data)} new items to cache.")
            existing_data.update(new_data)
            self.save_pickle(existing_data, f_name)
        else:
            print(f"❌No new data to cache for '{entity_type}'.")

    def _cache_taxon_info(self, **kwargs):
        taxids = kwargs.get("taxids", [])
        if not taxids:
            return None

        f_name = "gmmad2_taxon_info.pkl"
        existing_data = self.load_pickle(f_name) or {}
        taxids_to_query = [tid for tid in taxids if tid not in existing_data]
        if not taxids_to_query:
            print("✅All requested taxids are already in the cache.")
            return None

        taxon_info = self.ncbi_service.query_taxon_info_from_biothings(taxids=taxids_to_query)
        filtered_taxon_info, taxid_notfound = self.parser_helper.filter_taxon_info(taxon_info)

        if taxid_notfound:
            print(
                f"-> Found {len(taxid_notfound)} taxids that were not found: {taxid_notfound[:5]}"
            )
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
        print("\n⚙️---Updating cache with re-mapped taxids---")

        notfound_f_name = "gmmad2_taxon_info_notfound.pkl"
        main_cache_f_name = "gmmad2_taxon_info.pkl"
        taxid_notfound = self.load_pickle(notfound_f_name) or []
        main_taxon_info_cache = self.load_pickle(main_cache_f_name) or {}

        old_taxids_to_process = [tid for tid in taxid_notfound if tid not in main_taxon_info_cache]
        if not old_taxids_to_process:
            print("✅All 'not found' taxids are already in the cache.")
            return

        taxid_map = self.ncbi_service.get_taxid_mapping_from_ncbi_merged_dmp()
        new_taxid_map = self.ncbi_service.get_current_taxid_mapping(
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
            new_taxon_info = self.ncbi_service.query_taxon_info_from_biothings(
                taxids=new_taxids_to_query
            )
            valid_new_taxon_info, _ = self.parser_helper.filter_taxon_info(new_taxon_info)

            reverse_map = {v: k for k, v in new_taxid_map.items()}
            for new_id, new_info in valid_new_taxon_info.items():
                old_id = reverse_map.get(new_id)
                if old_id:
                    new_taxon_info_to_add[old_id] = new_info
                new_taxon_info_to_add[new_id] = new_info

        if not info_to_add_from_remap and not new_taxon_info_to_add:
            print("❌No new information to add to the cache.")
            return

        print(
            f"▶️Adding {len(info_to_add_from_remap) + len(new_taxon_info_to_add)} new/remapped entries to the main cache."
        )
        main_taxon_info_cache.update(info_to_add_from_remap)
        main_taxon_info_cache.update(new_taxon_info_to_add)
        self.save_pickle(main_taxon_info_cache, main_cache_f_name)

    def update_taxon_info_with_descriptions(self):
        """
        Loads the main taxon info cache, finds entries missing a description,
        queries for them, and saves the updated data back to the same file.
        """
        print("\n⚙️---Updating taxon info with NCIT descriptions---")

        main_cache_f_name = "gmmad2_taxon_info.pkl"
        desc_cache_f_name = "gmmad2_ncit_descriptions.pkl"

        taxon_info_cache = self.load_pickle(main_cache_f_name)
        if not taxon_info_cache:
            print(f"❌Main taxon info cache not found at '{main_cache_f_name}'. Cannot update.")
            return

        print(f"✅Loaded '{main_cache_f_name}' with {len(taxon_info_cache)} keys.")

        desc_cache = self.load_pickle(desc_cache_f_name) or {}
        names_need_desc = set()
        for info in taxon_info_cache.values():
            if "description" not in info or not info.get("description"):
                name = info["scientific_name"]
                names_need_desc.add(name)

        if not names_need_desc:
            print("✅All taxon entries already have descriptions.")
            return

        names_to_query = sorted([name for name in names_need_desc if name not in desc_cache])

        if names_to_query:
            ncit_service = NCItTaxonomyService()
            new_descriptions = ncit_service.run_async_query_ncit_taxon_descriptions(
                taxon_names=names_to_query
            )
            if new_descriptions:
                print(f"-> Found {len(new_descriptions)} new descriptions from API.")
                desc_cache.update(new_descriptions)
                self.save_pickle(desc_cache, desc_cache_f_name)
        else:
            print("✅All needed descriptions were already in the local description cache.")

        update_count = 0
        for taxid, info in taxon_info_cache.items():
            if ("description" not in info or not info.get("description")) and info.get(
                "scientific_name"
            ):
                name = info["scientific_name"]
                if name in desc_cache:
                    desc_dict = desc_cache[name]
                    taxon_info_cache[taxid]["description"] = (
                        desc_dict.get("description", "")
                        if isinstance(desc_dict, dict)
                        else str(desc_dict)
                    )
                    if isinstance(desc_dict, dict) and "xrefs" in desc_dict:
                        taxon_info_cache[taxid]["xrefs"] = {
                            **taxon_info_cache[taxid].get("xrefs", {}),
                            **desc_dict["xrefs"],
                        }
                    update_count += 1
        if update_count > 0:
            print(
                f"⚙️Updated {update_count} entries in the main taxon info cache with descriptions."
            )
            print(f"-> Saving '{main_cache_f_name}' with {len(taxon_info_cache)} keys.")
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
            print("✅All requested pubchem_cids are already in the cache.")
            return

        pubchem_service = PubChemService()
        pubchem_desc = pubchem_service.run_async_query_pug_pubchem_descriptions(cids=cids_to_query)

        if not pubchem_desc:
            print("-> API query did not return any new descriptions.")
            return
        print(f"-> Received {len(pubchem_desc)} new PubChem descriptions to cache.")
        existing_data.update(pubchem_desc)

        self.save_pickle(existing_data, f_name)

    def _cache_pubchem_mw(self, **kwargs):
        print("\n▶️---Caching PubChem Metabolite Physical Properties---")
        cids = kwargs.get("cids", [])
        if not cids:
            return

        existing_data = self.load_pickle("gmmad2_pubchem_mw.pkl") or {}
        cids_to_query = [cid for cid in cids if cid not in existing_data]
        if not cids_to_query:
            print("▶️All requested pubchem_cids are already in the cache.")
            return

        chem_property_service = ChemPropertyServices()
        new_chem_properties = chem_property_service.query_mw_and_xlogp_from_biothings(
            cids=cids_to_query
        )
        if not new_chem_properties:
            print("-> API query did not return any new properties.")
            return
        print(f"-> Received {len(new_chem_properties)} new PubChem properties to cache.")
        existing_data.update(new_chem_properties)
        self.save_pickle(existing_data, "gmmad2_pubchem_mw.pkl")

    def _cache_bigg_mapping(self, **kwargs):
        print("\n▶️---Caching BiGG Metabolite Mapping---")
        in_f = kwargs.get("downloads", "bigg_models_metabolites.txt")
        if not os.path.exists(in_f):
            print(f"❌BiGG mapping file '{in_f}' does not exist. Skipping caching.")
            return

        bigg_parser = BiGGMapper()
        bigg_map = bigg_parser.get_bigg_metabolite_mapping(in_f=in_f)

        if not bigg_map:
            print("❌No BiGG mapping data found.")
            return

        self.save_pickle(bigg_map, "gmmad2_bigg_metabolite_mapping.pkl")
        print(f"✅BiGG mapping cached with {len(bigg_map)} entries.")

    def _cache_uniprot_info(self, **kwargs):
        print("\n▶️---Caching UniProt Information---")
        uniprot_ids = kwargs.get("uniprots", [])
        if not uniprot_ids:
            return
        existing_data = self.load_pickle(self.COMBINED_GENE_CACHE_F_NAME) or {}
        uids_to_query = [uid for uid in uniprot_ids if uid not in existing_data]
        if not uids_to_query:
            print("✅All requested uniprot ids are already in the cache.")
            return

        uniprot_service = UniProtService()
        new_protein_info = uniprot_service.run_async_query_uniprot_names_and_functions(
            uniprot_ids=uids_to_query
        )

        if not new_protein_info:
            print("-> API query did not return any new descriptions.")
            return
        print(f"-> Received {len(new_protein_info)} new Uniprot entires to cache.")
        existing_data.update(new_protein_info)
        self.save_pickle(existing_data, self.COMBINED_GENE_CACHE_F_NAME)

    def _cache_gene_info(self, **kwargs):
        print("\n▶️---Caching Gene Information---")
        gene_ids = kwargs.get("gene_ids", [])
        if not gene_ids:
            return

        existing_data = self.load_pickle(self.COMBINED_GENE_CACHE_F_NAME) or {}
        gids_to_query = [gid for gid in gene_ids if gid not in existing_data]
        if not gids_to_query:
            print("✅All requested gene ids are already in the cache.")
            return

        go_service = GeneOntologyService()
        new_gene_info = go_service.query_gene_name_from_biothings(gene_ids=gids_to_query)

        if not new_gene_info:
            print("-> API query did not return any new descriptions.")
            return
        print(f"-> Received {len(new_gene_info)} new Gene descriptions to cache.")
        existing_data.update(new_gene_info)
        self.save_pickle(existing_data, self.COMBINED_GENE_CACHE_F_NAME)

    def _cache_pubmed_metadata(self, **kwargs):
        print("\n▶️---Caching PubMed Metadata---")
        pmids = kwargs.get("pmids", [])
        if not pmids:
            return

        existing_data = self.load_pickle("gmmad2_pubmed_metadata.pkl") or {}
        pmids_to_query = [pmid for pmid in pmids if pmid not in existing_data]
        if not pmids_to_query:
            print("✅All requested PMIDs are already in the cache.")
            return

        pubmed_service = PubMedService()
        new_pubmed_metadata = pubmed_service.query_pubmed_metadata(pmids=pmids_to_query)

        if not new_pubmed_metadata:
            print("-> API query did not return any new metadata.")
            return
        print(f"-> Received {len(new_pubmed_metadata)} new PubMed entries to cache.")
        existing_data.update(new_pubmed_metadata)
        self.save_pickle(existing_data, "gmmad2_pubmed_metadata.pkl")
