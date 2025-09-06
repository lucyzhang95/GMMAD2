from typing import Any, Dict

from tqdm.auto import tqdm

from utils.cache_manager import CacheHelper, RecordHelper

from .record_deduplicator import StreamDeduplicator


class RecordCacheManager(CacheHelper):
    """Manages the creation of a combined cache from individual relationship files."""

    COMBINED_FILENAME = f"gmmad2_parsed_records"

    def __init__(self, cache_dir=None):
        super().__init__(cache_dir)
        self.record_helper = RecordHelper()
        print("RecordCacheManager initialized.")
        print("=" * 50)

    def _standardize_record(self, record: Dict[str, Any]) -> Dict[str, Any]:
        """Standardizes a record by cleaning up empty values and formatting xrefs into a list."""
        clean_record = record.copy()

        for entity_key in ["subject", "object"]:
            if entity_key in clean_record and isinstance(clean_record[entity_key], dict):
                entity = clean_record[entity_key].copy()
                if "xrefs" in entity and isinstance(entity["xrefs"], dict):
                    entity["xrefs"] = [v for k, v in entity["xrefs"].items() if v]
                clean_record[entity_key] = entity

        return clean_record

    def create_deduplicated_jsonl(self, data_loader):
        """
        Process all relationship data, deduplicate by _id, and export to single JSONL file.
        Each record is one line in the JSONL output file.
        """
        print("\nCreating deduplicated GMMAD2 records as JSONL...")
        print("=" * 50)

        deduplicator = StreamDeduplicator()
        total_processed = 0

        print("\n>>> Processing microbe-disease associations...")
        with tqdm(desc="Processing microbe-disease") as pbar:
            for record in data_loader.load_microbe_disease_data():
                deduplicator.process_record(record)
                pbar.update(1)
                total_processed += 1

        print("\n>>> Processing microbe-metabolite associations...")
        with tqdm(desc="Processing microbe-metabolite") as pbar:
            for record in data_loader.load_microbe_metabolite_data():
                deduplicator.process_record(record)
                pbar.update(1)
                total_processed += 1

        print("\n>>> Processing metabolite-gene associations...")
        with tqdm(desc="Processing metabolite-gene") as pbar:
            for record in data_loader.load_metabolite_gene_data():
                deduplicator.process_record(record)
                pbar.update(1)
                total_processed += 1

        deduplicated_records = list(deduplicator.get_results())
        standardized_records = [self._standardize_record(record) for record in deduplicated_records]

        duplicates_removed = total_processed - len(standardized_records)

        print(f"\n[DONE] Processed {total_processed:,} total records")
        print(f"-> Removed {duplicates_removed:,} duplicates")
        print(f"-> Final unique records: {len(standardized_records):,}")

        jsonl_path = self.record_helper.save_jsonl(standardized_records, self.COMBINED_FILENAME)

        print("\n[DONE] JSONL export complete.")
        print("Finished.\n")

        return jsonl_path
