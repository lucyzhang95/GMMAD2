import json
from pathlib import Path
from typing import Any, Dict

from tqdm.auto import tqdm

from utils.cache_manager import CacheHelper, RecordHelper

from .record_deduplicator import MemoryEfficientDeduplicator, StreamDeduplicator


class RecordCacheManager(CacheHelper):
    """Manages the creation of a combined cache from individual relationship files."""

    COMBINED_FILENAME = "gmmad2_parsed_records"

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

    def create_deduplicated_jsonl_streamed(self, data_loader, use_memory_efficient=True):
        """
        Process all relationship data, deduplicate, and stream directly to JSONL file.

        Args:
            data_loader: The data loader instance
            use_memory_efficient: If True, uses disk-based deduplication for very large datasets
        """
        print(
            f"\nCreating deduplicated GMMAD2 records as JSONL (streaming, memory_efficient={use_memory_efficient})..."
        )
        print("=" * 50)

        if use_memory_efficient:
            deduplicator = MemoryEfficientDeduplicator(temp_dir=str(self.cache_dir / "temp_dedup"))
            print("Using memory-efficient deduplicator (disk-based)")
        else:
            deduplicator = StreamDeduplicator()
            print("Using in-memory deduplicator")

        output_path = Path(self.record_helper.cache_dir) / f"{self.COMBINED_FILENAME}.jsonl"

        total_processed = 0
        records_written = 0

        with open(output_path, "w", encoding="utf-8") as f:

            data_sources = [
                ("microbe-disease", data_loader.load_microbe_disease_data),
                ("microbe-metabolite", data_loader.load_microbe_metabolite_data),
                ("metabolite-gene", data_loader.load_metabolite_gene_data),
            ]

            for source_name, data_method in data_sources:
                print(f"\n>>> Processing {source_name} associations...")

                with tqdm(desc=f"Processing {source_name}") as pbar:
                    for record in data_method():
                        total_processed += 1

                        if use_memory_efficient:
                            processed_record = deduplicator.process_record(record)
                            if processed_record:
                                standardized_record = self._standardize_record(processed_record)
                                json_line = json.dumps(standardized_record, ensure_ascii=False)
                                f.write(json_line + "\n")
                                records_written += 1
                        else:
                            if not deduplicator.is_duplicate(record):
                                processed_record = deduplicator.process_record(record)
                                standardized_record = self._standardize_record(processed_record)
                                json_line = json.dumps(standardized_record, ensure_ascii=False)
                                f.write(json_line + "\n")
                                records_written += 1
                            else:
                                deduplicator.process_record(record)

                        pbar.update(1)

            if use_memory_efficient:
                print("\n>>> Writing remaining deduplicated records...")
                records_written = 0
                f.seek(0)
                f.truncate()

                with tqdm(desc="Writing final records") as pbar:
                    for record in deduplicator.get_all_records():
                        standardized_record = self._standardize_record(record)
                        json_line = json.dumps(standardized_record, ensure_ascii=False)
                        f.write(json_line + "\n")
                        records_written += 1
                        pbar.update(1)

                deduplicator.cleanup()

        duplicates_removed = total_processed - records_written

        print(f"\n[DONE] Processed {total_processed:,} total records")
        print(f"-> Removed {duplicates_removed:,} duplicates")
        print(f"-> Final unique records: {records_written:,}")
        print(f"-> Output file: {output_path}")

        print("\n[DONE] Streaming JSONL export complete.")
        print("Finished.\n")

        return str(output_path)

    def create_deduplicated_jsonl_simple_streaming(self, data_loader):
        """
        Simplified streaming approach - processes records one by one without deduplication.
        Use this if you're confident there are no duplicates or want maximum speed.
        """
        print("\nCreating GMMAD2 records as JSONL (simple streaming, no deduplication)...")
        print("=" * 50)

        output_path = (
            Path(self.record_helper.cache_dir) / f"{self.COMBINED_FILENAME}_no_dedup.jsonl"
        )
        total_processed = 0

        with open(output_path, "w", encoding="utf-8") as f:
            data_sources = [
                ("microbe-disease", data_loader.load_microbe_disease_data),
                ("microbe-metabolite", data_loader.load_microbe_metabolite_data),
                ("metabolite-gene", data_loader.load_metabolite_gene_data),
            ]

            for source_name, data_method in data_sources:
                print(f"\n>>> Processing {source_name} associations...")

                with tqdm(desc=f"Processing {source_name}") as pbar:
                    for record in data_method():
                        standardized_record = self._standardize_record(record)
                        json_line = json.dumps(standardized_record, ensure_ascii=False)
                        f.write(json_line + "\n")
                        total_processed += 1
                        pbar.update(1)

        print(f"\n[DONE] Processed {total_processed:,} total records")
        print(f"-> Output file: {output_path}")
        print("\n[DONE] Simple streaming export complete.")

        return str(output_path)

    def create_deduplicated_jsonl(self, data_loader):
        """
        Original method - loads all records into memory.
        Use create_deduplicated_jsonl_streaming() for better memory efficiency.
        """
        return self.create_deduplicated_jsonl_streamed(data_loader, use_memory_efficient=False)
