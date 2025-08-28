from tqdm.auto import tqdm

from utils.cache_manager import CacheHelper

from .record_deduplicator import StreamDeduplicator


class RecordCacheManager(CacheHelper):
    """Manages the creation of a combined cache from individual relationship files."""

    COMBINED_FILENAME = "gmmad2_parsed_records"

    def __init__(self, cache_dir=None):
        super().__init__(cache_dir)
        print("RecordCacheManager initialized.")

    def cache_combined_associations(self, data_loader):
        """
        Uses the DataLoader to parse each relationship and combines them into one file.
        The final format is a dictionary grouped by relationship type.
        """
        print("\nCreating full GMMAD2 association record cache...")

        combined_data = {}
        raw_counts = {}

        print("\n▶️ Processing microbe-disease associations...")
        deduplicator_md = StreamDeduplicator()
        record_iterator_md = data_loader.load_microbe_disease_data()

        with tqdm(desc="Deduplicating microbe-disease") as pbar:
            count = 0
            for record in record_iterator_md:
                deduplicator_md.process_record(record)
                pbar.update(1)
                count += 1
            raw_counts["microbe-disease"] = count
        combined_data["microbe-disease"] = list(deduplicator_md.get_results())

        print("\n▶️ Processing microbe-metabolite associations...")
        deduplicator_mm = StreamDeduplicator()
        record_iterator_mm = data_loader.load_microbe_metabolite_data()
        with tqdm(desc="Deduplicating microbe-metabolite") as pbar:
            count = 0
            for record in record_iterator_mm:
                deduplicator_mm.process_record(record)
                pbar.update(1)
                count += 1
            raw_counts["microbe-metabolite"] = count
        combined_data["microbe-metabolite"] = list(deduplicator_mm.get_results())

        print("\n▶️ Processing metabolite-gene associations...")
        deduplicator_mg = StreamDeduplicator()
        record_iterator_mg = data_loader.load_metabolite_gene_data()
        with tqdm(desc="Deduplicating metabolite-gene") as pbar:
            count = 0
            for record in record_iterator_mg:
                deduplicator_mg.process_record(record)
                pbar.update(1)
                count += 1
            raw_counts["metabolite-gene"] = count
        combined_data["metabolite-gene"] = list(deduplicator_mg.get_results())

        print("\n🎉 Deduplication complete. Finalizing cache...")
        for key, dedup_list in combined_data.items():
            raw_count = raw_counts.get(key, 0)
            removed = raw_count - len(dedup_list)
            print(
                f"-> {key}: Processed {raw_count:,} records, removed {removed:,} duplicates ({len(dedup_list):,} unique records)."
            )

        self.save_pickle(combined_data, f"{self.COMBINED_FILENAME}.pkl")
        self.save_json(combined_data, f"{self.COMBINED_FILENAME}.json")
        total_records = sum(len(v) for v in combined_data.values())
        print(
            f"\n🎉 Full GMMAD2 association record cache with {total_records:,} records saved at {self.cache_dir}."
        )
        print("Finished.\n")
        return combined_data
