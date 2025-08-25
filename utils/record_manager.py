from utils.cache_manager import CacheHelper


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

        combined_data = {
            "microbe-disease": list(data_loader.load_microbe_disease_data()),
            "microbe-metabolite": list(data_loader.load_microbe_metabolite_data()),
            "metabolite-gene": list(data_loader.load_metabolite_gene_data()),
        }

        self.save_pickle(combined_data, f"{self.COMBINED_FILENAME}.pkl")
        self.save_json(combined_data, f"{self.COMBINED_FILENAME}.json")
        total_records = sum(len(v) for v in combined_data.values())
        print(
            f"🎉Full GMMAD2 association record cache with {total_records} records saved at {self.cache_dir}."
        )
        print("Finished.\n")
        return combined_data
