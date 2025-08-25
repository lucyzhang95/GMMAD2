import os

from utils.cache_manager import CacheManager
from utils.ontology_mapper import BiGGMapper
from utils.reader import CSVReader


class DataCachePipeline:
    def __init__(self, cache_dir="cache", downloads_dir="downloads"):
        self.downloads_dir = downloads_dir
        self.cache_dir = cache_dir
        self.midi_path = os.path.join(self.downloads_dir, "disease_species.csv")
        self.mime_path = os.path.join(self.downloads_dir, "micro_metabolic.csv")
        self.mege_path = os.path.join(self.downloads_dir, "meta_gene_net.csv")
        self.csv_reader = CSVReader()
        self.cache_manager = CacheManager(cache_dir=self.cache_dir)
        print("Data Cache Pipeline initialized...")

    def _cache_midi_taxon_info(self):
        taxids = [
            line[5] for line in self.csv_reader.line_generator_for_microbe_disease(self.midi_path)
        ]
        self.cache_manager.cache_entity("taxon_info", taxids=taxids)

    def _cache_mime_taxon_info(self):
        taxids = [
            line[9]
            if (line[9] and line[9] != "not available")
            else line[16]
            if (line[16] and line[16] != "not available")
            else None
            for line in self.csv_reader.line_generator(self.mime_path)
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
            for line in self.csv_reader.line_generator(self.mime_path)
            if line[6] and line[6] != "not available"
        ]
        self.cache_manager.cache_entity("pubchem_description", cids=pubchem_cids)

    def _cache_mege_pubchem_descriptions(self):
        pubchem_cids = [
            line[3]
            for line in self.csv_reader.line_generator(self.mege_path)
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
            for line in self.csv_reader.line_generator(self.mege_path)
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
            for line in self.csv_reader.line_generator(self.mege_path)
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
            for line in self.csv_reader.line_generator(self.mime_path)
            if line[6] and line[6] != "not available"
        ]
        self.cache_manager.cache_entity("pubchem_mw", cids=pubchem_cids)

    def _cache_mege_chem_properties(self):
        pubchem_cids = [
            line[3]
            for line in self.csv_reader.line_generator(self.mege_path)
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
        bigg_mapper = BiGGMapper()
        bigg_map = bigg_mapper.get_bigg_metabolite_mapping(
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
