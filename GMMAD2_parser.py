import itertools
import logging
import os
import uuid
from typing import Iterator

from utils.cache_manager import CacheHelper
from utils.cache_pipeline import DataCachePipeline
from utils.parser_helper import ParserHelper
from utils.reader import CSVReader
from utils.record_manager import RecordCacheManager


class GMMAD2Parser(CacheHelper):
    """
    Parses and merges microbe-disease, microbe-metabolite,
    metabolite-gene associations.
    """

    def __init__(
        self,
        csv_reader: CSVReader,
        cache_pipeline: DataCachePipeline,  # require CacheManager and CacheHelper
        parser_helpers: ParserHelper,
        cache_dir="cache",
        downloads_dir="downloads",
    ):
        super().__init__(cache_dir)  # inherit from CacheHelper
        self.csv_reader = csv_reader
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
            "description": f"{line[18]}[MESH]",
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
            else f"uuid:{str(uuid.uuid4())}",  # use uuid if no primary id but has a compound name
            "name": name.lower(),
            "synonyms": pubchem_desc.get(cid, {}).get("synonyms", []),
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
            "predicate": "biolink:OrganismalEntityAsAModelOfDiseaseAssociation",
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
            "infores:metabolomics": "ECO:0001230",  # mass spec + manual
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
        for line in self.csv_reader.line_generator_for_microbe_disease(self.midi_path):
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
        for line in self.csv_reader.line_generator(self.mime_path):
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
        for line in self.csv_reader.line_generator(self.mege_path):
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


class DataLoader:
    def __init__(self, cache_dir="cache", downloads_dir="downloads"):
        csv_reader = CSVReader()
        cache_pipeline = DataCachePipeline()
        parser_helpers = ParserHelper()

        self.parser = GMMAD2Parser(
            csv_reader=csv_reader,
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
