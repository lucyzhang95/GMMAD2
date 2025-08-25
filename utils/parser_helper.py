from typing import List


class ParserHelper:
    @staticmethod
    def filter_taxon_info(taxon_info: list) -> tuple[dict, list]:
        """
        Filters the taxon_info dictionary to remove not found entries.

        :param taxon_info: Output of query_taxon_info_from_biothings.
        :return:
        """
        taxon_info_filtered = {t["query"]: t for t in taxon_info if "notfound" not in t.keys()}
        taxid_notfound = set(t["query"] for t in taxon_info if "notfound" in t.keys())
        return taxon_info_filtered, list(taxid_notfound)

    @staticmethod
    def fetch_taxon_names_from_taxon_info(taxon_info: dict) -> list[str]:
        """Extracts biothings names from the taxon_info dictionary."""
        taxon_names = set()
        for _, taxon in taxon_info.items():
            if "scientific_name" in taxon:
                taxon_names.add(taxon["scientific_name"].lower())
        return list(taxon_names)

    @staticmethod
    def get_organism_type(node) -> str:
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

    @staticmethod
    def get_suffix(identifier: str) -> str:
        return identifier.split(":", 1)[1].strip() if ":" in identifier else identifier.strip()

    @staticmethod
    def assign_primary_metabolite_id(id_hierarchy: List[tuple]) -> tuple[str, dict] or None:
        def classify_kegg(val):
            if val.startswith("C"):
                return "KEGG.COMPOUND"
            if val.startswith("G"):
                return "KEGG.GLYCAN"
            if val.startswith("D"):
                return "KEGG.DRUG"
            return "KEGG"

        xrefs = {}
        primary_id = None
        for val, prefix in id_hierarchy:
            if not val or str(val).strip().lower() == "not available":
                continue
            if prefix == "KEGG" and isinstance(val, str):
                prefix = classify_kegg(val)

            curie = f"{prefix}:{val}" if prefix else str(val)
            key_map = {
                "PUBCHEM.COMPOUND": "pubchem_cid",
                "DRUGBANK": "drugbank",
                "HMDB": "hmdb",
                "SMILES": "smiles",
                "BIGG.METABOLITE": "bigg",
                "KEGG.COMPOUND": "kegg",
                "KEGG.GLYCAN": "kegg",
                "KEGG.DRUG": "kegg",
            }
            key = key_map.get(prefix)
            if not key:
                continue

            if primary_id is None:
                primary_id = curie
            xrefs.setdefault(key, curie)

        return primary_id, xrefs

    @staticmethod
    def assign_primary_gene_id(line):
        ensembl_id = line[13]
        ncbi_id = line[14]
        hgnc_id = line[15]
        uniprotkb = line[16]

        # (line, prefix) pairs for ID hierarchy
        id_hierarchy = [
            (uniprotkb, "UniProtKB"),
            (ensembl_id, "ENSEMBL"),
            (ncbi_id, "NCBIGene"),
            (hgnc_id, "HGNC"),
        ]

        xrefs = {}
        primary_id = None

        for val, prefix in id_hierarchy:
            if not val or val.strip().lower() == "not available":
                continue

            curie = f"{prefix}:{val}" if prefix else val
            if prefix == "NCBIGene":
                key = "entrezgene"
            else:
                key = prefix.lower()

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
