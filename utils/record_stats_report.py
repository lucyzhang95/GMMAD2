import json
import os
import random
from collections import Counter, defaultdict
from datetime import datetime
from statistics import mean, median
from typing import Any, Dict, List

from utils.cache_manager import CacheHelper


class RecordStatsReporter(CacheHelper):
    """Generates comprehensive statistics for GMMAD2 parsed records."""

    def __init__(self, cache_dir="../cache", report_dir="../reports"):
        super().__init__(cache_dir)
        self.report_dir = report_dir
        os.makedirs(self.report_dir, exist_ok=True)
        print("RecordStatsReporter initialized.")

    def _safe_get_numeric_stats(self, values: List[Any]):
        """Safely calculate numeric statistics for a list of values."""
        numeric_values = []
        for v in values:
            if v is not None:
                try:
                    if isinstance(v, dict):
                        if "average" in v:
                            numeric_values.append(float(v["average"]))
                        elif "monoisotopic" in v:
                            numeric_values.append(float(v["monoisotopic"]))
                    else:
                        numeric_values.append(float(v))
                except (ValueError, TypeError):
                    continue

        if not numeric_values:
            return {"count": 0, "min": None, "max": None, "mean": None, "median": None}

        return {
            "count": len(numeric_values),
            "min": min(numeric_values),
            "max": max(numeric_values),
            "mean": round(mean(numeric_values), 3),
            "median": round(median(numeric_values), 3),
        }

    def _extract_curie_prefix(self, curie: str) -> str:
        """Extract prefix from CURIE format ID."""
        if ":" in curie:
            return curie.split(":", 1)[0]
        return "unknown"

    def _count_xrefs(self, node: Dict[str, Any]) -> Dict[str, int]:
        """Count xref types in a node."""
        xrefs = node.get("xrefs", {})
        if not xrefs:
            return {}

        xref_counts = {}
        for key, value in xrefs.items():
            if value:  # only count non-empty xrefs
                if isinstance(value, list):
                    xref_counts[key] = len(value)
                else:
                    xref_counts[key] = 1
        return xref_counts

    def _collect_unique_xrefs(self, node: Dict[str, Any]) -> Dict[str, set]:
        """Collect unique xref values by type from a node."""
        xrefs = node.get("xrefs", {})
        if not xrefs:
            return {}

        unique_xrefs = {}
        for key, value in xrefs.items():
            if value:  # only non-empty xrefs
                if isinstance(value, list):
                    unique_xrefs[key] = set(value)
                else:
                    unique_xrefs[key] = {value}
        return unique_xrefs

    def _analyze_molecular_weight(self, nodes: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze molecular weight data from nodes."""
        mw_data = {"average": [], "monoisotopic": []}

        for node in nodes:
            mw = node.get("molecular_weight", {})
            if isinstance(mw, dict):
                if "average_molecular_weight" in mw and mw["average_molecular_weight"] is not None:
                    mw_data["average"].append(mw["average_molecular_weight"])
                if (
                    "monoisotopic_molecular_weight" in mw
                    and mw["monoisotopic_molecular_weight"] is not None
                ):
                    mw_data["monoisotopic"].append(mw["monoisotopic_molecular_weight"])

        return {
            "average": self._safe_get_numeric_stats(mw_data["average"]),
            "monoisotopic": self._safe_get_numeric_stats(mw_data["monoisotopic"]),
        }

    def generate_record_stats(self):
        """Generate comprehensive statistics for GMMAD2 records."""
        print("\n⚙️ Generating GMMAD2 record statistics...")

        combined_data = self.load_pickle("gmmad2_parsed_records.pkl")

        if not combined_data:
            print(
                "❌ No data found. Please ensure gmmad2_parsed_records.pkl exists in cache directory."
            )
            return {}

        stats = {
            "data_source": "GMMAD2",
            "total_relationships": 3,
            "relationship_types": ["microbe-disease", "microbe-metabolite", "metabolite-gene"],
            "metadata": {
                "report_version": "1.0",
                "analysis_date": datetime.now().isoformat(),
            },
        }

        # total record counts
        total_records = 0
        relationship_counts = {}
        for rel_type, records in combined_data.items():
            count = len(records)
            relationship_counts[rel_type] = count
            total_records += count

        stats["metadata"]["total_records_analyzed"] = total_records
        stats["metadata"]["total_records_analyzed_by_relationship"] = relationship_counts

        all_records = []
        all_evidence_types = []
        all_subject_curies = []
        all_object_curies = []
        all_subject_xrefs = defaultdict(int)
        all_object_xrefs = defaultdict(int)

        # for each relationship type
        relationship_stats = {}

        for rel_type, records in combined_data.items():
            rel_stats = self._analyze_relationship(rel_type, records)
            relationship_stats[rel_type] = rel_stats

            all_records.extend(records)
            all_evidence_types.extend(rel_stats.get("evidence_types", {}).keys())
            all_subject_curies.extend(rel_stats.get("subject_curie_stats", {}).keys())
            all_object_curies.extend(rel_stats.get("object_curie_stats", {}).keys())

            # collect xrefs
            for xref_type, count in rel_stats.get("subject_xref_stats", {}).items():
                all_subject_xrefs[xref_type] += count
            for xref_type, count in rel_stats.get("object_xref_stats", {}).items():
                all_object_xrefs[xref_type] += count

        stats["relationship_analysis"] = relationship_stats

        # evidence types
        evidence_counter = Counter()
        for records in combined_data.values():
            for record in records:
                evidence_type = record.get("association", {}).get("evidence_type")
                if evidence_type:
                    if isinstance(evidence_type, list):
                        evidence_counter.update(evidence_type)
                    else:
                        evidence_counter[evidence_type] += 1

        stats["overall_evidence_types"] = dict(evidence_counter)

        # _id duplication check
        all_record_ids = []
        for records in combined_data.values():
            for record in records:
                all_record_ids.append(record.get("_id"))

        id_counter = Counter(all_record_ids)
        duplicates = {id_val: count for id_val, count in id_counter.items() if count > 1}

        count_groups = defaultdict(list)
        for id_val, count in id_counter.items():
            if count > 1:
                count_groups[count].append(id_val)

        # randomly select 3 IDs from each count group
        sampled_ids_by_count = {}
        for count, id_list in count_groups.items():
            sample_size = min(3, len(id_list))
            sampled_ids_by_count[count] = random.sample(id_list, sample_size)

        stats["overall_id_analysis"] = {
            "total_ids": len(all_record_ids),
            "unique_ids": len(id_counter),
            "duplicate_count": len(duplicates),
            "duplicates": sampled_ids_by_count,
            "duplicates_distribution": dict(Counter(k for k in id_counter.values() if k > 1)),
        }

        # CURIE analysis
        overall_subject_curies = []
        overall_object_curies = []
        for records in combined_data.values():
            for record in records:
                subject_id = record.get("subject", {}).get("id", "")
                object_id = record.get("object", {}).get("id", "")
                if subject_id:
                    overall_subject_curies.append(self._extract_curie_prefix(subject_id))
                if object_id:
                    overall_object_curies.append(self._extract_curie_prefix(object_id))

        stats["overall_subject_curie_stats"] = dict(Counter(overall_subject_curies))
        stats["overall_unique_subject_curie_count"] = len(set(overall_subject_curies))
        stats["overall_object_curie_stats"] = dict(Counter(overall_object_curies))
        stats["overall_unique_object_curie_count"] = len(set(overall_object_curies))

        # description analysis
        subject_desc_count = 0
        object_desc_count = 0
        for records in combined_data.values():
            for record in records:
                if record.get("subject", {}).get("description"):
                    subject_desc_count += 1
                if record.get("object", {}).get("description"):
                    object_desc_count += 1

        stats["overall_description_stats"] = {
            "subject_descriptions": subject_desc_count,
            "subject_description_percentage": round((subject_desc_count / total_records) * 100, 2),
            "object_descriptions": object_desc_count,
            "object_description_percentage": round((object_desc_count / total_records) * 100, 2),
        }

        overall_subject_xrefs = defaultdict(int)
        overall_object_xrefs = defaultdict(int)
        overall_subject_unique_xrefs = defaultdict(set)
        overall_object_unique_xrefs = defaultdict(set)

        for records in combined_data.values():
            for record in records:
                subject_xrefs = self._count_xrefs(record.get("subject", {}))
                object_xrefs = self._count_xrefs(record.get("object", {}))
                subject_unique_xrefs = self._collect_unique_xrefs(record.get("subject", {}))
                object_unique_xrefs = self._collect_unique_xrefs(record.get("object", {}))

                for xref_type, count in subject_xrefs.items():
                    overall_subject_xrefs[xref_type] += count
                for xref_type, count in object_xrefs.items():
                    overall_object_xrefs[xref_type] += count

                for xref_type, unique_values in subject_unique_xrefs.items():
                    overall_subject_unique_xrefs[xref_type].update(unique_values)
                for xref_type, unique_values in object_unique_xrefs.items():
                    overall_object_unique_xrefs[xref_type].update(unique_values)

        # Convert sets to counts for JSON serialization
        subject_unique_counts = {k: len(v) for k, v in overall_subject_unique_xrefs.items()}
        object_unique_counts = {k: len(v) for k, v in overall_object_unique_xrefs.items()}

        stats["overall_xref_stats"] = {
            "subject_xrefs": dict(overall_subject_xrefs),
            "object_xrefs": dict(overall_object_xrefs),
            "subject_unique_xrefs": subject_unique_counts,
            "object_unique_xrefs": object_unique_counts,
        }

        # publication analysis
        overall_publication_count = 0
        overall_pmid_count = 0
        for records in combined_data.values():
            for record in records:
                pub = record.get("association", {}).get("publications", {})
                pmid = record.get("association", {}).get("publications", {}).get("pmid")
                if pub:
                    overall_publication_count += 1
                if pmid:
                    overall_pmid_count += 1

        stats["overall_publication_stats"] = {
            "records_with_publication": overall_publication_count,
            "records_with_pmid": overall_pmid_count,
        }

        print(f"Generated comprehensive statistics for {total_records} total records.")
        return stats

    def _analyze_relationship(self, rel_type: str, records: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze statistics for a specific relationship type."""
        print(f"-> Analyzing {rel_type} relationship ({len(records)} records)...")

        rel_stats = {
            "record_count": len(records),
            "evidence_types": {},
            "subject_curie_stats": {},
            "object_curie_stats": {},
            "subject_description_count": 0,
            "object_description_count": 0,
            "subject_xref_stats": defaultdict(int),
            "object_xref_stats": defaultdict(int),
            "subject_unique_xref_stats": defaultdict(set),
            "object_unique_xref_stats": defaultdict(set),
            "publication_stats": {"records_with_pmid": 0},
        }

        # basic stats
        evidence_counter = Counter()
        subject_curie_counter = Counter()
        object_curie_counter = Counter()

        for record in records:
            # evidence types
            evidence_type = record.get("association", {}).get("evidence_type")
            if evidence_type:
                if isinstance(evidence_type, list):
                    evidence_counter.update(evidence_type)
                else:
                    evidence_counter[evidence_type] += 1

            # subject and object CURIE prefixes
            subject_id = record.get("subject", {}).get("id", "")
            object_id = record.get("object", {}).get("id", "")

            if subject_id:
                subject_curie_counter[self._extract_curie_prefix(subject_id)] += 1
            if object_id:
                object_curie_counter[self._extract_curie_prefix(object_id)] += 1

            # descriptions
            if record.get("subject", {}).get("description"):
                rel_stats["subject_description_count"] += 1
            if record.get("object", {}).get("description"):
                rel_stats["object_description_count"] += 1

            subject_xrefs = self._count_xrefs(record.get("subject", {}))
            object_xrefs = self._count_xrefs(record.get("object", {}))
            subject_unique_xrefs = self._collect_unique_xrefs(record.get("subject", {}))
            object_unique_xrefs = self._collect_unique_xrefs(record.get("object", {}))

            for xref_type, count in subject_xrefs.items():
                rel_stats["subject_xref_stats"][xref_type] += count
            for xref_type, count in object_xrefs.items():
                rel_stats["object_xref_stats"][xref_type] += count

            for xref_type, unique_values in subject_unique_xrefs.items():
                rel_stats["subject_unique_xref_stats"][xref_type].update(unique_values)
            for xref_type, unique_values in object_unique_xrefs.items():
                rel_stats["object_unique_xref_stats"][xref_type].update(unique_values)

            # publications
            pmid = record.get("association", {}).get("publications", {}).get("pmid")
            if pmid:
                rel_stats["publication_stats"]["records_with_pmid"] += 1

        rel_stats["evidence_types"] = dict(evidence_counter)
        rel_stats["subject_curie_stats"] = dict(subject_curie_counter)
        rel_stats["object_curie_stats"] = dict(object_curie_counter)
        rel_stats["subject_xref_stats"] = dict(rel_stats["subject_xref_stats"])
        rel_stats["object_xref_stats"] = dict(rel_stats["object_xref_stats"])

        rel_stats["subject_unique_xref_stats"] = {
            k: len(v) for k, v in rel_stats["subject_unique_xref_stats"].items()
        }
        rel_stats["object_unique_xref_stats"] = {
            k: len(v) for k, v in rel_stats["object_unique_xref_stats"].items()
        }

        # relationship-specific analysis
        if rel_type == "microbe-disease":
            rel_stats.update(self._analyze_midi_specific(records))
        elif rel_type == "microbe-metabolite":
            rel_stats.update(self._analyze_mime_specific(records))
        elif rel_type == "metabolite-gene":
            rel_stats.update(self._analyze_mege_specific(records))

        return rel_stats

    def _analyze_midi_specific(self, records: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze MIDI-specific statistics."""
        organism_types = []
        ranks = []
        qualifiers = []
        disease_sample_sizes = []
        healthy_sample_sizes = []
        disease_abundances = []
        healthy_abundances = []

        for record in records:
            subject = record.get("subject", {})
            association = record.get("association", {})

            # subject organism type and rank
            organism_type = subject.get("organism_type")
            if organism_type:
                organism_types.append(organism_type)

            rank = subject.get("rank")
            if rank:
                ranks.append(rank)

            # association data
            qualifier = association.get("qualifier")
            if qualifier:
                qualifiers.append(qualifier)

            disease_size = association.get("disease_sample_size")
            if disease_size:
                disease_sample_sizes.append(disease_size)

            healthy_size = association.get("healthy_sample_size")
            if healthy_size:
                healthy_sample_sizes.append(healthy_size)

            disease_mean = association.get("disease_abundance_mean")
            if disease_mean:
                disease_abundances.append(disease_mean)

            healthy_mean = association.get("healthy_abundance_mean")
            if healthy_mean:
                healthy_abundances.append(healthy_mean)

        return {
            "midi_organism_types": dict(Counter(organism_types)),
            "midi_ranks": dict(Counter(ranks)),
            "midi_qualifiers": dict(Counter(qualifiers)),
            "midi_disease_sample_size_stats": self._safe_get_numeric_stats(disease_sample_sizes),
            "midi_healthy_sample_size_stats": self._safe_get_numeric_stats(healthy_sample_sizes),
            "midi_disease_abundance_stats": self._safe_get_numeric_stats(disease_abundances),
            "midi_healthy_abundance_stats": self._safe_get_numeric_stats(healthy_abundances),
        }

    def _analyze_mime_specific(self, records: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze MIME-specific statistics."""
        habitats = []
        object_nodes = []
        xlogp_values = []

        for record in records:
            association = record.get("association", {})
            object_node = record.get("object", {})
            object_nodes.append(object_node)

            # association habitat
            habitat = association.get("association_habitat", [])
            if habitat:
                habitats.extend(habitat)

            # xLogP
            xlogp = object_node.get("xlogp")
            if xlogp is not None:
                xlogp_values.append(xlogp)

        return {
            "mime_association_habitats": dict(Counter(habitats)),
            "mime_object_molecular_weight_stats": self._analyze_molecular_weight(object_nodes),
            "mime_object_xlogp_stats": self._safe_get_numeric_stats(xlogp_values),
        }

    def _analyze_mege_specific(self, records: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze MEGE-specific statistics."""
        habitats = []
        subject_nodes = []
        xlogp_values = []
        publication_data = {"records_with_publications": 0}

        for record in records:
            association = record.get("association", {})
            subject_node = record.get("subject", {})
            subject_nodes.append(subject_node)

            # association habitat
            habitat = association.get("association_habitat", [])
            if habitat:
                habitats.extend(habitat)

            # xLogP
            xlogp = subject_node.get("xlogp")
            if xlogp is not None:
                xlogp_values.append(xlogp)

            # publication
            publications = association.get("publications", {})
            if publications and publications.get("pmid"):
                publication_data["records_with_publications"] += 1

        return {
            "mege_association_habitats": dict(Counter(habitats)),
            "mege_subject_molecular_weight_stats": self._analyze_molecular_weight(subject_nodes),
            "mege_subject_xlogp_stats": self._safe_get_numeric_stats(xlogp_values),
            "mege_publication_stats": publication_data,
        }

    def save_stats_report(self, stats) -> str:
        """Save the statistics report to JSON file."""
        report_path = os.path.join(self.report_dir, f"record_stats_{datetime.now().isoformat()}.json")
        with open(report_path, "w") as f:
            json.dump(stats, f, indent=2, sort_keys=True)

        print(f"💾 Statistics report saved to: {report_path}")
        return report_path

    def run_full_analysis(self):
        """Run complete statistical analysis and save report."""
        print("▶️ Starting comprehensive GMMAD2 record analysis...")

        stats = self.generate_record_stats()
        if stats:
            self.save_stats_report(stats)
            print("🎉 Analysis complete!")
        else:
            print("❌ Analysis failed - no data found.")

        return stats

    def export_all_duplicated_records(self) -> str:
        """Export all duplicated records to a JSON file."""
        print("\n▶️ Exporting all duplicated records...")

        combined_data = self.load_pickle("gmmad2_parsed_records.pkl")
        if not combined_data:
            print(
                "❌ No data found. Please ensure gmmad2_parsed_records.pkl exists in cache directory."
            )
            return ""

        all_records_with_type = []
        for rel_type, records in combined_data.items():
            for record in records:
                record_with_type = record.copy()
                record_with_type["_relationship_type"] = rel_type
                all_records_with_type.append(record_with_type)

        id_to_records = defaultdict(list)
        for record in all_records_with_type:
            record_id = record.get("_id")
            if record_id:
                id_to_records[record_id].append(record)

        # filter only duplicates (count > 1)
        duplicated_records = {}
        total_duplicate_records = 0

        for record_id, records in id_to_records.items():
            if len(records) > 1:
                duplicated_records[record_id] = {
                    "duplicate_count": len(records),
                    "records": records,
                }
                total_duplicate_records += len(records)

        export_data = {
            "metadata": {
                "export_date": datetime.now().isoformat(),
                "total_duplicate_ids": len(duplicated_records),
                "total_duplicate_records": total_duplicate_records,
                "description": "All records with duplicate _id values",
            },
            "duplicated_records": duplicated_records,
        }

        report_path = os.path.join(self.report_dir, f"all_duplicated_records_{datetime.now().isoformat()}.json")
        with open(report_path, "w") as f:
            json.dump(export_data, f, indent=2, sort_keys=True)

        print(f"💾 All duplicated records exported to: {report_path}")
        print(
            f"✅ Found {len(duplicated_records)} duplicate IDs with {total_duplicate_records} total records"
        )

        return report_path

    def export_sampled_duplicated_records(self) -> str:
        """Export randomly sampled duplicated records (3 examples per duplication count)."""
        print("\n▶️ Exporting sampled duplicated records...")

        combined_data = self.load_pickle("gmmad2_parsed_records.pkl")
        if not combined_data:
            print(
                "❌ No data found. Please ensure gmmad2_parsed_records.pkl exists in cache directory."
            )
            return ""

        all_records_with_type = []
        for rel_type, records in combined_data.items():
            for record in records:
                record_with_type = record.copy()
                record_with_type["_relationship_type"] = rel_type
                all_records_with_type.append(record_with_type)

        # duplicated records and group by count
        id_to_records = defaultdict(list)
        for record in all_records_with_type:
            record_id = record.get("_id")
            if record_id:
                id_to_records[record_id].append(record)

        # group duplicates by count
        count_groups = defaultdict(list)
        for record_id, records in id_to_records.items():
            if len(records) > 1:  # only duplicates
                count_groups[len(records)].append({"id": record_id, "records": records})

        # 3 examples from each count group
        sampled_duplicates = {}
        total_sampled_records = 0

        for count, duplicate_list in count_groups.items():
            sample_size = min(3, len(duplicate_list))
            sampled_items = random.sample(duplicate_list, sample_size)

            sampled_duplicates[f"count_{count}"] = {
                "duplicate_count": count,
                "total_ids_with_this_count": len(duplicate_list),
                "sampled_examples": {},
            }

            for item in sampled_items:
                sampled_duplicates[f"count_{count}"]["sampled_examples"][item["id"]] = {
                    "duplicate_count": count,
                    "records": item["records"],
                }
                total_sampled_records += len(item["records"])

        export_data = {
            "metadata": {
                "export_date": datetime.now().isoformat(),
                "sampling_strategy": "3 random examples per duplication count",
                "total_duplicate_count_groups": len(count_groups),
                "total_sampled_ids": sum(
                    len(group["sampled_examples"]) for group in sampled_duplicates.values()
                ),
                "total_sampled_records": total_sampled_records,
                "description": "Randomly sampled duplicate records (3 examples per duplication count)",
            },
            "count_distribution": {
                f"count_{count}": len(duplicate_list)
                for count, duplicate_list in count_groups.items()
            },
            "sampled_duplicates": sampled_duplicates,
        }

        report_path = os.path.join(self.report_dir, f"sampled_duplicated_records_{datetime.now().isoformat()}.json")
        with open(report_path, "w") as f:
            json.dump(export_data, f, indent=2, sort_keys=True)

        print(f"💾 Sampled duplicated records exported to: {report_path}")
        print(
            f"📊 Sampled {len(sampled_duplicates)} count groups with {total_sampled_records} total records"
        )

        for _, data in sampled_duplicates.items():
            count = data["duplicate_count"]
            total_with_count = data["total_ids_with_this_count"]
            sampled_count = len(data["sampled_examples"])
            print(f"-> {count}x duplicates: {sampled_count}/{total_with_count} IDs sampled")

        return report_path

    def export_both_duplicate_reports(self) -> Dict[str, str]:
        """Export both complete and sampled duplicate reports."""
        print("\n📋 Exporting duplicate record reports...")

        paths = {
            "all_duplicates": self.export_all_duplicated_records(),
            "sampled_duplicates": self.export_sampled_duplicated_records(),
        }

        print("🎉 Duplicate export complete!")
        return paths


if __name__ == "__main__":
    reporter = RecordStatsReporter()
    stats_report = reporter.run_full_analysis()

    paths = reporter.export_both_duplicate_reports()
    print(f"All duplicates: {paths['all_duplicates']}")
    print(f"Sampled duplicates: {paths['sampled_duplicates']}")
