import glob
import json
import os
import random
import sys
from collections import Counter, defaultdict
from datetime import datetime
from statistics import mean, median
from typing import Any, Dict, Iterator, List, Optional, Set


class GMMAD2RecordStatsReporter:
    """Generates comprehensive statistics for GMMAD2 parsed records using streaming approach."""

    def __init__(self, report_dir="../reports"):
        self.report_dir = report_dir
        os.makedirs(self.report_dir, exist_ok=True)
        print("GMMAD2RecordStatsReporter initialized.")

        self.all_subject_properties = set()
        self.all_object_properties = set()
        self.all_association_properties = set()

        self.overall_counters = {
            "_id": Counter(),
            "subject_id": Counter(),
            "object_id": Counter(),
            "subject_description": Counter(),
            "object_description": Counter(),
            "subject_xrefs": Counter(),
            "object_xrefs": Counter(),
            "association_publication": Counter(),
            "relationship_types": Counter(),
        }

        self.duplicate_records = defaultdict(list)

        self.relationship_property_counters = defaultdict(
            lambda: {
                "subject_properties": defaultdict(Counter),
                "object_properties": defaultdict(Counter),
                "association_properties": defaultdict(Counter),
                "record_count": 0,
            }
        )

    def _find_latest_jsonl_file(self, base_path: str = "gmmad2_parsed_records") -> Optional[str]:
        """Find the most recent timestamped JSONL file."""
        pattern = os.path.join("..", "records", f"{base_path}_*.jsonl")
        files = glob.glob(pattern)

        if not files:
            print(f"No JSONL files found matching pattern: {pattern}")
            return None

        latest_file = max(files, key=os.path.getmtime)
        print(f"Found latest JSONL file: {latest_file}")
        return latest_file

    def _load_jsonl_records(self, file_path: str) -> Iterator[Dict[str, Any]]:
        """Generator to load records one by one from JSONL file."""
        print(f"Streaming records from: {file_path}")

        record_count = 0
        with open(file_path, "r", encoding="utf-8") as f:
            for line_num, line in enumerate(f, 1):
                if line.strip():
                    try:
                        record = json.loads(line)
                        record_count += 1
                        if record_count % 10000 == 0:
                            print(f"-> Processed {record_count} records...")
                        yield record
                    except json.JSONDecodeError as e:
                        print(f"!!! Error parsing line {line_num}: {e}")
                        continue

        print(f"Total records processed: {record_count}")

    def _extract_all_properties(self, data: Dict[str, Any], prefix: str = "") -> Set[str]:
        """Recursively extract all property paths from a dictionary."""
        properties = set()

        if not isinstance(data, dict):
            return properties

        for key, value in data.items():
            current_path = f"{prefix}.{key}" if prefix else key
            properties.add(current_path)

            if isinstance(value, dict):
                nested_props = self._extract_all_properties(value, current_path)
                properties.update(nested_props)
            elif isinstance(value, list) and value and isinstance(value[0], dict):
                nested_props = self._extract_all_properties(value[0], current_path)
                properties.update(nested_props)

        return properties

    def _get_nested_value(self, data: Dict[str, Any], path: str) -> Any:
        """Get value from nested dictionary using dot notation path."""
        keys = path.split(".")
        current = data

        for key in keys:
            if isinstance(current, dict) and key in current:
                current = current[key]
            else:
                return None

        return current

    def _extract_curie_prefix(self, curie: str) -> str:
        """Extract prefix from CURIE format ID."""
        if isinstance(curie, str) and ":" in curie:
            return curie.split(":", 1)[0]
        return "unknown"

    def _safe_get_numeric_stats(self, values: List[Any]) -> Dict[str, Any]:
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

    def _process_single_record(self, record: Dict[str, Any]) -> None:
        """Process a single record and update streaming counters."""
        record_id = record.get("_id")
        if record_id:
            self.overall_counters["_id"][record_id] += 1
            if self.overall_counters["_id"][record_id] > 1:
                self.duplicate_records[record_id].append(record)

        # subject.id and object.id
        subject_id = record.get("subject", {}).get("id")
        if subject_id:
            self.overall_counters["subject_id"][subject_id] += 1

        object_id = record.get("object", {}).get("id")
        if object_id:
            self.overall_counters["object_id"][object_id] += 1

        # descriptions
        subject_desc = record.get("subject", {}).get("description")
        if subject_desc:
            self.overall_counters["subject_description"][subject_desc] += 1

        object_desc = record.get("object", {}).get("description")
        if object_desc:
            self.overall_counters["object_description"][object_desc] += 1

        # xrefs
        subject_xrefs = record.get("subject", {}).get("xrefs", [])
        if subject_xrefs and isinstance(subject_xrefs, list):
            for xref in subject_xrefs:
                if xref:
                    self.overall_counters["subject_xrefs"][str(xref)] += 1

        object_xrefs = record.get("object", {}).get("xrefs", [])
        if object_xrefs and isinstance(object_xrefs, list):
            for xref in object_xrefs:
                if xref:
                    self.overall_counters["object_xrefs"][str(xref)] += 1

        # association.publication
        pub = record.get("association", {}).get("publication")
        if pub:
            self.overall_counters["association_publication"][str(pub)] += 1

        # relationship type
        rel_type = record.get("association", {}).get("category")
        if rel_type:
            self.overall_counters["relationship_types"][rel_type] += 1
            self._process_relationship_record(record, rel_type)

    def _process_relationship_record(self, record: Dict[str, Any], rel_type: str) -> None:
        """Process record for relationship-specific statistics."""
        rel_stats = self.relationship_property_counters[rel_type]
        rel_stats["record_count"] += 1

        # subject properties
        subject = record.get("subject", {})
        if subject:
            subject_props = self._extract_all_properties(subject)
            self.all_subject_properties.update(subject_props)

            for prop in subject_props:
                value = self._get_nested_value(subject, prop)
                if value is not None:
                    if isinstance(value, list):
                        for v in value:
                            if v is not None:
                                rel_stats["subject_properties"][prop][str(v)] += 1
                    else:
                        rel_stats["subject_properties"][prop][str(value)] += 1

        # object properties
        obj = record.get("object", {})
        if obj:
            object_props = self._extract_all_properties(obj)
            self.all_object_properties.update(object_props)

            for prop in object_props:
                value = self._get_nested_value(obj, prop)
                if value is not None:
                    if isinstance(value, list):
                        for v in value:
                            if v is not None:
                                rel_stats["object_properties"][prop][str(v)] += 1
                    else:
                        rel_stats["object_properties"][prop][str(value)] += 1

        # association properties
        assoc = record.get("association", {})
        if assoc:
            assoc_props = self._extract_all_properties(assoc)
            self.all_association_properties.update(assoc_props)

            for prop in assoc_props:
                value = self._get_nested_value(assoc, prop)
                if value is not None:
                    if isinstance(value, list):
                        for v in value:
                            if v is not None:
                                rel_stats["association_properties"][prop][str(v)] += 1
                    else:
                        rel_stats["association_properties"][prop][str(value)] += 1

    def _collect_all_duplicates_in_second_pass(self, file_path: str) -> None:
        """Second pass to collect all duplicate records (including first occurrences)."""
        print("Second pass: collecting all duplicate records...")

        duplicate_ids = set(
            id_val for id_val, count in self.overall_counters["_id"].items() if count > 1
        )

        if not duplicate_ids:
            print("No duplicates found, skipping second pass")
            return

        for record in self._load_jsonl_records(file_path):
            record_id = record.get("_id")
            if record_id in duplicate_ids:
                if len(self.duplicate_records[record_id]) < self.overall_counters["_id"][record_id]:
                    self.duplicate_records[record_id].append(record)

    def _analyze_duplicates(self) -> Dict[str, Any]:
        """Analyze duplicate records and export samples."""
        print("Analyzing duplicate records...")

        duplicate_count_distribution = Counter()
        total_duplicate_records = 0

        for _record_id, count in self.overall_counters["_id"].items():
            if count > 1:
                duplicate_count_distribution[count] += 1
                total_duplicate_records += count

        sampled_duplicates = {}
        for count, _num_groups in duplicate_count_distribution.items():
            ids_with_count = [
                rid for rid, rcount in self.overall_counters["_id"].items() if rcount == count
            ]

            sample_size = min(3, len(ids_with_count))
            sampled_ids = random.sample(ids_with_count, sample_size)

            sampled_duplicates[count] = {}
            for record_id in sampled_ids:
                if record_id in self.duplicate_records:
                    sampled_duplicates[count][record_id] = self.duplicate_records[record_id]

        self._export_duplicate_records(self.duplicate_records, sampled_duplicates)

        duplicate_stats = {
            "total_duplicate_groups": len(duplicate_count_distribution),
            "total_duplicate_records": total_duplicate_records,
            "duplication_count_distribution": dict(duplicate_count_distribution),
            "max_duplicates_for_single_id": max(duplicate_count_distribution.keys())
            if duplicate_count_distribution
            else 0,
            "sampled_duplicates_by_count": {
                count: list(group_data.keys()) for count, group_data in sampled_duplicates.items()
            },
        }

        return duplicate_stats

    def _export_duplicate_records(
        self,
        all_duplicates: Dict[str, List[Dict]],
        sampled_duplicates: Dict[int, Dict[str, List[Dict]]],
    ) -> str:
        """Export duplicate records to JSON files."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        sampled_export_path = os.path.join(
            self.report_dir, f"gmmad2_sampled_duplicate_records_{timestamp}.json"
        )
        sampled_export_data = {
            "metadata": {
                "export_date": datetime.now().isoformat(),
                "description": "Randomly sampled GMMAD2 duplicate records (max 3 per duplication count)",
                "total_sampled_groups": sum(len(groups) for groups in sampled_duplicates.values()),
            },
            "duplicates_by_count": sampled_duplicates,
        }

        with open(sampled_export_path, "w", encoding="utf-8") as f:
            json.dump(sampled_export_data, f, indent=2, sort_keys=True, default=str)

        if len(all_duplicates) < 10000:
            all_export_path = os.path.join(
                self.report_dir, f"gmmad2_all_duplicate_records_{timestamp}.json"
            )
            all_export_data = {
                "metadata": {
                    "export_date": datetime.now().isoformat(),
                    "description": "All GMMAD2 duplicate records found",
                    "total_duplicate_groups": len(all_duplicates),
                    "total_duplicate_records": sum(len(rlist) for rlist in all_duplicates.values()),
                },
                "duplicate_records": all_duplicates,
            }

            with open(all_export_path, "w", encoding="utf-8") as f:
                json.dump(all_export_data, f, indent=2, sort_keys=True, default=str)

            print(f"All duplicate records exported to: {all_export_path}")
        else:
            print(f"Skipping full duplicate export (too many duplicates: {len(all_duplicates)})")

        print(f"Sampled duplicate records exported to: {sampled_export_path}")
        return sampled_export_path

    def _compile_overall_statistics(self) -> Dict[str, Any]:
        """Compile overall statistics from streaming counters."""
        total_records = sum(self.overall_counters["relationship_types"].values())

        stats = {
            "_id": {
                "raw_count": sum(self.overall_counters["_id"].values()),
                "unique_count": len(self.overall_counters["_id"]),
                "duplicate_count": sum(self.overall_counters["_id"].values())
                - len(self.overall_counters["_id"]),
                "duplicate_analysis": self._analyze_duplicates(),
            },
            "subject_id": {
                "raw_count": sum(self.overall_counters["subject_id"].values()),
                "unique_count": len(self.overall_counters["subject_id"]),
            },
            "object_id": {
                "raw_count": sum(self.overall_counters["object_id"].values()),
                "unique_count": len(self.overall_counters["object_id"]),
            },
            "subject_description": {
                "raw_count": sum(self.overall_counters["subject_description"].values()),
                "record_count": len(self.overall_counters["subject_description"]),
                "unique_count": len(self.overall_counters["subject_description"]),
                "percentage": round(
                    (len(self.overall_counters["subject_description"]) / total_records) * 100, 2
                )
                if total_records > 0
                else 0,
                "sample_values": list(self.overall_counters["subject_description"].keys())[:1],
            },
            "object_description": {
                "raw_count": sum(self.overall_counters["object_description"].values()),
                "record_count": len(self.overall_counters["object_description"]),
                "unique_count": len(self.overall_counters["object_description"]),
                "percentage": round(
                    (len(self.overall_counters["object_description"]) / total_records) * 100, 2
                )
                if total_records > 0
                else 0,
                "sample_values": list(self.overall_counters["object_description"].keys())[:1],
            },
            "subject_xrefs": {
                "raw_count": sum(self.overall_counters["subject_xrefs"].values()),
                "unique_count": len(self.overall_counters["subject_xrefs"]),
                "sample_values": list(self.overall_counters["subject_xrefs"].keys())[:5],
            },
            "object_xrefs": {
                "raw_count": sum(self.overall_counters["object_xrefs"].values()),
                "unique_count": len(self.overall_counters["object_xrefs"]),
                "sample_values": list(self.overall_counters["object_xrefs"].keys())[:5],
            },
            "association_publication": {
                "raw_count": sum(self.overall_counters["association_publication"].values()),
                "record_count": len(self.overall_counters["association_publication"]),
                "unique_count": len(self.overall_counters["association_publication"]),
                "percentage": round(
                    (len(self.overall_counters["association_publication"]) / total_records) * 100, 2
                )
                if total_records > 0
                else 0,
                "sample_values": list(self.overall_counters["association_publication"].keys())[:1],
            },
        }

        return stats

    def _compile_relationship_statistics(self) -> Dict[str, Any]:
        """Compile relationship-specific statistics from streaming counters."""
        relationship_stats = {}

        for rel_type, rel_data in self.relationship_property_counters.items():
            stats = {
                "record_count": rel_data["record_count"],
                "subject_properties": {},
                "object_properties": {},
                "association_properties": {},
            }

            # subject properties
            for prop, counter in rel_data["subject_properties"].items():
                stats["subject_properties"][prop] = {
                    "raw_count": sum(counter.values()),
                    "record_count": len(counter),
                    "unique_count": len(counter),
                    "percentage": round((len(counter) / rel_data["record_count"]) * 100, 2)
                    if rel_data["record_count"] > 0
                    else 0,
                    "sample_values": list(counter.keys())[:1],
                }

            # object properties
            for prop, counter in rel_data["object_properties"].items():
                stats["object_properties"][prop] = {
                    "raw_count": sum(counter.values()),
                    "record_count": len(counter),
                    "unique_count": len(counter),
                    "percentage": round((len(counter) / rel_data["record_count"]) * 100, 2)
                    if rel_data["record_count"] > 0
                    else 0,
                    "sample_values": list(counter.keys())[:1],
                }

            # association properties
            for prop, counter in rel_data["association_properties"].items():
                stats["association_properties"][prop] = {
                    "raw_count": sum(counter.values()),
                    "record_count": len(counter),
                    "unique_count": len(counter),
                    "percentage": round((len(counter) / rel_data["record_count"]) * 100, 2)
                    if rel_data["record_count"] > 0
                    else 0,
                    "sample_values": list(counter.keys())[:1],
                }

            relationship_stats[rel_type] = stats

        return relationship_stats

    def generate_comprehensive_stats(self, file_path: Optional[str] = None) -> Dict[str, Any]:
        """Generate comprehensive statistics using streaming approach."""
        if file_path is None:
            file_path = self._find_latest_jsonl_file()
            if file_path is None:
                return {}

        print(f"Starting streaming analysis of: {file_path}")

        # build counters
        print("First pass: analyzing all records...")
        for record in self._load_jsonl_records(file_path):
            self._process_single_record(record)

        # collect complete duplicate records
        if any(count > 1 for count in self.overall_counters["_id"].values()):
            self._collect_all_duplicates_in_second_pass(file_path)

        overall_stats = self._compile_overall_statistics()
        relationship_stats = self._compile_relationship_statistics()

        total_records = sum(self.overall_counters["relationship_types"].values())

        stats_report = {
            "metadata": {
                "data_source": "GMMAD2",
                "analysis_date": datetime.now().isoformat(),
                "source_file": file_path,
                "total_records": total_records,
                "relationship_types": list(self.overall_counters["relationship_types"].keys()),
                "relationship_type_counts": dict(self.overall_counters["relationship_types"]),
            },
            "overall_statistics": overall_stats,
            "relationship_type_analysis": relationship_stats,
            "property_summary": {
                "all_subject_properties": sorted(list(self.all_subject_properties)),
                "all_object_properties": sorted(list(self.all_object_properties)),
                "all_association_properties": sorted(list(self.all_association_properties)),
            },
        }

        return stats_report

    def save_stats_report(
        self, stats: Dict[str, Any], output_filename: Optional[str] = None
    ) -> str:
        """Save the statistics report to JSON file."""
        if output_filename is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_filename = f"gmmad2_comprehensive_stats_{timestamp}.json"

        report_path = os.path.join(self.report_dir, output_filename)

        with open(report_path, "w", encoding="utf-8") as f:
            json.dump(stats, f, indent=2, sort_keys=True, default=str)

        print(f"Statistics report saved to: {report_path}")
        return report_path

    def run_full_analysis(self, file_path: Optional[str] = None) -> Dict[str, Any]:
        """Run the complete analysis and save the report."""
        print("Starting full GMMAD2 record analysis...")

        stats = self.generate_comprehensive_stats(file_path)

        if stats:
            report_path = self.save_stats_report(stats)
            print(f"Analysis complete! Report saved to: {report_path}")

            metadata = stats.get("metadata", {})
            duplicate_info = (
                stats.get("overall_statistics", {}).get("_id", {}).get("duplicate_analysis", {})
            )

            print("\nSummary:")
            print(f"Total records analyzed: {metadata.get('total_records', 0)}")
            print(
                f"Unique _id count: {stats.get('overall_statistics', {}).get('_id', {}).get('unique_count', 0)}"
            )
            print(f"Duplicate records found: {duplicate_info.get('total_duplicate_records', 0)}")
            print(f"Duplicate groups found: {duplicate_info.get('total_duplicate_groups', 0)}")

            if duplicate_info.get("duplication_count_distribution"):
                print("Duplication distribution:")
                for count, groups in sorted(
                    duplicate_info.get("duplication_count_distribution", {}).items()
                ):
                    print(f"  {count} duplicates: {groups} groups")

            print(f"Relationship types found: {len(metadata.get('relationship_types', []))}")
            for rel_type, count in metadata.get("relationship_type_counts", {}).items():
                print(f"- {rel_type}: {count}")
        else:
            print("Analysis failed - no data found.")

        return stats


if __name__ == "__main__":

    file_path = sys.argv[1] if len(sys.argv) > 1 else None

    reporter = GMMAD2RecordStatsReporter()
    results = reporter.run_full_analysis(file_path)
