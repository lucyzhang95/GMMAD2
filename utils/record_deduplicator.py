import json
from typing import Dict, Iterator


def _merge_xrefs(target_xrefs: Dict, source_xrefs: Dict) -> None:
    """Merges source xrefs into target xrefs in-place."""
    for key, source_value in source_xrefs.items():
        if key not in target_xrefs:
            target_xrefs[key] = source_value
            continue

        target_value = target_xrefs[key]

        combined = []
        if isinstance(target_value, list):
            combined.extend(target_value)
        else:
            combined.append(target_value)

        if isinstance(source_value, list):
            combined.extend(source_value)
        else:
            combined.append(source_value)

        unique_items = list(dict.fromkeys(combined))

        if len(unique_items) == 1:
            target_xrefs[key] = unique_items[0]
        else:
            target_xrefs[key] = unique_items


def _merge_publications(target_assoc: Dict, source_assoc: Dict) -> None:
    """Merges publication info into a list in the target association node."""
    source_pubs = source_assoc.get("publications")
    if not source_pubs:
        return

    if "publications" not in target_assoc:
        target_assoc["publications"] = []
    elif not isinstance(target_assoc["publications"], list):
        target_assoc["publications"] = [target_assoc["publications"]]

    existing_pmids = {
        pub.get("pmid") for pub in target_assoc["publications"] if isinstance(pub, dict)
    }

    if isinstance(source_pubs, dict) and source_pubs.get("pmid") not in existing_pmids:
        target_assoc["publications"].append(source_pubs)
    elif isinstance(source_pubs, list):
        for pub in source_pubs:
            if isinstance(pub, dict) and pub.get("pmid") not in existing_pmids:
                target_assoc["publications"].append(pub)
                existing_pmids.add(pub.get("pmid"))

    if len(target_assoc["publications"]) == 1:
        target_assoc["publications"] = target_assoc["publications"][0]


def _create_fingerprint(record: Dict) -> str:
    """
    Creates a unique fingerprint for a record, ignoring fields
    that can be merged (xrefs, publications, _id).
    """
    subject = record.get("subject", {}).copy()
    obj = record.get("object", {}).copy()
    association = record.get("association", {}).copy()

    subject.pop("xrefs", None)
    obj.pop("xrefs", None)
    association.pop("publications", None)

    subject_fp = json.dumps(subject, sort_keys=True)
    obj_fp = json.dumps(obj, sort_keys=True)
    assoc_fp = json.dumps(association, sort_keys=True)

    return f"{subject_fp}|{obj_fp}|{assoc_fp}"


class StreamDeduplicator:
    """
    Processes a stream of records to remove duplicates and merge records
    based on predefined rules, optimized for low memory usage.
    """

    def __init__(self):
        self.processed_records: Dict[str, Dict] = {}

    def process_record(self, record: Dict) -> None:
        """
        Processes a single record, merging it with existing records if it's a
        duplicate or adding it to the collection if it's unique.
        """
        fingerprint = _create_fingerprint(record)

        if fingerprint not in self.processed_records:
            self.processed_records[fingerprint] = record
        else:
            existing_record = self.processed_records[fingerprint]

            # merge subject xrefs
            if "xrefs" in record.get("subject", {}):
                _merge_xrefs(
                    existing_record.setdefault("subject", {}).setdefault("xrefs", {}),
                    record["subject"]["xrefs"],
                )

            # merge object xrefs
            if "xrefs" in record.get("object", {}):
                _merge_xrefs(
                    existing_record.setdefault("object", {}).setdefault("xrefs", {}),
                    record["object"]["xrefs"],
                )

            #  merge publications
            if "publications" in record.get("association", {}):
                _merge_publications(
                    existing_record.setdefault("association", {}),
                    record["association"],
                )

    def get_results(self) -> Iterator[Dict]:
        """
        Returns an iterator for the final, deduplicated and merged records.
        """
        yield from self.processed_records.values()
