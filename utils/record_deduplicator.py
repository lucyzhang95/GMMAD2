import hashlib
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


def _merge_original_names(target_node: Dict, source_node: Dict) -> None:
    """Merges original_name fields into a list, avoiding duplicates."""
    source_name = source_node.get("original_name")
    if not source_name:
        return

    if "original_name" not in target_node:
        target_node["original_name"] = source_name
        return

    target_name = target_node["original_name"]

    if isinstance(target_name, list):
        names_list = target_name.copy()
    else:
        names_list = [target_name]

    if isinstance(source_name, list):
        names_list.extend(source_name)
    else:
        names_list.append(source_name)

    unique_names = list(dict.fromkeys(names_list))

    if len(unique_names) == 1:
        target_node["original_name"] = unique_names[0]
    else:
        target_node["original_name"] = unique_names


def _merge_evidence(target_assoc: Dict, source_assoc: Dict) -> None:
    """
    Merges evidence fields, prioritizing non-ECO:0000000 values.
    If both records have different non-ECO:0000000 evidence, creates a list.
    """
    source_evidence = source_assoc.get("evidence")
    if not source_evidence:
        return

    target_evidence = target_assoc.get("evidence")

    if not target_evidence:
        target_assoc["evidence"] = source_evidence
        return

    def is_default_evidence(evidence):
        return evidence == "ECO:0000000"

    target_list = target_evidence if isinstance(target_evidence, list) else [target_evidence]
    source_list = source_evidence if isinstance(source_evidence, list) else [source_evidence]

    all_evidence = target_list + source_list
    non_default_evidence = [e for e in all_evidence if not is_default_evidence(e)]
    unique_evidence = list(dict.fromkeys(non_default_evidence))

    if not unique_evidence:
        target_assoc["evidence"] = "ECO:0000000"
    elif len(unique_evidence) == 1:
        target_assoc["evidence"] = unique_evidence[0]
    else:
        target_assoc["evidence"] = unique_evidence


def _merge_primary_knowledge_source(target_assoc: Dict, source_assoc: Dict) -> None:
    """Merges primary_knowledge_source fields into a list if the rest of the fields are the same."""
    source_pks = source_assoc.get("primary_knowledge_source")
    if not source_pks:
        return

    target_pks = target_assoc.get("primary_knowledge_source")

    if not target_pks:
        target_assoc["primary_knowledge_source"] = source_pks
        return

    target_list = target_pks if isinstance(target_pks, list) else [target_pks]
    source_list = source_pks if isinstance(source_pks, list) else [source_pks]

    combined = target_list + source_list
    unique_pks = list(dict.fromkeys(combined))

    if len(unique_pks) == 1:
        target_assoc["primary_knowledge_source"] = unique_pks[0]
    else:
        target_assoc["primary_knowledge_source"] = unique_pks


def _create_fingerprint(record: Dict) -> str:
    """
    Creates a unique fingerprint for a record, ignoring only the specific fields
    that can be merged. Records must be identical in all other fields to be
    considered duplicates.

    Mergeable fields:
    - subject: xrefs, original_name
    - object: xrefs, original_name
    - association: publications, primary_knowledge_source, evidence
    """
    subject = record.get("subject", {}).copy()
    obj = record.get("object", {}).copy()
    association = record.get("association", {}).copy()

    subject.pop("xrefs", None)
    subject.pop("original_name", None)
    obj.pop("xrefs", None)
    obj.pop("original_name", None)
    association.pop("publications", None)
    association.pop("primary_knowledge_source", None)
    association.pop("evidence", None)

    subject_fp = json.dumps(subject, sort_keys=True)
    obj_fp = json.dumps(obj, sort_keys=True)
    assoc_fp = json.dumps(association, sort_keys=True)

    fingerprint_str = f"{subject_fp}|{obj_fp}|{assoc_fp}"

    return hashlib.md5(fingerprint_str.encode()).hexdigest()


def _merge_records(existing_record: Dict, new_record: Dict) -> None:
    """Merge new_record into existing_record."""
    if "xrefs" in new_record.get("subject", {}):
        _merge_xrefs(
            existing_record.setdefault("subject", {}).setdefault("xrefs", {}),
            new_record["subject"]["xrefs"],
        )

    if "original_name" in new_record.get("subject", {}):
        _merge_original_names(
            existing_record.setdefault("subject", {}),
            new_record["subject"],
        )

    if "xrefs" in new_record.get("object", {}):
        _merge_xrefs(
            existing_record.setdefault("object", {}).setdefault("xrefs", {}),
            new_record["object"]["xrefs"],
        )

    if "original_name" in new_record.get("object", {}):
        _merge_original_names(
            existing_record.setdefault("object", {}),
            new_record["object"],
        )

    if new_record.get("association"):
        existing_assoc = existing_record.setdefault("association", {})
        source_assoc = new_record["association"]

        if "publications" in source_assoc:
            _merge_publications(existing_assoc, source_assoc)

        if "primary_knowledge_source" in source_assoc:
            _merge_primary_knowledge_source(existing_assoc, source_assoc)

        if "evidence" in source_assoc:
            _merge_evidence(existing_assoc, source_assoc)


class StreamDeduplicator:
    """
    In-memory deduplicator.
    Processes a stream of records to remove duplicates and merge records based on the above rules.
    """

    def __init__(self):
        self.processed_records: Dict[str, Dict] = {}

    def is_duplicate(self, record: Dict) -> bool:
        """Check if the record is a duplicate without processing it."""
        fingerprint = _create_fingerprint(record)
        return fingerprint in self.processed_records

    def process_record(self, record: Dict) -> Dict:
        """
        Processes a single record, merging it with existing records if it's a
        duplicate or adding it to the collection if it's unique.
        Returns the processed record.
        """
        fingerprint = _create_fingerprint(record)

        if fingerprint not in self.processed_records:
            self.processed_records[fingerprint] = record.copy()
        else:
            existing_record = self.processed_records[fingerprint]
            _merge_records(existing_record, record)

        return self.processed_records[fingerprint]

    def get_results(self) -> Iterator[Dict]:
        """Returns an iterator for the final, deduplicated and merged records."""
        yield from self.processed_records.values()
