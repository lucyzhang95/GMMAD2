import csv
import os
from typing import Dict, Iterator, List

import pandas as pd


def line_generator(in_file: str | os.PathLike, delimiter=",", skip_header=False) -> Iterator[list]:
    with open(in_file, "r") as in_f:
        reader = csv.reader(in_f, delimiter=delimiter)
        if skip_header:
            next(reader)
        else:
            pass

        for line in reader:
            yield line


def check_line_fields(data: List[List]) -> Dict[int, int]:
    misaligned_lines = {}
    if not data:
        raise ValueError("Data is empty!")

    reference_length = len(data[0])

    for i, line in enumerate(data):
        if len(line) != reference_length:
            misaligned_lines[i] = len(line)

    return misaligned_lines


def get_columns(df: pd.DataFrame) -> list:
    return df.columns.tolist()


def check_missing_data(df: pd.DataFrame, choices=None):
    if choices is None:
        choices = ["unknown", "missing", "none", "n/a", "null", "not available", ""]
    lookup = {str(c).strip().lower() for c in choices}

    found_missing = {}
    for col in df.columns:
        is_na = df[col].isna()
        is_custom_missing = pd.Series([False] * len(df), index=df.index)
        if df[col].dtype == "object":
            is_custom_missing = df[col].str.strip().str.lower().isin(lookup)
        combined_mask = is_na | is_custom_missing
        if combined_mask.any():
            missing_vals = set(df.loc[combined_mask, col])
            found_missing[col] = missing_vals

    return found_missing
