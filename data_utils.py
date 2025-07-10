import pandas as pd


def get_columns(df: pd.DataFrame) -> list:
    return df.columns.tolist()


def check_missing_data(df: pd.DataFrame, choices=None):
    default_choices = ["unknown", "missing", "none", "n/a", "null", "not available", ""]
    lookup = {str(c).strip().lower() for c in (choices if choices is not None else default_choices)}

    found_missing = {}
    for col in get_columns(df):
        missing_vals = set()
        if df[col].isnull().any():
            missing_vals.add(None)
        for v in df[col].unique():
            if pd.isnull(v):
                continue
            if str(v).strip().lower() in lookup:
                missing_vals.update(v)

        if missing_vals:
            found_missing[col] = missing_vals

    return found_missing
