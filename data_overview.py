import os
import pandas as pd
import re


def col_value_contains_digits(value):
    return re.fullmatch(r"\d+", str(value)) is not None


def csv_data_overview(in_file, file_name):
    path = os.getcwd()
    infile_path = os.path.join(path, in_file, file_name)
    assert infile_path, f"The file {infile_path} does not exist."

    df = pd.read_csv(infile_path, low_memory=False)

    # Get the data types of each column
    print("\nData types for each column:")
    print(df.dtypes)

    columns = {idx: col_name for idx, col_name in enumerate(df.columns)}
    print(f"column names: {columns}")

    # Check for unique values in each column
    print("\nUnique values in each column:")
    for column in df.columns:
        unique_values = df[column].unique()
        print(f"Column: {column}")
        print(f"Number of unique values: {len(unique_values)}")
        digit_value = [val for val in unique_values if col_value_contains_digits(val)]
        other_value = [val for val in unique_values if not col_value_contains_digits(val)]
        if len(unique_values) > 10:
            print(f"Unique values contain digits: {digit_value[:9]}")
            print(f"Unique values no digits: {other_value[:9]}")
        else:
            print(f"Unique values contain digits: {unique_values}")
            print(f"Unique values no digits: {other_value}")

        print("-" * 50)


if __name__ == "__main__":
    # csv_data_overview("data", "micro_metabolic.csv")
    # disease_meta.csv contains predictions, so better not to process it now, but can use for comparison later
    # csv_data_overview("data", "disease_meta.csv")
    csv_data_overview("data", "meta_gene_net.csv")
