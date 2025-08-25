import csv
import os
from typing import Iterator


class CSVReader:
    """Lightweight generator for large CSV/TXT tables."""

    @staticmethod
    def line_generator(
        in_file: str | os.PathLike, delimiter=",", skip_header=True
    ) -> Iterator[list]:
        """Generates lines from a CSV file, yielding each line as a list of strings
        This function opens the specified CSV file, skips the header row, and yields each following line as a list of strings.

        :param skip_header:
        :param in_file: The path to the CSV file.
        :param delimiter: The character used to separate values in the CSV file (default is comma).
        :return: An iterator that yields each line of the CSV file as a list of strings.
        """
        with open(in_file, "r") as in_f:
            reader = csv.reader(in_f, delimiter=delimiter)
            if skip_header:
                next(reader)
            else:
                pass

            for line in reader:
                yield line

    @staticmethod
    def line_generator_for_microbe_disease(
        in_file: str | os.PathLike, skip_header=True
    ) -> Iterator[list[str]]:
        """Yield each CSV line as a list of exactly 24 fields, specially handling disease_species.csv
        rejoining misaligned columns with commas for the 'disease' and 'disease_info' columns.
        """
        EXPECTED_COUNT = 24
        HEAD_COUNT = 2
        TAIL_COUNT = 5
        FIXED_COUNT = 15

        with open(in_file, "r", encoding="utf-8") as f:
            if skip_header:
                next(f)
            else:
                pass

            for line in f:
                parts = [part.strip() for part in line.split(",")]
                if len(parts) == EXPECTED_COUNT:
                    yield parts
                else:
                    alteration_idx = next(
                        i for i, p in enumerate(parts) if p in ("Increase", "Decrease")
                    )
                    head = parts[:HEAD_COUNT]
                    fixed_start = alteration_idx - (FIXED_COUNT - 1)
                    disease = ",".join(parts[2:fixed_start])
                    fixed = parts[fixed_start : alteration_idx + 1]
                    disease_info = ",".join(parts[alteration_idx + 1 : len(parts) - TAIL_COUNT])
                    tail = parts[-TAIL_COUNT:]
                    new_line = head + [disease] + fixed + [disease_info] + tail
                    assert (
                        len(new_line) == EXPECTED_COUNT
                    ), f"Expected {EXPECTED_COUNT} cols, got {len(new_line)}"
                    yield new_line
