import os
from typing import Dict

from .reader import CSVReader


class BiGGMapper:
    """BiGG metabolite mapping helper."""

    def __init__(self):
        self.csv_reader = CSVReader()

    def get_bigg_metabolite_mapping(
        self, in_f: str | os.PathLike, delimiter: str = "\t", skip_header: bool = True
    ) -> Dict[str, str]:
        bigg_map: Dict[str, str] = {}
        for fields in self.csv_reader.line_generator(
            in_f, delimiter=delimiter, skip_header=skip_header
        ):
            bigg_id = fields[1].strip()
            metab_key = fields[2].strip().lower()
            bigg_map[metab_key] = bigg_id

        return bigg_map
