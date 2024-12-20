"""Preprocess the results from a pipeline and extract information from them

This takes the .tar.gz file(s) made by running the `heatmaps` and `geo_heatmap`
pipelines and extracts useful information from them.
"""

import logging
import os
import re
from functools import reduce
from pathlib import Path
from typing import Callable

import pandas as pd

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)


def find_re(regex: str) -> Callable:
    def matcher(x: str):
        return bool(re.match(regex, x))

    return matcher


def negate(fn: Callable):
    def wrapped(*args, **kwargs):
        return not fn(*args, **kwargs)

    return wrapped


strip_flags = negate(find_re(".*flag.*"))

REQUIRED_FILES = {
    "data/deas": {"action": "merge", "filters": [strip_flags], "id_col": "sample"},
    "data/geo": {
        "action": "merge",
        "filters": [strip_flags, find_re(".*\.dea\.csv")],
        "id_col": "gene_id",
    },
}


def remove_suffixes(path: Path):
    """Remove all the suffixes from a path, returning the "clean" filename.

    Example:
        remove_suffixes(Path("path/to/file.txt.gz")) == "file"
    """
    return path.name.split(".")[0]


def merge_deas(files: list[Path], merge_col: str = "sample") -> pd.DataFrame:
    data = {}
    for file in files:
        log.info(f"Reading in {file}...")
        data[remove_suffixes(Path(file))] = pd.read_csv(file)
    # The files have all the same structure: a col with 'sample' and one with
    # 'ranking'. We must rename the 'ranking' col with the name of the file
    # and then do a many-way merge
    log.info("Merging data...")
    renamed = {k: v.rename(columns={"ranking": k}) for k, v in data.items()}
    merged = reduce(lambda x, y: pd.merge(x, y, on=merge_col), renamed.values())

    return merged


def main(args):
    log.info(f"Reading in {args.input_dir}")

    if args.slug == "auto":
        slug = remove_suffixes(args.input_tarball)
    else:
        slug = args.slug
    log.info(f"Slug is '{slug}'.")

    log.info("Looking for files...")
    files = args.input_dir.iterdir()
    names = [x.name for x in files]

    all_found = True
    for internal_path in REQUIRED_FILES.keys():
        if internal_path in names:
            log.info(f"Found {internal_path}.")
        else:
            log.error(f"{internal_path} does not exist!")
            all_found = False
    if not all_found:
        raise ValueError("Could not find some required files in archive.")

    for key, value in REQUIRED_FILES.items():
        if value["action"] == "merge":
            log.info(f"Merging contents of {key}")
            files_to_merge = [x for x in names if x.startswith(key) and x != key]
            for fn in value["filters"]:
                files_to_merge = filter(fn, files_to_merge)
            merged = merge_deas(files_to_merge, value.get("id_col", "sample"))

            target = args.output_dir / f"{slug}_{remove_suffixes(Path(key))}.csv"
            log.info(f"Saving to {target}...")

            merged.to_csv(target, index=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "input_dir",
        type=Path,
        help="Path to the directory with all DEAs (usually data/)",
    )
    parser.add_argument(
        "output_dir", type=Path, help="Path to a dir to save the outputs in"
    )
    parser.add_argument(
        "--slug",
        type=str,
        default="auto",
        help="Override the slug give to the output files",
    )

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    main(args)
