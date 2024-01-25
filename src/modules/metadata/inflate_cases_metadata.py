import pandas as pd
from pathlib import Path
from typing import Callable

from sys import stdout, stderr

def eprint(*args, **kwargs):
    print(*args, file=stderr, **kwargs)

def clean_ids(ids: list[str]) -> list[str]:
    for id in ids:
        if not id.startswith("TCGA") and len(id) != 15:
            yield id
            continue
        yield id[:12]

def split_path(path: str, default_col: Callable = lambda: "samples") -> (Path, str):
    """For compactdness I use the /path/to/file.csv@column_with_ids shorthand

    This function parses it to a Path and a string with the column, if any.
    """

    parts = path.split("@")
    if len(parts) == 1:
        return (Path(path), default_col())
    
    if len(parts) == 2:
        return (Path(parts[0]), parts[1])

    raise ValueError(f"Path {path} could not be parsed as valid - too many @ symbols.")


def yield_delim(path: Path) -> str:
    ext = path.suffix
    match ext:
        case '.tsv':
            return '\t'
        case '.csv':
            return ','
        case _:
            eprint(f"Warning: unknown input type '{ext}'. Fallback to ','")
            return ','

def check_col(data, col):
    if col not in data.columns:
        raise ValueError(f"{col} not in columns: {data.columns}")

def main(case_path, sample_path, case_col, sample_col):
    cases = pd.read_csv(case_path, sep=yield_delim(case_path), encoding_errors='replace')
    samples = pd.read_csv(sample_path, sep=yield_delim(sample_path), encoding_errors='replace')

    check_col(cases, case_col)
    check_col(samples, sample_col)

    eprint("Inflating case metadata to sample metadata")

    samples["case_id"] = list(clean_ids(samples[sample_col].to_list()))

    data = samples.merge(cases, left_on="case_id", right_on=case_col, how="outer", validate="many_to_one")

    data = data.drop(columns=["case_id"])

    return data


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    
    parser.add_argument("case_metadata", type=str)
    parser.add_argument("sample_metadata", type=str)

    args = parser.parse_args()
    
    case_path, case_col = split_path(args.case_metadata)
    sample_path, sample_col = split_path(args.sample_metadata)

    data = main(case_path=case_path, sample_path=sample_path, case_col=case_col, sample_col=sample_col)

    data.to_csv(stdout, index=False)
