import json
import os
import sys
import tempfile
from functools import partial, reduce
from pathlib import Path
from subprocess import run

import pandas as pd

# This is essentially copy-paste the code from select_and_run.py


def replace(string, pattern, replacement):
    string = str(string)
    return string.replace(str(pattern), str(replacement))


def run_wrapper(
    keyvalue,
    input_matrix_path,
    input_metadata_path,
    output_dir,
    delimiter,
    case_only,
    control_only,
):
    key, value = keyvalue
    set_meta = partial(
        replace,
        pattern="<meta>",
        replacement=input_metadata_path.expanduser().absolute(),
    )
    print(f"Processing {key}.")

    assert not (
        case_only and control_only
    ), "Cannot set both case_only and control_only"

    # Make the "case" file
    if not control_only:
        print(f"Making input file {key}_case")
        args = [
            "metasplit",
        ]
        args.extend([set_meta(x) for x in value["case"]])
        args.extend(
            [
                input_matrix_path,
                output_dir / f"{key}_case",
                "--ignore_missing",
                "--input_delimiter",
                delimiter,
                "--always_include",
                "sample",
            ]
        )
        args = [str(x) for x in args]
        print(f"Executing {' '.join(args)}")
        run(args, check=True)
        case = pd.read_csv(output_dir / f"{key}_case")

    # Make the "control" file
    if not case_only:
        print(f"Making input file {key}_control")
        args = [
            "metasplit",
        ]
        args.extend([set_meta(x) for x in value["control"]])
        args.extend(
            [
                input_matrix_path,
                output_dir / f"{key}_control",
                "--ignore_missing",
                "--input_delimiter",
                delimiter,
                "--always_include",
                "sample",
            ]
        )
        args = [str(x) for x in args]
        print(f"Executing {' '.join(args)}")
        run(args, check=True)
        control = pd.read_csv(output_dir / f"{key}_control")

    # Merge case and control together, if we need to
    if case_only:
        matrix = case
    elif control_only:
        matrix = control
    else:
        matrix = pd.merge(case, control, on="sample", how="outer")

    # Calculate row means
    # It's OK to rely on "numeric_only" as all cols except one should be
    # numeric, and they are all relevant
    matrix[key] = matrix.mean(axis=1, numeric_only=True)
    matrix = matrix.filter(items=["sample", key], axis="columns")

    # Delete the useless input files
    if not control_only:
        os.remove(output_dir / f"{key}_case")
    if not case_only:
        os.remove(output_dir / f"{key}_control")

    return matrix


def main(
    queries: dict,
    input_matrix_path: Path,
    input_metadata_path: Path,
    output_path: Path,
    delimiter=",",
    case_only=False,
    control_only=False,
):
    with tempfile.TemporaryDirectory() as tmp:
        run = partial(
            run_wrapper,
            input_matrix_path=input_matrix_path,
            input_metadata_path=input_metadata_path,
            output_dir=Path(tmp),
            delimiter=delimiter,
            case_only=case_only,
            control_only=control_only,
        )
        matrices = list(map(run, queries.items()))

    result: pd.DataFrame = reduce(
        lambda x, y: pd.merge(x, y, how="outer", on="sample"), matrices
    )

    result.to_csv(output_path, index=False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "queries_file", type=Path, help="JSON file with the queries to launch"
    )
    parser.add_argument(
        "input_matrix", type=Path, help="Input (big) expression matrix to subset"
    )
    parser.add_argument(
        "input_metadata", type=Path, help="Input metadata matrix to use to subset"
    )
    parser.add_argument("output_file", type=Path, help="Output .csv file path")
    parser.add_argument("--delimiter", default=",", help="Delimiter for the input")
    parser.add_argument(
        "--case-only",
        action="store_true",
        help="Calculate expression using only case samples?",
    )
    parser.add_argument(
        "--control-only",
        action="store_true",
        help="Calculate expression using only control samples?",
    )

    args = parser.parse_args()

    if args.control_only and args.case_only:
        print("ERROR: Cannot specify both --case-only and --control-only.")
        sys.exit(1)

    with args.queries_file.open("r") as stream:
        queries = json.load(stream)

    main(
        queries=queries,
        input_matrix_path=args.input_matrix,
        input_metadata_path=args.input_metadata,
        output_path=args.output_file,
        delimiter=args.delimiter,
        case_only=args.case_only,
        control_only=args.control_only,
    )
