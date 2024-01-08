from pathlib import Path
import json
from subprocess import run
from functools import partial
import os

import multiprocessing as mp

def replace(string, pattern, replacement):
    string = str(string)
    return string.replace(str(pattern), str(replacement))

def run_wrapper(keyvalue, input_matrix_path, input_metadata_path, output_dir, delimiter):
    key, value = keyvalue
    set_meta = partial(replace, pattern = "<meta>", replacement = input_metadata_path.expanduser().absolute())
    print(f"Processing {key}.")
    # Make the "case" file
    args = ["metasplit",]
    args.extend([set_meta(x) for x in value["case"]])
    args.extend([input_matrix_path, output_dir / f"{key}_case", "--ignore_missing", "--input_delimiter", delimiter, "--always_include", "sample"])
    args = [str(x) for x in args]
    run(args, check=True)

    # Make the "control" file
    args = ["metasplit",]
    args.extend([set_meta(x) for x in value["control"]])
    args.extend([input_matrix_path, output_dir / f"{key}_control", "--ignore_missing", "--input_delimiter", delimiter, "--always_include", "sample"])
    args = [str(x) for x in args]
    run(args, check=True)

    # Now we can run run_deseq.R
    dea_args = ["generanker", output_dir / f"{key}_case", output_dir / f"{key}_control", "--output-file", output_dir / f"{key}_deseq.csv", "deseq_shrinkage", "--id-col", "sample"]
    run(dea_args)

    # Delete the useless input files
    os.remove(output_dir / f"{key}_case")
    os.remove(output_dir / f"{key}_control")


def main(
        queries = dict,
        input_matrix_path = Path,
        input_metadata_path = Path,
        output_dir = Path,
        delimiter = ",",
        cpus = None,
):
    run = partial(run_wrapper, input_matrix_path = input_matrix_path, input_metadata_path=input_metadata_path, output_dir=output_dir, delimiter=delimiter)
    print("Spawning pool of workers...")
    with mp.Pool(cpus or mp.cpu_count()) as pool:
        pool.map(run, queries.items())


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("queries_file", type=Path, help="JSON file with the queries to launch")
    parser.add_argument("input_matrix", type=Path, help="Input (big) expression matrix to subset")
    parser.add_argument("input_metadata", type=Path, help="Input metadata matrix to use to subset")
    parser.add_argument("output_dir", type=Path, help="Output directory to save files to")
    parser.add_argument("--delimiter", default = ",", help="Delimiter for the input")
    parser.add_argument("--cpus", type=int, help="Number of CPUS to use. If unspecified, runs with one thread per available core.")

    args = parser.parse_args()

    with args.queries_file.open("r") as stream:
        queries = json.load(stream)
    
    main(
        queries=queries,
        input_matrix_path=args.input_matrix,
        input_metadata_path=args.input_metadata,
        output_dir=args.output_dir,
        delimiter=args.delimiter,
        cpus=args.cpus
    )


