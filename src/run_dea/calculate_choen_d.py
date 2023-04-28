#!/usr/bin/env python

from pathlib import Path
from statistics import variance as var, mean
import json
import logging
from logging import StreamHandler

import pandas as pd
from tqdm import tqdm
from colorama import Fore, Style, Back

## >>>> Logging setup
class ColorFormatter(logging.Formatter):
    # Change this dictionary to suit your coloring needs!
    COLORS = {
        "WARNING": Fore.YELLOW,
        "ERROR": Fore.RED,
        "DEBUG": Style.BRIGHT + Fore.MAGENTA,
        "INFO": Fore.GREEN,
        "CRITICAL": Style.BRIGHT + Fore.RED,
    }

    def format(self, record):
        reset = Fore.RESET + Back.RESET + Style.NORMAL
        color = self.COLORS.get(record.levelname, "")
        if color:
            record.name = Style.BRIGHT + Fore.CYAN + record.name + reset
            if record.levelname != "INFO":
                record.msg = color + record.msg + reset
            record.levelname = color + record.levelname + reset
        return logging.Formatter.format(self, record)


# Harry parses stuff now!
log = logging.getLogger("Harry")  # Keep this at the module level name
log.setLevel(logging.DEBUG)
log.propagate = False
# Keep this at DEBUG - set levels in handlers themselves

format = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
console_formatter = ColorFormatter(format)

stream_h = StreamHandler()
stream_h.setFormatter(console_formatter)
stream_h.setLevel(logging.DEBUG)

log.addHandler(stream_h)
## <<<< Logging setup

def choen_d(case: list[float], control: list[float]) -> float:
    """Calculate Cohen's d metric from two series of numerical values

    Args:
        case (list[float]): The 'case' or 'test' values to compare to.
        control (list[float]): The 'control' or 'background" values to compare with.

    Returns:
        float: The calculated Cohen's D.
    """
    assert any(case.isna()) is False, "Case must not have missing values"
    assert any(control.isna()) is False, "Control must not have missing values"

    assert len(case) > 1, "Case must have at least two elements"
    assert len(control) > 1, "Control must have at least two elements"

    # Calculate pooled variance
    pooled_var = (
        ((len(case) - 1) * var(case) + (len(control) - 1) * var(control))
        / (len(case) + len(control) - 2)
    ) ** 0.5

    if pooled_var == 0:
        # The pooled variance is 0. All numbers from case and control are 
        # indentical. Thus, the result is 0
        return 0

    return (mean(case) - mean(control)) / pooled_var


def process_chunk(chunk: pd.DataFrame, samples: dict) -> pd.DataFrame:
    def _process_df_row(row: pd.Series) -> pd.Series:
        result = {}
        for comparison, comp_samples in samples.items():
            result[comparison] = choen_d(
                case = row[comp_samples['case']],
                control = row[comp_samples['control']]
            )
        
        result = pd.Series(result)

        return result
    
    processed_chunk = {}
    for label, row in chunk.iterrows():
        processed_chunk[label] = _process_df_row(row)
    
    return pd.DataFrame(processed_chunk)


def main(
    metadata_path: Path, matrix_path: Path, matches_path: Path, outfile_path: Path, chunksize: int
) -> None:
    
    # Read the data that we can load into memory to memory
    log.info("Reading in metadata...")
    metadata: pd.DataFrame = pd.read_csv(metadata_path, index_col=0)

    metadata = metadata.dropna(subset="Sample")

    # Read a tiny bit of the big file and get the colnames
    reader = pd.read_csv(matrix_path, chunksize=1, sep = "\t")
    sample_ids = reader.get_chunk(1).columns.to_list()

    metadata = metadata[[x in sample_ids for x in metadata["Sample"]]]

    log.info("Reading in matches file...")
    with matches_path.open("r") as stream:
        matches = json.load(stream)

    # TODO: Some way to check if the metadata and data are compatible?

    # Extract the sample names of the case/control matches
    log.info("Generating sample selection lists...")
    samples = {}
    for comparison, match in matches.items():
        log.info(f"Generating {comparison}")
        # TODO: 'Sample' should be a var
        control_samples = metadata["Sample"][metadata[match["control_variable"]] == match["control_value"]]
        case_samples = metadata["Sample"][metadata[match["case_variable"]] == match["case_value"]]

        control_samples = control_samples.to_list()
        case_samples = case_samples.to_list()

        assert len(control_samples) != 0, f"{comparison} selected zero control samples."
        assert len(case_samples) != 0, f"{comparison} selected zero case samples."

        #assert any(pd.isna(x) for x in control_samples), f"{comparison} created NA control samples."
        #assert any([pd.isna(x) for x in case_samples]), f"{comparison} created NA case samples."

        samples[comparison] = {'control': control_samples, 'case': case_samples}

    log.info("Getting handler for matrix...")
    reader = pd.read_csv(matrix_path, chunksize=chunksize, sep="\t")

    log.info("Starting to calculate....")
    write_header = True # A hack to write the header just once
    for chunk in tqdm(reader):
        processed_chunk = process_chunk(chunk, samples=samples)
        
        processed_chunk.to_csv(
            outfile_path, header=write_header, mode = 'a'
        )
        write_header = False

    log.info("Finished writing chunks! Done.")
    
    return None




if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "meta",
        type=Path,
        help="A path to a metadata .csv file with as rows the sample names and as columns variables that can be used to subset the sample space with.",
    )
    parser.add_argument(
        "expression_matrix",
        type=Path,
        help="A path to a expression matrix .csv file with as rows the genes and as columns the log2(x + 1) expression of that gene in that sample.",
    )
    parser.add_argument(
        "matches",
        type=Path,
        help=r"A path to a json file with a dictionary of comparison_name : {'case_variable': var, 'case_value': value 'control_variable': variable, 'control_value': value} to query the metadata file with.",
    )

    parser.add_argument(
        "outfile",
        type=Path,
        help="A path to the output .csv file with cohen's d values",
    )

    parser.add_argument(
        "--chunksize",
        type=int,
        help="Number of lines to read for each parsed chunk",
        default=10,
    )

    args = parser.parse_args()

    main(
        metadata_path=args.meta,
        matrix_path=args.expression_matrix,
        matches_path=args.matches,
        outfile_path=args.outfile,
        chunksize=args.chunksize,
    )
