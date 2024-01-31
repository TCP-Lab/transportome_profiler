#!/usr/bin/env python

"""Metasample

This script performs stratified sampling on an input .csv file, reducing the
number of columns and rows.

To run, you must install `pandas` and have `xsv` (https://github.com/BurntSushi/xsv)
in your PATH.
"""

## --- LICENSE ---
# This is free and unencumbered software released into the public domain.
#
# Anyone is free to copy, modify, publish, use, compile, sell, or
# distribute this software, either in source code form or as a compiled
# binary, for any purpose, commercial or non-commercial, and by any
# means.
#
# In jurisdictions that recognize copyright laws, the author or authors
# of this software dedicate any and all copyright interest in the
# software to the public domain. We make this dedication for the benefit
# of the public at large and to the detriment of our heirs and
# successors. We intend this dedication to be an overt act of
# relinquishment in perpetuity of all present and future rights to this
# software under copyright law.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
#
# For more information, please refer to <https://unlicense.org>
## --- ---

from typing import TextIO, Optional
from math import ceil
import pandas as pd
import subprocess as sb
from random import sample
from pathlib import Path
import os
import tempfile

import logging

logging.basicConfig(level=logging.DEBUG)

log = logging.getLogger(__name__)


def sample_size_string_to_percentage(sample_size: str, total_size: int) -> int:
    """Convert a sample size string to a percentage.

    E.g. '30%' > 0.3, '200' (with a 400 total size) > 0.5, etc...

    Fails if the computed size is greater than 1 (i.e. more than 100%)
    """
    if not sample_size.endswith("%"):
        sample_size = int(sample_size) / total_size
    else:
        sample_size = float(sample_size.strip("%")) / 100

    if sample_size > 1:
        raise ValueError(
            f"Computed sample size {sample_size * total_size} is greater than total size {total_size}"
        )

    return sample_size


class NumberCompressor:
    """Compresses sequences of numbers to call xsv in a more concise way

    If the selection was raw (e.g. 1,2,3,4,5) you might hit the cap of command
    length. Using a compressor makes it so that such long strings are compressed
    and cause such problems less frequently.
    """

    def __init__(self) -> None:
        self.compressed = []
        self.expected_number = None
        self.buffer = []

    def flush(self):
        if not self.buffer:
            return
        if len(self.buffer) == 1:
            self.compressed.append(f"{self.buffer[0]}")
        else:
            self.compressed.append(f"{min(self.buffer)}-{max(self.buffer)}")
        self.buffer = []

    def gobble(self, i: int) -> None:
        if self.expected_number is None or i == self.expected_number:
            self.expected_number = i + 1
            self.buffer.append(i)
            return

        self.flush()
        self.buffer.append(i)
        self.expected_number = i + 1


def exec(*args, **kwargs) -> str:
    res = sb.run(
        *args, **kwargs, encoding="UTF-8", capture_output=True, errors="replace"
    )

    if res.returncode != 0:
        raise OSError(f"Process exited with code {res.returncode}:\n{res.stderr}")

    return res.stdout.strip()


def xsv_select(
    file: Path,
    var: str,
    delim: str = ",",
    include_header: bool = False,
    output_file: Optional[Path] = None,
) -> list[str]:
    assert file.exists(), f"Cannot run xsv on file {file} that does not exist"

    command = ["xsv", "select", "-d", delim, var, file.as_posix()]
    if output_file:
        command.extend(["-o", output_file.as_posix()])
    values: list[str] = exec(command).split("\n")
    if not include_header:
        values.pop(0)

    return values


def metasample(
    input: Path,
    output: Path,
    metadata: TextIO,
    meta_rowname_var: str,
    metavars: list[str],
    sample_size: str,
    row_sample_size: str,
    always_include: str,
):
    log.debug("Reading input header")
    with input.open("r+") as stream:
        input_header = stream.readline().split(",")
    log.debug("Getting input size")
    with input.open("r+") as stream:
        input_size = sum(1 for _ in stream)

    log.debug("Parsing variables")
    sample_size = sample_size_string_to_percentage(sample_size, len(input_header))
    # This is retarded but bear with me
    row_sample_size = sample_size_string_to_percentage(row_sample_size, input_size)
    row_sample_size = ceil(row_sample_size * input_size)

    log.debug("Reading metadata")
    metadata = pd.read_csv(metadata, encoding_errors="replace")
    log.debug(f"Metadata cols: {metadata.columns.to_list()}")

    if not all([x in metadata.columns for x in metavars]):
        # TODO: Make better error msg
        raise ValueError("Some selected metadata variables were not found.")

    if not meta_rowname_var in metadata.columns:
        raise ValueError("Row name column not fonud in metadata")

    log.debug("Grouping metadata")

    def get_labels(x):
        return x.index.to_list()

    metadata = metadata.set_index(meta_rowname_var)
    groups = metadata.groupby(metavars).apply(get_labels)

    selected = []
    log.debug("Sampling groups")
    for group in groups:
        new_group = sample(group, ceil(len(group) * sample_size))
        log.debug(
            f"Sampled group. Was {len(group)}, now is {len(new_group)} ({((len(new_group) - len(group))/len(group))*100:.2f})"
        )
        selected.extend(new_group)

    always_include = always_include.split(",")
    selected.extend(always_include)

    log.debug("Compressing selection")
    compressor = NumberCompressor()
    for i, item in enumerate(input_header, start=1):
        if item in selected:
            compressor.gobble(i)

    compressor.flush()

    with tempfile.NamedTemporaryFile() as temp:
        log.debug("Executing xsv - selecting columns")
        compressed_str = ",".join(compressor.compressed)
        xsv_select(input, compressed_str, include_header=True, output_file=Path(temp.name))
        log.debug("Executing xsv - sampling rows")
        exec(
            [
                "xsv",
                "sample",
                "-o",
                f"{output.as_posix()}",
                f"{row_sample_size}",
                f"{temp.name}",
            ]
        )

    print("Done!")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "input_matrix", help="A path to a csv table to sample", type=Path
    )
    parser.add_argument(
        "output_matrix", help="A path to the output csv table to create", type=Path
    )
    parser.add_argument("metadata", help="A path to the input metadata to consider")
    parser.add_argument(
        "row_names_var",
        help="The name of the variable that holds the row names in the input",
    )
    parser.add_argument(
        "--metavars",
        help="A comma-separated list of variables in the metadata to consider when sampling",
    )
    parser.add_argument(
        "--always-include",
        help="A comma-separated list of column names in the input to always include",
    )
    parser.add_argument(
        "col_sample_size",
        help="Either an integer of the number of columns in the output or a percentage (with a percent sign) of the initial columns to sample (round up).",
    )
    parser.add_argument(
        "row_sample_size", help="Same as col_sample_size but for row number."
    )

    args = parser.parse_args()

    if args.metavars:
        metavars = args.metavars.split(",")
    else:
        with args.input_matrix.open("r+") as stream:
            header = next(stream).split(",")
        metavars = [x for x in header if x != args.row_names_var]

    metasample(
        input=args.input_matrix,
        output=args.output_matrix,
        metadata=args.metadata,
        metavars=metavars,
        meta_rowname_var=args.row_names_var,
        sample_size=args.col_sample_size,
        row_sample_size=args.row_sample_size,
        always_include=args.always_include,
    )
