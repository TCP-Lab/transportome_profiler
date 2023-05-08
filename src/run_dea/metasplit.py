from pathlib import Path
from dataclasses import dataclass
from enum import Enum
from typing import Optional

import subprocess as sb
import re

FILE_VAR_REGEX = re.compile(r"^(.*?)@(.*?)\?(.*?)$")
SELECTION_REGEX = re.compile(r"\[(.*?)(?:,(.*?))*\]")


class ReturnCodeError(Exception):
    pass


class MissingHeaderError(Exception):
    pass


class NoSelectionError(Exception):
    pass


class SelectionSign(Enum):
    POSITIVE = "+"
    NEGATIVE = "-"

@dataclass
class Selection:
    id_variable: str
    """The variable to select resulting IDs from"""
    filter_variable: str
    """The variable to filter by"""
    sign: SelectionSign
    filter_values: list[str]
    """The values to filter by"""


class MetaPath:
    def __init__(self, meta_string: str) -> None:
        matches = FILE_VAR_REGEX.search(meta_string)
        if not matches:
            raise ValueError("Could not match input string. Is it malformed?")
        # This has:
        # Group 1: The file path
        # Group 2: The variable to select
        # Group 3: all the rest of the query (except the first ?)

        self.original: str = meta_string

        self.file: Path = Path(matches.group(1)).expanduser().resolve()
        assert self.file.exists(), f"File {self.file} not found."
        self.selection_var: str = matches.group(2)

        selections = []
        logic_strings = matches.group(3).split("?")
        for l_str in logic_strings:
            if "!=" in l_str:
                sign = SelectionSign.NEGATIVE
                pieces = l_str.split("!=")
            elif "=" in l_str:
                sign = SelectionSign.POSITIVE
                pieces = l_str.split("=")
            else:
                raise ValueError(f"Invalid selection string {l_str}")

            if match := SELECTION_REGEX.match(pieces[1]):
                # If this matches, then the selection is of the form
                # [aaa,bbb,...] and we need to unpack it
                values = [x for x in match.groups() if x is not None]
            else:
                # If that did not match, then we consider this as just one
                # selection
                values = [pieces[1]]

            selections.append(
                Selection(
                    id_variable=pieces[0],
                    sign=sign,
                    filter_values=values,
                )
            )

        self.selections: list[Selection] = selections

        ## Validity checks
        headers = exec(["xsv", "headers", "-j", self.file]).split("\n")
        assert (
            self.variable in headers
        ), f"Header {self.variable} not in metadata headers"

    def __str__(self) -> str:
        return f"{type(self).__name__} object :: file {self.file} selecting {self.selection_var} on {self.variable} with {self.values}"


def exec(*args, **kwargs) -> str:
    res = sb.run(*args, **kwargs, encoding="UTF-8", capture_output=True)

    if res.returncode != 0:
        raise ReturnCodeError(
            f"Process exited with code {res.returncode}:\n{res.stderr}"
        )

    return res.stdout.strip()


def xsv_select(
    file: Path,
    var: str,
    delim: str = ",",
    include_header: bool = False,
    output_file: Optional[Path] = None,
) -> list(str):
    assert file.exists(), f"Cannot run xsv on file {file} that does not exist"

    command = ["xsv", "select", "-d", delim, var, file]
    if output_file:
        command.extend(["-o", output_file])
    values: list[str] = exec(command).split("\n")
    if not include_header:
        values.pop(0)

    return values

def indexes_of(list: list[str], selection: list[str]) -> list[int]:
    return [i for i, x in enumerate(list) if x in selection]


def main(
    metadata: list[MetaPath],
    input_file: Path,
    output_file: Path,
    ignore_missing: bool = False,
    input_delimiter: str = ",",
) -> None:
    assert input_file.exists(), f"Input csv {input_file} does not exist."

    cols = []
    # We can now select the columns of interest
    for meta in metadata:
        print(f"Running on {meta.file.name}...")
        # Get the IDs that we need to select
        identifiers = xsv_select(meta.file, meta.selection_var or "0", delim=input_delimiter)
        indexes = []
        for selection in meta.logical_variables:
            positive = []
            negative = []
            # For every selection we need to run a positive-negative selection
            # on the variables
            selection_values = xsv_select(meta.file, selection.id_variable)

            for filter in selection.filter_values:
                if selection.sign == SelectionSign.POSITIVE:
                    positive.extend[indexes_of(selection_values, filter.value)]
                else:
                    negative.extend[indexes_of(selection_values, filter.value)]
                # Now we do a positive - negative
                indexes.extend(list(set(positive) - set(negative)))

        print(f"Found {len(indexes)} for metadata {meta}")

        cols.extend(identifiers[indexes])

    if not cols:
        raise NoSelectionError("Nothing was selected by the metadata directives")

    target_headers = exec(
        ["xsv", "headers", "-j", "-d", input_delimiter, input_file]
    ).split("\n")
    print(len(target_headers))

    if ignore_missing:
        result = [x for x in result if x in target_headers]

        if not result:
            raise NoSelectionError("Nothing survived after removing missing headers")
    else:
        if any([x not in target_headers for x in result]):
            raise MissingHeaderError(
                "Some metadata selected headers are not in the subsetted matrix"
            )

    print(f"Selecting {len(result)} results...")
    selection_str = ",".join([f'"{x}"' for x in result])
    exec(
        [
            "xsv",
            "select",
            "-d",
            input_delimiter,
            selection_str,
            input_file,
            "-o",
            output_file,
        ]
    )

    print("Done!")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("metadata", type=str, nargs="+", help="Metadata to use")
    parser.add_argument("input_csv", type=Path, help="csv to process")
    parser.add_argument("output_csv", type=Path, help="output csv")
    parser.add_argument(
        "--ignore_missing",
        action="store_true",
        help="If set, ignore IDs in the metadata not in the target csv to subset",
    )
    parser.add_argument(
        "--input_delimiter",
        type=str,
        default=",",
        help="The delimiter to use in the input file",
    )

    args = parser.parse_args()

    main(
        metadata=[MetaPath(x) for x in args.metadata],
        input_file=args.input_csv,
        output_file=args.output_csv,
        ignore_missing=args.ignore_missing,
        input_delimiter=args.input_delimiter,
    )
