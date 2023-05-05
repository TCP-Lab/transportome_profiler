from pathlib import Path
from typing import Optional

import subprocess as sb
import re

META_STRING_REGEX = re.compile(r"(.*?)(?:\?(.*?))*?@(.*?)=(.*?)$")
VALUE_REGEX = re.compile(r"\\")

class ReturnCodeError(Exception):
    pass

class MissingHeaderError(Exception):
    pass

class NoSelectionError(Exception):
    pass

class MetaPath:
    def __init__(self, meta_string: str) -> None:
        matches = META_STRING_REGEX.search(meta_string)

        self.original: str = meta_string
        self.file: Path = Path(matches.group(1)).expanduser().resolve()

        assert self.file.exists(), f"File {self.file} not found."

        self.variable: str = matches.group(3)
        self.selection_var: Optional[str] = matches.group(2) or None

        self.values: Optional[list[str]] = None

        value: Optional[str] = matches.group(4) or None
        # Split a possible value
        if value and value.startswith("[") and value.endswith("]"):
            value = value.strip("[]").split(",")
            self.values = value
        elif value:
            self.values = [value]
        else:
            self.values = None
        
        ## Validity checks
        headers = exec(["xsv", "headers", "-j", self.file]).split("\n")

        assert self.variable in headers, f"Header {self.variable} not in metadata headers"
    
    def __str__(self) -> str:
        return f"{type(self).__name__} object :: file {self.file} selecting {self.selection_var} on {self.variable} with {self.values}"

def exec(*args, **kwargs) -> str:
    res =  sb.run(*args, **kwargs, encoding="UTF-8", capture_output=True)

    if res.returncode != 0:
        raise ReturnCodeError(f"Process exited with code {res.returncode}:\n{res.stderr}")

    return res.stdout.strip()

def main(
    metadata: list[MetaPath],
    input_file: Path,
    output_file: Path,
    ignore_missing: bool = False,
    input_delimiter: str = ","
) -> None:
    
    assert input_file.exists(), f"Input csv {input_file} does not exist."

    cols = []
    for meta in metadata:
        print(f"Running on {meta.file.name}...")
        values = exec(["xsv", "select", meta.variable, meta.file]).split("\n")
        values.pop(0)
        identifiers = exec(["xsv", "select", meta.selection_var or 0, meta.file]).split("\n")
        identifiers.pop(0)

        result = [identifiers[i] for (i, x) in enumerate(values) if x in meta.values]

        print(f"Found {len(result)} for metadata {meta}")

        cols.extend(list(set(result)))

    
    if not cols:
        raise NoSelectionError("Nothing was selected by the metadata directives")

    target_headers = exec(["xsv", "headers", "-j", "-d", input_delimiter, input_file]).split("\n")
    print(len(target_headers))

    if ignore_missing:
        result = [x for x in result if x in target_headers]

        if not result:
            raise NoSelectionError("Nothing survived after removing missing headers")
    else:
        if any([x not in target_headers for x in result]):
            raise MissingHeaderError("Some metadata selected headers are not in the subsetted matrix")
    
    print(f"Selecting {len(result)} results...")
    selection_str = ",".join([f'"{x}"' for x in result])
    exec(["xsv", "select", "-d", input_delimiter, selection_str, input_file, "-o", output_file])

    print("Done!")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "metadata", type = str, nargs="+", help = "Metadata to use"
    )
    parser.add_argument(
        "input_csv", type= Path, help="csv to process"
    )
    parser.add_argument(
        "output_csv", type= Path, help="output csv"
    )
    parser.add_argument(
        "--ignore_missing", action="store_true", help="If set, ignore IDs in the metadata not in the target csv to subset"
    )
    parser.add_argument(
        "--input_delimiter",type = str, default=",", help="The delimiter to use in the input file"
    )

    args = parser.parse_args()

    main( metadata=[MetaPath(x) for x in args.metadata], input_file=args.input_csv, output_file=args.output_csv, ignore_missing = args.ignore_missing, input_delimiter=args.input_delimiter)