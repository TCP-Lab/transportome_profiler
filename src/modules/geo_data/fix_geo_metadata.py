"""
After hours of trying to find a fix for the horrible, horrible metadata
headers, I've chosen to fuck it, and do all of them manually.

This takes the original metadata file and the accession number and spits out
the fixed metadata.

Just for fun, I wanted to do this with as little state as I could, so it's
function fest!
"""
from pathlib import Path
import sys
import pandas as pd
from typing import Callable, Iterable, Optional, Any
import re


def eprint(*args, **kwargs):
    """Print to stderr, not stdout, ya fool"""
    print(*args, **kwargs, file=sys.stderr)


def lmap(*args, **kwargs):
    """A not lazy map which converts immediately to a list"""
    return list(map(*args, **kwargs))


def identity(original: Any) -> Any:
    """No-op. Take the input and return the output"""
    return original


def apply_to(from_column: str, to_column: str, callable: Callable) -> Callable:
    """Apply a function to a column, saving the output it another column"""

    def wrapped(data: pd.DataFrame):
        data[to_column] = callable(data[from_column])
        return data

    return wrapped


def apply_at(column: str, callable: Callable) -> Callable:
    """Apply a function to a column, replacing it

    A specific case of `apply_to`, but with the same columns.
    """
    return apply_to(column, column, callable)


def take_field(sep: str, position: int) -> Callable:
    """Split a string with `sep` and take the field at `position`"""

    def wrapped(data: Iterable):
        return lmap(lambda x: x.split(sep)[position], data)

    return wrapped


def chain(*args):
    """Chain multiple transformations in one call, applying them in order"""

    def wrapped(input):
        for arg in args:
            input = arg(input)
        return input

    return wrapped


def replace(_from: str, to: str) -> Callable:
    """Replace all of `_from` to `to` in all items of an iterable.

    The call to "map" could be taken out of the fuction, but it would get
    a bit too verbose for my liking.
    """

    def wrapped(data: Iterable):
        return lmap(lambda x: x.replace(_from, to), data)

    return wrapped


def follow_map(map: dict, default: Optional[str] = None) -> Callable:
    """Replace all instances in an iterable with new stuff

    Takes a hashmap with as keys the strings to replace and as values the
    strings to replace to, as in {'old': 'new'}.
    """

    def wrapped(data: Iterable):
        def try_convert(x):
            try:
                return map[x]
            except KeyError as e:
                if default:
                    return default
                raise e

        return lmap(try_convert, data)

    return wrapped


def follow_regmap(map: dict, default: Optional[str] = None) -> Callable:
    """Replace all instances that match a regex in an iterable with new stuff

    Like `follow_map`, but taking regex expressions as keys and replacing if
    any of them, in order, find a match.
    """

    def wrapped(data: Iterable):
        def try_convert(x):
            for key in map.keys():
                if re.search(key, x):
                    return map[key]
            raise ValueError(f"No regex key matches pattern: {x}")

        return lmap(try_convert, data)

    return wrapped


def set_to(column: str, fill: str) -> Callable:
    """Set a column to some value, usually a static string"""

    def wrapped(data: pd.DataFrame):
        data[column] = fill
        return data

    return wrapped


OPERATIONS: dict[str, Callable] = {
    "GSE22260": chain(
        apply_to(
            "characteristics_ch1",
            "status",
            follow_map(
                {
                    "tissue: Prostate cancer tissue": "case",
                    "tissue: Normal prostate tissue": "control",
                }
            ),
        ),
        set_to("type", "COAD"),
    ),
    "GSE29580": chain(
        apply_to(
            "characteristics_ch1",
            "status",
            follow_map(
                {
                    "tissue type: Bowel adenocarcinoma tissue": "case",
                    "tissue type: normal bowel tissue": "control",
                }
            ),
        ),
        set_to("type", "COAD"),
    ),
    "GSE121842": chain(
        apply_to(
            "source_name_ch1",
            "status",
            follow_regmap(
                {
                    "colorectal cancer tissue": "case",
                    "colorectal pericarcinomatous tissue": "control",
                }
            ),
        ),
        set_to("type", "COAD"),
    ),
    "GSE159857": chain(
        apply_to(
            "characteristics_ch1_5",
            "status",
            follow_regmap({"Normal$": "control", "Tumor$": "case"}),
        ),
        apply_to(
            "characteristics_ch1_2",
            "type",
            follow_regmap({"Adenoma$": "LUAD", "Squamus$": "LUSC"}),
        ),
    ),
}


def main(metadata: Path, accession: str):
    data = pd.read_csv(metadata)
    result = OPERATIONS[accession](data)
    result.to_csv(sys.stdout, index=False)


def stem_of(path: Path) -> str:
    while path.suffixes:
        path = Path(path.stem)
    return str(path)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("metadata", type=Path)
    parser.add_argument(
        "--id",
        type=str,
        help="Accession ID. If not specified, uses the name of the metadata file.",
    )

    args = parser.parse_args()

    main(
        metadata=args.metadata, accession=args.id if args.id else stem_of(args.metadata)
    )
