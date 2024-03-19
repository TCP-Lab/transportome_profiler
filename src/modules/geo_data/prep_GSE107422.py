from pathlib import Path
import gzip
import tarfile
from io import BytesIO, StringIO
import pandas as pd
import sys
import re
from copy import copy
from tqdm import tqdm
from functools import reduce

def eprint(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)

def extract_refseq_id(original: str) -> str:
    matcher = re.compile(r"gi\|[0-9]+?\|ref\|(N[M|R]_[0-9]+?)\.[0-9]+\|")
    res = matcher.match(original)
    
    return res.group(1)

def main(args):
    
    files = {}
    with tarfile.open(args.raw, "r") as conn:
        members = conn.getmembers()
        for member in members:
            if not member.isfile():
                continue
            
            files[member.name] = BytesIO(conn.extractfile(member).read())

    eprint(f"Read {len(files)} files. Decompressing")
    
    data = []
    for k, v in files.items():
        eprint(f"Decompressing {k}")
        conn = gzip.decompress(v.read())
        data.append(pd.read_csv(StringIO(conn.decode("utf-8")), sep = "\t"))

    eprint("Finished loading csvs. Parsing IDs")
    # Each id is of format `gi|...|ref|<refseq_id>|`. We need the refseq ID
    def rename_ids(x):
        y = copy(x)
        y["unique_id"] = list(map(extract_refseq_id, x["unique_id"]))

        return y

    data = list(map(rename_ids, tqdm(data)))
    
    eprint("Merging...")
    merged = reduce(lambda x, y: pd.merge(x, y, on="unique_id", how="outer"), tqdm(data))

    out = StringIO()
    merged.to_csv(out, index=False)
    out.seek(0)
    print(out.read())


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("raw", help="Path to raw compressed data", type = Path)

    args = parser.parse_args()

    main(args)
