from pathlib import Path
import pandas as pd
import sys
import Levenshtein


def find_match(item, possibilities):
    scores = {p: Levenshtein.jaro_winkler(p, item) for p in possibilities}
    
    return max(scores, key=scores.get)


def eprint(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)

def main(args):
    with args.data.open("r") as stream:
        header = stream.readline().split(",")
        header = [x.strip().strip('"') for x in header]
        header = [x for x in header if x != "ensg"]
    
    meta = pd.read_csv(args.metadata)
    metadata = meta["sample_id"]

    # Now we can compare the two
    hits = {}
    for item in metadata:
        hits[item] = find_match(item, header)

    # We have the matches
    meta["sample_id"] = list(map(lambda x: hits[x], meta["sample_id"]))

    meta.to_csv(sys.stdout, index = False)
    

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("metadata", type=Path)
    parser.add_argument("data", type=Path)

    args = parser.parse_args()

    main(args)
