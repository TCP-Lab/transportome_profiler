"""Retrieve a series file and take out specific columns from it"""
import sys
import csv

def series_to_dict(text: str, strip_prefix: bool = True):
    keys = {}
    prefix = "!Series_" if strip_prefix else "!"
    for line in text.split("\n"):
        if not line.startswith("!"):
            continue
        line = line.strip()
        line = line[len(prefix):]
        pieces = line.split("\t")
        if len(pieces) < 2:
            continue
        id = pieces[0]
        pieces = pieces[1:] if len(pieces) > 2 else pieces[1]
        keys[id] = [x.strip('"') for x in pieces]

    return keys


def main(input_stream, output_stream, args):
    file = input_stream.read()
    data = series_to_dict(file, not args.keep_sample_prefix)
    
    if args.id_col not in data:
        raise ValueError(f"ID col {args.id_col} not in input file")

    id = data[args.id_col]
    metadata_cols = []
    for col in args.meta_cols:
        if col not in data:
            raise ValueError(f"Meta col {col} not in input file")
        if len(data[col]) != len(id):
            raise ValueError(f"Meta col {col} has bad length: expected {len(id)}, got {len(data[col])}")

        metadata_cols.append([col] + data[col])
    
    # Join together the values
    id = ["sample_id"] + id

    writer = csv.writer(output_stream)

    for tup in zip(id, *metadata_cols):
        writer.writerow(tup)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("id_col", help="Name of column with the sample IDs")
    parser.add_argument("--keep-sample-prefix", help="Keep the 'Sample_' prefix?")
    parser.add_argument("meta_cols", help="Column(s) to extract", nargs="+")

    args = parser.parse_args()

    input_stream = sys.stdin
    output_stream = sys.stdout

    main(input_stream= input_stream, output_stream=output_stream, args = args)

