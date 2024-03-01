import sys

BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/{}nnn/{}/matrix/{}_series_matrix.txt.gz"

import shutil
import urllib.request
from contextlib import closing
import io
from gzip import GzipFile

# Everything is an exception. Some series have many matrix files.
# This fixes these ambiguities.
# Fuck me.
EXCEPTIONS = {
    # This has both RNA-seq and miRNA, and we want the metadata for RNA-Seq
    "GSE121842": "GSE121842-GPL20795"
}

def eprint(*args, **kwargs):
    print(*args, **kwargs, file = sys.stderr)

def fetch_ftp(url):
    eprint(f"Fetching {url}...")
    stream = io.BytesIO()
    with closing(urllib.request.urlopen(url)) as r:
        shutil.copyfileobj(r, stream)

    stream.seek(0)
    return GzipFile(fileobj=stream).read()

def main(output_stream, args):
    truncated = args.accession[:-3]
    accession = args.accession if args.accession not in EXCEPTIONS else  EXCEPTIONS[args.accession]
    url = BASE_URL.format(
        truncated,
        args.accession,
        accession
    )

    data = fetch_ftp(url)

    output_stream.writelines(data.decode("utf-8"))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("accession", help="Accession ID of series to download")

    args = parser.parse_args()
    output_stream = sys.stdout

    main(output_stream, args)

