import sys

BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/{truncated}nnn/{acc}/matrix/{acc}_series_matrix.txt.gz"

import shutil
import urllib.request
from contextlib import closing
import io
from gzip import GzipFile

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
    url = BASE_URL.format(truncated=truncated, acc=args.accession)

    data = fetch_ftp(url)

    output_stream.writelines(data.decode("utf-8"))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("accession", help="Accession ID of series to download")

    args = parser.parse_args()
    output_stream = sys.stdout

    main(output_stream, args)

