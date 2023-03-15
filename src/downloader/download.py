
from pathlib import Path
import requests as rq
from tqdm import tqdm

def main(out_dir: Path):
    out_dir = out_dir.expanduser()

    with Path("/home/hedmad/Files/repos/TCGA_dea/src/tcga_study_ids.txt").open("r") as stream:
        lines = stream.readlines()
    
    lines = [x.strip() for x in lines]

    server = "https://gdc-hub.s3.us-east-1.amazonaws.com/download/"
    base_file = "{}.htseq_counts.tsv.gz" 
    for study in tqdm(lines):
        file = base_file.format(study)
        url = server + file

        bytes = rq.get(url=url).content

        with (out_dir / file).open("wb+") as stream:
            stream.write(bytes)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("output_dir", type=Path, help="Where the files will be saved")

    args = parser.parse_args()

    main(args.output_dir)
