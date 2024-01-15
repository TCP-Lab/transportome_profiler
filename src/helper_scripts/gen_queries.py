from pathlib import Path
import json
from copy import copy

TEMPLATE = {
    "case": [
        "<meta>@sample?_study=TCGA&_primary_site=<site>&_sample_type!=[Solid Tissue Normal,Control Analyte]"
    ],
    "control": [
        "<meta>@sample?_study=[TCGA,GTEX]&_primary_site=<site>&_sample_type=[Normal Tissue,Solid Tissue Normal]&_sample_type!=Cell Line"
    ],
}


def main(matches: Path, out: Path):
    with matches.open("r") as stream:
        dict_matches = json.load(stream)

    tcga_to_gtex = {}
    for key, value in dict_matches.items():
        rep_str = ",".join(value)
        if "," in rep_str:
            rep_str = f"[{rep_str}]"

        tcga_to_gtex[key] = {
            "case": [x.replace("<site>", rep_str) for x in TEMPLATE["case"]],
            "control": [x.replace("<site>", rep_str) for x in TEMPLATE["control"]],
        }

    with out.open("w+") as stream:
        json.dump(tcga_to_gtex, stream, indent=4)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "matches_file", type=Path, help="Path to the input matches file"
    )
    parser.add_argument("output_file", type=Path, help="Path to the output file")

    args = parser.parse_args()

    main(args.matches_file, args.output_file)
