
from pathlib import Path

import pandas as pd

# TODO: Write documentation

def main(gtex_meta_path: Path, expression_matrix: Path, tcga_meta_path: Path, output_path : Path) -> None:
    
    # Read a tiny bit of the big file and get the colnames
    reader = pd.read_csv(expression_matrix, chunksize=1, sep = "\t")
    sample_ids = reader.get_chunk(1).columns.to_list()

    # Make a dataframe with the tcga_ids to TSS Code
    tcga_ids = [x for x in sample_ids if x.startswith("TCGA")]
    tcga_ids_to_tss = pd.DataFrame({
        "Sample": tcga_ids,
        "TSS Code": [x.split("-")[1] for x in tcga_ids if x.startswith("TCGA")]
    })

    tcga_meta = pd.read_csv(tcga_meta_path, index_col=0)
    gtex_meta = pd.read_csv(gtex_meta_path, sep = "\t")

    tcga_meta = tcga_meta.merge(tcga_ids_to_tss, on = "TSS Code", how = "outer")
    
    merged_metadata = pd.concat([tcga_meta, gtex_meta])

    merged_metadata.to_csv(output_path, index = False)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("gtex_meta", help="GTEX metadata file (tsv)", type=Path)
    parser.add_argument("expression_matrix", help="Expression matrix to source the TCGA IDs from", type=Path)
    parser.add_argument("tcga_meta", help="TCGA metadata file (csv)", type=Path)
    parser.add_argument("output_path", help="The output path")

    args = parser.parse_args()

    main(gtex_meta_path=args.gtex_meta, expression_matrix=args.expression_matrix, tcga_meta_path=args.tcga_meta, output_path=args.output_path)

