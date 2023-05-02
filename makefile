# This makefile can be used to run the whole analysis and rebuild the paper.
# The external input files can be dowloaded by the `retrieve_sources` script.
#
# IMPORTANT :: Before running this, you have to link this project to a data
# 	directory with the `link` utility. Read the README for more information.
#
# The makefile does the following:
# 0. Setup a Python venv to use throughout the makefile;
# 1. Retrieve all the data files through `retrieve_sources`;
# 2. Generates the genesets from the MTP-DB with `geneset_maker` (Python)
# 3. Compiles the `cohen_calculator` executable with `cargo` (Rust);
# 4. Runs the executable with the TCGA/GTEX data plus some static data files
#      (see the `run_dea` README) to generate the Cohen's D matrix;
# 5. Runs the pre-ranked GSEA with the Cohen's D matrix and the genesets (R);
# 6. Generates the output plots with the `gsea_plotting_graphs.R` script (R).
# 7. Compiles the output .tex with `tectonic` from the .tex files in the
#      ./paper directory with the output plots from the analysis.

.PHONY = all venv make_genesets run_dea

data_dir = ./data/
expression_matrix = $(data_dir)/in/tcga_target_gtex
local_mtpdb = $(data_dir)/in/MTPDB.sqlite

all: $(data_dir)/cohen_d_matrix.csv

## --- 0 --- Make the python virtual environment
venv: venv/touchfile

venv/touchfile: requirements.txt
	test -d venv || virtualenv venv
	. venv/bin/activate; pip install -Ur requirements.txt
	touch venv/touchfile

## --- 1 ---
$(expression_matrix) $(local_mtpdb) &: ./src/retrieve_sources
	./src/retrieve_sources $(data_dir)

## --- 2 ---
# This makes way too many (variable) files, so I made the target phony.
make_genesets: \
		$(local_mtpdb) \
		./src/geneset_maker/make_genesets.py \
		./src/geneset_maker/basic_gene_lists.json

	. venv/bin/activate; python \
		./src/geneset_maker/make_genesets.py $(local_mtpdb) ./src/geneset_maker/basic_gene_lists.json \
		--prune_direction "topdown" \
		--verbose

## --- 3 --- Compile the cohen calculator
./src/run_dea/cohen_calculator/target/release/cohen_calculator: ./src/run_dea/cohen_calculator/src/main.rs
	cargo build --release --manifest-path ./src/run_dea/cohen_calculator/Cargo.toml

## --- 4 --- Calculate the cohen's D matrix
$(data_dir)/cohen_d_matrix.csv: $(expression_matrix) \
		./src/run_dea/metadata.csv \
		./src/run_dea/cohen_calculator/target/release/cohen_calculator \

	./src/run_dea/cohen_calculator/target/release/cohen_calculator \
		"$(expression_matrix)" ./src/run_dea/metadata.csv ./src/run_dea/tcga_to_gtex_matches.json \
		"$(data_dir)/cohen_d_matrix.csv"

./src/run_dea/metadata.csv: \
		./src/run_dea/tcga_metadata.csv \
		./src/run_dea/GTEX_phenotype.tsv \
		./src/run_dea/prep_metadata.csv \
		$(expression_matrix)
	
	. venv/bin/activate; python ./src/run_dea/prep_metadata.py \
		"./src/run_dea/GTEX_phenotype.tsv" "$(expression_matrix)" "./src/run_dea/tcga_metadata.csv" \
		"./src/run_dea/metadata.csv"
