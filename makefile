# This makefile can be used to run the whole analysis and rebuild the paper.
# The external input files can be dowloaded by the `retrieve_sources` script.
#
# IMPORTANT :: Before running this, you have to link this project to a data
# 	directory with the `link` utility. Read the README for more information.
#
# The makefile does the following:
# 0. Setup a Python venv to use throughout the makefile;
# 1. Retrieve all the data files;
# 2. Generates the genesets from the MTP-DB with `geneset_maker` (Python)
# 3. Split the GTEX files to single .csv files based on the metadata columns
# 4. Run the DEA with Deseq2 from the GTEX expression + GTEX metadata +
#    a TCGA file +   
# 5. Runs the pre-ranked GSEA with the Cohen's D matrix and the genesets (R);
# 6. Generates the output plots with the `gsea_plotting_graphs.R` script (R).
# 7. Compiles the output .tex with `tectonic` from the .tex files in the
#      ./paper directory with the output plots from the analysis.

.PHONY = all venv make_genesets run_dea retrieve_sources clean clean_input

TCGA_ID = LAML ACC BLCA LGG BRCA CESC CHOL LCML COAD  ESC GBM HNSC KICH KIRC KIRP LIHC LUAD LUSC DLBC MESO OV PAAD PCPG PRAD READ SARC STAD TGCT THYM UCS UCEC UVM

data_dir = ./data
expression_matrix = $(data_dir)/in/tcga_target_gtex
local_mtpdb = $(data_dir)/in/MTPDB.sqlite

all: "$(data_dir)/out/paper/main.pdf"

# Clean just the intermediary files - this is to re-run the analysis quickly
clean:
	rm -f $(data_dir)/cohen_d_matrix.csv
	rm -f $(data_dir)/metadata.csv
	find $(data_dir)/out/ -type f -prune ! -name '.*' -exec rm {} +
	rm -f ./venv/touchfile
	rm -rf ./venv

# Clean the input files. This is a harder clean than merely "clean", and to
# re-run the analysis you'll need to redownload the inputs
restart: clean
	rm -f $(expression_matrix)
	rm -f $(local_mtpdb)

# Clean every file made by make and remove the links made by ./link
scrub: restart
	rm -f $(data_dir)
	rm -f ./paper/resources/images/generated

## --- 0 --- Make the python virtual environment
venv: venv/touchfile

venv/touchfile: requirements.txt
	python3.11 -m venv venv
	. venv/bin/activate; pip install -Ur requirements.txt
	touch venv/touchfile

## --- 1 --- Retrieve sources
# Retrieve all the TCGA expression files from Xena - the $(@F) variables expands
# to the filename of the target (TCGA-xxxx in this case)
$(data_dir)/in/tcga_raw/TCGA-$(TCGA_ID):
	mkdir -p $(@D)
	wget -4 --inet4-only -O $@.gz https://gdc-hub.s3.us-east-1.amazonaws.com/download/$(@F).htseq_counts.tsv.gz
	gunzip $@.gz

$(data_dir)/in/gtex_raw_counts:
	mkdir -p $(@D)
	wget -4 --inet4-only -O $@.gz https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/gtex_gene_expected_count.gz 
	gunzip $@.gz

$(local_mtpdb):
	mkdir -p $(@D)
	wget -O $@.tar.gz https://github.com/CMA-Lab/MTP-DB/releases/download/0.23.17-beta/MTPDB_v0.23.17-beta.sqlite.tar.gz --inet4-only
	tar -xvzf $@.tar.gz -C $(@D)
	rm $@.tar.gz

## --- 4 --- Calculate the cohen's D matrix
$(data_dir)/cohen_d_matrix.csv: $(expression_matrix) \
		$(data_dir)/metadata.csv \
		./src/run_dea/cohen_calculator/target/release/cohen_calculator \

	./src/run_dea/cohen_calculator/target/release/cohen_calculator \
		"$(expression_matrix)" $(data_dir)/metadata.csv ./src/run_dea/tcga_to_gtex_matches.json \
		"$(data_dir)/cohen_d_matrix.csv"

$(data_dir)/metadata.csv: \
		venv/touchfile \
		./src/run_dea/tcga_metadata.csv \
		./src/run_dea/GTEX_phenotype.tsv \
		./src/run_dea/prep_metadata.py \
		$(expression_matrix)
	
	. venv/bin/activate; python ./src/run_dea/prep_metadata.py \
		"./src/run_dea/GTEX_phenotype.tsv" "$(expression_matrix)" "./src/run_dea/tcga_metadata.csv" \
		"$(data_dir)/metadata.csv"

## --- 2 ---
# The fact that it run is detected in the `run_gsea` step by the `all.txt` file
$(data_dir)/genesets/all.txt: \
		venv/touchfile \
		$(local_mtpdb) \
		./src/geneset_maker/make_genesets.py \
		./src/geneset_maker/basic_gene_lists.json

	mkdir -p $(@D)

	. venv/bin/activate; python \
		./src/geneset_maker/make_genesets.py $(local_mtpdb) ./src/geneset_maker/basic_gene_lists.json \
		$(data_dir)/genesets \
		--prune_direction "topdown" \
		--verbose

## --- 5 --- Run the pre-ranked GSEA
"$(data_dir)/out/enrichments/done.flag": \
		$(data_dir)/genesets/all.txt \
		./src/gsea_runner/run_gsea.R \
		$(data_dir)/cohen_d_matrix.csv

	mkdir -p $(@D)

	Rscript ./src/gsea_runner/run_gsea.R \
		"$(data_dir)/cohen_d_matrix.csv" "$(data_dir)/genesets" "$(data_dir)/out/enrichments"

	touch "$(data_dir)/out/enrichments/done.flag"

## --- 6 --- Generate plots
"$(data_dir)/out/figures/enrichments/done.flag": \
		"$(data_dir)/out/enrichments/done.flag" \
		$(data_dir)/genesets/all.txt \
		./src/gsea_runner/gsea_plotting_graphs.R

	mkdir -p $(@D)

	Rscript ./src/gsea_runner/gsea_plotting_graphs.R \
		$(data_dir)/out/enrichments/ \
		$(data_dir)/out/figures/enrichments
	
	touch "$(data_dir)/out/figures/enrichments/done.flag"

## --- 7 --- Generate paper
"$(data_dir)/out/paper/main.pdf": \
		"$(data_dir)/out/figures/enrichments/done.flag" \
		./paper

	mkdir -p $(@D)
	tectonic ./paper/main.tex -o $(data_dir)/out/paper/
