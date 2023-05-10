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

data_dir = ./data
expression_matrix = $(data_dir)/in/tcga_target_gtex
local_mtpdb = $(data_dir)/in/MTPDB.sqlite

all: $(data_dir)/out/figures/enrichments/done.flag

# Clean just the intermediary files - this is to re-run the analysis quickly
clean:
	rm -f $(data_dir)/deas/*
	rm -f $(data_dir)/enrichments/*
	rm -rf $(data_dir)/genesets
	find $(data_dir)/out/ -type f -prune ! -name '.*' -exec rm {} +
	rm -f ./env/touchfile
	rm -rf ./env

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
env: env/touchfile

env/touchfile: requirements.txt
	python3.11 -m venv env
	. env/bin/activate; pip install -Ur requirements.txt
	touch env/touchfile

## --- 1 --- Retrieve sources
# Retrieve all the expression files from xena
$(data_dir)/in/tcga_target_gtex:
	mkdir -p $(@D)
	wget -4 --inet4-only -O $@.gz https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_gene_expected_count.gz 
	gunzip -f $@.gz

$(data_dir)/in/selected_metadata:
	mkdir -p $(@D)
	wget -4 --inet4-only -O $@.tsv.gz https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz
	gunzip -f $@.tsv.gz
	xsv slice $@.tsv > $@
	rm $@.tsv

$(local_mtpdb):
	mkdir -p $(@D)
	wget -O $@.tar.gz https://github.com/CMA-Lab/MTP-DB/releases/download/0.23.17-beta/MTPDB_v0.23.17-beta.sqlite.tar.gz --inet4-only
	tar -xvzf $@.tar.gz -C $(@D)
	mv $(data_dir)/in/MTPDB_v0.23.17-beta.sqlite $@
	rm $@.tar.gz

## --- --- Calculate the DEA files
$(data_dir)/deas/flag.txt: \
	env/touchfile \
	$(data_dir)/in/tcga_target_gtex \
	$(data_dir)/in/selected_metadata \
	./src/run_dea/select_and_run.py \
	./src/run_dea/run_deseq.R \
	./src/run_dea/tcga_gtex_queries.json

	mkdir -p $(@D)

	. env/bin/activate; python \
		./src/run_dea/select_and_run.py \
		./src/run_dea/tcga_gtex_queries.json \
		$(data_dir)/in/tcga_target_gtex \
		$(data_dir)/in/selected_metadata \
		$(@D) \
		./src/run_dea/run_deseq.R \
		--delimiter '\t'

	touch $@

## --- 2 ---
# The fact that it run is detected in the `run_gsea` step by the `all.txt` file
$(data_dir)/genesets/all.txt: \
		env/touchfile \
		$(local_mtpdb) \
		./src/geneset_maker/make_genesets.py \
		./src/geneset_maker/basic_gene_lists.json

	mkdir -p $(@D)

	. env/bin/activate; python \
		./src/geneset_maker/make_genesets.py $(local_mtpdb) ./src/geneset_maker/basic_gene_lists.json \
		$(data_dir)/genesets \
		--prune_direction "topdown" \
		--verbose

## --- 5 --- Run the pre-ranked GSEA
$(data_dir)/out/enrichments/done.flag: \
		$(data_dir)/genesets/all.txt \
		./src/gsea_runner/run_gsea.R \
		$(data_dir)/deas/flag.txt

	mkdir -p $(@D)

	Rscript ./src/gsea_runner/run_gsea.R \
		"$(data_dir)/deas/" "$(data_dir)/genesets" "$(data_dir)/out/enrichments"

	touch $(data_dir)/out/enrichments/done.flag

## --- 6 --- Generate plots
$(data_dir)/out/figures/enrichments/done.flag: \
		$(data_dir)/out/enrichments/done.flag \
		$(data_dir)/genesets/all.txt \
		./src/gsea_runner/gsea_plotting_graphs.R

	mkdir -p $(@D)

	Rscript ./src/gsea_runner/gsea_plotting_graphs.R \
		$(data_dir)/out/enrichments/ \
		$(data_dir)/out/figures/enrichments
	
	touch "$(data_dir)/out/figures/enrichments/done.flag"

