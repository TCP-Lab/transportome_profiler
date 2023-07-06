# This makefile can be used to run the whole analysis and rebuild the paper.
# The external input files can be dowloaded by the `retrieve_sources` script.
#
# IMPORTANT :: Before running this, you have to link this project to a data
# 	directory with the `link` utility. Read the README for more information.
#
# The makefile does the following:
# - Setup a Python venv to use throughout the makefile;
# - Retrieve all the data files;
# - Generates the genesets from the MTP-DB with `geneset_maker` (Python)
# - Split the GTEX files to single .csv files based on the metadata columns and
#   run the DEA with Deseq2 from the splitted files.
# - Runs the pre-ranked GSEA with the DEG tables from DeSeq2 and the genesets (R);
# - Generates the output plots with the `gsea_plotting_graphs.R` script (R).

# The phony targets are as follows:
# - all -> run the whole analysis
# - env -> regenerate the python virtual environment
# - clean -> remove the output files
# - restart -> clean and remove the input files and the virtual env
# - scrub -> clean and remove the ./data/ link

data_dir = ./data
expression_matrix = $(data_dir)/in/tcga_target_gtex
local_mtpdb = $(data_dir)/in/MTPDB.sqlite
split_threads = 3

PHONY += clean
# Clean just the intermediary files - this is to re-run the analysis quickly
clean:
	rm -f $(data_dir)/deas/*
	rm -f $(data_dir)/enrichments/*
	rm -rf $(data_dir)/genesets
	find $(data_dir)/out/ -type f -prune ! -name '.*' -exec rm {} +

PHONY += restart
# Clean the input files. This is a harder clean than merely "clean", and to
# re-run the analysis you'll need to redownload the inputs
restart: clean
	rm -f $(expression_matrix)
	rm -f $(local_mtpdb)
	rm -rf ./env

PHONY += scrub
# Clean every file made by make and remove the links made by ./link
# the ending / removes the actual linked directory, the one without the /
# removes the soft link
scrub: restart
	rm -rf $(data_dir)/
	rm -f $(data_dir)

PHONY += env
## >> Make the python virtual environment
env: env/touchfile

PHONY += env/touchfile
env/touchfile: requirements.txt
	python3.11 -m venv env
	. env/bin/activate; pip install -Ur requirements.txt
	touch env/touchfile

## >> Retrieve sources
$(data_dir)/in/tcga_target_gtex:
	mkdir -p $(@D)
	wget -4 --inet4-only -O $@.gz https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_gene_expected_count.gz
	gunzip -f $@.gz

$(data_dir)/in/selected_metadata:
	mkdir -p $(@D)
	wget -4 --inet4-only -O $@.tsv.gz https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGTEX_phenotype.txt.gz
	gunzip -f $@.tsv.gz
	xsv fmt -d '\t' $@.tsv > $@
	rm $@.tsv

# TODO: once we update the releases to output a static tag, change this to
# point at the latest release.
# TODO2: once we publish, change this to point at a static tag again.
$(local_mtpdb):
	mkdir -p $(@D)
	wget -O $@.tar.gz https://github.com/CMA-Lab/MTP-DB/releases/latest/download/MTPDB.sqlite.gz --inet4-only
	tar -xvzf $@.tar.gz -C $(@D)
	mv $(data_dir)/in/MTPDB_v0.23.17-beta.sqlite $@
	rm $@.tar.gz

$(data_dir)/in/ensg_data.csv : ./src/gsea_runner/ensg_data.csv.gz
	mkdir -p $(@D)
	cp ./src/gsea_runner/ensg_data.csv.gz $@.gz
	gunzip $@.gz

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
		--delimiter '\t' \
		--cpus 2

	touch $@

## --- 2 ---
# The fact that it run is detected in the `run_gsea` step by the `all.txt` file
$(data_dir)/genesets/all.txt: \
		env/touchfile \
		$(local_mtpdb) \
		./src/geneset_maker/make_genesets.py \
		./src/geneset_maker/basic_gene_lists.json

	rm -rf $(@D)
	mkdir -p $(@D)

	. env/bin/activate; python \
		./src/geneset_maker/make_genesets.py $(local_mtpdb) ./src/geneset_maker/basic_gene_lists.json \
		$(data_dir)/genesets \
		--prune_direction "topdown" \
		--prune_similarity 0.9 \
		--verbose

## --- 5 --- Run the pre-ranked GSEA
$(data_dir)/out/enrichments/done.flag: \
		$(data_dir)/genesets/all.txt \
		./src/gsea_runner/run_gsea.R \
		$(data_dir)/deas/flag.txt \
		$(data_dir)/in/ensg_data.csv

	mkdir -p $(@D)

	Rscript ./src/gsea_runner/run_gsea.R \
		"$(data_dir)/deas/" "$(data_dir)/genesets" "$(@D)" \
		--ensg-hugo-data $(data_dir)/in/ensg_data.csv

	touch $@

$(data_dir)/out/absolute_enrichments/done.flag: \
		$(data_dir)/genesets/all.txt \
		./src/gsea_runner/run_gsea.R \
		$(data_dir)/deas/flag.txt \
		$(data_dir)/in/ensg_data.csv

	mkdir -p $(@D)

	Rscript ./src/gsea_runner/run_gsea.R \
		"$(data_dir)/deas/" "$(data_dir)/genesets" "$(@D)" \
		--ensg-hugo-data $(data_dir)/in/ensg_data.csv \
		--absolute

	touch $@

## --- 6 --- Generate plots
/tmp/genesets_tree/tree.txt: $(data_dir)/genesets/all.txt
	mkdir -p $(@D)
	tree $(data_dir)/genesets -df | head -n -2 > $@

ALL += $(data_dir)/out/figures/enrichments/done.flag
$(data_dir)/out/figures/enrichments/done.flag: \
		$(data_dir)/out/enrichments/done.flag \
		$(data_dir)/genesets/all.txt \
		./src/plotting/gsea_plotting_graphs.R \
		./src/plotting/general_heatmap.R

	mkdir -p $(@D)

	Rscript ./src/plotting/gsea_plotting_graphs.R \
		$(data_dir)/out/enrichments/ \
		$(@D)

	touch $@

ALL += $(data_dir)/out/figures/enrichments/pancan_heatmap.png
$(data_dir)/out/figures/enrichments/pancan_heatmap.png: \
		$(data_dir)/out/enrichments/done.flag \
		./src/plotting/general_heatmap.R \
		/tmp/genesets_tree/tree.txt

	Rscript ./src/plotting/general_heatmap.R \
		$(data_dir)/out/enrichments/ \
		/tmp/genesets_tree/tree.txt \
		$@ \
		--height 15

## This is duplicated from above but i'm not sure how to merge them.
ALL += $(data_dir)/out/figures/absolute_enrichments/done.flag
$(data_dir)/out/figures/absolute_enrichments/done.flag: \
		$(data_dir)/out/absolute_enrichments/done.flag \
		$(data_dir)/genesets/all.txt \
		./src/plotting/gsea_plotting_graphs.R \
		./src/plotting/general_heatmap.R

	mkdir -p $(@D)

	Rscript ./src/plotting/gsea_plotting_graphs.R \
		$(data_dir)/out/absolute_enrichments/ \
		$(@D)

	touch $@

ALL +=$(data_dir)/out/figures/absolute_enrichments/pancan_heatmap.png
$(data_dir)/out/figures/absolute_enrichments/pancan_heatmap.png: \
		$(data_dir)/out/absolute_enrichments/done.flag \
		./src/plotting/general_heatmap.R \
		/tmp/genesets_tree/tree.txt

	Rscript ./src/plotting/general_heatmap.R \
		$(data_dir)/out/absolute_enrichments/ \
		/tmp/genesets_tree/tree.txt \
		$@ \
		--height 15

ALL +=$(data_dir)/out/figures/combined_heatmap.png
$(data_dir)/out/figures/combined_heatmap.png: \
		$(data_dir)/out/absolute_enrichments/done.flag \
		$(data_dir)/out/enrichments/done.flag \
		./src/plotting/fused_general_heatmap.R \
		/tmp/genesets_tree/tree.txt

	Rscript --no-save --no-restore --verbose ./src/plotting/fused_general_heatmap.R \
		$(data_dir)/out/enrichments/ \
		$(data_dir)/out/absolute_enrichments/ \
		/tmp/genesets_tree/tree.txt \
		$@ \
		--height 15


build_paper_command = cd ./paper/src/ && latexmk -lualatex -f -quiet -gg -synctex=1 -interaction=nonstopmode -file-line-error main.tex

## --- 7 --- Build the paper
$(data_dir)/out/paper.pdf: \
		$(data_dir)/out/figures/enrichments/done.flag \
		$(data_dir)/out/figures/enrichments/pancan_heatmap.pdf \
		./paper/src

	$(build_paper_command)

	mv ./paper/src/main.pdf $@

paper: $(data_dir)/out/paper.pdf

thin_paper:
	$(build_paper_command)


PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all
