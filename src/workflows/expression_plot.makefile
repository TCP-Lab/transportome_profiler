#? Generate the plot with the expressed/not expressed information for genesets

SHELL = /bin/bash


OPTS=./data/in/config/heatmaps_runtime_options.json

PRUNE_SIMILARITY ?= $(shell cat $(OPTS) | jq -r '.prune_similarity')
PRUNE_DIRECTION ?= $(shell cat $(OPTS) | jq -r '.prune_direction')

# Shorthands
mods = ./src/modules
rexec = Rscript --no-save --no-restore --verbose

## -- Generic rules
## --- Decompress sources
./data/%: ./data/in/%.gz
	gunzip -cfv $< > $@

%.csv: %.tsv
	xsv fmt -d '\t' $< > $@

# If we have the input data as-is in the data/in folder, but we need it in
# data/ we can just copy it.
./data/%: ./data/in/%
	cp $< $@

## --- Generate the genesets from the MTPDB
./data/genesets.json ./data/genesets_repr.txt &: \
		./data/MTPDB.sqlite \
		$(mods)/make_genesets.py \
		./data/in/config/gene_lists/no_experimental_ions.json

	# Min recurse set size is 0 otherwise pruning is aleatory!!
	python $(mods)/make_genesets.py ./data/MTPDB.sqlite ./data/in/config/gene_lists/no_experimental_ions.json \
		./data/genesets.json ./data/genesets_repr.txt \
		--prune_direction $(PRUNE_DIRECTION) \
		--prune_similarity $(PRUNE_SIMILARITY) \
		--min_recurse_set_size 0 \
		--verbose

## --- Calculate the expressed/not expressed matrix based on tumor type (TCGA/GTEX)
./data/expression_means.csv: \
	./data/expression_matrix_tpm.csv \
	./data/expression_matrix_metadata.csv \
	$(mods)/calc_expression_means.py \
	./data/in/config/DEA_queries/dea_queries.json

	mkdir -p $(@D)

	python $(mods)/calc_expression_means.py \
		./data/in/config/DEA_queries/dea_queries.json \
		./data/expression_matrix_tpm.csv \
		./data/expression_matrix_metadata.csv \
		$(@)

## --- Calculate the expressed/not expressed matrix based on tumor type (TCGA only)
./data/expression_means_TCGA.csv: \
	./data/expression_matrix_tpm.csv \
	./data/expression_matrix_metadata.csv \
	$(mods)/calc_expression_means.py \
	./data/in/config/DEA_queries/dea_queries.json

	mkdir -p $(@D)

	python $(mods)/calc_expression_means.py \
		./data/in/config/DEA_queries/dea_queries.json \
		./data/expression_matrix_tpm.csv \
		./data/expression_matrix_metadata.csv \
		$(@) \
		--case-only

## --- Calculate the expressed/not expressed matrix based on tumor type (GTEX only)
./data/expression_means_GTEX.csv: \
	./data/expression_matrix_tpm.csv \
	./data/expression_matrix_metadata.csv \
	$(mods)/calc_expression_means.py \
	./data/in/config/DEA_queries/dea_queries.json

	mkdir -p $(@D)

	python $(mods)/calc_expression_means.py \
		./data/in/config/DEA_queries/dea_queries.json \
		./data/expression_matrix_tpm.csv \
		./data/expression_matrix_metadata.csv \
		$(@) \
		--control-only

ALL += ./data/out/figures/expression_means.png
./data/out/figures/expression_means.png: \
	./data/expression_means.csv \
	./src/modules/plotting/plot_expression_means.R \
	./data/genesets.json \
	./data/genesets_repr.txt

	mkdir -p $(@D)

	$(rexec) ./src/modules/plotting/plot_expression_means.R $< \
		./data/genesets.json ./data/genesets_repr.txt $@ \
		--res 400 --height 15 --expression_threshold 1 \
		--renames ./data/in/config/tcga_renames_sample.json \
		--extra_title "TCGA + GTEX data"

ALL += ./data/out/figures/expression_means_TCGA_only.png
./data/out/figures/expression_means_TCGA_only.png: \
	./data/expression_means_TCGA.csv \
	./src/modules/plotting/plot_expression_means.R \
	./data/genesets.json \
	./data/genesets_repr.txt

	$(rexec) ./src/modules/plotting/plot_expression_means.R $< \
		./data/genesets.json ./data/genesets_repr.txt $@ \
		--res 400 --height 15 --expression_threshold 1 \
		--renames ./data/in/config/tcga_renames_sample.json \
		--extra_title "TCGA only"

ALL += ./data/out/figures/expression_means_GTEX_only.png
./data/out/figures/expression_means_GTEX_only.png: \
	./data/expression_means_GTEX.csv \
	./src/modules/plotting/plot_expression_means.R \
	./data/genesets.json \
	./data/genesets_repr.txt

	$(rexec) ./src/modules/plotting/plot_expression_means.R $< \
		./data/genesets.json ./data/genesets_repr.txt $@ \
		--res 400 --height 15 --expression_threshold 1 \
		--renames ./data/in/config/tcga_renames_sample.json \
		--extra_title "GTEX only"

ALL += ./data/out/avg_expression.csv
./data/out/avg_expression.csv: \
	./data/expression_means_GTEX.csv \
	./src/modules/get_avg_expression.R \
	./data/ensg_data.csv

	$(rexec) ./src/modules/get_avg_expression.R $< \
		./data/ensg_data.csv \
		$@ 0

## --- Calculate the UPSET plots
./data/out/figures/%_expression_upset.png: \
	./data/%_single_geneset.json \
	./data/expression_means.csv \
	./src/modules/plotting/plot_expression_upset.R

	$(rexec) ./src/modules/plotting/plot_expression_upset.R \
		./data/expression_means.csv $< $@ \
		--expression_threshold 1 \
		--renames ./data/in/config/tcga_renames_sample.json \
		--extra_title "Upset of expressed genes - ${*}" \
		--png

./data/%_single_geneset.json: ./data/genesets.json
	jq 'to_entries[] | select(.value.name == "$*") | .value.data' $< > $@

ALL += ./data/out/figures/whole_transportome_expression_upset.png \
	   ./data/out/figures/channels_expression_upset.png \
	   ./data/out/figures/transporters_expression_upset.png

PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all
