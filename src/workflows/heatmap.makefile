#? Generate the general heatmap for all transporters
#?
#? This takes in the whole expression matrix, splits it, runs ranking
#? against the splits, runs GSEA against the ranks with genesets from
#? the MTP-DB and produces the large pancancer heatmaps.

SHELL = /bin/bash

.ONESHELL:

OPTS=./data/in/config/heatmaps_runtime_options.json
# Number of threads to use to parallelize the ranking process
N_THREADS ?= $(shell cat $(OPTS) | jq -r '.threads')
# Method to use
RANK_METHOD ?= $(shell cat $(OPTS) | jq -r '.rank_method')
PRUNE_SIMILARITY ?= $(shell cat $(OPTS) | jq -r '.prune_similarity')
PRUNE_DIRECTION ?= $(shell cat $(OPTS) | jq -r '.prune_direction')
ALPHA_THRESHOLD ?= $(shell cat $(OPTS) | jq -r '.alpha_threshold')

# Option switches
SAVE_EXTRA_PLOTS ?= $(shell cat $(OPTS) | jq -r '.save_extra_plots')
ifeq ($(SAVE_EXTRA_PLOTS), true)
_gsea_runtime_flags += "--save-plots"
endif

RUN_UNWEIGHTED ?= $(shell cat $(OPTS) | jq -r '.run_unweighted')
ifeq ($(RUN_UNWEIGHTED), true)
_gsea_runtime_flags += "--unweighted"
endif

CLUSTER_COLS ?= $(shell cat $(OPTS) | jq -r '.cluster_heatmap_cols')
ifeq ($(CLUSTER_COLS), false)
_heatmap_plot_flags += "--no_cluster"
endif

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

## --- Calculate the ranking files from the expression matrix
./data/deas/flag.txt: \
	./data/expression_matrix.csv \
	./data/expression_matrix_metadata.csv \
	$(mods)/ranking/select_and_run.py \
	./data/in/config/DEA_queries/dea_queries.json

	mkdir -p $(@D)

	python $(mods)/ranking/select_and_run.py \
		./data/in/config/DEA_queries/dea_queries.json \
		./data/expression_matrix.csv \
		./data/expression_matrix_metadata.csv \
		$(@D) \
		--cpus $(N_THREADS) \
		--method $(RANK_METHOD)

	touch $@

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

## --- Run the pre-ranked GSEA
./data/out/enrichments/done.flag: \
		./data/genesets.json \
		$(mods)/run_gsea.R \
		./data/deas/flag.txt \
		./data/ensg_data.csv

	mkdir -p $(@D)

	$(rexec) $(mods)/run_gsea.R \
		"./data/deas/" "./data/genesets.json" "$(@D)" \
		--ensg-hugo-data ./data/ensg_data.csv \
		$(_gsea_runtime_flags)

	touch $@

# --- Run the pre-ranked GSEA, but keep absolute enrichments only
./data/out/absolute_enrichments/done.flag: \
		./data/genesets.json \
		$(mods)/run_gsea.R \
		./data/deas/flag.txt \
		./data/ensg_data.csv

	mkdir -p $(@D)

	$(rexec) $(mods)/run_gsea.R \
		"./data/deas/" "./data/genesets.json" "$(@D)" \
		--ensg-hugo-data ./data/ensg_data.csv \
		--absolute \
		$(_gsea_runtime_flags)

	touch $@

## --- Make the large heatmap with all the results
ALL +=./data/out/figures/deregulation_heatmap.png
./data/out/figures/deregulation_heatmap.png: \
		./data/out/absolute_enrichments/done.flag \
		./data/out/enrichments/done.flag \
		$(mods)/plotting/plot_large_heatmap.R \
		./data/genesets_repr.txt \
		./data/genesets.json

	mkdir -p $(@D)

	$(rexec) $(mods)/plotting/plot_large_heatmap.R \
		./data/out/enrichments/ \
		./data/genesets.json \
		./data/genesets_repr.txt \
		$@ \
		--dots_gsea_results ./data/out/absolute_enrichments/ \
		--alpha 0.20 \
		--extra_title "alpha 0.20, $(RANK_METHOD)" \
		--height 15 \
		$(_heatmap_plot_flags)

## --- Make the output summary plots
data/filter_genes.txt: data/genesets.json
	cat $< | jq -r '.[] | select(.name == "whole_transportome").data | @csv' > $@

data/merged_deas.csv: ./data/deas/flag.txt
	mkdir -p ${@D}
	# This renames the header slot 'ranking' to the name of the file
	# This magical $ thing is to replace the .csv to nothing
	# the pattern is {variable/pattern/replacement}
	for name in ./data/deas/*_deseq.csv; do \
		replacement=$$(basename "$${name/_deseq.csv}"); \
		sed "1s/ranking/$${replacement}/" "$${name}" > $${name/.csv/}.renames.csv; \
	done
	# Move one random file to be the "base" that we can reduce into
	mv $$(find "./data/deas/" -name "*.renames.csv" -print -quit) ./data/base
	# The xsv select is there to remove the second "sample" col that "join" retains.
	for item in ./data/deas/*.renames.csv; do \
		xsv join sample ./data/base sample $${item} --full | xsv select '!sample[1]' > ./data/tmp;\
		mv ./data/tmp ./data/base; \
	done
	mv ./data/base $@
	# Since the .renames files may conflict with the old pipeline, I
	# just delete them here instead of changing the previous steps.
	# Sorry!
	rm ./data/deas/*.renames.csv

## ---- Shared dysregulation plots ---

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

ALL +=./data/suppressed_merged_deas.csv
./data/suppressed_merged_deas.csv: \
	./data/merged_deas.csv \
	./data/expression_means.csv \
	${mods}/suppress_not_expressed.R
	${rexec} ${mods}/suppress_not_expressed.R ./data/merged_deas.csv ./data/expression_means.csv $@ \
		--expression-threshold 0

## "Top dysregulation" plots.
# The names here need to follow a pattern, since I deconstruct the file names
# to set variables:
# top_disregulation_thr_{{threshold}}_set_{{geneset}}.png
ALL +=./data/out/figures/top_disregulation_thr_1_set_whole_transportome.png
ALL +=./data/out/figures/top_disregulation_thr_1.5_set_whole_transportome.png
ALL +=./data/out/figures/top_disregulation_thr_2_set_whole_transportome.png
ALL +=./data/out/figures/top_disregulation_thr_1_set_channels.png
ALL +=./data/out/figures/top_disregulation_thr_1.5_set_channels.png
ALL +=./data/out/figures/top_disregulation_thr_2_set_channels.png
ALL +=./data/out/figures/top_disregulation_thr_1_set_transporters.png
ALL +=./data/out/figures/top_disregulation_thr_1.5_set_transporters.png
ALL +=./data/out/figures/top_disregulation_thr_2_set_transporters.png
./data/out/figures/top_disregulation_thr_%.png: \
		./data/suppressed_merged_deas.csv \
		./data/genesets.json \
		${mods}/plotting/plot_shared_dysregulation.R \
		./data/ensg_data.csv
	THR=$$(echo '$@' | rg '.*top_disregulation_thr_([0-9,.]+)_set_(.+?).png' -or '$$1')
	SET=$$(echo '$@' | rg '.*top_disregulation_thr_([0-9,.]+)_set_(.+?).png' -or '$$2')
	echo "$${THR}"
	echo "$${SET}"
	mkdir -p ${@D}
	cat ./data/genesets.json | jq -r ".[] | select(.name == \"$${SET}\").data | @csv" > /tmp/selected_genes.csv
	${rexec} ${mods}/plotting/plot_shared_dysregulation.R $@ $< data/ensg_data.csv \
		--selected_genes /tmp/selected_genes.csv --static_threshold $${THR} --renames data/in/config/tcga_renames.json \
		--png --res 300

ALL +=./data/out/figures/upset.png
./data/out/figures/upset.png: \
		./data/merged_deas.csv \
		./data/filter_genes.txt \
		${mods}/plotting/plot_general_upset.R \
		./data/ensg_data.csv
	mkdir -p ${@D}
	${rexec} ${mods}/plotting/plot_general_upset.R $@ $< \
		--selected_genes data/filter_genes.txt --png --res 400

ALL +=./data/out/figures/correlation.png
./data/out/figures/correlation.png: \
		./data/genesets.json \
		./data/genesets_repr.txt
	mkdir -p ${@D}
	${rexec} ${mods}/plotting/plot_genesets_correlation.R \
		./data/genesets.json ./data/genesets_repr.txt \
		$@ --png --res 400

PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all
