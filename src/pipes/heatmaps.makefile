#? Generate the general heatmap for all transporters
#?
#? This takes in the whole expression matrix, splits it, runs ranking
#? against the splits, runs GSEA against the ranks with genesets from
#? the MTP-DB and produces the large pancancer heatmaps.


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
		./data/in/config/gene_lists/basic.json

	python $(mods)/make_genesets.py ./data/MTPDB.sqlite ./data/in/config/gene_lists/basic.json \
		./data/genesets.json ./data/genesets_repr.txt \
		--prune_direction $(PRUNE_DIRECTION) \
		--prune_similarity $(PRUNE_SIMILARITY) \
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
		--extra_title "alpha 0.20, metric $(RANK_METHOD)" \
		--height 15 \
		$(_heatmap_plot_flags)

PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all
