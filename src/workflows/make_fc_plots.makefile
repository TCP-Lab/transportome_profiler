rexec = Rscript --no-save --no-restore --verbose
mods = ./src/modules
N_THREADS = 1

method = "norm_fold_change"

./data/%: ./data/in/%.gz
	gunzip -cfv $< > $@

%.csv: %.tsv
	xsv fmt -d '\t' $< > $@


## --- Calculate the ranking files from the expression matrix
./data/fc_deas/flag.txt: \
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
		--method $(method)

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

## --- Make the fold change plot
ALL += data/out/figures/fold_change_heatmap.png
./data/out/figures/fold_change_heatmap.png: \
	data/fc_deas/flag.txt \
	$(mods)/plotting/plot_fc.R\
	data/genesets.json
	
	cat data/genesets.json | \
		jq -r '.. | select(.["name"]? == "whole_transportome") | .data[]' > \
		/tmp/filter_genes.txt
	
	mkdir -p $(@D)
	$(rexec) $(mods)/plotting/plot_fc.R \
		data/fc_deas $@ --png \
		--filter_genes /tmp/filter_genes.txt \
		--extra_title $(method) \
		--height 12

PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all
