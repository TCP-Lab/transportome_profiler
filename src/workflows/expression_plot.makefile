#? Generate the plot with the expressed/not expressed information for genesets

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

ALL += ./data/out/figures/expression_means.png
./data/out/figures/expression_means.png: \
	./data/expression_means.csv \
	./src/modules/plotting/plot_expression_means.R \
	./data/genesets.json \
	./data/genesets_repr.txt

	$(rexec) ./src/modules/plotting/plot_expression_means.R $< \
		./data/genesets.json ./data/genesets_repr.txt $@ \
		--res 400 --height 15 --expression_threshold 0 \
		--extra_title "TCGA + GTEX data"

ALL += ./data/out/figures/expression_means_TCGA_only.png
./data/out/figures/expression_means_TCGA_only.png: \
	./data/expression_means_TCGA.csv \
	./src/modules/plotting/plot_expression_means.R \
	./data/genesets.json \
	./data/genesets_repr.txt

	$(rexec) ./src/modules/plotting/plot_expression_means.R $< \
		./data/genesets.json ./data/genesets_repr.txt $@ \
		--res 400 --height 15 --expression_threshold 0 \
		--extra_title "TCGA only"

PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all
