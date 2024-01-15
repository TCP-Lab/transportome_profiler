#? Generate the general heatmap for all transporters

# Number of threads to use to parallelize the ranking process
N_THREADS = 3

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
	./data/in/dea_queries.json

	mkdir -p $(@D)

	python $(mods)/ranking/select_and_run.py \
		./data/in/dea_queries.json \
		./data/expression_matrix.csv \
		./data/expression_matrix_metadata.csv \
		$(@D) \
		--cpus $(N_THREADS)

	touch $@

## --- Generate the genesets from the MTPDB
./data/genesets.json ./data/genesets_repr.txt: \
		./data/MTPDB.sqlite \
		$(mods)/make_genesets.py \
		./data/in/basic_gene_lists.json

	python $(mods)/make_genesets.py ./data/MTPDB.sqlite ./data/in/basic_gene_lists.json \
		./data/genesets.json ./data/genesets_repr.txt \
		--prune_direction "topdown" \
		--prune_similarity 0.9 \
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
		--save-plots

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
		--absolute

	touch $@

## --- Make the large heatmap with all the results
ALL +=./data/out/figures/deregulation_heatmap.png
./data/out/figures/deregulation_heatmap.png: \
		./data/out/absolute_enrichments/done.flag \
		./data/out/enrichments/done.flag \
		$(mods)/plotting/plot_large_heatmap.R \
		./data/genesets_repr.txt \
		./data/genesets.json

	$(rexec) $(mods)/plotting/plot_large_heatmap.R \
		./data/out/enrichments/ \
		./data/genesets.json \
		./data/genesets_repr.txt \
		$@ \
		--dots_gsea_results ./data/out/absolute_enrichments/ \
		--height 15

PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all
