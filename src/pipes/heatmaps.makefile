#? Generate the general heatmap for all transporters

# Number of threads to use to parallelize the ranking process
N_THREADS = 3

# Shorthands
mods = ./src/modules
rexec = Rscript --no-save --no-restore --verbose

## --- Decompress sources
./data/tcga_target_gtex: ./data/in/expression_matrix.csv.gz
	gunzip -cf $< | xsv fmt -d '\t' > $@

./data/selected_metadata: ./data/in/expression_matrix_metadata.tsv.gz
	gunzip -cf $< | xsv fmt -d '\t' > $@

./data/MTPDB.sqlite: ./data/in/MTPDB.sqlite.gz
	gunzip -cf $< > $@

./data/ensg_data.csv : ./data/in/ensg_data.csv.gz
	gunzip -cf $< > $@

## --- Calculate the ranking files from the expression matrix
./data/deas/flag.txt: \
	./data/tcga_target_gtex \
	./data/selected_metadata \
	$(mods)/ranking/select_and_run.py \
	./data/in/dea_queries.json

	mkdir -p $(@D)

	python $(mods)/ranking/select_and_run.py \
		./data/in/dea_queries.json \
		./data/tcga_target_gtex \
		./data/selected_metadata \
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
		--ensg-hugo-data ./data/ensg_data.csv

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

## --- Generate the output GSEA plots from the gsea enrichments
ALL += ./data/out/figures/enrichments/done.flag
./data/out/figures/enrichments/done.flag: \
		./data/out/enrichments/done.flag \
		./data/genesets.json \
		$(mods)/plotting/gsea_plotting_graphs.R \
		$(mods)/plotting/general_heatmap.R

	mkdir -p $(@D)

	$(rexec) $(mods)/plotting/gsea_plotting_graphs.R \
		./data/out/enrichments/ ./data/genesets.json \
		$(@D)

	touch $@

## --- Make the large heatmap with all the results
ALL += ./data/out/figures/enrichments/pancan_heatmap.png
./data/out/figures/enrichments/pancan_heatmap.png: \
		./data/out/enrichments/done.flag \
		$(mods)/plotting/general_heatmap.R \
		./data/genesets_repr.txt

	$(rexec) $(mods)/plotting/general_heatmap.R \
		./data/out/enrichments/ \
		./data/genesets_repr.txt \
		./data/genesets.json \
		$@ \
		--height 15

## This is duplicated from above but i'm not sure how to merge them.
# They are the same plots but with the absolute enrichments.
ALL += ./data/out/figures/absolute_enrichments/done.flag
./data/out/figures/absolute_enrichments/done.flag: \
		./data/out/absolute_enrichments/done.flag \
		./data/genesets.json \
		$(mods)/plotting/gsea_plotting_graphs.R \
		$(mods)/plotting/general_heatmap.R

	mkdir -p $(@D)

	$(rexec) $(mods)/plotting/gsea_plotting_graphs.R \
		./data/out/absolute_enrichments/ ./data/genesets.json \
		$(@D)

	touch $@

ALL +=./data/out/figures/absolute_enrichments/pancan_heatmap.png
./data/out/figures/absolute_enrichments/pancan_heatmap.png: \
		./data/out/absolute_enrichments/done.flag \
		$(mods)/plotting/general_heatmap.R \
		./data/genesets_repr.txt

	$(rexec) $(mods)/plotting/general_heatmap.R \
		./data/out/absolute_enrichments/ \
		./data/genesets_repr.txt \
		./data/genesets.json \
		$@ \
		--height 15

ALL +=./data/out/figures/combined_heatmap.png
./data/out/figures/combined_heatmap.png: \
		./data/out/absolute_enrichments/done.flag \
		./data/out/enrichments/done.flag \
		$(mods)/plotting/fused_general_heatmap.R \
		./data/genesets_repr.txt

	$(rexec) $(mods)/plotting/fused_general_heatmap.R \
		./data/out/enrichments/ \
		./data/out/absolute_enrichments/ \
		./data/genesets_repr.txt \
		./data/genesets.json
		$@ \
		--height 15

PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all

