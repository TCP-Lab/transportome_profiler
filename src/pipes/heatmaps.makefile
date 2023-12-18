#? Generate the general heatmap for all transporters
#?
#? This makefile can be used to run the whole analysis and rebuild the paper.
#? The external input files can be dowloaded by the `retrieve_sources` script.

local_mtpdb = ./data/MTPDB.sqlite
split_threads = 3
mods = ./src/modules

PHONY += env
## >> Make the python virtual environment
env: env/touchfile

PHONY += env/touchfile
env/touchfile: ./src/requirements.txt
	python3.11 -m venv env
	. env/bin/activate; pip install -Ur ./src/requirements.txt
	touch env/touchfile

## >> Retrieve sources
./data/tcga_target_gtex: ./data/in/expression_matrix.csv.gz
	gunzip -cf $< > $@

./data/selected_metadata: ./data/in/expression_matrix_metadata.tsv.gz
	gunzip -cf $< | xsv fmt -d '\t' > $@

# TODO: once we update the releases to output a static tag, change this to
# point at the latest release.
# TODO2: once we publish, change this to point at a static tag again.
$(local_mtpdb): ./data/in/MTPDB.sqlite.gz
	gunzip -cf $< > $@

./data/ensg_data.csv : ./data/in/ensg_data.csv.gz
	gunzip -cf $< > $@

## --- --- Calculate the DEA files
./data/deas/flag.txt: \
	env/touchfile \
	./data/tcga_target_gtex \
	./data/selected_metadata \
	$(mods)/run_dea/select_and_run.py \
	$(mods)/run_dea/run_deseq.R \
	./data/in/dea_queries.json

	mkdir -p $(@D)

	. env/bin/activate; python \
		$(mods)/run_dea/select_and_run.py \
		./data/in/dea_queries.json \
		./data/tcga_target_gtex \
		./data/selected_metadata \
		$(@D) \
		$(mods)/run_dea/run_deseq.R \
		--delimiter '\t' \
		--cpus 2

	touch $@

## --- 2 ---
# The fact that it run is detected in the `run_gsea` step by the `all.txt` file
./data/genesets/all.txt: \
		env/touchfile \
		$(local_mtpdb) \
		$(mods)/geneset_maker/make_genesets.py \
		$(mods)/geneset_maker/basic_gene_lists.json

	rm -rf $(@D)
	mkdir -p $(@D)

	. env/bin/activate; python \
		$(mods)/geneset_maker/make_genesets.py $(local_mtpdb) $(mods)/geneset_maker/basic_gene_lists.json \
		./data/genesets \
		--prune_direction "topdown" \
		--prune_similarity 0.9 \
		--verbose

## --- 5 --- Run the pre-ranked GSEA
./data/out/enrichments/done.flag: \
		./data/genesets/all.txt \
		$(mods)/gsea_runner/run_gsea.R \
		./data/deas/flag.txt \
		./data/ensg_data.csv

	mkdir -p $(@D)

	Rscript $(mods)/gsea_runner/run_gsea.R \
		"./data/deas/" "$(data_dir)/genesets" "$(@D)" \
		--ensg-hugo-data ./data/ensg_data.csv

	touch $@

./data/out/absolute_enrichments/done.flag: \
		./data/genesets/all.txt \
		$(mods)/gsea_runner/run_gsea.R \
		./data/deas/flag.txt \
		./data/ensg_data.csv

	mkdir -p $(@D)

	Rscript $(mods)/gsea_runner/run_gsea.R \
		"./data/deas/" "$(data_dir)/genesets" "$(@D)" \
		--ensg-hugo-data ./data/ensg_data.csv \
		--absolute

	touch $@

## --- 6 --- Generate plots
/tmp/genesets_tree/tree.txt: ./data/genesets/all.txt
	mkdir -p $(@D)
	tree ./data/genesets -df | head -n -2 > $@

ALL += ./data/out/figures/enrichments/done.flag
./data/out/figures/enrichments/done.flag: \
		./data/out/enrichments/done.flag \
		./data/genesets/all.txt \
		$(mods)/plotting/gsea_plotting_graphs.R \
		$(mods)/plotting/general_heatmap.R

	mkdir -p $(@D)

	Rscript $(mods)/plotting/gsea_plotting_graphs.R \
		./data/out/enrichments/ \
		$(@D)

	touch $@

ALL += ./data/out/figures/enrichments/pancan_heatmap.png
./data/out/figures/enrichments/pancan_heatmap.png: \
		./data/out/enrichments/done.flag \
		$(mods)/plotting/general_heatmap.R \
		/tmp/genesets_tree/tree.txt

	Rscript $(mods)/plotting/general_heatmap.R \
		./data/out/enrichments/ \
		/tmp/genesets_tree/tree.txt \
		$@ \
		--height 15

## This is duplicated from above but i'm not sure how to merge them.
ALL += ./data/out/figures/absolute_enrichments/done.flag
./data/out/figures/absolute_enrichments/done.flag: \
		./data/out/absolute_enrichments/done.flag \
		./data/genesets/all.txt \
		$(mods)/plotting/gsea_plotting_graphs.R \
		$(mods)/plotting/general_heatmap.R

	mkdir -p $(@D)

	Rscript $(mods)/plotting/gsea_plotting_graphs.R \
		./data/out/absolute_enrichments/ \
		$(@D)

	touch $@

ALL +=./data/out/figures/absolute_enrichments/pancan_heatmap.png
./data/out/figures/absolute_enrichments/pancan_heatmap.png: \
		./data/out/absolute_enrichments/done.flag \
		$(mods)/plotting/general_heatmap.R \
		/tmp/genesets_tree/tree.txt

	Rscript --no-save --no-restore --verbose $(mods)/plotting/general_heatmap.R \
		./data/out/absolute_enrichments/ \
		/tmp/genesets_tree/tree.txt \
		$@ \
		--height 15

ALL +=./data/out/figures/combined_heatmap.png
./data/out/figures/combined_heatmap.png: \
		./data/out/absolute_enrichments/done.flag \
		./data/out/enrichments/done.flag \
		$(mods)/plotting/fused_general_heatmap.R \
		/tmp/genesets_tree/tree.txt

	Rscript --no-save --no-restore --verbose $(mods)/plotting/fused_general_heatmap.R \
		./data/out/enrichments/ \
		./data/out/absolute_enrichments/ \
		/tmp/genesets_tree/tree.txt \
		$@ \
		--height 15

PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all

