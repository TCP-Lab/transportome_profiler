#? Generate a heatmap from the GEO samples
#?
#? The main task of this pipeline is to preprocess the data that comes in
#? in a bajillion different formats.

.ONESHELL:

OPTS=./data/in/config/heatmaps_runtime_options.json

# Number of threads to use to parallelize the ranking process
N_THREADS ?= $(shell cat $(OPTS) | jq -r '.threads')
# Method to use
RANK_METHOD ?= $(shell cat $(OPTS) | jq -r '.rank_method')
PRUNE_SIMILARITY ?= $(shell cat $(OPTS) | jq -r '.prune_similarity')
PRUNE_DIRECTION ?= $(shell cat $(OPTS) | jq -r '.prune_direction')
ALPHA_THRESHOLD ?= $(shell cat $(OPTS) | jq -r '.alpha_threshold')

# Shorthands
mods = ./src/modules
rexec = Rscript --no-save --no-restore --verbose

## -- Generic rules
## --- Decompress sources
./data/%: ./data/in/%.gz
	mkdir -p $(@D)
	gunzip -cfv $< > $@

%.csv: %.tsv
	xsv fmt -d '\t' $< > $@

data/%: data/in/%
	cp $< $@

data/geo/%.meta.csv: data/geo/%.rawmeta.csv
	python src/modules/geo_data/fix_geo_metadata.py $< > $@

# GSE159857 is special: it has both LUAD and LUSC samples.
# This rule splits it into the two parts
data/geo/GSE159857_LUAD.counts.csv data/geo/GSE159857_LUAD.meta.csv \
	data/geo/GSE159857_LUSC.counts.csv data/geo/GSE159857_LUSC.meta.csv: \
	data/geo/GSE159857.meta.csv data/geo/GSE159857.counts.csv
	metasplit "$<@ena_run?type=LUAD" data/geo/GSE159857.counts.csv data/geo/GSE159857_LUAD.counts.csv --always_include gene_id
	metasplit "$<@ena_run?type=LUSC" data/geo/GSE159857.counts.csv data/geo/GSE159857_LUSC.counts.csv --always_include gene_id

	cat $< | xsv search LUAD > data/geo/GSE159857_LUAD.meta.csv
	cat $< | xsv search LUSC > data/geo/GSE159857_LUSC.meta.csv

data/geo/%.dea.csv: data/geo/%.meta.csv data/geo/%.counts.csv
	# Run metasplit on the data to separate case and controls
	metasplit "$<@ena_run?status=case" $(word 2, $^) $@.case.unlogged --always_include gene_id
	metasplit "$<@ena_run?status=control" $(word 2, $^) $@.control.unlogged --always_include gene_id
	cat $@.case.unlogged | src/helper_scripts/log_values > $@.case
	cat $@.control.unlogged | src/helper_scripts/log_values > $@.control
	generanker --output-file $@ --id-col gene_id $@.case $@.control $(RANK_METHOD)
	#rm $@.case $@.control $@.case.unlogged $@.control.unlogged

# GSE254461 is specialer: it has a bunch of different tumor types
# This rule splits it into the two parts
data/geo/GSE254461_BRCA.counts.csv data/geo/GSE254461_BRCA.meta.csv \
	data/geo/GSE254461_KI.counts.csv data/geo/GSE254461_KI.meta.csv \
	data/geo/GSE254461_LIV.counts.csv data/geo/GSE254461_LIV.meta.csv \
	data/geo/GSE254461_COAD.counts.csv data/geo/GSE254461_COAD.meta.csv: \
	data/geo/GSE254461.meta.csv data/geo/GSE254461.counts.csv
	for id in BRCA KI LIV COAD; do
		metasplit "$<@ena_run?type=$$id" data/geo/GSE254461.counts.csv data/geo/GSE254461_$${id}.counts.csv --always_include gene_id
		cat $< | xsv search $$id > data/geo/GSE254461_$${id}.meta.csv
	done

## --- Generate the genesets from the MTPDB
./data/genesets.json ./data/genesets_repr.txt: \
		./data/MTPDB.sqlite \
		$(mods)/make_genesets.py \
		./data/in/config/gene_lists/no_experimental_ions.json

	python $(mods)/make_genesets.py ./data/MTPDB.sqlite ./data/in/config/gene_lists/no_experimental_ions.json \
		./data/genesets.json ./data/genesets_repr.txt \
		--prune_direction $(PRUNE_DIRECTION) \
		--prune_similarity $(PRUNE_SIMILARITY) \
		--min_recurse_set_size 0 \
		--verbose

## -- Run the pre-ranked GSEA
./data/out/geo_enrichments/%.gsea.csv: data/geo/%.dea.csv ./data/genesets.json \
		./data/ensg_data.csv
	mkdir -p $(@D)
	$(rexec) $(mods)/run_gsea_once.R \
		$< "./data/genesets.json" $@ \
		--ensg-hugo-data data/ensg_data.csv \
		$(_gsea_runtime_flags)

./data/out/absolute_geo_enrichments/%.gsea.csv: data/geo/%.dea.csv ./data/genesets.json \
		./data/ensg_data.csv
	mkdir -p $(@D)
	$(rexec) $(mods)/run_gsea_once.R \
		$< "./data/genesets.json" $@ \
		--ensg-hugo-data data/ensg_data.csv \
		--absolute \
		$(_gsea_runtime_flags)

# A list of all GSEs that we have.
GEO = GSE22260 GSE29580 GSE121842 GSE159857_LUAD GSE159857_LUSC GSE60052 \
	  GSE254461_BRCA GSE254461_KI GSE254461_LIV GSE254461_COAD GSE103001 \
	  GSE119224 GSE151352 GSE159260 GSE234927
# Make the requirements for this aggregative rule
gseas = $(addprefix data/out/geo_enrichments/,$(addsuffix .gsea.csv,$(GEO)))
abs_gseas = $(addprefix data/out/absolute_geo_enrichments/, $(addsuffix .gsea.csv,$(GEO)))

## --- Make the large heatmap with all the results
ALL +=./data/out/figures/geo_deregulation_heatmap.png
./data/out/figures/geo_deregulation_heatmap.png: \
		$(gseas) \
		$(abs_gseas) \
		$(mods)/plotting/plot_large_heatmap.R \
		./data/genesets_repr.txt \
		./data/genesets.json

	mkdir -p $(@D)

	$(rexec) $(mods)/plotting/plot_large_heatmap.R \
		./data/out/geo_enrichments/ \
		./data/genesets.json \
		./data/genesets_repr.txt \
		$@ \
		--dots_gsea_results ./data/out/absolute_geo_enrichments/ \
		--alpha $(ALPHA_THRESHOLD) \
		--extra_title "alpha $(ALPHA_THRESHOLD), metric $(RANK_METHOD)" \
		--height 15 \
		--renames data/in/config/geo_renames.json \
		$(_heatmap_plot_flags)


deas = $(addprefix data/geo/, $(addsuffix .dea.csv,$(GEO)))

## --- Make the output summary plots
data/filter_genes.txt: data/genesets.json
	cat $< | jq -r '.[] | select(.name == "whole_transportome").data | @csv' > $@

data/geo_merged_deas.csv: $(deas)
	mkdir -p ${@D}
	# This renames the header slot 'ranking' to the name of the file
	# This magical $ thing is to replace the .csv to nothing
	# the pattern is {variable/pattern/replacement}
	for name in ./data/geo/*.dea.csv; do \
		replacement=$$(basename "$${name/.dea.csv}"); \
		sed "1s/ranking/$${replacement}/" "$${name}" > $${name/.csv/}.renames.csv; \
	done
	# Move one random file to be the "base" that we can reduce into
	mv $$(find "./data/geo/" -name "*.renames.csv" -print -quit) ./data/base
	# The xsv select is there to remove the second "gene_id" col that "join" retains.
	for item in ./data/geo/*.renames.csv; do \
		xsv join gene_id ./data/base gene_id $${item} --full | xsv select '!gene_id[1]' > ./data/tmp;\
		mv ./data/tmp ./data/base; \
	done
	mv ./data/base $@
	# Since the .renames files may conflict with the old pipeline, I
	# just delete them here instead of changing the previous steps.
	# Sorry!
	rm ./data/geo/*.renames.csv

ALL +=./data/out/figures/top_geo_disregulation.png
./data/out/figures/top_geo_disregulation.png: \
		./data/geo_merged_deas.csv \
		./data/filter_genes.txt \
		${mods}/plotting/plot_shared_dysregulation.R \
		./data/ensg_data.csv
	mkdir -p ${@D}
	${rexec} ${mods}/plotting/plot_shared_dysregulation.R $@ $< data/ensg_data.csv \
		--id_col gene_id \
		--selected_genes data/filter_genes.txt --png --res 400

ALL +=./data/out/figures/geo_upset.png
./data/out/figures/geo_upset.png: \
		./data/geo_merged_deas.csv \
		./data/filter_genes.txt \
		${mods}/plotting/plot_general_upset.R \
		./data/ensg_data.csv
	mkdir -p ${@D}
	${rexec} ${mods}/plotting/plot_general_upset.R $@ $< \
		--id_col gene_id \
		--selected_genes data/filter_genes.txt --png --res 400

PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all

#.PRECIOUS = *.gsea.csv
