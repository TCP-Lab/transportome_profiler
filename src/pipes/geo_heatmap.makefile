#? Generate a heatmap from the GEO samples
#?
#? The main task of this pipeline is to preprocess the data that comes in
#? in a bajillion different formats.

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
	gunzip -cfv $< > $@

%.csv: %.tsv
	xsv fmt -d '\t' $< > $@

%.csv: %.xls
	xls2csv -x $< -c $@

# GSE121842
GEO += data/geo/GSE121842.csv
data/geo/GSE121842.csv: data/GSE121842.csv
	mkdir -p $(@D)
	cat $< | panid "GeneID:hgnc_symbol>ensg:ensg" | sponge | xsv select !GeneType | \
		xsv search -s ensg "ENSG" | sponge | src/helper_scripts/log_values | sponge > $@

# GSE107422
#GEO += data/geo/GSE107422.csv
#data/geo/GSE107422.csv: data/in/GSE107422_RAW.tar
#	mkdir -p $(@D)
#	# This is the worst file I've ever seen.
#	# The tar has a bunch of .tsv (but named .txt) gzip-compressed files
#	# with the colnames: unique_id	<name_of_sample>
#	# where "unique_id" is a frankenstein code with (seemingly) this pattern:
#	# gi|<some number>|ref|<refseq ID?>
#	# ...
#	# Who did this?
#	python src/modules/geo_data/prep_GSE107422.py $< | sponge | \
#		panid "unique_id:refseq_rna_id>ensg:ensg" | sponge | \
#		xsv search -s ensg "ENSG" > $@

#GEO += data/geo/GSE201284.csv
#data/geo/GSE201284.csv: data/GSE201284.csv
#	mkdir -p $(@D)
#	cat $< | panid "gene_id:hgnc_symbol>ensg:ensg" | sponge | \
#		xsv search -s ensg "ENSG" > $@

data/geo/GSE159857.csv: data/GSE159857.csv
	mkdir -p $(@D)
	cat $< | panid "GeneSymbol:hgnc_symbol>ensg:ensg" | sponge | \
		xsv search -s ensg "ENSG" > $@

data/geo/%.metadata.series: data/geo/%.csv
	python src/modules/geo_data/get_series.py `basename $< .csv` > $@

COORDS = data/in/config/series_coordinates.json
data/geo/%.metadata.raw : data/geo/%.metadata.series $(COORDS)
	cat $< | python src/modules/geo_data/meta_from_series.py \
		$$(jq ".$$(basename $@ .metadata).id? // .id" $(COORDS) -r) \
		$$(jq ".$$(basename $@ .metadata).var? // .var" $(COORDS) -r) > $@

data/geo/%.metadata: data/geo/%.metadata.raw
	python src/modules/geo_data/fix_geo_metadata.py $< > $@

GEO += data/geo/GSE159857_LUAD.csv
GEO += data/geo/GSE159857_LUSC.csv
data/geo/GSE159857_LUAD.csv data/geo/GSE159857_LUAD.metadata \
	data/geo/GSE159857_LUSC.csv data/geo/GSE159857_LUSC.metadata: \
	data/geo/GSE159857.metadata data/geo/GSE159857.csv
	metasplit "$<@sample_id?type=LUAD" data/geo/GSE159857.csv data/geo/GSE159857_LUAD.csv --always_include ensg
	metasplit "$<@sample_id?type=LUSC" data/geo/GSE159857.csv data/geo/GSE159857_LUSC.csv --always_include ensg

	cat $< | xsv search LUAD > data/geo/GSE159857_LUAD.metadata
	cat $< | xsv search LUSC > data/geo/GSE159857_LUSC.metadata

data/geo/%.dea.csv: data/geo/%.metadata data/geo/%.csv
	# Run metasplit on the data to separate case and controls
	metasplit "$<@sample_id?var_0=case" $(word 2, $^) $@.case --always_include ensg
	metasplit "$<@sample_id?var_0=control" $(word 2, $^) $@.control --always_include ensg
	generanker --id-col ensg $@.case $@.control $(RANK_METHOD) > $@
	rm $@.case $@.control

## --- Generate the genesets from the MTPDB
./data/genesets.json ./data/genesets_repr.txt: \
		./data/MTPDB.sqlite \
		$(mods)/make_genesets.py \
		./data/in/config/gene_lists/basic.json

	python $(mods)/make_genesets.py ./data/MTPDB.sqlite ./data/in/config/gene_lists/basic.json \
		./data/genesets.json ./data/genesets_repr.txt \
		--prune_direction $(PRUNE_DIRECTION) \
		--prune_similarity $(PRUNE_SIMILARITY) \
		--verbose

## -- Run the pre-ranked GSEA
./data/out/geo_enrichments/%.gsea.csv: data/geo/%.dea.csv ./data/genesets.json
	mkdir -p $(@D)
	$(rexec) $(mods)/run_gsea_once.R \
		$< "./data/genesets.json" $@ \
		--ensg-hugo-data data/in/ensg_data.csv \
		$(_gsea_runtime_flags)

./data/out/absolute_geo_enrichments/%.gsea.csv: data/geo/%.dea.csv ./data/genesets.json
	mkdir -p $(@D)
	$(rexec) $(mods)/run_gsea_once.R \
		$< "./data/genesets.json" $@ \
		--ensg-hugo-data data/in/ensg_data.csv \
		--absolute \
		$(_gsea_runtime_flags)

# Make the requirements for this aggregative rule
gseas = $(addprefix data/out/geo_enrichments/, $(notdir $(GEO:.csv=.gsea.csv)))
abs_gseas = $(addprefix data/out/absolute_geo_enrichments/, $(notdir $(GEO:.csv=.gsea.csv)))

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
		--alpha 0.20 \
		--extra_title "alpha 0.20, metric $(RANK_METHOD)" \
		--height 15 \
		$(_heatmap_plot_flags)

PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all
