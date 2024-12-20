#? Generate summary plots and other bits and bobs
#?
#? This takes the output of running the heatmaps pipeline over and over
#? with all the different ranking methods and generates the final plots

.ONESHELL:

# Shorthands
mods = ./src/modules
rexec = Rscript --no-save --no-restore --verbose
target = data/

%.csv: %.tsv
	xsv fmt -d '\t' $< > $@

data/extracted_results/%_deas.csv data/extracted_results/%_geo.csv &:
	mkdir -p ${@D}
	python src/preprocess_result.py data $(@D)

ALL_DEAS = $(addprefix data/extracted_results/,$(addsuffix _deas.csv,${ALL_RES}) $(addsuffix _geo.csv,${ALL_RES}))

ALL = ${ALL_DEAS}

# We have to choose just one tarball to extract the genesets from
# This is pretty arbitrary
data/genesets.json: data/in/results/deseq_shrinkage.tar
	tar -xvf $< $@ > $@

data/filter_genes.txt: data/genesets.json
	cat $< | jq -r '.[] | select(.name == "whole_transportome").data | @csv' > $@

ALL += data/out/plots/shared_dysregulation.png
data/out/plots/shared_dysregulation.png: data/extracted_results/deseq_shrinkage_deas.csv \
		data/filter_genes.txt ${mods}/plotting/plot_shared_dysregulation.R \
		data/in/ensg_data.csv
	mkdir -p ${@D}
	${rexec} ${mods}/plotting/plot_shared_dysregulation.R $@ $< data/in/ensg_data.csv \
		--selected_genes data/filter_genes.txt --png --res 400

ALL += data/out/plots/continuous_congruency_colorectal.png
data/out/plots/continuous_congruency_colorectal.png: ${ALL_DEAS} data/filter_genes.txt \
		${mods}/plotting/plot_ranks.R
	mkdir -$ ${@D}
	${rexec} ${mods}/plotting/plot_ranks.R

PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all
