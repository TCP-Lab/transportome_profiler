#? Generate summary plots and other bits and bobs
#?
#? This takes the output of running the heatmaps pipeline over and over
#? with all the different ranking methods and generates the final plots

.ONESHELL:

# Shorthands
mods = ./src/modules
rexec = Rscript --no-save --no-restore --verbose

%.csv: %.tsv
	xsv fmt -d '\t' $< > $@

ALL_RES = norm_bws_test deseq_shrinkage norm_cohen_d norm_fold_change \
		  norm_s2n_ratio

data/extracted_results/%_deas.csv data/extracted_results/%_geo.csv &: data/in/results/%.tar.gz
	mkdir -p ${@D}
	python src/preprocess_result.py $< $(@D)

ALL = $(addprefix data/extracted_results/,$(addsuffix _deas.csv,${ALL_RES}) $(addsuffix _geo.csv,${ALL_RES}))

PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all
