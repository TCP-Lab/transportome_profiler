rexec = Rscript --no-save --no-restore --verbose
mods = ./src/modules
N_THREADS = 1

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
		--method "norm_fold_change"

	touch $@

## --- Make the fold change plot
ALL += data/out/figures/fold_change_heatmap.png
./data/out/figures/fold_change_heatmap.png: \
	data/fc_deas/flag.txt

	mkdir -p $(@D)
	$(rexec) $(mods)/plotting/plot_fc.R \
		data/fc_deas $@ --png \
		--height 16

PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all
