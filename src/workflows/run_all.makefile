# This makefile runs ALL the analyses

RANK_METHODS = fold_change deseq_shrinkage cohen_d norm_cohen_d norm_fold_change s2n_ratio norm_s2n_ratio bws_test norm_bws_test

ALL_OUTPUTS += $(addsuffix .tar.gz,$(addprefix packaged_output_,${RANK_METHODS})) extra_plots.tar.gz

packaged_output_%.tar.gz:
	# This sed call just takes out the ranking method string
	kerblam run heatmap -- \
	-e RANK_METHOD=$(shell echo $@ | sed -n "s/packaged_output_\(\S*\).tar.gz/\1/p")
	kerblam run geo_heatmap -- \
	-e RANK_METHOD=$(shell echo $@ | sed -n "s/packaged_output_\(\S*\).tar.gz/\1/p")

	# Copy the large tables so they get checked out
	cp ./data/merged_deas.csv ./data/out/merged_deas.csv
	cp ./data/geo_merged_deas.csv ./data/out/geo_merged_deas.csv
	kerblam data pack --output-only $@
	# We're done - let's clean up so we can start fresh
	kerblam data clean -y

extra_plots.tar.gz:
	kerblam run expression_plot
	kerblam data pack --output-only $@
	kerblam data clean -y

ALL += data/out/transportome_profiler_results.tar
data/out/transportome_profiler_results.tar: $(ALL_OUTPUTS)
	tar -cvf $@ $^

PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all
