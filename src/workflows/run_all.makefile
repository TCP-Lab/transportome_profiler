# This makefile runs ALL the analyses

RANK_METHODS = fold_change deseq_shrinkage cohen_d norm_cohen_d norm_fold_change s2n_ratio norm_s2n_ratio bws_test norm_bws_test

ALL += $(addsuffix .tar.gz,$(addprefix packaged_output_,${RANK_METHODS}))

packaged_output_%.tar.gz:
	RANK_METHOD=$(shell echo $@ | sed -n "s/packaged_output_\(\S*\).tar.gz/\1/p") \
	kerblam run heatmaps -l
	
	RANK_METHOD=$(shell echo $@ | sed -n "s/packaged_output_\(\S*\).tar.gz/\1/p") \
	kerblam run geo_heatmap -l
	
	chown $$USER -R ./data
	tar cf $@ ./data/out
	kerblam data clean

PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all
