.DELETE_ON_ERROR:

#? Generate the TCGA metadata for the clinical stratification pipeline

rexec = Rscript --no-save --no-restore --verbose

# Generic decompression and conversion rules
./data/%: ./data/in/%.gz
	gunzip -cfv $< > $@

%.csv: %.tsv
	xsv fmt -d '\t' $< > $@

./data/tcga_clinical_metadata.csv: \
	./src/modules/metadata/preprocess_tcga_clinical_metadata.R \
	./data/raw_tcga_clinical_metadata.csv
	
	mkdir -p $(@D)
	$(rexec) $< ./data/raw_tcga_clinical_metadata.csv > $@

./data/raw_tcga_clinical_metadata.csv: \
	./src/modules/metadata/retrieve_tcga_clinical_metadata.py

	python $< 'ALL' > $@

./data/out/tcga_missing_metadata.csv: \
	./src/modules/metadata/inspect_tcga_metadata.R \
	./data/tcga_clinical_metadata.csv

	$(rexec) $< ./data/tcga_clinical_metadata.csv > $@

ALL+=./data/harmonized_metadata.csv
./data/harmonized_metadata.csv: \
	./src/modules/metadata/inflate_cases_metadata.py \
	./data/tcga_clinical_metadata.csv \
	./data/expression_matrix_metadata.csv

	python $< \
		./data/tcga_clinical_metadata.csv@submitter_id \
		./data/expression_matrix_metadata.tsv@sample | \
		xsv select '!gender' > $@

PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all

