.DELETE_ON_ERROR:

rexec = Rscript --no-save --no-restore --verbose

ALL+=./data/tcga_clinical_metadata.csv
./data/tcga_clinical_metadata.csv: \
	./src/modules/preprocess_tcga_clinical_metadata.R \
	./data/raw_tcga_clinical_metadata.csv
	
	mkdir -p $(@D)
	$(rexec) $< ./data/raw_tcga_clinical_metadata.csv > $@

./data/raw_tcga_clinical_metadata.csv: \
	./src/modules/retrieve_tcga_clinical_metadata.py

	python $< 'ALL' > $@

ALL+=./data/out/tcga_missing_metadata.csv
./data/out/tcga_missing_metadata.csv: \
	./src/modules/inspect_tcga_metadata.R \
	./data/tcga_clinical_metadata.csv

	$(rexec) $< ./data/tcga_clinical_metadata.csv > $@

PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all

