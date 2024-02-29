#? Generate a heatmap from the GEO samples
#? 
#? The main task of this pipeline is to preprocess the data that comes in
#? in a bajillion different formats.

# Shorthands
mods = ./src/modules
rexec = Rscript --no-save --no-restore --verbose

## -- Generic rules
## --- Decompress sources
./data/%: ./data/in/%.gz
	gunzip -cfv $< > $@

%.csv: %.tsv
	xsv fmt -d '\t' $< > $@

# GSE121842
ALL += ./data/geo/GSE121842_clean.csv
data/geo/GSE121842_clean.csv: data/GSE121842.csv
	mkdir -p $(@D)
	cat $< | panid "GeneID:hgnc_symbol>ensg:ensg" | sponge | xsv select !GeneType | \
		xsv search -s ensg "ENSG" > $@
	
# GSE107422
ALL += data/geo/GSE107422_clean.csv
data/geo/GSE107422_clean.csv: data/in/GSE107422_RAW.tar
	mkdir -p $(@D)
	# This is the worst file I've ever seen.
	# The tar has a bunch of .tsv (but named .txt) gzip-compressed files
	# with the colnames: unique_id	<name_of_sample>
	# where "unique_id" is a frankenstein code with (seemingly) this pattern:
	# gi|<some number>|ref|<refseq ID?>
	# ...
	# Who did this?
	python src/modules/geo_data/prep_GSE107422.py $< | sponge | \
		panid "unique_id:refseq_rna_id>ensg:ensg" | sponge | \
		xsv search -s ensg "ENSG" > $@

PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all
