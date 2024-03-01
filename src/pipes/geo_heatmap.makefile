#? Generate a heatmap from the GEO samples
#? 
#? The main task of this pipeline is to preprocess the data that comes in
#? in a bajillion different formats.

DELETE_ON_ERROR:

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
		xsv search -s ensg "ENSG" > $@
	
# GSE107422
GEO += data/geo/GSE107422.csv
data/geo/GSE107422.csv: data/in/GSE107422_RAW.tar
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

GEO += data/geo/GSE201284.csv
data/geo/GSE201284.csv: data/GSE201284.csv
	mkdir -p $(@D)
	cat $< | panid "gene_id:hgnc_symbol>ensg:ensg" | sponge | \
		xsv search -s ensg "ENSG" > $@

GEO += data/geo/GSE159857.csv
data/geo/GSE159857.csv: data/GSE159857.csv
	mkdir -p $(@D)
	cat $< | panid "GeneSymbol:hgnc_symbol>ensg:ensg" | sponge | \
		xsv search -s ensg "ENSG" > $@

data/geo/%.metadata.series: data/geo/%.csv
	python src/modules/geo_data/get_series.py `basename $< .csv` > $@

data/geo/%.metadata.unaligned: data/geo/%.metadata.series data/in/series_coordinates.json
	cat $< | python src/modules/geo_data/meta_from_series.py \
		$$(jq ".$$(basename $@ .metadata).id? // .id" data/in/series_coordinates.json -r) \
		$$(jq ".$$(basename $@ .metadata).var? // .var" data/in/series_coordinates.json -r) > $@	

data/geo/%.metadata: data/geo/%.metadata.unaligned data/geo/%.csv
	# The recipy order is important, first the metadata, then the csv
	python src/modules/geo_data/heuristically_align_metadata.py $^ > $@

ALL += data/all_geo_data.csv
data/all_geo_data.csv: $(GEO)
	./src/helper_scripts/fuse_csv $^ > $@

ALL += data/all_geo_data.metadata
data/all_geo_data.metadata: $(GEO:.csv=.metadata)
	./src/helper_scripts/fuse_csv $^ > $@


PHONY += all
all: $(ALL)

.PHONY = $(PHONY)
.DEFAULT_GOAL := all
