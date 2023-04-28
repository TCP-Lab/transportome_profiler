 
.PHONY = all analyze clean retrieve make_genesets run_dea

data_dir = ./data/ 

all:

retrieve: ./src/retrieve_sources
	./src/retrieve_sources $(data_dir)

# The data dowloaded by retrieve
$(data_dir)/in/tcga_target_gtex $(data_dir)/in/MTPDB.sqlite &: retrieve

# This makes way too many (variable) files, so I made the target phony.
make_genesets: \
		$(data_dir)/in/MTPDB.sqlite \
		./src/geneset_maker/make_genesets.py \
		./src/geneset_maker/basic_gene_lists.json

	python \
		./src/geneset_maker/make_genesets.py $(data_dir)/in/MTPDB.sqlite ./src/geneset_maker/basic_gene_lists.json \
		--prune_direction "topdown" \
		--verbose

run_dea: $(data_dir)/in/tcga_target_gtex \
		./src/run_dea/GTEX_phenotype \
		./src/run_dea/deseq_run.R 

	Rscript ./src/run_dea/deseq_run.R "$(data_dir)/in/tcga_target_gtex" "$(data_dir)" "./src/run_dea/GTEX_phenotype"

run_gsea: 