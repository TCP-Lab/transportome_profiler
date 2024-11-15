#!/usr/bin/env bash

echo -e "n\ttissue\tstatus"
tar -xOf data/in/results/deseq_shrinkage.tar deseq_shrinkage.log | \
    grep "gene_ranker\.ranker: Loaded a" | cut -d ' ' -f 7,15 | \
    sed 's/ /\t/' | \
    sed 's/data\/deas\///' | sed 's/Head_n_Neck/HeadAndNeck/' | sed 's/_cancer//' | sed 's/_/\t/'
tar -xOf data/in/results/deseq_shrinkage.tar geo_deseq_shrinkage.log | \
    grep "gene_ranker\.ranker: Loaded a" | cut -d ' ' -f 7,15 | \
    sed 's/ /\t/' | \
    sed 's/data\/geo\///' | sed 's/.dea.csv./\t/'
