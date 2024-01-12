#!/usr/bin/env bash

#? Make the test expression data from the TCGA/GTEX dataset

echo "Extracting sources"
gunzip -cv ./data/in/expression_matrix.tsv.gz > ./data/expression_matrix.tsv
gunzip -cv ./data/in/expression_matrix_metadata.tsv.gz > ./data/expression_matrix_metadata.tsv

echo "Reformatting as .csv"
xsv fmt -d '\t' ./data/expression_matrix_metadata.tsv > ./data/expression_matrix_metadata.csv
xsv fmt -d '\t' ./data/expression_matrix.tsv > ./data/expression_matrix.csv

echo "Running metasample"
./src/metasample.py --metavars _primary_site --always-include sample \
    ./data/expression_matrix.csv ./data/test_expression_matrix.csv \
    ./data/expression_matrix_metadata.csv sample \
    "25%" 5000

echo "Compressing and saving as input file"
gzip -cv ./data/test_expression_matrix.csv > ./data/in/test_expression_matrix.csv.gz
