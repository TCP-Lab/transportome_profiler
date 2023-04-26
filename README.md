# Transportome Profiler

This repository contains the code for the analysis on the expression profile of the transportome in Cancer based on the [MTP-DB](https://github.com/CMA-Lab/MTP-DB).

This is a two-step process. The database is queried for information by the script in `src/geneset_maker`. The algorithm generates gene sets ready for use by GSEA. We then generate a series of DEG tables with `src/run_dea` based on the comparison of gene expression from TCGA (cancer) and GTEx (healthy) tissues. Finally, GSEA is called by `src/gsea_runner` in a pre-ranked manner on all the DEG tables with all of the genesets, making enrichment tables.
The enrichment tables are then processed by a script in `src/gsea_runner` to make enrichment plots.
