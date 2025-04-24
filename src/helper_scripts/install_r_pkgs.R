#!/usr/bin/env Rscript

install.packages("pak", repos = sprintf(
  "https://r-lib.github.io/p/pak/stable/%s/%s/%s",
  .Platform$pkgType,
  R.Version()$os,
  R.Version()$arch
))

packages <- c(
    "tidyverse",
    "BiocManager",
    "DESeq2",
    "fgsea",
    "argparser",
    "biomaRt",
    "assertthat",
    "uuid",
    "grid",
    "reshape2",
    "rjson",
    "archive",
    "glmGamPoi",
    "extrafont",
    "ComplexUpset",
    "devtools",
    "TCP-Lab/compare_ranks"
)

pak::pkg_install(packages)
