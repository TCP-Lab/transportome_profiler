#!/usr/bin/env Rscript

options(repos = "http://cloud.r-project.org/")

install <- function(pkg) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
        install.packages(pkg)
    } else {
        cat(paste0("Skipping ", pkg, " as it is already installed.\n"))
    }
}

install_bioc <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        BiocManager::install(pkg)
    } else {
        cat(paste0("Skipping ", pkg, " as it is already installed.\n"))
    }
}

install("tidyverse")
install("BiocManager")
install_bioc("DESeq2")
install_bioc("fgsea")
install("argparser")
install_bioc("biomaRt")
install("ggraph")
install("assertthat")
install("uuid")
install("grid")
install("reshape2")
install("rjson")
install_bioc("glmGamPoi")
install("extrafont")
