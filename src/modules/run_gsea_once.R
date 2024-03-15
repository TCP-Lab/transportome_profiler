#!/usr/bin/env Rscript

# This script runs GSEA on all the DEG tables with all the genesets
# and saves the resulting deg tables to an output folder.

set.seed(42)

# If you are running this from RStudio, you can skip this >>>>>>>>>>>>>>>>>>
if (sys.nframe() == 0L) {
  # Parsing arguments
  requireNamespace("argparser")

  parser <- argparser::arg_parser("Run GSEA on one DEG table")

  parser |>
    argparser::add_argument(
      "input_deg_table", help="DEG table to run on.", type="character"
    ) |>
    argparser::add_argument(
      "input_genesets_json", help = "JSON file with genesets tree",
      type = "character"
    ) |>
    argparser::add_argument(
      "output_deg_table", help = "Output file with GSEA results",
      type = "character"
    ) |>
    argparser::add_argument(
      "--ensg-hugo-data", help = "A .csv file with the 'ensembl_gene_id', 'hgnc_symbol', 'gene_biotype' columns. If unspecified, downloads from ENSEMBL",
      type = "character", default = NA
    ) |>
    argparser::add_argument(
      "--low-memory", help = "If specified, saves only output tables, skipping plots, and using less memory.",
      flag = TRUE, type = "logical"
    ) |>
    argparser::add_argument(
      "--absolute", help = "If specified, runs GSEA by sorting on the absolute values instead of the real values.",
      flag = TRUE, type = "logical"
    ) |>
    argparser::add_argument(
      "--unweighted", help = "If specified, runs GSEA in an unweighted manner",
      flag = TRUE, type = "logical"
    ) |>
    argparser::add_argument(
      "--save-plots", help = "If specified, also save GSEA plots alongside tables.",
      flag = TRUE, type = "logical"
    ) -> parser

  args <- argparser::parse_args(parser)
}
# To here <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

suppressMessages({
  requireNamespace("biomaRt")
  requireNamespace("fgsea")
  requireNamespace("rjson")
  library(tidyverse, quietly = TRUE)
})

# This suppresses the messages from readr::read_table
options(readr.num_columns = 0)

#' This takes a node json from `Bonsai` and loads it as gene sets
#'
#' @param file The node list file
#' @param biomart_data A data.frame with at least the "ensembl_gene_id" and
#'   "hgnc_symbol" columns, with the correspondence from ENSG to symbol.
#'   Such a table can be retrieved from Biomart with biomaRt.
#' @returns A list of vectors, each with the gene symbols of that gene set

# Each entry in the json looks like this:
# "{id}": {"name": "{name}", "data": [data], "parent": "{id}"}
load_genesets_from_node_json <- function(file, biomart_data) {
  data <- rjson::fromJSON(readr::read_file(file))

  # We need a list of id: data
  data <- lapply(data, \(x) {x[["data"]]})

  return(data)
}

purge_ensg_versions <- function(data, id_col = "gene_id") {
  data |> mutate("{id_col}" := str_remove(.data[[id_col]], "\\.[0-9]+$")) -> data
  data |> distinct(.data[[id_col]], .keep_all = TRUE) -> data
  data
}


#' Make a frame ready for GSEA from a DEG file, made by BioTEA or DESeq2.
#'
#' Some compatibility is needed to parse DESeq2 files, as they have no
#' "SYMBOL" column.
#'
#' Uses biomart data to retain only coding genes.
#'
#' @param deg_file (full) Path to the DEG file that needs to be extracted, as .csv.
#' @param biomart_data A data.frame with at least the "ensembl_gene_id" and
#'   "gene_biotype" columns.
#'   Such a table can be retrieved from Biomart with biomaRt.
#' @param absolute Boolean. If TRUE, will compute the absolute values of the
#'   statistic instead of the actual ones.
#' @param id_col ID of the column with the row/gene ids
#' @param rank_col ID of the column with the rankings
#'
#' @returns A named vector of gene_names : statistic, ready for fgsea::fgsea
extract_ranks <- function(deg_file, biomart_data, id_col="gene_id", rank_col="ranking", absolute = FALSE) {
  data <- read_csv(deg_file, show_col_types = FALSE)

  if (! id_col %in% colnames(data)) {
    stop(paste0("Cannot find ID column ", id_col, " in input data frame"))
  }
  if (! rank_col %in% colnames(data)) {
    stop(paste0("Cannot find ID column ", rank_col, " in input data frame"))
  }

  named_vec <- data[[rank_col]]
  names(named_vec) <- purge_ensg_versions(data, id_col = id_col)[[id_col]]

  coding <- biomart_data$ensembl_gene_id[biomart_data$gene_biotype == "protein_coding"]

  named_vec <- named_vec[names(named_vec) %in% coding]

  if (any(is.na(named_vec))) {
    warning("Some values in named vector are NA. Setting to 0")
    named_vec[is.na(named_vec)] <- 0
  }

  if (absolute) {
    named_vec <- abs(named_vec)
  }

  named_vec
}


#' Run GSEA with many genesets on some data
#'
#' This is a tiny wrapper for future compatibility in case we need to change
#' from fgsea.
#'
#' @param genesets A list of genesets to check, with each geneset a vector of
#'   gene names.
#' @param ranks A named list of gene names: statistic to use as ranked list for gsea.
#'
#' @returns A table with GSEA results. See fgsea:fgsea for details.
run_gsea <- function(genesets, ranks, unweighted = FALSE) {
  result <- fgsea::fgsea(
    pathways = genesets,
    stats = ranks,
    gseaParam = if (unweighted) {0} else {1}
  )

  result
}

#' (run and) Plot GSEA
#'
#' @param genesets A list of genesets to check, with each geneset a vector of
#'   gene names.
#' @param ranks A named list of gene names: statistic to use as ranked list for gsea.
#'
#' @returns A list of geneset: ggplot object with the generated GSEA plots.
plot_gsea <- function(genesets, ranks) {
  results <- list()
  for (i in seq_along(genesets)) {
    p <- fgsea::plotEnrichment(genesets[[i]], ranks)
    results[[names(genesets)[i]]] <- p
  }

  return(results)
}


run_one_gsea <- function(input_file_path, genesets_path, biomart_data, output_path, absolute = FALSE, unweighted = FALSE) {
  cat("Loading genesets...\n")
  genesets <- load_genesets_from_node_json(genesets_path)
  cat(paste0("Loaded ", length(genesets), " genesets.\n"))

  cat(paste0("Running GSEA on ", input_file_path, "\n"))
  ranks <- extract_ranks(input_file_path, biomart_data, absolute = absolute)

  result <- run_gsea(genesets, ranks, unweighted)

  cat(paste0("Saving data to ", output_path, ".csv"), "\n")
  save_result(result, dirname(output_path), basename(output_path))
}


#' Save a result from `run_all_gsea` to a file.
#'
#' @param result An item from a list made by `run_all_gsea`.
#' @param out_dir The output directory.
#' @param name The name to give to the output file
#' @param plot If `result` is a plot, pass plot = TRUE to treat it as one.
#'
#' @returns NULL
save_result <- function(result, out_dir, name, plot = FALSE) {
  out_path <- file.path(out_dir, name)

  if (plot) {
    for (i in seq_along(result)) {
      pdf(paste0(out_path, names(result)[i], ".pdf"), width = 12, height = 8)
      print(result[[i]])
      dev.off()
    }
    return(TRUE)
  }

  write_csv(result, out_path)
  return(TRUE)
}

if (FALSE){ # LOCAL DEBUGGING RUN ONLY

embl <- biomaRt::useEnsembl(biomart = "genes")
hs.embl <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl", mart = embl)
ensg_data <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
  mart = hs.embl
)

write_csv(ensg_data, "./data/ensg_data.csv")

results <- run_all_gsea(
  "/home/hedmad/Files/data/transportome_profiler/deas",
  "/home/hedmad/Files/data/transportome_profiler/genesets",
  ensg_data
)

save_results(results, out_dir = "/home/hedmad/Files/data/transportome_profiler/out/enrichments", skip_plots = TRUE)

} # ----------------------------------------------------------------------------

# If you are running this from RStudio, you can skip this >>>>>>>>>>>>>>>>>>
if (sys.nframe() == 0L) {

  if (is.na(args$ensg_hugo_data)) {
    embl <- biomaRt::useEnsembl(biomart = "genes")
    hs.embl <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl", mart = embl)
    ensg_data <- biomaRt::getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
      mart = hs.embl
    )
  } else {
    ensg_data <- read.csv(args$ensg_hugo_data, header = TRUE)
  }

  run_one_gsea(
    args$input_deg_table,
    args$input_genesets_json,
    ensg_data,
    output_path = args$output_deg_table,
    absolute = args$absolute,
    unweighted = args$unweighted
  )
}
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
