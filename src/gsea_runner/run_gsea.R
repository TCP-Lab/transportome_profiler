#!/usr/bin/env Rscript

# This script runs GSEA on all the DEG tables with all the genesets
# and saves the resulting deg tables to an output folder.

# If you are running this from RStudio, you can skip this >>>>>>>>>>>>>>>>>>
if (sys.nframe() == 0L) {
  # Parsing arguments
  requireNamespace("argparser")

  parser <- argparser::arg_parser("Run GSEA on DEG tables")

  parser |>
    argparser::add_argument(
      "input_expr_matrix", help="Expression matrix to run GSEA on.", type="character"
    ) |>
    argparser::add_argument(
      "input_genesets_folder", help = "Folder with input genesets as .txt files",
      type = "character"
    ) |>
    argparser::add_argument(
      "output_dir", help = "Output directory",
      type = "character"
    ) |>
    argparser::add_argument(
      "--low-memory", help = "If specified, saves only output tables, skipping plots, and using less memory.",
      flag = TRUE, type = "logical"
    ) |>
    argparser::add_argument(
      "--save-plots", help = "If specified, also save GSEA plots alongside tables.",
      flag = TRUE, type = "logical"
    ) -> parser

  args <- argparser::parse_args(parser)
}
# To here <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

requireNamespace("biomaRt")
requireNamespace("fgsea")
library(tidyverse, quietly = TRUE)

# This suppresses the messages from readr::read_table
options(readr.num_columns = 0)

load_genesets <- function(folder, biomart_data) {
  #' Load all genesets from FOLDER.
  #'
  #' Genesets should have been made by make_genesets.py
  #'
  #' @param folder The folder to find files in, all in plain one-record-per-line
  #'   .txt format.
  #' @param biomart_data A data.frame with at least the "ensembl_gene_id" and
  #'   "hgnc_symbol" columns, with the correspondence from ENSG to symbol.
  #'   Such a table can be retrieved from Biomart with biomaRt.
  #' @returns A list of vectors, each with the gene symbols of that gene set

  files <- list.files(folder, full.names = TRUE, recursive = TRUE)
  
  # Remove the "all.txt" file
  files <- files[! endsWith(files, "all.txt")]

  data <- list()
  for (file in files) {
    # We need the "clean" file id,
    # e.g. '/a/b/whole_transportome/data.txt' -> '/whole_transportome'
    file |> 
      str_remove("\\/data\\.txt$") |>
      str_remove(paste0("^", folder)) -> id
    
    data[[ id ]] <- read_table(file, col_names = "ensg")[["ensg"]]
  }

  filter_values <- reduce(data, c)
  filter_values <- unique(filter_values)

  # Convert from ENSG to gene symbol
  biomart_data |> select(all_of(c("ensembl_gene_id", "hgnc_symbol"))) -> biomart_data
  ensg_to_symbol <- function(ensgs) {
    symbols <- biomart_data$hgnc_symbol[biomart_data$ensembl_gene_id %in% ensgs]

    return(symbols)
  }

  data <- lapply(data, ensg_to_symbol)

  return(data)
}

extract_ranks <- function(deg_file, biomart_data) {
  #' Make a frame ready for GSEA from a DEG file, made by the cohen's D calculator.
  #' 
  #' This will filter out the non-coding genes as well as add the gene symbols
  #' from BioMart instead of ENSGs.
  #'
  #' @param deg_file (full) Path to the DEG file that needs to be extracted, as .csv.
  #' @param biomart_data A data.frame with at least the "ensembl_gene_id" and
  #'   "gene_biotype" columns.
  #'   Such a table can be retrieved from Biomart with biomaRt.
  #'
  #' @returns A list of named vector of gene_names : statistic, ready for fgsea::fgsea
  data <- read_csv(deg_file, show_col_types = FALSE)

  # Assert that we find a gene_id col
  stopifnot("gene_id" %in% colnames(data))

  # Get rid of the ENSG version -- this causes some rows to collide.
  # We need to get rid of the collisions in the original data
  # - Compute the would-be genes
  data["gene_id"] <- gsub("\\.[0-9]+?$", "", data[["gene_id"]], perl = TRUE)
  data <- na.omit(data)
  
  data |> column_to_rownames("gene_id") -> data
  
  # Avoid duplicated row-names
  biomart_data |> distinct(ensembl_gene_id, .keep_all = TRUE) -> biomart_data
  biomart_data |> column_to_rownames("ensembl_gene_id") -> biomart_data
  
  # Keep only ENSGs that have a gene symbol AND are coding
  biomart_data |>
    dplyr::filter(hgnc_symbol != "") |>
    dplyr::filter(gene_biotype == "protein_coding") -> biomart_data
  
  # Drop the lines in the larger data that do not have a match in the biomart data
  data <- data[rownames(data) %in% rownames(biomart_data), ]
  
  # Make a static gene_symbols that we can slap onto every vector
  gene_symbols <- biomart_data[row.names(data), "hgnc_symbol"]
  
  result <- as.list(data)
  
  for (i in seq_along(result)) {
    names(result[[i]]) <- gene_symbols
  }

  return(result)
}


#' Run GSEA with many genesets on some data
#'
#' This is a tiny wrapper for future compatibility in case we need to change
#' from fgsea.
#'
#' @param genesets A list of genesets to check, with each geneset a vector of
#'   gene names.
#' @param ranks A list of named vector of gene names: statistic to use as ranked list for gsea.
#'
#' @returns A table with GSEA results. See fgsea:fgsea for details.
run_gsea <- function(genesets, ranks) {
  result <- fgsea::fgsea(
    pathways = genesets,
    stats = ranks,
    nPermSimple = 1000,
    nproc=8
  )
  
  ## TODO: FIXME!
  # For some reason we get NAs if we calculate cohen's D that we did not use to
  # get with the t-stat. We are not sure why this is.
  # See the issue: https://github.com/CMA-Lab/transportome_profiler/issues/1
  
  # We band-aid fix this by setting all pval and padj that are NA to 1
  # and all NES and log2err to 0.
  
  result[is.na(result$pval), "pval"] <- 1
  result[is.na(result$padj), "padj"] <- 1
  result[is.na(result$NES), "NES"] <- 0
  result[is.na(result$log2err), "log2err"] <- 0

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


#' Run GSEA on all DEG tables in a folder, with all genesets from another folder.
#'
#' @param input_rank_matrix (Full) path to the input rank matrix made by 
#'   `cohen_calculator`. Loaded by `extract_ranks`.
#' @param output_dir The output directory to save output files to. If NA, does
#'   not save files, and instead returns a list of results.
#' @param genesets_folder_patdh (Full) path to the folder with genesets, as .txt
#'   files with one gene id per row.
#' @param biomart_data A data.frame with at least the "ensembl_gene_id",
#'   "hgnc_symbol" and "gene_biotype" columns.
#'   Such a table can be retrieved from Biomart with biomaRt.
#'
#' @returns A list of values with file names as names and GSEA results as values.
run_all_gsea <- function(input_rank_matrix, genesets_folder_path, biomart_data, output_dir = NA) {
  cat("Loading genesets...\n")
  genesets <- load_genesets(genesets_folder_path, biomart_data = biomart_data)
  
  cat("Loading ranks...\n")
  ranks <- extract_ranks(input_rank_matrix, biomart_data)

  results <- list()
  for (i in seq_along(ranks)) {
    current_name <- names(ranks)[i]
    cat(paste0("Running GSEA on ", current_name, "\n"))

    if (! is.na(output_dir)) {
      result <- run_gsea(genesets, ranks[[i]])

      cat(paste0("Saving data to ", paste0(file.path(output_dir, fnames(ranks)[i]), ".csv"), "\n"))
      save_result(result, output_dir, current_name)
    } else {
      results[[current_name]] <- run_gsea(genesets, ranks[[i]])
      results[[paste0("plot_", current_name)]] <- plot_gsea(genesets, ranks[[i]])
    }
  }

  if (length(results) > 0) {
    return(results)
  }
  return(NULL)
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

#' Save all results from `run_all_gsea` to a series of files.
#'
#' Detects plots to save from their name in the list, starting with `plot`.
#' Filenames are taken from the corresponding list names.
#'
#' @param results The list generated by `run_all_gsea`.
#' @param out_dir The output directory to save to.
#' @param skip_plots If TRUE, does not save plots to output directory.
#'
#' @returns NULL
save_results <- function(results, out_dir, skip_plots = FALSE) {
  wrap <- function(x, name) {
    cat(paste0("Saving ", name, "..."))

    is_plot <- startsWith(name, "plot")
    if (is_plot & skip_plots) {
      cat(" .. Skipped\n")
      return()
    }
    save_result(x, out_dir, name, is_plot)
    cat(".. OK\n")
  }
  # I can't make it work with sapply so, get a for loop
  for (i in seq_along(results)) {
    wrap(results[[i]], names(results)[i])
  }
}

if (FALSE){ # LOCAL DEBUGGING RUN ONLY

embl <- biomaRt::useEnsembl(biomart = "genes")
hs.embl <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl", mart = embl)
ensg_data <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
  mart = hs.embl
)

results <- run_all_gsea(
  "/home/hedmad/Files/data/transportome_profiler/cohen_d_matrix.csv",
  "/home/hedmad/Files/data/transportome_profiler/genesets/",
  ensg_data
)

save_results(results, out_dir = "/home/hedmad/Files/data/mtpdb/gsea_output/", skip_plots = TRUE)

} # ----------------------------------------------------------------------------

# If you are running this from RStudio, you can skip this >>>>>>>>>>>>>>>>>>
if (sys.nframe() == 0L) {

  if (args$low_memory && args$save_plots) {
    cat("WARNING: Low memory mode. Cannot save plots!")
  }

  embl <- biomaRt::useEnsembl(biomart = "genes")
  hs.embl <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl", mart = embl)
  ensg_data <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
    mart = hs.embl
  )

  if (args$low_memory) {
    run_all_gsea(
      args$input_expr_matrix,
      args$input_genesets_folder,
      ensg_data,
      output_dir = args$output_dir
    )
  } else {
    results <- run_all_gsea(
      args$input_expr_matrix,
      args$input_genesets_folder,
      ensg_data
    )

    save_results(results, out_dir = args$output_dir, skip_plots = !args$save_plots)
  }
}
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
