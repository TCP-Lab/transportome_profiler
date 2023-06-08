# Run DESeq2 on two .csv files, one with the case and one with the controls

## >>> Skip this if you are running this from RStudio
# Parsing arguments
requireNamespace("argparser")

parser <- argparser::arg_parser("Run GSEA on DEG tables")

parser |>
  argparser::add_argument(
    "case_csv_file", help=".csv file with the 'case' expression values", type="character"
  ) |>
  argparser::add_argument(
    "control_csv_file", help = ".csv file with the 'control' expression values",
    type = "character"
  ) |>
  argparser::add_argument(
    "output_csv", help = "Path to the output csv file, including its filename",
    type = "character"
  ) |>
  argparser::add_argument(
    "--case_rownames_col", help = "Name of the column with rownames in the case file",
    type = "character", default = "gene_id"
  ) |>
  argparser::add_argument(
    "--control_rownames_col", help = "Name of the column with rownames in the control file",
    type = "character", default = "gene_id"
  ) -> parser

args <- argparser::parse_args(parser)

# <<<<<<< To here <<<<<<<<<<<<<<<<<,

# Functions definitions
library(tidyverse)
requireNamespace("DESeq2")
requireNamespace("assertthat")

main <- function(case_path, control_path, output_path, case_rownames_col, control_rownames_col) {
  cat("Loading case data...\n")
  read_csv(case_path, show_col_types = FALSE) |> column_to_rownames(case_rownames_col) -> case_data
  cat(paste0(
    "Loaded a ", nrow(case_data), " rows by ", ncol(case_data),
    " columns case expression matrix from '", case_path, "'.\n"
  ))
  cat("Loading control data...\n")
  read_csv(control_path, show_col_types = FALSE) |> column_to_rownames(control_rownames_col) -> control_data 
  cat(paste0(
    "Loaded a ", nrow(control_data), " rows by ", ncol(control_data),
    " columns control expression matrix from '", control_path, "'.\n"
  ))

  cat("Running deseq...\n")
  result <- run_deseq(case_data, control_data)
  
  cat("Saving result...\n")
  write.csv(result, file = output_path, row.names=FALSE)
}

rowmerge <- function(x, y, ...) {
  res <- merge(x, y, by='row.names', ...)
  row.names(res) <- res$Row.names
  res$Row.names <- NULL
  
  res
}

run_deseq <- function(case_frame, control_frame) {
  
  # Make a dummy metadata frame
  metadata <- data.frame(
    row.names = c(colnames(case_frame), colnames((control_frame))),
    status = factor(c(rep("case", ncol(case_frame)), rep("control", ncol(control_frame))), levels = c("control", "case"))
  )
  
  # Merge the case and control frames
  data <- rowmerge(case_frame, control_frame, all = FALSE) # inner merge
  
  # Assert that the order of the cols and rows is the same
  data <- data[, row.names(metadata)]

  # The data must be un-logged
  data <- round((2 ** data) -1)

  DESeq2::DESeqDataSetFromMatrix(
    data, metadata, ~ status
  ) |>
    DESeq2::DESeq() -> dsq_obj

  DESeq2::results(dsq_obj, name = "status_case_vs_control") -> dsq_res
  
  cat("Deseq finished running!\n")
  as.data.frame(dsq_res) |> rownames_to_column("gene_id")
}

main(
  case_path = args$case_csv_file, control_path = args$control_csv_file,
  output_path = args$output_csv,
  case_rownames_col = args$case_rownames_col,
  control_rownames_col = args$control_rownames_col
)

