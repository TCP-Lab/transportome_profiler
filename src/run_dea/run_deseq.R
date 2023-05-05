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
  read_csv(case_path) |> column_to_rownames(case_rownames_col) -> case_data
  read_csv(control_path) |> column_to_rownames(control_rownames_col) -> control_data 

  assertthat::assert_that(is.numeric(case_data), msg = "The case matrix must be numeric")
  assertthat::assert_that(is.numeric(control_data), msg = "The control matrix must be numeric")
  
  result <- run_deseq(case_data, control_data)
  
  write_csv(result, output_path)
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
    status = c(rep("case", ncol(case_frame)), rep("control", ncol(control_frame)))
  )
  
  # Merge the case and control frames
  data <- rowmerge(case_frame, control_frame, all = FALSE) # inner merge
  
  # Assert that the order of the cols and rows is the same
  data <- data[, row.names(metadata)]
  
  DESeq2::DESeqDataSetFromMatrix(
    data, metadata, ~ status
  ) |>
    DESeq2::DESeq(dsq_data) |>
    DESeq2::results(contrast = c("status", "case", "control")) -> dsq_res
  
  dsq_res
}

main(
  case_path = args$case_csv_file, control_path = args$control_csv_file,
  output_path = args$output_path,
  case_rownames_col = args$case_rownames_col,
  control_rownames_col = args$control_rownames_col
)

