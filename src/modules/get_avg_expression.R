library(tidyverse)

# This is not really configurable but it's so small

args <- commandArgs(trailingOnly = TRUE)

expr_matrix_path <- args[1]
ensg_data_path <- args[2]
output_data_path <- args[3]
expression_threshold <- args[4] |> as.numeric()

purge_ensg_versions <- function(data, id_col = "gene_id") {
  data |> mutate("{id_col}" := str_remove(.data[[id_col]], "\\.[0-9]+$")) -> data
  data |> distinct(.data[[id_col]], .keep_all = TRUE) -> data
  data
}

print(paste0("Reading in matrix ", expr_matrix_path))
expr_matrix <- read_csv(expr_matrix_path) |> purge_ensg_versions(id_col = "sample")

print(paste0("Reading in ENSG data ", ensg_data_path))
ensg_data <- read_csv(ensg_data_path)

# Filter coding only

coding_genes <- ensg_data |>
  filter(gene_biotype == "protein_coding") |>
  select(ensembl_gene_id) |>
  unlist()

expr_matrix_coding <- expr_matrix |> filter(sample %in% coding_genes)

print(head(expr_matrix_coding))

means <- expr_matrix_coding |>
  select(!sample) |>
  mutate(across(everything(), \(x) {
    x > expression_threshold
  })) |>
  summarise(across(everything(), \(x) {
    sum(x, na.rm = TRUE)
  })) |>
  mutate(across(everything(), \(x) {
    x / nrow(expr_matrix_coding)
  }))

print(paste0("Writing output means ", output_data_path))
write_csv(means, output_data_path)
