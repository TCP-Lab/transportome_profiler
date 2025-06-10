#!/usr/bin/env Rscript
options(warn = 1)

if (!exists("LOCAL_DEBUG")) {
    # Parsing arguments
    requireNamespace("argparser")

    parser <- argparser::arg_parser("Plot a large dotplot with number of expressed genes per genesets.")

    parser |>
        argparser::add_argument(
            "file_to_suppress",
            help = ".csv file with genes as rows and samples as columns to suppress", type = "character"
        ) |>
        argparser::add_argument(
            "expression_values",
            help = ".csv file with genes as rows and samples as columns to source expression values from. Must coincide with the shape of the file to suppress, and should have the same column and row names", type = "character"
        ) |>
        argparser::add_argument(
            "output_file",
            help = "The path to the output file",
            type = "character"
        ) |>
        argparser::add_argument(
            "--suppress-id-col",
            help = "The name of the column with row IDs in the file to suppress",
            default = "sample",
            type = "character"
        ) |>
        argparser::add_argument(
            "--expression-id-col",
            help = "The name of the column with row IDs in the file with expression values",
            default = "sample",
            type = "character"
        ) |>
        argparser::add_argument(
            "--expression-threshold",
            help = "Expression threshold to filter 'expressed' vs 'non expressed' with, > abs(thr).",
            default = 0, type = "numerical"
        ) -> parser

    args <- argparser::parse_args(parser)
}

suppressMessages({
    options(tidyverse.quiet = TRUE)
    library(tidyverse)
    library(assertthat)
})

suppress_data <- function(data, expression_data, threshold) {
    # Check if the shape, colnames and rownames of the data are the same.
    assert_that(all(dim(data) == dim(expression_data)), msg = "Dimensions of the two objects are not the same")

    # Since the dimensions are the same, we can just check if A %in% B, without having
    # to check if B %in% A since both are the same size.
    assert_that(all(row.names(data) %in% row.names(expression_data)))
    assert_that(all(colnames(data) %in% colnames(expression_data)))

    # Everything is in place - we can start suppressing the cells that need it
    # First, make sure that the columns and rows are in the same order
    expression_data <- expression_data[row.names(data), colnames(data)]

    suppress_mask <- expression_data |> mutate(across(everything(), \(x) {
        res <- rep(1, length(x))
        res[x <= threshold] <- NA
        res
    }))

    res <- mapply(`*`, data, suppress_mask)
    # We lost the rownames, and the process made the data a matrix
    res <- as.data.frame(res)
    row.names(res) <- row.names(data)

    res
}

{
    # Unit tests
    # Note - at last one row in each col must be retained or
    # identical() complains that the all-NA col is now "logical" instead of
    # "numeric" as it should
    tdata <- data.frame(
        row.names = c("A", "B", "C", "D"),
        col_a = c(1, 2, 3, 4),
        col_b = c(1, 2, 3, 4),
        col_c = c(1, 2, 3, 4)
    )
    edata <- data.frame(
        row.names = c("A", "B", "D", "C"), # C and D are swapped
        col_a = c(10, 10, 10, 10),
        col_c = c(10, 0.5, 10, 0.5), # col c and b are swapped
        col_b = c(10, 0.1, 0.3, 0.9)
    )
    res <- data.frame(
        row.names = c("A", "B", "C", "D"),
        col_a = c(1, 2, 3, 4),
        col_b = c(1, NA, NA, NA),
        col_c = c(1, NA, NA, 4)
    )
    assert_that(identical(suppress_data(tdata, edata, 1), res))
}

main <- function(args) {
    data <- read_csv(args$file_to_suppress) |> column_to_rownames(args$suppress_id_col)
    expression_data <- read_csv(args$expression_values) |> column_to_rownames(args$expression_id_col)

    result <- suppress_data(data, expression_data, args$expression_threshold) |> rownames_to_column(args$suppress_id_col)

    result |> write_csv(args$output_file)
}

main(args)
