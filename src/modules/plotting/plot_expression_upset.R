#!/usr/bin/env Rscript
options(warn = 1)

if (! exists("LOCAL_DEBUG")) {
    # Parsing arguments
    requireNamespace("argparser")

    parser <- argparser::arg_parser("Plot a large dotplot with number of expressed genes per genesets.")

    parser |>
        argparser::add_argument(
            "input_expression_means", help=".csv file to read expression means with.", type="character"
        ) |>
        argparser::add_argument(
            "set", help="Shallow JSON file to read a set of genes from, of the type ['gene_1', 'gene_2', ...].", type="character"
        ) |>
        argparser::add_argument(
            "output_file", help = "The path to the output file",
            type = "character"
        ) |>
        argparser::add_argument(
            "--expression_threshold", help = "Expression threshold to filter 'expressed' vs 'non expressed' with.",
            default = 0, type = "numerical"
        ) |>
        argparser::add_argument(
            "--extra_title", help = "Extra title to add to the figure",
            type = "character", default = NULL
        ) |>
        argparser::add_argument(
            "--png", help = "Save PNG plots instead of PDF.",
            flag = TRUE
        ) |>
        argparser::add_argument(
            "--res", help = "Resolution of plot, in pixels per inch.",
            default= 400, type = "numerical"
        ) |>
        argparser::add_argument(
            "--width", help = "Plot width, in inches.",
            default = 10, type = "numerical"
        ) |>
        argparser::add_argument(
            "--renames", help = "JSON file with renames to apply to sample names when plotting",
            default = NA, type = "character"
        ) |>
        argparser::add_argument(
            "--height", help = "Plot height, in inches.",
            default = 10, type = "numerical"
        ) -> parser

    args <- argparser::parse_args(parser)
}

suppressMessages({
    options(tidyverse.quiet = TRUE)
    library(tidyverse)
    library(ComplexUpset)
    requireNamespace("UpSetR")
    requireNamespace("stringi")
    requireNamespace("jsonlite")
})

purge_ensg_versions <- function(data, id_col = "gene_id") {
    data |> mutate("{id_col}" := str_remove(.data[[id_col]], "\\.[0-9]+$")) -> data
    data |> distinct(.data[[id_col]], .keep_all = TRUE) -> data
    data
}


gen_plot_data <- function(
        input_data,
        geneset,
        expressed_threshold = 0,
        id_col = "sample"
) {
    data <- input_data |> purge_ensg_versions(id_col=id_col) |> column_to_rownames(id_col)

    # First, we need to bin the data
    data <- data |> mutate(across(where(is.numeric), \(x) {x > expressed_threshold}))

    # Now, for all the types of tumors (columns != to id_col) we need to check if the
    # genes in the geneset are expressed or not.
    print(length(geneset))
    plot_data <- list()
    for (name in colnames(data)) {
        plot_data[[name]] <- {
            # We take out just this col, and just the rows in the geneset
            subs <- data[geneset, name, drop=FALSE]
            # We can now just take whatever is set to 1 here, and discard the rest
            names <- rownames(subs)[subs[[name]]]
            names
        }
    }

    UpSetR::fromList(plot_data)
}

main <- function(
    input_means_path,
    input_set_path,
    output_file_path,
    expression_threshold,
    extra_title = NA,
    renames = NA,
    res = 300,
    width = 1600,
    height = 600,
    png = FALSE
) {
    input_means <- read_csv(input_means_path)
    input_set <- jsonlite::read_json(input_set_path) |> unlist()

    data <- gen_plot_data(input_means, input_set, expressed_threshold = expression_threshold)

    p <- upset(
        data,
        names(data),
        n_intersections = 15,
        base_annotations=list(
            'Intersection size'=intersection_size(
                text_mapping=aes(label=paste0(
                    !!upset_text_percentage(),
                    '\n',
                    '(',
                    !!get_size_mode('exclusive_intersection'),
                    ')'
                ))
            )
        )
    ) + ggtitle(extra_title)

    print(p)
}


if (exists("LOCAL_DEBUG")) {
    main(
        input_means_path = "data/expression_means.csv",
        input_set_path = "data/test_transp.json",
        output_file_path = "data/out/figures/whole_transportome_test_upset.png",
        expression_threshold = 0,
        extra_title = "Transporters",
        renames = NA,
        res = 300,
        width = 1600,
        height = 600,
        png = TRUE
    )
} else {
    main(
        input_expression_means = args$input_expression_means,
        input_tree = args$input_tree,
        genesets_file = args$genesets,
        out_file = args$output_file,
        no_cluster = args$no_cluster,
        extra_title = args$extra_title,
        save_png = TRUE,
        png_res = args$res,
        plot_width = args$width,
        plot_height = args$height,
        renames = args$renames,
        expressed_threshold = args$expression_threshold
    )
}
