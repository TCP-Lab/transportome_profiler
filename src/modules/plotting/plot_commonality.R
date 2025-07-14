#!/usr/bin/env Rscript
options(warn = 1)

if (! exists("LOCAL_DEBUG")) {
    # Parsing arguments
    requireNamespace("argparser")

    parser <- argparser::arg_parser("Plot commonality plots")

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
            default = 7, type = "numerical"
        ) |>
        argparser::add_argument(
            "--renames", help = "JSON file with renames to apply to sample names when plotting",
            default = NA, type = "character"
        ) |>
        argparser::add_argument(
            "--height", help = "Plot height, in inches.",
            default = 3, type = "numerical"
        ) -> parser

    args <- argparser::parse_args(parser)
}

suppressMessages({
    options(tidyverse.quiet = TRUE)
    library(tidyverse)
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

    freqs <- UpSetR::fromList(plot_data) |> rowSums() |> table()
    
    tibble(
        freq = as.vector(freqs),
        nums = seq_along(freqs)
    )
}

main <- function(
    input_means_path,
    input_set_path,
    output_file_path,
    expression_threshold,
    extra_title = NA,
    renames = NA,
    res = 300,
    width = 10,
    height = 10,
    png = FALSE
) {
    input_means <- read_csv(input_means_path)
    input_set <- jsonlite::read_json(input_set_path) |> unlist()

    if (!is.na(renames)) {
        renames <- jsonlite::read_json(renames)
    } else {
        renames <- NA
    }

    plot_data <- gen_plot_data(input_means, input_set, expressed_threshold = expression_threshold)
    
    plot_title <- if (!is.na(extra_title)) {
        paste0("Commonality plot - ", extra_title)
    } else {
        "Commonality plot"
    }

    apply_renames <- function(x) {
        if (is.na(renames)) {
            return(x)
        }

        if (!x %in% names(renames)) {
            return(x)
        }

        renames[[x]]
    }
    
    # Pretty dynamic labels
    label_is_high <- plot_data$freq > max(plot_data$freq) * 0.9
    # This is not the actual length in say, pixels, but if the labels are short,
    # you can assume each character takes about the same space on screen,
    # So you just multiply by some constant the number of characters in a string.
    label_length <- str_length(plot_data$freq)

    p <- ggplot(plot_data, aes(x = nums, y = freq)) +
        geom_bar(stat = "identity", fill = "orange") +
        geom_text(
            aes(
                label = freq,
                y = freq + ifelse(label_is_high, -0.05, 0.05) * label_length * max(freq)
            ),
            angle = 90,
            vjust = 0.5, # This centers the labels *horizontally* (90 angle)
            hjust = ifelse(label_is_high, 0, 1)
        ) +
        scale_x_continuous(
            breaks = seq(1, 19, by = 1),
            labels = seq(1, 19, by = 1),
        ) +
        theme_minimal() +
        ggtitle(plot_title) +
        guides(size = "none") +
        xlab("Number of tissue types gene is expressed in") +
        ylab("Number of genes") +
        theme(panel.grid.minor.x = element_blank())
    
    if (png) {
        png(output_file_path, width = width, height = height, res = res, units = "in")
    } else {
        pdf(output_file_path, width = width, height = height)
    }
    print(p)
    dev.off()
}


if (exists("LOCAL_DEBUG")) {
    main(
        input_means_path = "data/expression_means.csv",
        input_set_path = "data/test_channels.json",
        output_file_path = "data/out/figures/whole_transportome_test_upset.png",
        expression_threshold = 0,
        extra_title = "Transporters",
            renames = "data/in/config/tcga_renames.json",
        res = 300,
        width = 7,
        height = 4,
        png = TRUE
    )
} else {
    main(
        input_means_path = args$input_expression_means,
        input_set_path = args$set,
        output_file_path = args$output_file,
        expression_threshold = args$expression_threshold,
        extra_title = args$extra_title,
        renames = args$renames,
        res = args$res,
        width = args$width,
        height = args$height,
        png = args$png
    )
}

