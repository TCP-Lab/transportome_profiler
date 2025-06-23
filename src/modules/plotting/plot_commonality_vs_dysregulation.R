options(warn = 1)

suppressMessages({
    options(tidyverse.quiet = TRUE)
    library(tidyverse)
    library(assertthat)
    library(reshape2)
    library(patchwork)
    extrafont::loadfonts()
    requireNamespace("gridExtra")
})

{
    # This removes the version from an ensembl ID
    strip_ensg_version <- function(x) {
        str_split_1(x, "\\.")[1]
    }

    assert_that(are_equal(strip_ensg_version("ENSG0000.12"), "ENSG0000"))
    assert_that(are_equal(strip_ensg_version("ENSG0000.1"), "ENSG0000"))
}


if (!exists("LOCAL_DEBUG")) {
    # Parsing arguments
    requireNamespace("argparser")
    
    parser <- argparser::arg_parser("Plot a shared dysregulation plot")
    
    parser |>
        argparser::add_argument(
            "output_file",
            help = "Path to the output file",
            type = "character"
        ) |>
        argparser::add_argument(
            "input_results",
            help = "Input results to parse", type = "character"
        ) |>
        argparser::add_argument(
            "ensg_to_hugo",
            help = "Map from ENSGs to HUGO symbols",
            type = "character",
            default = NULL,
        ) |>
        argparser::add_argument(
            "--selected_genes",
            help = "File with comma-separated list of ENSGs to select before plotting",
            type = "character",
            default = NULL,
        ) |>
        argparser::add_argument(
            "--id_col",
            help = "Name of the column holding ENSG IDs",
            type = "character",
            default = "sample",
        ) |>
        argparser::add_argument(
            "--absolute",
            help = "Plot aboslute dysregulation instead of signed one",
            flag = TRUE,
        ) |>
        argparser::add_argument(
            "--extra-title",
            help = "Extra title string to include in the title",
            default = NA,
        ) |>
        argparser::add_argument(
            "--res",
            help = "Resolution of plot, in pixels per inch.",
            default = 400, type = "numerical"
        ) |>
        argparser::add_argument(
            "--png",
            help = "Save png plot instead of pdf",
            flag = TRUE
        ) |>
        argparser::add_argument(
            "--width",
            help = "Plot width, in inches.",
            default = 9, type = "numerical"
        ) |>
        argparser::add_argument(
            "--height",
            help = "Plot height, in inches.",
            default = 6, type = "numerical"
        ) -> parser
    
    args <- argparser::parse_args(parser)
} else {
    args <- list()
    args$input_results <- "./data/suppressed_merged_deas.csv"
    args$ensg_to_hugo <- "./data/ensg_data.csv"
    args$selected_genes <- "./data/filter_genes.txt"
    args$absolute <- FALSE
    args$id_col <- "sample"
    args$extra_title <- "some_extra_title"
    args$png <- TRUE
}

smart_mean <- function(x, absolute = FALSE) {
    new_x <- c()
    for (i in x) {
        if (is.na(i)) {
            next
        }
        
        if (!is.numeric(i)) {
            next
        }
        
        new_x <- c(new_x, i)
    }
    
    if (length(new_x) == 0) {
        return(NA)
    }
    
    if (absolute) {
        mean(abs(new_x))
    } else {
        mean(new_x)
    }
}

generate_plot_data <- function(data, absolute = FALSE) {
    # We need a frame with ENSG, # of NA and mean of row
    
    pdata <- data |> rowwise() |> mutate(
        na_num = sum(is.na(across(where(is.numeric)))),
        mean = smart_mean(across(where(is.numeric)), absolute = absolute),
        ensg = sample,
        .keep = "none"
    )

    pdata
}


main <- function(args) {
    data <- read_csv(args$input_results)
    ensg_to_hugo <- read_csv(args$ensg_to_hugo)
    selected_genes <- if (!is.null(args$selected_genes)) {
        read_file(args$selected_genes) |>
            str_split_1(",") |> str_remove_all("\"") |> str_remove_all("\\n")
    } else {
        NULL
    }
    
    data[[args$id_col]] <- map(data[[args$id_col]], strip_ensg_version) |> unlist()
    
    # The # of NAs maximum to plot is -1 (for the id col)
    max_number_of_nas <- ncol(data) - 1
    plot_title <- if (!is.na(args$extra_title)) {
        paste0("Commonality vs disregulation plot - ", args$extra_title)
    } else {
        "Commonality vs disregulation plot"
    }
    
    if (args$absolute) {
        plot_title <- paste0(plot_title, " (Absolute)")
    }
    
    # Drop all rows where the ID col is NA
    data <- data[!is.na(data[[args$id_col]]), ]
    
    if (!is.null(selected_genes)) {
        col <- args$id_col
        data |> filter(data[[col]] %in% selected_genes) -> data
        print(dim(data))
    }
    
    plot_data <- generate_plot_data(data, absolute = args$absolute)
    # We swap the number of NAs to the non-number of NAs
    plot_data <- plot_data |> mutate(
        non_na_num = max_number_of_nas - na_num
    ) |> filter(non_na_num != 0) # Some rows have all NAs, so they make no sense to be plotted 
    
    p <- ggplot(plot_data, aes(x = non_na_num, y = mean)) +
        geom_boxplot(
            aes(x = factor(non_na_num)),
            outlier.alpha = 0.5,
            outlier.colour = "darkred",
            fill = "orange"
        ) +
        theme_minimal() +
        xlab("Number of tissues gene is expressed in") +
        ylab("Algebraic mean value of metric") +
        ggtitle(plot_title)
    
    if (args$png) {
        png(filename = args$output_file, width=args$width, height = args$height, units = "in", res=args$res)
    } else {
        pdf(file = args$output_file, width=args$width, height = args$height)
    }
    print(p)
    dev.off()
}

main(args)
