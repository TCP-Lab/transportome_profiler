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
            "--extra-title",
            help = "Extra title string to include in the title",
            default = NA,
        ) |>
        argparser::add_argument(
            "--res",
            help = "Resolution of plot, in pixels per inch.",
            default = 300, type = "numerical"
        ) |>
        argparser::add_argument(
            "--png",
            help = "Save png plot instead of pdf",
            flag = TRUE
        ) |>
        argparser::add_argument(
            "--width",
            help = "Plot width, in inches.",
            default = 7, type = "numerical"
        ) |>
        argparser::add_argument(
            "--height",
            help = "Plot height, in inches.",
            default = 3.5, type = "numerical"
        ) -> parser
    
    args <- argparser::parse_args(parser)
} else {
    args <- list()
    args$input_results <- "./data/suppressed_merged_deas.csv"
    args$selected_genes <- "./data/filter_genes_small.txt"
    args$output_file <- "./data/out/figures/test_commonality_vs_disregulation.png"
    args$extra_title <- "some_extra_title"
    args$simple <- TRUE
    args$width <- 8
    args$height <- 8
    args$res <- 300
    args$png <- TRUE
}

generate_plot_data <- function(data) {
    # We need a frame with ENSG, # of NA and mean of row
    
    pdata <- data
    pdata$na_num <- pdata |>
        select(where(is.numeric)) |>
        apply(1, \(x) {sum(is.na(x))})
    
    pdata <- pdata |> select(sample, na_num)
    
    pdata <- pdata |> group_by(na_num) |>
        mutate(count = n())
    
    pdata <- pdata |> select(count, na_num) |> unique()
    
    pdata
}

create_plot <- function(plot_data, title = NA) {
    
    # Pretty dynamic labels
    label_is_high <- plot_data$count > max(plot_data$count) * 0.9
    # This is not the actual length in say, pixels, but if the labels are short,
    # you can assume each character takes about the same space on screen,
    # So you just multiply by some constant the number of characters in a string.
    label_length <- str_length(plot_data$count)
    
    p <- ggplot(plot_data, aes(x = non_na_num, y = count)) +
        geom_bar(stat = "identity", fill = "orange") +
        geom_text(
            aes(
                label = count,
                y = count + ifelse(label_is_high, -0.05, 0.05) * label_length * max(count)
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
        ggtitle(title) +
        guides(size = "none") +
        xlab("Number of tissue types gene is expressed in") +
        ylab("Number of genes") +
        theme(panel.grid.minor.x = element_blank())
    
    p
}

main <- function(args) {
    data <- read_csv(args$input_results)
    
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
        paste0("Commonality plot - ", args$extra_title)
    } else {
        "Commonality plot"
    }

    # Drop all rows where the ID col is NA
    data <- data[!is.na(data[[args$id_col]]), ]
    
    if (!is.null(selected_genes)) {
        cat("Filtering to genes of interest...\n")
        col <- args$id_col
        data |> filter(data[[col]] %in% selected_genes) -> data
    }
    
    plot_data <- generate_plot_data(data)
    # We swap the number of NAs to the non-number of NAs
    plot_data <- plot_data |> mutate(
        non_na_num = max_number_of_nas - na_num
    ) |> filter(non_na_num != 0) # Some rows have all NAs, so they make no sense to be plotted 

    p <- create_plot(plot_data, title = plot_title)
    
    cat(paste0("Saving plot to ", args$output_file, "\n"))
    
    if (args$png) {
        png(filename = args$output_file, width=args$width, height = args$height, units = "in", res=args$res)
    } else {
        pdf(file = args$output_file, width=args$width, height = args$height)
    }
    print(p)
    dev.off()
}

main(args)

