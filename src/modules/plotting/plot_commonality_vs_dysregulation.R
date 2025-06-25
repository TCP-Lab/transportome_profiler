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
            default = 16, type = "numerical"
        ) |>
        argparser::add_argument(
            "--height",
            help = "Plot height, in inches.",
            default = 25, type = "numerical"
        ) -> parser
    
    args <- argparser::parse_args(parser)
} else {
    args <- list()
    args$input_results <- "./data/suppressed_merged_deas.csv"
    args$ensg_to_hugo <- "./data/ensg_data.csv"
    args$selected_genes <- "./data/filter_genes.txt"
    args$absolute <- TRUE
    args$id_col <- "sample"
    args$extra_title <- "some_extra_title"
    args$png <- TRUE
}

generate_plot_data <- function(data, absolute = FALSE) {
    # We need a frame with ENSG, # of NA and mean of row
    
    pdata <- data
    pdata$na_num <- pdata |> select(where(is.numeric)) |> apply(1, \(x) {sum(is.na(x))})

    pdata |> melt(id.vars = c("sample", "na_num"))
}

create_plot <- function(plot_data, plot_y_limits, title) {
    plot_data <- plot_data |> group_by(non_na_num) |> mutate(count = sum(!is.na(value)))
    p <- ggplot(plot_data, aes(x = non_na_num, y = value)) +
        geom_hline(yintercept = 0, linewidth = 1, color = "gray") +
        geom_jitter(
            width = 0.25,
            alpha = 0.25,
            color = "blue",
            na.rm = TRUE
        ) +
        scale_color_gradient(low = "black", high = "red", name = "St. Dev.") +
        geom_label(
            aes(x = non_na_num, y = plot_y_limits[2] + 3, label = count),
            angle = 90
        ) +
        geom_point(
            aes(x = non_na_num, y = plot_y_limits[2] + 1, size = count),
            color = "red"
        ) +
        theme_minimal() +
        ggtitle(title) +
        scale_x_continuous(
            breaks = seq(1, 19, by = 1)
        ) +
        scale_y_continuous(
            breaks = seq(plot_y_limits[1], plot_y_limits[2], by = 2),
            limits = c(plot_y_limits[1], plot_y_limits[2] + 4)
        ) +
        guides(size = "none") +
        theme(panel.grid.minor.x = element_blank(), legend.position = "right")
    
    p
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
    
    all_plots <- list()

    for (value in unique(plot_data$variable)) {
        this_plot_data <- plot_data[plot_data$variable == value, ]
        plot_y_limits <- c(floor(min(this_plot_data$value, na.rm = TRUE)), ceiling(max(this_plot_data$value, na.rm = TRUE)))
        all_plots[[value]] <- create_plot(this_plot_data, plot_y_limits = plot_y_limits, title = value)
    }
    
    p <- {all_plots |> reduce(\(x, y) {x + y})} + plot_layout(ncol = 3)
    p <- p + xlab("Number of tissues gene is expressed in") +
        ylab("Mean value of metric") +
        ggtitle(plot_title)
    
    print(p)
    
    if (args$png) {
        png(filename = args$output_file, width=args$width, height = args$height, units = "in", res=args$res)
    } else {
        pdf(file = args$output_file, width=args$width, height = args$height)
    }
    print(p)
    dev.off()
}

main(args)

