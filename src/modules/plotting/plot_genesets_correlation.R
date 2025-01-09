#!/usr/bin/env Rscript
options(warn = 1)

if (! exists("LOCAL_DEBUG")) {
    # Parsing arguments
    requireNamespace("argparser")
    
    parser <- argparser::arg_parser("Plot a large correlation matrix with correlations between genesets.")
    
    parser |>
        argparser::add_argument(
            "genesets", help="The JSON file with the genesets", type="character"
        ) |>
        argparser::add_argument(
            "input_tree", help="A txt file with the folder structure of the tree to be parsed", type="character"
        ) |>
        argparser::add_argument(
            "output_file", help = "The path to the output file",
            type = "character"
        ) |>
        argparser::add_argument(
            "--extra_title", help = "Extra title to add to the figure",
            type = "character", default = NULL
        ) |>
        argparser::add_argument(
            "--res", help = "Resolution of plot, in pixels per inch.",
            default= 400, type = "numerical"
        ) |>
        argparser::add_argument(
            "--png", help = "Save as PNG istead of PDF.",
            flag = TRUE
        ) |>
        argparser::add_argument(
            "--width", help = "Plot width, in inches.",
            default = 15, type = "numerical"
        ) |>
        argparser::add_argument(
            "--height", help = "Plot height, in inches.",
            default = 15, type = "numerical"
        ) -> parser
    
    args <- argparser::parse_args(parser)
}

suppressMessages({
    options(tidyverse.quiet = TRUE)
    library(tidyverse)
    requireNamespace("stringi")
    requireNamespace("reshape2")
    extrafont::loadfonts()
    library(patchwork)
})

parse_tree_labels <- function(tree, genesets) {
    #' Parse the raw tree from the file to an usable labels dataframe
    #'
    #' @param tree A vector of lines where each line is from the tree output
    #' @param genesets The JSON geneset data
    
    result <- data.frame(original = tree, order = rev(seq_along(tree)))
    
    remove_backbone <- function(x) {
        # Remove all backbone chars from an input
        str_remove_all(x, "[└─│├]") |> str_trim()
    }
    result$id <- remove_backbone(result$original)
    
    result$backbone <- str_split_i(result$original, "─ ", 1) |> paste0("─ ")
    result$backbone[1] <- ""
    
    result$rev_backbone <- stringi::stri_reverse(result$backbone) |>
        str_replace_all("├", "┤") |>
        str_replace_all("└", "┘")
    
    result$reverse <- paste0(result$id, result$reverse_backbone)
    
    result$pretty <- paste0(
        result$backbone,
        sapply(result$id, \(id) {genesets[[id]]$name})
    )
    result$rev_pretty <- paste0(
        sapply(result$id, \(id) {genesets[[id]]$name}),
        result$rev_backbone
    )
    
    result
}

# This is the same as make_genesets.py, the Sorens-Dice coefficient.
calc_correlation <- function(set_1, set_2) {
    set_1 <- unique(set_1)
    set_2 <- unique(set_2)
    
    2 * (length(intersect(set_1, set_2)) / (length(set_1) + length(set_2))) 
}

gen_plot_data <- function(
        genesets, tree
) {
    # Generate the two-way correlations of the genesets
    combinations <- expand.grid(x = names(genesets), y = names(genesets))
    combinations$value <- apply(combinations, 1, function(x) {
        calc_correlation(genesets[[x[1]]]$data, genesets[[x[2]]]$data)
    })
    
    combinations
}

tile_size <- function(x, min_size = 0.5,  max_size = 1) {
    size <- scales::rescale(x, from = c(0, 1), to = c(min_size, max_size))
    size[x == 0] <- 0
    size
}

create_correlation_heatmap <- function(
        plot_data,
        labels,
        extra_title = NULL
) {
    fig_title <- if (!is.null(extra_title)) {
        paste0("Genesets Congruency - ", extra_title)
    } else {
        "Geneset Congruency"
    }
    
    get_label <- function(x) {
        sapply(x, \(x) {
            labels$rev_pretty[labels$id == x]
    })
    }
    
    p <- ggplot(plot_data, aes(fill = value, x = factor(x, levels = rev(labels$id[labels$order])), y = factor(y, levels = rev(labels$id[labels$order])))) +
        geom_tile(aes(width = tile_size(value, min_size = 0.3), height = tile_size(value, min_size = 0.3))) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1)) +
        ggtitle(fig_title) +
        ylab("Gene Set") + xlab("Gene Set") +
        scale_fill_gradientn( # not a typo - it is really called gradientn
            "Congruency",
            colours = c("purple", "skyblue", "lightgoldenrod", "darkorange")
        ) +
        scale_x_discrete(labels = get_label) +
        scale_y_discrete(labels = get_label, limits = rev) +
        theme(
            text = element_text(family = "FiraCode Nerd Font", size = 10),
            panel.grid = element_blank(),
            panel.grid.major = element_line(colour = scales::alpha("black", 0.05))
        )
    
    p
}

main <- function(
        input_tree,
        genesets_file,
        out_file,
        extra_title = NULL,
        save_png = FALSE,
        png_res = 300,
        plot_width = 10,
        plot_height = 6
) {
    cat("Reading in genesets\n")
    genesets <- jsonlite::fromJSON(read_file(genesets_file))
    tree <- read_lines(input_tree)
    
    plot_data <- gen_plot_data(
        tree = tree,
        genesets = genesets
    )
    
    labels <- parse_tree_labels(tree, genesets)
    large_plot <- create_correlation_heatmap(plot_data, labels, extra_title)
    
    # Save plot to output
    if (is.null(out_file)) {
        print(large_plot)
        invisible()
    } else {
        if (save_png) {
            png(filename = out_file, width = plot_width, height = plot_height, units = "in", res = png_res)
        } else {
            pdf(file = out_file, width = plot_width, height = plot_height)
        }
        print(large_plot)
        dev.off()
    }
}

if (exists("LOCAL_DEBUG")) {
    main(
        input_tree = "~/Files/repos/tprof/data/genesets_repr.txt",
        genesets_file = "~/Files/repos/tprof/data/genesets.json",
        out_file = NULL,
        save_png = TRUE,
        png_res = 500,
        plot_height = 10
    )
} else {
    main(
        input_tree = args$input_tree,
        genesets_file = args$genesets,
        out_file = args$output_file,
        save_png = args$png,
        png_res = args$res,
        plot_width = args$width,
        plot_height = args$height
    )
}
