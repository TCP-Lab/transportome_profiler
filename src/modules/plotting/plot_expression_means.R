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
            "--expression_threshold", help = "Expression threshold to filter 'expressed' vs 'non expressed' with.",
            default = 0, type = "numerical"
        ) |>
        argparser::add_argument(
            "--no_cluster", help = "Skip clustering of x axis labels",
            flag=TRUE, type = "logical"
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
    library(reshape2)
    requireNamespace("stringi")
    requireNamespace("reshape2")
    extrafont::loadfonts()
    library(patchwork)
})

purge_ensg_versions <- function(data, id_col = "gene_id") {
    data |> mutate("{id_col}" := str_remove(.data[[id_col]], "\\.[0-9]+$")) -> data
    data |> distinct(.data[[id_col]], .keep_all = TRUE) -> data
    data
}

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

gen_plot_data <- function(
        input_tree,
        input_data,
        genesets,
        renames = NULL,
        expressed_threshold = 0,
        id_col = "sample",
        cluster_x_axis = FALSE
) {
    cat(paste0(
        "Generating plot data from input", input_tree,
        ", reading mean expression data from ", input_data,
        "\n"
    ))
    tree <- read_lines(input_tree)
    labels <- parse_tree_labels(tree, genesets)
    
    data <- read_csv(input_data) |> purge_ensg_versions(id_col=id_col) |> column_to_rownames(id_col)
    
    # For the heatmap we will need a binary matrix with yes/no expression
    # so I binarize the data here
    data <- data |> mutate(across(where(is.numeric), \(x) {x > expressed_threshold}))
    
    # Compute the information about the genesets
    plot_data <- list()
    for (geneset in names(genesets)) {
        plot_data[[geneset]] <- sapply(names(data), \(x) {
            # This is a vector of only the TRUE genes for this tumor type
            genes <- row.names(data)[data[[x]]]
            geneset_genes <- genesets[[geneset]]$data
            
            sum(geneset_genes %in% genes) / length(geneset_genes)
        })
    }
    plot_data <- data.frame(plot_data) |>
        rownames_to_column(id_col) |> melt(id.vars = id_col, variable.name = "pathway")
    
    # The transformations have fucked up the pathway IDs, so I fix them here
    plot_data$pathway <- str_replace_all(plot_data$pathway, "\\.", "-") |> 
        str_replace_all("^X", "") # This is to fix R adding an X if the name starts with a number
    
    # Add the size of the genesets
    plot_data$size <- sapply(plot_data$pathway, \(x) {
        length(genesets[[x]]$data)  
    })
    
    # Set the order of the column samples based on hclust
    # - We need to make a matrix from the input
    # - This is fine to be ran on just one type of data, based on what we want to
    #   sort by. I choose to run it on the relative data.
    if (cluster_x_axis) {
        clust_data <- reshape2::dcast(plot_data, id_col ~ pathway, value.var = "NES")
        clust_data <- clust_data |> column_to_rownames(id_col)
        
        clust <- hclust(dist(clust_data), method = "ward.D")
        
        # Set the order of the labels. The actual values will be set in the plot
        plot_data$fac_id <- factor(plot_data[[id_col]], levels = clust$labels[clust$order])
    } else {
        # If not clustered, resort to alphabetical clustering
        plot_data$fac_id <- factor(plot_data[[id_col]], levels = sort(unique(plot_data[[id_col]])))
    }
    # If we need to rename the factors, do it now
    plot_data$fac_id <- plyr::revalue(plot_data$fac_id, renames)
    
    # Add the information about the labels
    plot_data$label <- sapply(plot_data$pathway, \(x) {
        labels[labels$id == x, "rev_pretty"]
    })
    plot_data$order <- sapply(plot_data$pathway, \(x) {
        labels[labels$id == x, "order"]
    })
    
    plot_data$fac_pathway <- factor(plot_data$pathway,  levels = labels$id[labels$order])
    
    plot_data
}




create_large_heatmap <- function(
        plot_data,
        extra_title = NULL
) {
    fig_title <- if (!is.null(extra_title)) {
        paste0("Geneset Expression Overview - ", extra_title)
    } else {
        "Geneset Expression Overview"
    }
    
    p <- ggplot(plot_data, aes(fill = value, x = fac_id, y = fac_pathway)) +
        geom_tile() +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        ggtitle(fig_title) +
        scale_fill_gradientn( # not a typo - it is really called gradientn
            "Expression\nRate",
            colours = rev(c("purple", "violet", "lightblue", "white")),
            breaks = c(0, 0.25, 0.50, 0.75, 1),
            labels = c("0 %", "25 %", "50 %", "75 %", "100 %"),
            limits = c(0, 1)
        ) +
        ylab("Gene Set") + xlab("Cohort") +
        scale_y_discrete(breaks = plot_data$fac_pathway, labels = plot_data$label) +
        theme(text = element_text(family = "FiraCode Nerd Font", size = 10))
    
    p
}

create_set_size_plot <- function(plot_data) {
    # Select just one category - they all share the same set size
    filtered_plot_data <- plot_data |> filter(fac_id == levels(plot_data$fac_id)[1])
    filtered_plot_data$fac_id <- "Set\nsize"
    filtered_plot_data$fac_id <- filtered_plot_data$fac_id |> as.vector() |> as.factor()
    
    p <- filtered_plot_data |>
        ggplot(aes(fill = size, x = fac_id, y = fac_pathway)) +
        geom_tile() +
        geom_text(aes(label = size), family = "FiraCode Nerd Font", size = 3) +
        theme_minimal() +
        scale_fill_gradient(low="#d7ffd4", high = "#10c400", guide = "colourbar") +
        scale_y_discrete(breaks = filtered_plot_data$fac_pathway, labels = NULL) +
        theme(
            text = element_text(family = "FiraCode Nerd Font", size = 10),
            panel.grid = element_blank()
        ) +
        ylab(NULL) + xlab(NULL) +
        guides(fill = guide_legend(title = "Set size"))
    
    p
}




main <- function(
        input_expression_means,
        input_tree,
        genesets_file,
        out_file,
        no_cluster = FALSE,
        extra_title = NULL,
        save_png = FALSE,
        png_res = 300,
        plot_width = 10,
        plot_height = 6,
        expressed_threshold = 0,
        renames = NA
) {
    cat(paste0("Reading in genesets from ", genesets_file, "\n"))
    genesets <- jsonlite::fromJSON(read_file(genesets_file))
    
    if (! is.na(renames)) {
        cat("Reading in renames\n")
        # This gives a named list. We need a named vector
        renames <- jsonlite::read_json(renames)
        nms <- names(renames)
        renames <- unlist(renames)
        names(renames) <- nms
    } else {
        renames <- NULL
    }
    
    pdata <- gen_plot_data(
        input_tree,
        input_expression_means,
        genesets,
        expressed_threshold = expressed_threshold,
        renames = renames
    )
    
    large_plot <- create_large_heatmap(pdata, extra_title = extra_title)
    
    # Add the set size points
    set_size_plot <- create_set_size_plot(pdata)
    
    large_plot <- large_plot + plot_spacer() + set_size_plot + plot_layout(
        ncol = 3, nrow = 1,
        widths = c(15, -1.1, 2),
        guides = "collect"
    ) + theme(legend.position = "right")
    
    # Save plot to output
    if (is.null(out_file)) {
        print(large_plot)
        return(invisible())
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
        input_dir = "~/Files/repos/tprof/data/out/enrichments",
        input_tree = "~/Files/repos/tprof/data/genesets_repr.txt",
        genesets_file = "~/Files/repos/tprof/data/genesets.json",
        input_dot_dir = "~/Files/repos/tprof/data/out/absolute_enrichments",
        out_file = NULL,
        save_png = TRUE,
        alpha = 0.20,
        png_res = 500,
        plot_height = 10
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

