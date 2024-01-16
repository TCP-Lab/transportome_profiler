 #!/usr/bin/env Rscript

options(warn = 1)

if (! exists("LOCAL_DEBUG")) {
  # Parsing arguments
  requireNamespace("argparser")
  
  parser <- argparser::arg_parser("Plot a large heatmat with relative and absolute GSEA output.")
  
  parser |>
    argparser::add_argument(
      "input_gsea_results", help="Folder with GSEA output .csv files to read.", type="character"
    ) |>
    argparser::add_argument(
      "--dots_gsea_results", help="Folder with GSEA output .csv files to read to display as dots.", type="character",
      default = NULL
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
      "--no_cluster", help = "Skip clustering of x axis labels",
      flag=TRUE, type = "logical"
    ) |>
    argparser::add_argument(
      "--alpha", help = "Alpha level to plot at",
      default=0.05, type = "numerical"
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
      "--height", help = "Plot height, in inches.",
      default = 10, type = "numerical"
    ) -> parser
  
  args <- argparser::parse_args(parser)
}

suppressMessages({
  options(tidyverse.quiet = TRUE)
  library(tidyverse)
  requireNamespace("stringi")
  requireNamespace("reshape2")
  extrafont::loadfonts()
})

parse_tree_labels <- function(tree, genesets) {
  #' Parse the raw tree from the file to an usable labels dataframe
  #' 
  #' @param tree A vector of lines where each line is from the tree output
  #' @param genesets The JSON geneset data
  
  result <- data.frame(original = tree, order = length(tree):1)
  
  remove_backbone <- function (x) {
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
    input_dir,
    genesets,
    alpha = 0.05,
    cluster_x_axis = FALSE
  ) {
  tree <- read_lines(input_tree)
  labels <- parse_tree_labels(tree, genesets)
  
  input_files <- list.files(input_dir, ".csv$", full.names = TRUE)
  enrichments <- lapply(input_files, function(x) {read.csv(x)})
  
  # Get the 'pathway' var to look like the paths in the labels
  # this means getting rid of the /whole_transportome leading bit
  enrichments <- lapply(enrichments, \(frame) {
    frame$pathway_name <- sapply(frame$pathway, \(id) {genesets[[id]]$name}, simplify = TRUE)
    frame$label <- sapply(frame$pathway, \(id) {labels[labels$id == id, "rev_pretty"]})
    frame
  })
  
  # For the heatmap we will need a melted matrix. So I first combine all
  # the enrichment frames
  enrichments <- lapply(seq_along(enrichments), function(i) {
    frame <- enrichments[[i]]
    frame$id <- str_remove_all(input_files, ".csv")[i] |> str_split_i("/", -1) |>
      str_remove_all("_deseq") |> str_replace_all("_", " ") # Further clean the inputs
    
    frame
  })
  
  # We can now join all the frames together
  plot_data <- reduce(enrichments, rbind)
  
  # Set the alpha values manually
  plot_data$alpha_from_padj <- 0.40
  plot_data$alpha_from_padj[plot_data$padj < alpha] <- 1
  
  # Set the order of the column samples based on hclust
  # - We need to make a matrix from the input
  # - This is fine to be ran on just one type of data, based on what we want to
  #   sort by. I choose to run it on the relative data.
  if (cluster_x_axis) {
    clust_data <- reshape2::dcast(plot_data, id ~ pathway, value.var = "NES")
    clust_data |> column_to_rownames("id") -> clust_data
    
    clust <- hclust( dist( clust_data ), method = "ward.D")
    
    # Set the order of the labels. The actual values will be set in the plot
    plot_data$fac_id <- factor(plot_data$id, levels = clust$labels[clust$order])
  } else {
    # If not clustered, resort to alphabetical clustering
    plot_data$fac_id <- factor(plot_data$id, levels = plot_data$id[order(plot_data$pathway)])
  }
  plot_data$fac_pathway <- factor(plot_data$pathway,  levels = labels$id[labels$order])
  
  plot_data$show <- FALSE
  plot_data$show[plot_data$padj < alpha] <- TRUE
  
  plot_data
}

create_large_heatmap <- function(
    plot_data
) {
  p <- ggplot(plot_data, aes(fill = NES, x = fac_id, y = fac_pathway, alpha = alpha_from_padj)) +
    geom_tile() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle("Deregulation Overview") +
    ylab("Gene Set") + xlab("Cohort") +
    scale_fill_gradientn( # not a typo - it is really called gradientn
      "NES",
      colours = c("purple", "gray", "gray", "darkorange")
    ) + 
    scale_alpha_identity("alpha_from_padj") +
    scale_y_discrete(breaks = plot_data$fac_pathway, labels = plot_data$label) +
    theme(text = element_text(family = "FiraCode Nerd Font", size = 10))
  
  p
}

add_dots <- function(
  original_plot,
  plot_data
) {
  original_plot +
    geom_point(data = plot_data, aes(x = fac_id, y = fac_pathway, alpha = show), shape = 8)
}


main <- function(
    input_dir,
    input_tree,
    genesets_file,
    out_file,
    input_dot_dir = NULL,
    no_cluster = FALSE,
    save_png = FALSE,
    png_res = 300,
    plot_width = 10,
    plot_height = 6,
    alpha = 0.20
) {
  genesets <- jsonlite::fromJSON(read_file(genesets_file))
  relative_plot_data <- gen_plot_data(
    input_tree = input_tree,
    input_dir = input_dir,
    genesets = genesets,
    cluster_x_axis = (! no_cluster)
  )
  large_plot <- create_large_heatmap(relative_plot_data)
  
  if (! is.null(input_dot_dir) ) {
    dot_plot_data <- gen_plot_data(input_tree = input_tree, input_dir = input_dot_dir, genesets = genesets)
    large_plot <- add_dots(large_plot, dot_plot_data)
  }
  
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
    input_dir = "~/Files/repos/transportome_profiler/data/out/enrichments",
    input_tree = "~/Files/repos/transportome_profiler/data/genesets_repr.txt",
    genesets_file = "~/Files/repos/transportome_profiler/data/genesets.json",
    out_file = NULL,
    save_png = TRUE,
    alpha = 0.20,
    png_res = 500,
    plot_height = 10
    )
} else {
  main(
    input_dir = args$input_gsea_results,
    input_tree = args$input_tree,
    genesets_file = args$genesets,
    out_file = args$output_file,
    input_dot_dir = args$dots_gsea_results,
    no_cluster = args$no_cluster,
    alpha = args$alpha,
    save_png = TRUE,
    png_res = args$res,
    plot_width = args$width,
    plot_height = args$height
  )
}

