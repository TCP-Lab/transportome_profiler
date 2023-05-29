 #!/usr/bin/env Rscript

if (sys.nframe() == 0L) {
  # Parsing arguments
  requireNamespace("argparser")
  
  parser <- argparser::arg_parser("Plot a large heatmat with GSEA output")
  
  parser |>
    argparser::add_argument(
      "input_gsea_results", help="Folder with GSEA output .csv files to read.", type="character"
    ) |>
    argparser::add_argument(
      "input_tree", help="A txt file with the folder structure of the tree to be parsed", type="character"
    ) |>
    argparser::add_argument(
      "output_file", help = "The path to the output file",
      type = "character"
    ) |>
    argparser::add_argument(
      "--res", help = "Resolution of plot, in pixels per inch.",
      default= 400, type = "logical"
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

RUN_LOCAL <- FALSE

options(tidyverse.quiet = TRUE)
library(tidyverse)
requireNamespace("stringi")

parse_tree_labels <- function(tree) {
  #' Parse the raw tree from the file to an usable labels dataframe
  #' 
  #' The dataframe has the following cols:
  #' - original: The original, raw tree
  #' - cropped: The raw tree, with just the paths from the whole_transportome dir
  #' - pretty: The prettified tree, without the full paths
  #' - backbone: Just the backbone of the tree, with no paths
  #' - rev_backbone: Just the backbone of the tree, but reversed.
  #' - rev_pretty: The prettified tree, but reversed (with the backbone on the
  #'   right)
  #' - paths: The full paths to the tree (without the backbone)
  #' - clean_paths: The relative paths to the tree (without the backbone)
  #' 
  #' @param tree A vector of lines where each line is from the tree output
  
  result <- data.frame(original = tree, order = length(tree):1)
  
  result$cropped <- str_split_i(result$original, "whole_transportome", 2)
  fix_crop <- result$cropped
  # The "whole_transportome" is reduced to "". This fixes it.
  fix_crop[fix_crop == ""] <- "Whole Transportome"
  
  remove_backbone <- function (x) {
    # Remove all backbone chars from an input
    str_remove_all(x, "[└─│├]") |> str_trim()
  }
  result$paths <- remove_backbone(result$original)
  
  result$backbone <- str_split_i(result$original, "─ ", 1) |> paste0("─ ")
  result$backbone[1] <- ""
  
  result$rev_backbone <- stringi::stri_reverse(result$backbone) |>
    str_replace_all("├", "┤") |>
    str_replace_all("└", "┘")
  
  result$pretty <- paste0(result$backbone, str_split_i(fix_crop, "/", -1))
  result$rev_pretty <- paste0(str_split_i(fix_crop, "/", -1), result$rev_backbone)
  
  result
}

if (RUN_LOCAL) {
  tree <- read_lines("/tmp/geneset_tree/tree.txt")
  labs <- parse_tree_labels(tree)
}

main <- function(
    input_dir,
    input_tree,
    out_file,
    save_png = FALSE,
    png_res = 300,
    plot_width = 10,
    plot_height = 6,
    alpha = 0.20
) {
  # Load and parse the tree labels
  tree <- read_lines(input_tree)
  labels <- parse_tree_labels(tree)
  
  input_files <- list.files(input_dir, ".csv", full.names = TRUE)
  enrichment_data <- lapply(input_files, \(x) {read.csv(x)})
  
  # Get the 'pathway' var to look like the paths in the labels
  # this means getting rid of the /whole_transportome leading bit
  enrichment_data <- lapply(enrichment_data, \(frame) {frame$clean_pathway <- str_split_i(frame$pathway, "whole_transportome", 2) ; frame})
  
  # For the heatmap we will need a melted matrix. So I first combine all
  # the enrichment frames
  enrichment_data <- lapply(seq_along(enrichment_data), function(i) {
    frame <- enrichment_data[[i]]
    frame$id <- str_remove_all(input_files, ".csv")[i] |> str_split_i("/", -1)
    
    frame
  })
  
  # We can now join all the frames together
  plot_data <- reduce(enrichment_data, rbind)
  
  plot_data$alpha_from_padj <- 0.40
  plot_data$alpha_from_padj[plot_data$padj < alpha] <- 1
  
  plot_data$fac_id <- factor(plot_data$id)
  print(labels)
  plot_data$fac_cpath <- factor(plot_data$clean_pathway,  levels = labels$cropped[labels$order])
  
  p <- ggplot(plot_data, aes(fill = NES, x = fac_id, y = fac_cpath, alpha = alpha_from_padj)) +
    geom_tile() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    ggtitle(paste0("Deregulation Overview (Alpha ", alpha, ")")) +
    ylab("Pathway") + xlab("Cohort") +
    scale_fill_gradient2("NES", low = "purple", mid = "darkgray", high = "darkorange") +
    scale_alpha_identity("alpha_from_padj") +
    scale_y_discrete(breaks = labels$cropped, labels = labels$rev_pretty) +
    theme(text = element_text(family = "FiraCode Nerd Font", size = 10))
  
  
  # Save plot to output
  if (save_png) {
    png(filename = out_file, width = plot_width, height = plot_height, units = "in", res = png_res)
  } else {
    pdf(file = out_file, width = plot_width, height = plot_height)
  }
  print(p)
  dev.off()
}

if (RUN_LOCAL) {
  main(
    input_dir = "/home/hedmad/Desktop/banana/out/enrichments",
    input_tree = "/tmp/geneset_tree/tree.txt",
    out_file = "/home/hedmad/Files/repos/transportome_profiler/data/out/figures/full_heatmap.png",
    save_png = TRUE,
    alpha = 0.20,
    png_res = 300,
    plot_height = 10
    )
}

main(
    input_dir = args$input_gsea_results,
    input_tree = args$input_tree,
    out_file = args$output_file,
    save_png = TRUE,
    png_res = args$res,
    plot_width = args$width,
    plot_height = args$height
)
