#!/usr/bin/env Rscript

#' This script converts GSEA output files from `fgsea` into flower plots.
#'
#' Flower plots are the ones with the graph in the circular fashion, with
#' colored dots as the leaves.

if (sys.nframe() == 0L) {
  # Parsing arguments
  requireNamespace("argparser")

  parser <- argparser::arg_parser("Plot Tree graphs from GSEA outputs")

  parser |>
    argparser::add_argument(
      "input_gsea_result", help=".csv file with GSEA output", type="character"
    ) |>
    argparser::add_argument(
      "genesets", help="JSON file with geneset information", type="character"
    ) |>
    argparser::add_argument(
      "output_path", help = "Output file to save",
      type = "character"
    ) |>
    argparser::add_argument(
      "--png", help = "If specified, saves plots as PNG. Use `--res` to set the resolution",
      flag = TRUE, type = "logical"
    ) |>
    argparser::add_argument(
      "--res", help = "Resolution of PNG plots, in pixels per inch.",
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

suppressMessages({
  library(tidyverse, quietly=TRUE)
  requireNamespace("RColorBrewer")
  requireNamespace("igraph")
  requireNamespace("uuid")
  requireNamespace("rjson")
  library(ggraph) # This needs to be a library() call
  library(grid)
  library(assertthat)
})


#' Read a series of .csv files from a directory
#'
#' This is for the output files from fgsea::fgsea in the `run_gsea.R` file
#'
#' @param input_dir The input directory to read from. Reads ALL files.
#'
#' @returns A list with file names as names and tibbles as values.
read_results <- function(input_dir) {
  files <- list.files(input_dir)

  # Remove flag files
  files <- files[! endsWith(files, ".flag")]

  # Remove non-csv files
  files <- files[endsWith(files, ".csv")]

  res <- list()
  for (i in seq_along(files)) {
    res[[files[i]]] <- read_csv(file.path(input_dir, files[i]), show_col_types = FALSE)
  }
  return(res)
}

#' Convert a result table to a graph structure for plotting
#'
#' The graph structure is derived from the list names, using `\` as node name
#' splitting character. E.g. `a\b\c` becomes `a -> b -> c` in the graph.
#'
#' @param result The data.frame with the GSEA results. Needs at least the
#'   "pathway", "NES" and "padj" columns, to add the data to the graph.
#' @param genesets The JSON genesets with the information on the graph structure
#' @param base_edges A data.frame with base edges
result_to_graph <- function(result, genesets) {
  # This fun gets a results list (w/o plots) and converts it to a dataframe
  # that can be used by ggraph and igraph
  genesets_data <- rjson::fromJSON(readr::read_file(genesets))
  # We can now build the frame of edges, where each row is a source -> sink edge.
  insert <- function(x, item) {
    x[[length(x) + 1]] <- item

    x
  }

  edges <- list()
  for (i in seq_along(genesets_data)) {
    parent <- genesets_data[[i]]$parent
    edges[[i]] <- c(parent, names(genesets_data)[i])
  }

  edges <- as.data.frame(do.call(rbind, edges))
  colnames(edges) <- c("source", "sink")

  edges |> distinct() -> edges

  # Now we have a frame of edges, so we can grab the node data from the results

  vertice_data <- list()
  vertices <- unique(unlist(edges))
  vertices <- vertices[!is.null(vertices)]
  for (i in seq_along(vertices)) {
    item <- vertices[i]

    NES <- result[result$pathway == item, "NES"]
    padj <- result[result$pathway == item, "padj"]
    if (is.null(NES)) {
      NES <- 0
    }
    if (is.null(padj)) {
      padj <- 1
    }
    item_data <- c(
      item,
      genesets_data[[item]]$name,
      NES,
      padj
    )

    assert_that(length(item_data) == 4)

    vertice_data[[i]] <- item_data

  }
  vertice_frame <- as.data.frame(do.call(rbind, vertice_data))

  colnames(vertice_frame) <- c("uuid", "human_label", "NES", "padj")
  # Convert to numbers
  vertice_frame |> mutate(NES = as.numeric(NES), padj = as.numeric(padj)) -> vertice_frame

  vertice_frame$NES[is.na(vertice_frame$NES)] <- 0

  return(igraph::graph_from_data_frame(edges, vertices = vertice_frame))
}

# ----

make_colours <- function(palette, values) {
  colour_fun <- colorRamp(palette)

  values <- (values-min(values))/(max(values)-min(values))

  col_values <- colour_fun(values)

  colours <- apply(col_values, 1, function(x) {
    x[is.na(x)] <- 0
    rgb(x[1], x[2], x[3], maxColorValue = 255)
  })
}

parse_name_to_label <- function(x) {
  sapply(x, \(value){
    if (startsWith(value, "carried_solute::")) {
      value <- str_remove(value, "carried_solute::")
    }

    replacements <- c(
      "outward", "inward", "voltage independent", "voltage gated", "not LG", "ligand gated",
      "voltage independent", "voltage gated", "not LG", "LG", ""
    )
    names(replacements) <- c(
      "direction::out", "direction::in", "is_voltage_gated::0", "is_voltage_gated::1",
      "is_ligand_gated::0", "is_ligand_gated::1", "is_voltage_gated::0.0", "is_voltage_gated::1.0",
      "is_ligand_gated::0.0", "is_ligand_gated::1.0", "whole_transportome"
    )

    if (value %in% names(replacements)) {
      value <- replacements[value]
    }

    value
  })
}

get_point_coords <- function(p) {
  ggp <- ggplot_build(p)

  return(ggp$data[[1]][, c("label", "x", "y")])

}

calculate_angle_from_pos <- function(pos_dataframe, specials = NULL) {
  pos_dataframe$angle <- apply(pos_dataframe[,c("x", "y")], 1, \(x) {atan(x[2] / x[1])})
  # Replace every NaN (like, 0/0) with angle = 0
  pos_dataframe[is.na(pos_dataframe)] <- 0

  # Change the positions to be slightly more outward
  # 1. we calculate the hypothenuse + the dodge value
  scaling_factor <- 0.01
  new_coords <- apply(pos_dataframe, 1, \(row) {
    name <- row["label"]; y <- as.numeric(row["y"]); angle <- as.numeric(row["angle"])
    # It's important that we
    x <- as.numeric(row["x"]);
    hypothenuse <- sqrt(x ** 2 + y ** 2) + (min(max(str_length(name), 5), 15) * scaling_factor)
    if (y < 0) {
      new_y <- hypothenuse * abs(sin(angle)) * - 1
    } else {
      new_y <- hypothenuse * abs(sin(angle))
    }

    if (x < 0) {
      new_x <- hypothenuse * abs(cos(angle)) * -1
    } else {
      new_x <- hypothenuse * abs(cos(angle))
    }

    return(c(new_x, new_y))
  })

  pos_dataframe$dodged_x <- new_coords[1,] # it is filled row-wise
  pos_dataframe$dodged_y <- new_coords[2,]

  # From radians to degrees
  pos_dataframe$angle <- pos_dataframe$angle * 180 / pi

  # Detect which labels do not lay on the outer circle, so that we can
  # label them differently
  hypothenuses <- sqrt(pos_dataframe$x ** 2 + pos_dataframe$y ** 2)

  pos_dataframe
}

plot_result <- function(result, genesets_file, title = "") {

  data_graph <- result_to_graph(result, genesets_file)
  colours <- make_colours(c("blue", "gray", "red"), as.numeric(igraph::vertex_attr(data_graph, "NES")))

  # "Plot" a graph with just the labels
  p <- ggraph(data_graph, layout='igraph', algorithm = "tree", circular = TRUE) +
    coord_fixed() +
    geom_node_text(
      aes(
        label = parse_name_to_label(igraph::vertex_attr(data_graph, "human_label"))
      )
    )

  plot_labels <- calculate_angle_from_pos(get_point_coords(p))

  expand_vec <- c(0.05, 0.05)

  # Now we have the angles, we can build the real plot
  pp <- ggraph(data_graph, layout='igraph', algorithm = "tree", circular = TRUE) +
    geom_edge_diagonal(aes(alpha = after_stat(index)), show.legend = FALSE) +
    coord_fixed() +
    scale_edge_colour_distiller(palette = "RdPu") +
    geom_node_point(
      aes(
        size = -log10(as.numeric(igraph::vertex_attr(data_graph, "padj"))),
        color = colours
      ),
      alpha = (as.numeric(as.numeric(igraph::vertex_attr(data_graph, "padj")) < 0.05) + 0.2 ),
      show.legend = setNames(c(FALSE, FALSE, FALSE), c("color", "size", "alpha"))
    ) +
    geom_node_text(
      aes(
        x = plot_labels$dodged_x,
        y = plot_labels$dodged_y,
        label = plot_labels$label
      ), angle = plot_labels$angle,
      linesize = 2.5
    ) +
    scale_color_manual(values = colours, limits = colours, guide = guide_legend(title = "NES")) +
    theme(legend.position = "bottom", panel.background = element_blank()) +
    # Give more space to the plot area so the lables are drawn properly
    scale_x_continuous(expand = expand_vec) + scale_y_continuous(expand = expand_vec) +
    ggtitle(title)

  return(pp)
}

if (sys.nframe() == 0L) {
  print(paste0("Plotting flower plot for ", args$input_gsea_result))

  args$input_gsea_results |>
    stringr::str_remove("\\.csv") |>
    stringr::str_split_i() -> fname
  p <- plot_result(
    read.csv(args$input_gsea_result),
    args$genesets,
    title = fname
  )

  if (args$png) {
    png(args$output_path, width = args$width, height = args$height, units = "in", res = args$res)
    print(p)
    dev.off()
  } else {
    pdf(args$output_path, width = args$width, height = args$height)
  }
}
