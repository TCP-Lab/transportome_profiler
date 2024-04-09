options(warn = 1)

if (!exists("LOCAL_DEBUG")) {
  # Parsing arguments
  requireNamespace("argparser")

  parser <- argparser::arg_parser("Plot a fold change heatmap")

  parser |>
    argparser::add_argument(
      "input_results_dir",
      help = "Folder with output .csv files to read.", type = "character"
    ) |>
    argparser::add_argument(
      "output_file",
      help = "The path to the output file",
      type = "character"
    ) |>
    argparser::add_argument(
      "--extra_title",
      help = "Extra title to add to the figure",
      type = "character", default = NULL
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
      default = 10, type = "numerical"
    ) |>
    argparser::add_argument(
      "--renames",
      help = "JSON file with renames to apply to sample names when plotting",
      default = NULL, type = "character"
    ) |>
    argparser::add_argument(
      "--filter_genes",
      help = "Text file with ENSGs to filter input with",
      default = NULL, type = "character"
    ) |>
    argparser::add_argument(
      "--height",
      help = "Plot height, in inches.",
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
  library(gplots)
})

gen_plot_data <- function(
    input_dir) {
  input_files <- list.files(input_dir, ".csv$", full.names = TRUE)
  fold_changes <- lapply(input_files, function(x) {
    read.csv(x)
  })

  fold_changes <- lapply(seq_along(fold_changes), function(i) {
    frame <- fold_changes[[i]]
    frame$id <- str_remove_all(input_files, ".csv")[i] |>
      str_split_i("/", -1) |>
      str_remove_all("_deseq") |>
      str_replace_all("_", " ") # Further clean the inputs

    frame
  })

  plot_data <- reduce(fold_changes, rbind)

  plot_data
}

strip_version <- function(x) {
  return(str_split_i(x, "\\.", 1))
}

main <- function(
    input_results_dir,
    output_file,
    extra_title = NULL,
    res = 300,
    width = 12,
    height = 8,
    save_png = FALSE,
    filter_genes = NA,
    renames = NULL) {
  plot_data <- gen_plot_data(input_results_dir)

  # Heatmap() needs a matrix, so we have to recast
  cat("Reading input frame...\n")
  plot_data <- reshape2::recast(plot_data, sample ~ id, measure.var = "ranking") |>
    column_to_rownames("sample")

  if (!is.na(filter_genes)) {
    cat("Filtering...\n")
    genes <- read.table(filter_genes)
    new_rows <- strip_version(row.names(plot_data))
    plot_data <- plot_data[!duplicated(new_rows), ]
    row.names(plot_data) <- new_rows[!duplicated(new_rows)]
    plot_data <- plot_data[unlist(genes), ]
  }

  cat("Dropping rownames...\n")
  rownames(plot_data) <- NULL
  cat("Casting to matrix...\n")
  tmp <- as.matrix(plot_data)
  cat("Setting NAs to zero...\n")
  tmp[is.na(tmp)] <- 0

  hmp_call <- function() {
    heatmap.2(
      tmp,
      trace = "none",
      col = "bluered",
      xlab = "Tumor Type",
      ylab = "Gene Expression",
      main = ifelse(
        is.na(extra_title),
        "Metric heatmap",
        paste0("Metric heatmap - ", extra_title)
      ),
      scale = "row",
      margins = c(10, 4),
      labRow = ""
    )
  }

  # Save plot to output
  if (is.null(output_file)) {
    hmp_call()
    return(invisible())
  } else {
    if (save_png) {
      png(filename = output_file, width = width, height = height, units = "in", res = res)
    } else {
      pdf(file = output_file, width = width, height = height)
    }
    hmp_call()
    dev.off()
  }
}

main(
  input_results_dir = args$input_results_dir,
  output_file = args$output_file,
  extra_title = args$extra_title,
  res = args$res,
  width = args$width,
  height = args$height,
  renames = args$renames,
  save_png = args$png,
  filter_genes = args$filter_genes
)
