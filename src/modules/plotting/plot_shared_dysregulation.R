options(warn = 1)

if (!exists("LOCAL_DEBUG")) {
  # Parsing arguments
  requireNamespace("argparser")
  
  parser <- argparser::arg_parser("Plot a shared dysregulation plot")
  
  parser |>
    argparser::add_argument(
      "input_results_dir",
      help = "Folder with output .tar.gz files to read.", type = "character"
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
      "--height",
      help = "Plot height, in inches.",
      default = 10, type = "numerical"
    ) |>
  argparser::add_argument(
    "--renames",
    help = "JSON file with renames to apply to sample names when plotting",
    default = NULL, type = "character"
  ) -> parser
  
  args <- argparser::parse_args(parser)
}

suppressMessages({
  options(tidyverse.quiet = TRUE)
  library(tidyverse)
  library(archive)
})

#' Create a named list of input tar files
find_input_paths <- function(input_dir) {
  x <- list.files(input_dir, pattern = ".*\\.tar\\.gz", full.names=TRUE)
  
  cat(paste0("Found ", length(x), " potential input files.\n"))
  
  res <- list()
  for (item in x) {
    clean_filename <- str_split_1(item, "/") |> tail(1) |> str_split_1("\\.") |> head(1) |> unlist()
    res[clean_filename] <- item
  }
  
  res
}

#' Get all the DEAs with ranking and all from a specific tarball.
#' 
#' Pass geo=TRUE to extract GEO deas instead
takeout_deas <- function(tarball, geo=FALSE) {
  path <- if (geo) {
    "TODO: Fill me up"
  } else {
    "data/deas/"
  }
  
}

find_input_paths("/home/hedmad/Files/repos/tprof/data/in/final")

