options(warn = 1)

if (!exists("LOCAL_DEBUG")) {
  # Parsing arguments
  requireNamespace("argparser")

  parser <- argparser::arg_parser("Plot a shared dysregulation plot")

  parser |>
    argparser::add_argument(
      "output_file",
      help = "The path to the output file",
      type = "character"
    ) |>
    argparser::add_argument(
      "input_results",
      help = "Folder with output .tar.gz files to read.", type = "character",
      nargs = Inf
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
    library(assertthat)
    library(ComplexUpset)
})

data <- read_csv("/home/hedmad/Files/repos/tprof/data/extracted_results/norm_fold_change_deas.csv")

strip_ensg_version <- function(x) {
    str_split_1(x, "\\.")[1]
}

assert_that(are_equal(strip_ensg_version("ENSG0000.12"), "ENSG0000"))

prep_data <- function(data, id_col = "sample") {
    data[[id_col]] <- map(data[[id_col]], strip_ensg_version)
    data |> column_to_rownames(id_col)
}

pdata <- prep_data(data)

extract_top_dysregulated <- function(data, n = 100) {
    types <- list()

    for (id in colnames(data)) {
        types[[id]] <- list()
        types[[id]][["down"]] <- sort_by(rownames(data), data[[id]])[1:n]
        types[[id]][["up"]] <- sort_by(rownames(data), decreasing=TRUE, data[[id]])[1:n]
    }

    types
}

top_dys <- extract_top_dysregulated(pdata)

transpose_dataframe <- function(x) {
    old_cols <- colnames(x)
    old_rows <- rownames(x)

    x_mtrx <- as.matrix(x)

    x_t_mtrx <- t(x_mtrx)

    x_t <- as.data.frame(x_t_mtrx)

    colnames(x_t) <- old_rows
    rownames(x_t) <- old_cols

    x_t
}

transpose_dataframe(data_frame(some = c(1, 2,3), stuff = c(5, 6, NA)))

#' This is the only function from UpSetR that I use, so I stole it.
#' It's 'fromList' originally
prep_for_upset <- function(input) {
    elements <- unique(unlist(input))
    data <- unlist(lapply(input, function(x) {
        x <- as.vector(match(elements, x))
    }))
    data[is.na(data)] <- as.integer(0)
    data[data != 0] <- as.integer(1)
    data <- data.frame(matrix(data, ncol = length(input), byrow = F))
    data <- data[which(rowSums(data) != 0), ]
    names(data) <- names(input)
    return(data)
}

plot_dysregulation <- function(data) {
    filtered_data_up <- list()
    filtered_data_down <- list()
    for (name in names(data)) {
        filtered_data_up[[name]] <- data[[name]][["up"]]
        filtered_data_down[[name]] <- data[[name]][["down"]]
    }

    prep_for_upset(filtered_data_up) -> dt_up
    prep_for_upset(filtered_data_down) -> dt_down

    dt_up$direction <- "up"
    dt_down$direction <- "down"

    dt <- rbind(dt_up, dt_down)

    p <- upset(
        dt, names(data),
        min_size = 5,
        base_annotations = list(
            'Intersection size' = intersection_size(
                mapping = aes(fill = direction),
                text=list(
                    vjust=-1,
                    hjust=-1,
                    angle=90
                )
            )
        ),
        set_sizes=FALSE
    )

    p <- p + ggtitle(paste0(
        "Overlap between top ",
        length(data[[1]][[1]]), " genes"
    ))

    p
}


plot_dysregulation(top_dys)
