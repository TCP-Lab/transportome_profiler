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
            help = "Path to large input deregulation table.", type = "character"
        ) |>
        argparser::add_argument(
            "--extra_title",
            help = "Extra title to add to the figure",
            type = "character", default = NULL
        ) |>
        argparser::add_argument(
            "--selected_genes",
            help = "Path to a txt list of genes to select among all",
            type = "character", default = NULL
        ) |>
        argparser::add_argument(
            "--id_col",
            help = "Name of the column with ENSG IDs",
            type = "character", default = "sample"
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

strip_ensg_version <- function(x) {
    str_split_1(x, "\\.")[1]
}

assert_that(are_equal(strip_ensg_version("ENSG0000.12"), "ENSG0000"))

prep_data <- function(data, id_col = "sample") {
    data[[id_col]] <- map(data[[id_col]], strip_ensg_version)
    data |> column_to_rownames(id_col)
}

extract_top_dysregulated <- function(data, n = 100) {
    types <- list()

    for (id in colnames(data)) {
        types[[id]] <- list()
        types[[id]][["down"]] <- sort_by(rownames(data), data[[id]])[1:n]
        types[[id]][["up"]] <- sort_by(rownames(data), decreasing = TRUE, data[[id]])[1:n]
    }

    types
}

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
            "Intersection size" = intersection_size(
                mapping = aes(fill = direction),
                text = list(
                    vjust = 0.5,
                    hjust = -1,
                    angle = 90
                )
            )
        ),
        set_sizes = FALSE,
        matrix = (
            intersection_matrix() + scale_y_discrete(
                # This does not work, see ComplexUpset issue 206
                labels = str_replace_all(names(data), "\\_", " "),
            )
        )
    )

    p <- p + ggtitle(paste0(
        "Overlap between top ",
        length(data[[1]][[1]]), " genes"
    ))

    p
}

purge_ensg_versions <- function(data, id_col = "gene_id") {
    data |> mutate("{id_col}" := str_remove(.data[[id_col]], "\\.[0-9]+$")) -> data
    data |> distinct(.data[[id_col]], .keep_all = TRUE) -> data
    data
}


main <- function(input_results,
                 output_file,
                 selected_genes = NULL,
                 extra_title = NULL,
                 res = 300,
                 width = 12,
                 height = 8,
                 save_png = FALSE,
                 id_col = "sample") {
    data <- read_csv(input_results, show_col_types = FALSE) |> purge_ensg_versions(id_col)
    selected_genes <- if (!is.null(selected_genes)) {
        read_file(selected_genes) |>
            str_split_1(",") |>
            str_remove_all("\"") |>
            str_remove_all("\\n")
    } else {
        NULL
    }

    if (!is.null(selected_genes)) {
        # This does not work since we need to strip the ensg version!!
        data |> filter(!!as.symbol(id_col) %in% selected_genes) -> data
    }

    pdata <- prep_data(data, id_col)
    top_dys <- extract_top_dysregulated(pdata)
    p <- plot_dysregulation(top_dys)

    # Save plot to output
    if (is.null(output_file)) {
        print(p)
    } else {
        if (save_png) {
            png(filename = output_file, width = width, height = height, units = "in", res = res)
        } else {
            pdf(file = output_file, width = width, height = height)
        }
        print(p)
        graphics.off()
    }
}


if (!exists("LOCAL_DEBUG")) {
    main(
        args$input_results,
        args$output_file,
        args$selected_genes,
        args$extra_title,
        args$res,
        args$width,
        args$height,
        args$png,
        args$id_col
    )
} else {
    main(
        "data/geo_merged_deas.csv",
        NULL,
        id_col = "gene_id",
        selected_genes = "data/filter_genes.txt"
    ) |> print()
}
