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
    library(assertthat)
    library(ComplexUpset)
    library(reshape2)
    library(patchwork)
    requireNamespace("gridExtra")
})

data <- read_csv(
  "/home/hedmad/Files/repos/tprof/data/extracted_results/deseq_shrinkage_deas.csv",
  show_col_types = FALSE
)
names <- read_file("/home/hedmad/Files/repos/tprof/names.txt") |>
    str_split_1(",") |> str_remove_all("\"") |> str_remove_all("\\n")
ensg_data <- read_csv(
    "/home/hedmad/Files/repos/tprof/data/ensg_data.csv"
)

strip_ensg_version <- function(x) {
    str_split_1(x, "\\.")[1]
}

assert_that(are_equal(strip_ensg_version("ENSG0000.12"), "ENSG0000"))
assert_that(are_equal(strip_ensg_version("ENSG0000.1"), "ENSG0000"))

prep_data <- function(data, id_col = "sample") {
    data[[id_col]] <- map(data[[id_col]], strip_ensg_version)
    data |> column_to_rownames(id_col)
}

pdata <- prep_data(data)

extract_top_dysregulated <- function(data, n = 100, thr = 1, gene_filter = NULL) {
    types <- list()

    if (!is.null(gene_filter)) {
        cat("Filtering...\n")
        data <- data[row.names(data) %in% gene_filter,]
    }

    for (id in colnames(data)) {
        types[[id]] <- list()

        tum_data <- data |> rownames_to_column("names") |> select(c(!!id, names)) |>
            rename(value = !!id)

        tum_data <- tum_data |> filter(abs(value) > thr)

        tum_up <- tum_data |> filter(value > 0)
        tum_down <- tum_data |> filter(value < 0)

        types[[id]][["down"]] <- sort_by(tum_down$names, tum_down$value) |> head(n = n)
        types[[id]][["up"]] <- sort_by(tum_up$names, decreasing=TRUE, tum_up$value) |> head(n = n)
    }

    types
}

top_dys <- extract_top_dysregulated(pdata, thr=1, gene_filter = names, n=Inf)

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

assert_that(are_equal(
    transpose_dataframe(data.frame(some = c(1, 2, 3), stuff = c(5, 6, NA), row.names=c("A", "B", "C"))),
    data.frame(A=c(1,5), B=c(2,6), C=c(3,NA), row.names= c("some", "stuff"))
))

#' This is the only function from UpSetR that I use, so I stole it
#' I had to edit it in the end, so good thing that I stole it.
#' It's 'fromList' originally
#'
#' I edited this so much that now I probably do it backwards: i first set
#' the values to 0 and 1, and then I update them with real values.
#' Perhaps there's a better way...
prep_for_upset <- function(input, values) {
    elements <- unique(unlist(input))
    data <- unlist(lapply(input, function(x) {
        x <- as.vector(match(elements, x))
    }))
    data[is.na(data)] <- as.integer(0)
    data[data != 0] <- as.integer(1)
    data <- data.frame(matrix(data, ncol = length(input), byrow = F))
    row.names(data) <- elements
    data <- data[which(rowSums(data) != 0), ]
    colnames(data) <- names(input)

    # Convert the 1s to actual numbers
    for (col in colnames(data)) {
        if (! col %in% names(values)) {
            next
        }

        rep_positions <- data[[col]] == 1
        data[[col]][rep_positions] <- values[row.names(data), col][rep_positions]
    }

    return(data)
}

assert_that(
    are_equal(
        prep_for_upset(
            list(
                test_tumor = c("A", "B", "C"),
                other_tumor = c("A", "D")
            ),
            data.frame(
                test_tumor = c("A" = 10, "B" =12.2, "C"=11, "D"=-25),
                other_tumor = c("A" = 2, "B" =6, "C"=100, "D"=-20)
            )
        ),
        data.frame(
            test_tumor = c(10,12.2,11,0),
            other_tumor=c(2,0,0,-20),
            row.names=c("A", "B", "C", "D")
        )
    )
)

gen_plot_data <- function(data, values) {
    filtered_data_up <- list()
    filtered_data_down <- list()
    for (name in names(data)) {
        filtered_data_up[[name]] <- data[[name]][["up"]]
        filtered_data_down[[name]] <- data[[name]][["down"]]
    }

    # This should not work. It only works bc there are no genes that are up
    # in some tt and down in another.

    prep_for_upset(filtered_data_up, values) -> dt_up
    prep_for_upset(filtered_data_down, values) -> dt_down

    dt_up$direction <- "up"
    dt_down$direction <- "down"

    dt_up |> rownames_to_column("name") -> dt_up
    dt_down |> rownames_to_column("name") -> dt_down

    dt <- rbind(dt_up, dt_down)
    assert_that(are_equal(nrow(dt), nrow(dt_up) + nrow(dt_down)))

    dt
}

plot_dysregulation <- function(data) {
    dt <- gen_plot_data(data)

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

select_initial <- function(x) {
    str_split_1(x, "_")[1]
}

tt_renames <- list(
    "Head_n_Neck_cancer_deseq" = "Head and Neck"
)

rename_tumor_type <- function(x) {
    print(x)
    if (x %in% names(tt_renames)) {
        return(tt_renames[[x]])
    }
    return(select_initial(x))
}


plot_shared_genes <- function(
        data,
        data_renames,
        values,
        types_renames_fn = rename_tumor_type,
        n_genes = 50
    ) {
    dt <- gen_plot_data(data, values) |> tibble() |> arrange(name)

    dt <- merge(
        dt,
        data_renames[, c("ensembl_gene_id", "hgnc_symbol")],
        by.x = "name", by.y = "ensembl_gene_id",
        all.x = TRUE, all.y = FALSE
    )

    binary_dt <- dt |> mutate(
        across(where(is.numeric),
               \(x) {as.numeric(as.logical(abs(x)))}
    ))

    bar_dt <- binary_dt |> mutate(row_value = rowSums(select_if(binary_dt, is.numeric)))

    bar_order <- bar_dt |> group_by(name) |> summarise(total = sum(row_value))

    bar_dt <- merge(bar_dt, bar_order, by="name")

    # TODO: This is bad - if there are genes with =/= ensg but the same hgnc,
    # this causes a collision. But it's impossible (?) to revalue the labels
    # even with scale_x_discrete (or at least I could not get it to work).
    # I plot hgnc_symbols directly and hope for the best

    # I need to double the number of genes here since there are two rows
    # per gene in bar_dt, so this way we keep both
    bar_dt <- bar_dt |> top_n(n_genes * 2, total)

    bar_plot <- ggplot(bar_dt, aes(x=reorder(hgnc_symbol, total), y=row_value)) +
        geom_bar(stat = "identity", aes(fill = direction), position = "stack") +
        scale_fill_manual(values = c("purple", "darkorange")) +
        theme_minimal() +
        theme(
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "none",
            legend.title = element_blank()
        ) +
        #scale_y_reverse() +
        scale_x_discrete(position = "top") +
        coord_flip()

    bar_x_values <- layer_scales(bar_plot)$x$range$range

    dot_dt <- dt |> filter(hgnc_symbol %in% bar_x_values) |>
        select(!c("name", "direction"))

    dot_dt_clust <- select(dot_dt, !c("hgnc_symbol"))

    dot_dt_clust.col <- as.dendrogram(hclust(dist(dot_dt_clust), method = "ward.D2"))

    # ordering based on clustering
    col.ord <- order.dendrogram(dot_dt_clust.col)

    dot_dt <- dot_dt |> melt(id.vars = c("hgnc_symbol")) |>
        filter(value != 0)

    dot_dt$hgnc_symbol <- factor(dot_dt$hgnc_symbol, levels=bar_x_values)
    dot_dt$variable <- factor(dot_dt$variable, levels=colnames(dot_dt_clust)[col.ord])

    levels(dot_dt$variable) <- sapply(levels(dot_dt$variable), types_renames_fn)

    dot_plot <- ggplot(
        dot_dt,
        aes(
            y=hgnc_symbol,
            x=variable,
            fill=value
        )) +
        scale_fill_gradient2(
            low = "purple", mid="white", high="darkorange",
            name = "Score"
        ) +
        geom_tile() +
        theme_minimal() +
        theme(
            axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            #axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "left",
        )

    combined <- dot_plot +
        ggtitle(paste0("Top ", min(n_genes, length(bar_x_values)), " dysregulated genes shared across tumor types")) +
        bar_plot +
        plot_layout(widths = c(2, 0.5))

    print(combined)
}

plot_shared_genes(top_dys, ensg_data, pdata)
