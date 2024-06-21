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
        
        types[[id]][["down"]] <- sort_by(tum_up$names, tum_up$value) |> head(n = n)
        types[[id]][["up"]] <- sort_by(tum_down$names, decreasing=TRUE, tum_down$value) |> head(n = n)
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
prep_for_upset <- function(input) {
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
    return(data)
}

assert_that(
    are_equal(
        prep_for_upset(
            list(
                test_tumor = c("A", "B", "C"),
                other_tumor = c("A", "D")
            )
        ),
        data.frame(
            test_tumor = c(1,1,1,0),
            other_tumor=c(1,0,0,1),
            row.names=c("A", "B", "C", "D")
        )
    )
)

gen_plot_data <- function(data) {
    filtered_data_up <- list()
    filtered_data_down <- list()
    for (name in names(data)) {
        filtered_data_up[[name]] <- data[[name]][["up"]]
        filtered_data_down[[name]] <- data[[name]][["down"]]
    }
    
    # This should not work. It only works bc there are no genes that are up
    # in some tt and down in another.
    
    prep_for_upset(filtered_data_up) -> dt_up
    prep_for_upset(filtered_data_down) -> dt_down
    
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

plot_shared_genes <- function(data, data_renames) {
    dt <- gen_plot_data(data) |> tibble() |> arrange(name) 
    
    dt <- dt |> mutate(row_value = rowSums(select_if(dt, is.numeric)))
    
    bar_order <- dt |> group_by(name) |> summarise(total = sum(row_value))
    
    dt <- merge(dt, bar_order, by="name")
    
    dt <- merge(
        dt,
        data_renames[, c("ensembl_gene_id", "hgnc_symbol")],
        by.x = "name", by.y = "ensembl_gene_id",
        all.x = TRUE, all.y = FALSE
    )
    print(tibble(dt))
    
    # Keep only the top n genes
    
    # TODO: This is bad - if there are genes with =/= ensg but the same hgnc,
    # this causes a collision. But it's impossible (?) to revalue the labels
    # even with scale_x_discrete (or at least I could not get it to work).
    # I plot hgnc_symbols directly and hope for the best
    dt <- dt |> filter(total > 13)
    
    bar_plot <- ggplot(dt, aes(x=reorder(hgnc_symbol, total), y=-row_value)) +
        geom_bar(stat = "identity", aes(fill = direction), position = "stack") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
        coord_flip()
    
    dot_plot 
}

plot_shared_genes(top_dys, ensg_data)
