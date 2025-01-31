options(warn = 1)

suppressMessages({
    options(tidyverse.quiet = TRUE)
    library(tidyverse)
    library(assertthat)
    library(ComplexUpset)
    library(reshape2)
    library(patchwork)
    requireNamespace("gridExtra")
})

{
    # This removes the version from an ensembl ID
    strip_ensg_version <- function(x) {
        str_split_1(x, "\\.")[1]
    }

    assert_that(are_equal(strip_ensg_version("ENSG0000.12"), "ENSG0000"))
    assert_that(are_equal(strip_ensg_version("ENSG0000.1"), "ENSG0000"))
}

#' Perform basic data preparation before starting processing
#'
#' This strips the (eventual) ensembl id version from the data IDs and
#' moves the ID to row names both to use native R and to standardize the
#' name of the column in further processing.
prep_data <- function(data, id_col = "sample") {
    data[[id_col]] <- map(data[[id_col]], strip_ensg_version)
    data |> column_to_rownames(id_col)
}

#' Extract the top dysregulated genes from the data
#'
#' This can apply three types of filters:
#' - One, it filters based on the score, keeping genes strictly above or equal `thr`
#'   (for upregulated) or below or equal `-thr` (for downregulated)
#' - Two, it further narrows down the number of genes keeping only the top-N
#'   in the list, sorted by score (or |score|).
#' - Three, it keeps only genes that are specified in `gene_filter`, discarding
#'   all the rest.
#'
#' To turn off the filters, set `n` to `Inf` (the default), `gene_filter` to
#' `NULL` (the default) and `quantile` to `100`.
extract_top_dysregulated <- function(data, quantile = 5, n = Inf, gene_filter = NULL) {
    types <- list()

    if (!is.null(gene_filter)) {
        cat("Filtering...\n")
        data <- data[row.names(data) %in% gene_filter,]
    }

    thr <- quantile(unlist(data), probs = c(quantile/100, (100-quantile)/100), na.rm = TRUE)
    lower_thr <- thr[1]
    upper_thr <- thr[2]

    for (id in colnames(data)) {
        types[[id]] <- list()

        tum_data <- data |> rownames_to_column("names") |> select(c(!!id, names)) |>
            rename(value = !!id)

        # I use >= so that we at least have 1 gene.
        tum_data <- tum_data |> filter(value >= upper_thr | value <= lower_thr)

        # I filter on 0, but there are no zeroes, due to the thr applied from 
        # before! (well, unless the thr is 0, but...)
        tum_up <- tum_data |> filter(value > 0)
        tum_down <- tum_data |> filter(value <= 0)

        types[[id]][["down"]] <- sort_by(tum_down$names, tum_down$value) |> head(n = n)
        types[[id]][["up"]] <- sort_by(tum_up$names, decreasing=TRUE, tum_up$value) |> head(n = n)
    }

    types
}


{
    #' This does `t(data_frame(...))` but properly keeps the rownames and colnames
    #' put, which `t()` does not do (since it casts the DF to a matrix).
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
        transpose_dataframe(
            data.frame(some = c(1, 2, 3), stuff = c(5, 6, NA), row.names=c("A", "B", "C"))
        ),
        data.frame(A=c(1,5), B=c(2,6), C=c(3,NA), row.names= c("some", "stuff"))
    ))
}

{
    #' This is the only function from UpSetR that I use, so I stole it
    #' I had to edit it in the end, so good thing that I stole it.
    #' It's 'fromList' originally, converting lists of IDs to a matrix of logical
    #' values.
    #'
    #' I edited this so much that now I probably do it backwards: i first set
    #' the values to 0 and 1, and then I update them with real values.
    #' Perhaps there's a better way?
    #'
    #' This takes a list of lists (such as that produced by top_dys) with up and
    #' down regulated genes, and the original matrix with real values, and
    #' creates a matrix suitable for upset plots.
    #'
    #' See the assert statements below
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

    update_to_real <- function(upset_data, values) {

        # TODO: This looks a lot like matrix multiplication to me.
        # The upset_data is a matrix, and so is values.
        # We could just multiply the two and obtain this same effect.
        # Maybe look into it?
        for (col in colnames(upset_data)) {
            if (! col %in% names(values)) {
                next
            }

            rep_positions <- upset_data[[col]] == 1
            upset_data[[col]][rep_positions] <- values[row.names(upset_data), col][rep_positions]
        }

        upset_data
    }

    assert_that(
        are_equal(
            prep_for_upset(
                list(
                    test_tumor = c("A", "B", "C"),
                    other_tumor = c("A", "D")
                )
            ) |> update_to_real(
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
}


#' This does the heavy lifting of converting the prepped data to
#' data suitable for plotting. It applies a series of transformations
#' and re-binds up- and down-regulated genes together to allow for
#' ggplot to plot them.
#'
#' If values is not `NULL`, uses `update_to_real` to update the upset data
#' with the values in `values`. Otherwise, leaves them binarized.
gen_plot_data <- function(data, values = NULL) {
    filtered_data_up <- list()
    filtered_data_down <- list()
    for (name in names(data)) {
        filtered_data_up[[name]] <- data[[name]][["up"]]
        filtered_data_down[[name]] <- data[[name]][["down"]]
    }

    # TODO: This should not work. It only works bc there are no genes that are up
    # in some tt and down in another.

    if (! is.null(values)) {
        prep_for_upset(filtered_data_up) |> update_to_real(values) -> dt_up
        prep_for_upset(filtered_data_down) |> update_to_real(values) -> dt_down
    } else {
        prep_for_upset(filtered_data_up) -> dt_up
        prep_for_upset(filtered_data_down) -> dt_down
    }
    # If some genes are in "up" but not in "down", add them in, and vice-versa


    dt_up$direction <- "up"
    dt_down$direction <- "down"

    dt_up |> rownames_to_column("name") -> dt_up
    dt_down |> rownames_to_column("name") -> dt_down

    dt <- rbind(dt_up, dt_down)
    assert_that(are_equal(nrow(dt), nrow(dt_up) + nrow(dt_down)))

    dt
}

plot_dysregulation <- function(data, min_intersection_size = 5) {
    dt <- gen_plot_data(data)

    suppressWarnings({
        p <- upset(
            dt, names(data),
            min_size = min_intersection_size,
            base_annotations = list(
                'Intersection size' = intersection_size(
                    mapping = aes(fill = direction)
                )
            ),
            set_sizes=FALSE
        )
    })


    p <- p + ggtitle(paste0(
        "Overlap between top ",
        length(data[[1]][[1]]), " genes"
    ))

    p
}

{
    strip_tick <- function(x) {
        str_replace_all(x, "_", " ")
    }

    assert_that(are_equal(strip_tick("some_nice_text"), "some nice text"))


    make_tumor_type_renamer <- function(tt_map) {
        wrap <- function(x) {
            if (x %in% names(tt_map)) {
                return(tt_map[[x]])
            }
            return(strip_tick(x))
        }

        wrap
    }

    assert_that(are_equal(
        make_tumor_type_renamer(
            list("original" = "new")
        )("original"),
        "new"
    ))
    assert_that(are_equal(
        make_tumor_type_renamer(
            list("original"= "new")
        )("wow_how_nice"),
        strip_tick("wow_how_nice")
    ))
}

#' Plot the shared genes map
#'
#' This is probably better handled by ComplexHeatmap, but it grew organically
#' so I'm not changing it at this time.
plot_shared_genes <- function(
        data,
        data_renames,
        values,
        types_renames_fn = NULL,
        n_genes = 50
    ) {
    dt <- gen_plot_data(data, values) |> tibble() |> arrange(name)
    gene_order <- rowSums(abs(values)) |> sort(decreasing = TRUE)

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
    
    get_score_of <- function(name) {
        gene_order[name]
    }

    bar_totals <- bar_dt |> group_by(name) |> summarise(total = sum(abs(row_value)))
    bar_sums <- bar_dt |> group_by(name) |> reframe(order = get_score_of(name))

    bar_dt <- list(bar_dt, bar_sums, bar_totals) |>
        purrr::reduce(\(x, y) {merge(x, y, by="name", all=TRUE)}) |>
        distinct()
    # Here we define the order of both plots
    # Sort first by "total", then break ties by "order"
    bar_dt <- bar_dt[order(bar_dt$total, bar_dt$order, decreasing = TRUE), ]

    # TODO: This is bad - if there are genes with =/= ensg but the same hgnc,
    # this causes a collision. But it's impossible (?) to revalue the labels
    # even with scale_x_discrete (or at least I could not get it to work).
    # I plot hgnc_symbols directly and hope for the best
    
    # I need to keep the first n_genes unique genes. Some may repeat, so I
    # can't just bar_dt[1:n_genes, ] directly
    selection <- unique(bar_dt$name)[1:n_genes]
    bar_dt <- bar_dt[bar_dt$name %in% selection, ]

    bar_plot <- ggplot(bar_dt, aes(x=factor(hgnc_symbol, levels=rev(unique(hgnc_symbol))), y=row_value)) +
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
        scale_x_discrete(position = "top", breaks=bar_dt$hgnc_symbol) +
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

    renames <- if (!is.null(types_renames_fn)) {
        sapply(levels(dot_dt$variable), types_renames_fn)
    } else {
        levels(dot_dt$variable)
    }

    max_value_color <- max(abs(dot_dt$value))

    dot_plot <- ggplot(
        dot_dt,
        aes(
            y=hgnc_symbol,
            x=variable,
            fill=value
        )) +
        scale_fill_gradient2(
            low = "purple", mid="white", high="darkorange",
            name = "Score",
            limits = c(-max_value_color, max_value_color),
            n.breaks = 5
        ) +
        geom_tile() +
        theme_minimal() +
        scale_x_discrete(labels = renames) +
        theme(
            axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
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

if (!exists("LOCAL_DEBUG")) {
    # Parsing arguments
    requireNamespace("argparser")
    
    parser <- argparser::arg_parser("Plot a shared dysregulation plot")
    
    parser |>
        argparser::add_argument(
            "output_file",
            help = "Path to the output file",
            type = "character"
        ) |>
        argparser::add_argument(
            "input_result",
            help = "Input results to parse", type = "character"
        ) |>
        argparser::add_argument(
            "ensg_to_hugo",
            help = "Map from ENSGs to HUGO symbols",
            type = "character",
            default = NULL,
        ) |>
        argparser::add_argument(
            "--selected_genes",
            help = "File with comma-separated list of ENSGs to select before plotting",
            type = "character",
            default = NULL,
        ) |>
        argparser::add_argument(
            "--id_col",
            help = "Name of the column holding ENSG IDs",
            type = "character",
            default = "sample",
        ) |>
        argparser::add_argument(
            "--quantile",
            help = "Quantile to cut the expression at to consider up- (>100-quantile) or down- (<quantile) regulated",
            type = "numerical",
            default = 5,
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
            default = 9, type = "numerical"
        ) |>
        argparser::add_argument(
            "--height",
            help = "Plot height, in inches.",
            default = 16, type = "numerical"
        ) |>
        argparser::add_argument(
            "--renames",
            help = "JSON file with renames to apply to sample names when plotting",
            default = NULL
        ) -> parser
    
    args <- argparser::parse_args(parser)
}

main <- function(args) {
    data <- read_csv(
        args$input_result,
        show_col_types = FALSE
    )

    selected_genes <- if (!is.null(args$selected_genes)) {
        read_file(args$selected_genes) |>
            str_split_1(",") |> str_remove_all("\"") |> str_remove_all("\\n")
    } else {
        NULL
    }

    if (!is.na(args$renames)) {
        tt_renames <- jsonlite::fromJSON(read_lines(args$renames))
    } else {
        # Don't judge me - it's the only rename I need for the TCGA.
        # Just so I don't have to make a new config file just for TCGA.
        # Don't judge me.
        tt_renames <- list(
            "Head_n_Neck_cancer_deseq" = "Head and Neck"
        )
    }
    
    ensg_data <- read_csv(args$ensg_to_hugo)
    
    pdata <- prep_data(data, args$id_col)
    
    top_dys <- extract_top_dysregulated(pdata, quantile=args$quantile, gene_filter = selected_genes)
    
    #plot_dysregulation(top_dys, 3)
    
    if (args$png) {
        png(filename = args$output_file, width=args$width, height = args$height, units = "in", res=args$res)
    } else {
        pdf(file = args$output_file, width=args$width, height = args$height)
    }
    plot_shared_genes(top_dys, ensg_data, pdata, types_renames_fn = make_tumor_type_renamer(tt_renames))
    dev.off()
}

main(args)
