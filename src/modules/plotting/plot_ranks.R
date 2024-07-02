library(compare.ranks)
library(tidyverse)

files <- list.files("data/extracted_results", pattern = "deas", full.names = TRUE)

data <- lapply(files, read_csv)

selected_tumor <- lapply(data, \(x) x[c("sample", "Colorectal_cancer_deseq")])
ranks <- lapply(selected_tumor, \(x) {x$sample[order(x$Colorectal_cancer_deseq)]} )

renames <- list(
    "deseq_shrinkage_deas" = "shrinkage",
    "norm_bws_test_deas" = "BWS",
    "norm_cohen_d_deas" = "Cohen's_D",
    "norm_fold_change_deas" = "FC",
    "norm_s2n_ratio_deas" = "S2N"
)
names(ranks) <- str_split_i(files, "\\/", 3) |> str_split_i("\\.", 1) |> sapply(\(x) renames[[x]])

set.seed(1)
ensg_of_interest <- read_file("data/filter_genes.txt") |>
    str_split_1(",") |> str_remove_all("\"") |> str_remove_all("\\n") |> sample(size = 250)

strip_ensg_version <- function(x) {
    str_split_1(x, "\\.")[1]
}

colors <- c("#FF5733", "#33FF57", "#3357FF", "#FF33A1", "#A133FF", "#33FFA1", "#A1FF33", "gray", "#FF5733", "#5733A1", "#A1A133", "#A133A1", "#FFA133", "#33A1FF", "#FF33A1", "#33FF33", "#FF33FF", "#33FF57", "#5733A1", "black")

selected_ranks <- lapply(ranks, \(x) {sapply(x, strip_ensg_version) |> Filter(f = \(x) x %in% ensg_of_interest)})

comparisons <- compare_two_way_ranks(selected_ranks)
names(comparisons) <- sapply(names(comparisons), \(x) if (x != "step_fraction") {str_replace_all(x, "_", " ")} else {x} )

p <- plot_continuous_congruency(comparisons) + scale_color_manual(values=colors)
print(p)

png("data/out/plots/continuous_congruency_colorectal.png", width=13, height=9, units="in", res=400)
print(p)
dev.off()

