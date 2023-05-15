
library(tidyverse)
requireNamespace("gplots")
requireNamespace("reshape2")

temp_dir <- "/home/hedmad/Files/data/transportome_profiler/tests/"
out_files <- list.files(temp_dir, full.names = TRUE)

# Temp vars
files <- out_files

TREE_LABELS <- ".
├── ./carried_solute::anion
├── ./carried_solute::Ca2+
├── ./carried_solute::cation
├── ./carried_solute::Cl-
├── ./carried_solute::HCO3-
├── ./carried_solute::Na+
├── ./pores
│   ├── ./pores/aquaporins
│   └── ./pores/channels
│       ├── ./pores/channels/carried_solute::anion
│       ├── ./pores/channels/carried_solute::Ca2+
│       │   ├── ./pores/channels/carried_solute::Ca2+/gating_mechanism::ligand
│       │   └── ./pores/channels/carried_solute::Ca2+/gating_mechanism::voltage
│       ├── ./pores/channels/carried_solute::cation
│       │   ├── ./pores/channels/carried_solute::cation/gating_mechanism::leakage
│       │   └── ./pores/channels/carried_solute::cation/gating_mechanism::voltage
│       ├── ./pores/channels/carried_solute::Cl-
│       ├── ./pores/channels/carried_solute::K+
│       │   ├── ./pores/channels/carried_solute::K+/gating_mechanism::leakage
│       │   ├── ./pores/channels/carried_solute::K+/gating_mechanism::ligand
│       │   └── ./pores/channels/carried_solute::K+/gating_mechanism::voltage
│       ├── ./pores/channels/carried_solute::Na+
│       │   ├── ./pores/channels/carried_solute::Na+/gating_mechanism::ligand
│       │   └── ./pores/channels/carried_solute::Na+/gating_mechanism::voltage
│       ├── ./pores/channels/gating_mechanism::ligand
│       │   ├── ./pores/channels/gating_mechanism::ligand/carried_solute::Cl-
│       │   └── ./pores/channels/gating_mechanism::ligand/carried_solute::HCO3-
│       └── ./pores/channels/gating_mechanism::voltage
│           ├── ./pores/channels/gating_mechanism::voltage/carried_solute::Ba2+
│           └── ./pores/channels/gating_mechanism::voltage/carried_solute::Cs+
└── ./transporters
    ├── ./transporters/atp_driven
    │   ├── ./transporters/atp_driven/ABC
    │   │   ├── ./transporters/atp_driven/ABC/carried_solute::lipid
    │   │   ├── ./transporters/atp_driven/ABC/direction::in
    │   │   └── ./transporters/atp_driven/ABC/direction::out
    │   ├── ./transporters/atp_driven/carried_solute::cation
    │   ├── ./transporters/atp_driven/carried_solute::lipid
    │   ├── ./transporters/atp_driven/carried_solute::phospholipid
    │   ├── ./transporters/atp_driven/direction::in
    │   │   ├── ./transporters/atp_driven/direction::in/carried_solute::H+
    │   │   └── ./transporters/atp_driven/direction::in/carried_solute::lipid
    │   ├── ./transporters/atp_driven/direction::out
    │   │   ├── ./transporters/atp_driven/direction::out/carried_solute::cation
    │   │   ├── ./transporters/atp_driven/direction::out/carried_solute::H+
    │   │   └── ./transporters/atp_driven/direction::out/carried_solute::phospholipid
    │   └── ./transporters/atp_driven/pumps
    │       ├── ./transporters/atp_driven/pumps/carried_solute::cation
    │       │   ├── ./transporters/atp_driven/pumps/carried_solute::cation/direction::in
    │       │   └── ./transporters/atp_driven/pumps/carried_solute::cation/direction::out
    │       ├── ./transporters/atp_driven/pumps/carried_solute::H+
    │       ├── ./transporters/atp_driven/pumps/carried_solute::phospholipid
    │       ├── ./transporters/atp_driven/pumps/direction::in
    │       └── ./transporters/atp_driven/pumps/direction::out
    ├── ./transporters/carried_solute::amine
    ├── ./transporters/carried_solute::anion
    │   ├── ./transporters/carried_solute::anion/direction::in
    │   └── ./transporters/carried_solute::anion/direction::out
    ├── ./transporters/carried_solute::A*P
    ├── ./transporters/carried_solute::Ca2+
    ├── ./transporters/carried_solute::cation
    │   ├── ./transporters/carried_solute::cation/direction::in
    │   └── ./transporters/carried_solute::cation/direction::out
    ├── ./transporters/carried_solute::H+
    │   ├── ./transporters/carried_solute::H+/direction::in
    │   └── ./transporters/carried_solute::H+/direction::out
    └── ./transporters/solute_carriers
        ├── ./transporters/solute_carriers/carried_solute::amine
        ├── ./transporters/solute_carriers/carried_solute::amino acid
        ├── ./transporters/solute_carriers/carried_solute::cation
        │   ├── ./transporters/solute_carriers/carried_solute::cation/net_charge::0.0
        │   └── ./transporters/solute_carriers/carried_solute::cation/port_type::symport
        ├── ./transporters/solute_carriers/carried_solute::Cl-
        ├── ./transporters/solute_carriers/carried_solute::his
        └── ./transporters/solute_carriers/carried_solute::Na+
            ├── ./transporters/solute_carriers/carried_solute::Na+/net_charge::0.0
            └── ./transporters/solute_carriers/carried_solute::Na+/port_type::symport
"

str_reverse <- function(x) {
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
}

treelabs_to_lab_matrix <- function(treelabs) {
  splitted <- str_split_1(treelabs, "\\n") |> head(-1) # remove last \n
  
  splitted_no_ident <- str_split_i(splitted, "\\.", 2)
  ident <- str_split_i(splitted, "\\.", 1) # keep just the tree on the left
  
  original_path <- paste0("/whole_transportome", splitted_no_ident)
  
  leaves <- str_split_i(original_path, "/", -1)
  
  # We can now compute an "inverted" tree graph
  inverted_ident <- str_reverse(ident)
  # The special characters ├ and └ need to be replaced manually
  inverted_ident |> str_replace_all("├", "┤") |> str_replace_all("└", "┘") -> inverted_ident
  
  result <- data.frame(
    original = splitted,
    idents = paste0(ident, leaves),
    paths = original_path,
    no_idents = splitted_no_ident,
    inverted = paste0(leaves, inverted_ident),
    order = 1:length(splitted) # to preserve the right order
  )
  
  print(result)
}

lab_matrix <- treelabs_to_lab_matrix(TREE_LABELS)

main <- function(files, lab_matrix) {
  data <- map(files, read_csv, show_col_types = FALSE, progress = FALSE)
  
  # Grab the file names
  files |> str_split_i("/", -1) |> str_split_i("\\.", 1) |> str_split_i("_", 1) -> file_names
  
  print(file_names)  
  # For a heatmap we need a matrix.
  # We can build it by combining all the NES values for each pathway
  
  # Set the order of the rownames to be standard
  row_order <- data[[1]][["pathway"]]
  heat_matrix <- matrix(nrow = length(row_order), ncol = length(files))
  colnames(heat_matrix) <- file_names
  rownames(heat_matrix) <- row_order
  for (i in seq_along(files)) {
    frame <- as.matrix(data[[i]])
    row.names(frame) <- frame[, "pathway"]
    frame <- frame[row_order, ]
    heat_matrix[, i] <- as.numeric(frame[, "NES"])
  }
  
  value_vars <- colnames(heat_matrix)
  
  # Fuse with label matrix
  heat_matrix <- heat_matrix[match(row.names(heat_matrix), lab_matrix$paths), ]
  
  heat_matrix[is.na(heat_matrix)] <- 0
  
  gplots::heatmap.2(
    heat_matrix,
    col = colorRampPalette(c("blue", "gray", "red"))(n = 100),
    dendrogram = "none",
    trace = "none",
    labRow = lab_matrix$idents,
    Rowv = lab_matrix$order
  )
}

main(out_files, lab_matrix)


main_ggplot <- function(files, lab_matrix) {
  data <- map(files, read_csv, show_col_types = FALSE, progress = FALSE)
  
  # Grab the file names
  files |> str_split_i("/", -1) |> str_split_i("\\.", 1) |> str_split_i("_", 1) -> file_names
  
  print(file_names)  
  # For a heatmap we need a matrix.
  # We can build it by combining all the NES values for each pathway
  
  # Set the order of the rownames to be standard
  row_order <- data[[1]][["pathway"]]
  heat_matrix <- matrix(nrow = length(row_order), ncol = length(files))
  colnames(heat_matrix) <- file_names
  rownames(heat_matrix) <- row_order
  for (i in seq_along(files)) {
    frame <- as.matrix(data[[i]])
    row.names(frame) <- frame[, "pathway"]
    frame <- frame[row_order, ]
    heat_matrix[, i] <- as.numeric(frame[, "NES"])
  }
  
  # calculate dendrogram of cols
  dendro_order <- hclust(dist(t(heat_matrix)), method = "ward.D2")$order
  print(dendro_order)
  
  value_vars <- colnames(heat_matrix)
  
  # Fuse with label matrix
  heat_matrix <- heat_matrix[match(row.names(heat_matrix), lab_matrix$paths), ]
  heat_matrix <- cbind(heat_matrix, lab_matrix)
  
  heat_matrix |> reshape2::melt(
    measure.vars = value_vars
  ) -> melted_data

  p <- ggplot(melted_data, aes(x = variable, y = reorder(inverted, order, decreasing = TRUE), fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "gray", high = "red") + 
    scale_x_discrete(limits = colnames(heat_matrix)[dendro_order]) +
    # Rotate labels
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle("A heatmap") + xlab("Tumor Type") + ylab("Enrichment Type")
  
  # legend tile
  
  print(p)
  
  invisible(p)
}

main_ggplot(out_files, lab_matrix)


