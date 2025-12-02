
library(dplyr)
library(ggplot2)
library(Seurat)
library(Signac)
library(patchwork)
library(SeuratDisk)

setwd("/mnt/research/Rama_lab/Aastha_work/ABC_atlas/data/processed")
# 1. Read the expression matrix
counts <- ReadMtx(
  mtx = "hippocampus_expr.mtx",
  features = "genes.tsv",
  cells = "barcodes.tsv",
  feature.column = 1
)

counts

# 2. Read metadata
meta <- read.csv("hippocampus_selected_classes_metadata.csv",
                 stringsAsFactors = FALSE)

# IMPORTANT: make sure the column you used in barcodes.tsv is here
# You used `selected_labels_in_order = obs_index[mask]` (cell labels),


# so we expect a 'cell_label' column in meta:
stopifnot("cell_label" %in% colnames(meta))

# 3. Align metadata rows to columns of counts
rownames(meta) <- meta$cell_label

# Reorder / subset rows to match the barcodes used in counts
meta <- meta[colnames(counts), ]

# Sanity check: rownames(meta) should equal colnames(counts)
stopifnot(all(rownames(meta) == colnames(counts)))

# 4. Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = counts,
  meta.data = meta
)

# 5. Set default identity to ABC class (if 'class' column exists)
if ("class" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj$cell_class <- seurat_obj@meta.data$class
  Idents(seurat_obj) <- "cell_class"
}

print(seurat_obj)
table(Idents(seurat_obj))

# 6. Save for later
saveRDS(seurat_obj, file = "seurat_hippocampus_subset.rds")
