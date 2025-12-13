library(Seurat)
library(Signac)
library(zellkonverter)
library(GenomeInfoDb)
library(GenomicRanges)
library(SeuratExtend)
library(ggrepel)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(patchwork)
library(EnsDb.Mmusculus.v79)
library(Matrix)
library(SingleCellExperiment)
library(dplyr)
library(knitr)


# Working directory and data ------------------------------------------------

setwd("/mnt/research/Rama_lab/Aastha_work/ABC_atlas/data/processed")

abc_hip <- readRDS("seurat_hippocampus_subset.rds")
str(abc_hip)
dim(abc_hip)
DefaultAssay(abc_hip) <- "RNA"

#1.  Quality control -----------------------------------------------------------

# Calculate % mitochondrial genes
abc_hip[["percent.mt"]] <- PercentageFeatureSet(abc_hip, pattern = "^mt-")

# Pre-filter violin plots
p1 <- VlnPlot(abc_hip,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3)
ggsave("prefilter_violinplot.png", p1, width = 10, height = 4, dpi = 300)

# Pre-filter scatter plots
p_scatter1 <- FeatureScatter(abc_hip, "nCount_RNA", "percent.mt")
p_scatter2 <- FeatureScatter(abc_hip, "nCount_RNA", "nFeature_RNA")
ggsave("prefilter_scatterplots.png",
       p_scatter1 + p_scatter2,
       width = 10, height = 4, dpi = 300)


#2. Filtering -----------------------------------------------------------------

abc_hip <- subset(
  abc_hip,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 11000 &
    percent.mt < 10
)

# Post-filter QC plots
p2 <- VlnPlot(abc_hip,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
              ncol = 3)
ggsave("postfilter_violinplot.png", p2, width = 10, height = 4, dpi = 300)

p_scatter3 <- FeatureScatter(abc_hip, "nCount_RNA", "percent.mt")
p_scatter4 <- FeatureScatter(abc_hip, "nCount_RNA", "nFeature_RNA")
ggsave("postfilter_scatterplots.png",
       p_scatter3 + p_scatter4,
       width = 10, height = 4, dpi = 300)


#3.  Normalization, and Feature selection --------------------------------------

abc_hip <- NormalizeData(abc_hip)
abc_hip <- FindVariableFeatures(abc_hip, selection.method = "vst", nfeatures = 5000)
abc_hip <- ScaleData(abc_hip)
abc_hip <- RunPCA(abc_hip, features = VariableFeatures(abc_hip))

# Elbow plot
p_elbow <- ElbowPlot(abc_hip)
ggsave("elbowplot.png", p_elbow, width = 6, height = 4, dpi = 300)


#4.  Clustering and UMAP -------------------------------------------------------

abc_hip <- FindNeighbors(abc_hip, dims = 1:15)
abc_hip <- FindClusters(abc_hip)
abc_hip <- RunUMAP(abc_hip, dims = 1:15)

# UMAP colored by subclass, class
umap_subclass <- DimPlot2(abc_hip, reduction = "umap", features = "subclass")
umap_class    <- DimPlot2(abc_hip, reduction = "umap", features = "class")
ggsave("dimplot_subclass.png", umap_subclass, width = 8, height = 6, dpi = 300)
ggsave("dimplot_class.png",    umap_class,    width = 8, height = 6, dpi = 300)

# Mixed metadata UMAP
DimPlot2(abc_hip, features = c("neurotransmitter", "class", "subclass"),
         theme = NoAxes())
ggsave("dimplot_mix.png", width = 10, height = 8, dpi = 300)

# Cluster-labelled UMAP
cluster_umap <- DimPlot(abc_hip, label = TRUE, repel = TRUE)
ggsave("dimplot_cluster.png", cluster_umap, width = 10, height = 8, dpi = 900)

# Save object for downstream work
#saveRDS(abc_hip, "ABC_hip.rds")

abc_hip <- readRDS("ABC_hip.rds")

#5.  Marker discovery ----------------------------------------------------------

hip_markers <- FindAllMarkers(
  abc_hip,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5
)

# Top markers per cluster
top3_markers <- hip_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 3)

heat_top3 <- DoHeatmap(abc_hip, features = top3_markers$gene) + NoLegend()
ggsave("hippocampus_top3markers_heatmap.png",
       heat_top3, width = 12, height = 10)


#6.  Differential expression: Glutamatergic neurons ---------------------------
Idents(abc_hip) <- "class"

glut_markers <- FindMarkers(
  abc_hip,
  ident.1 = "04 DG-IMN Glut",
  ident.2 = NULL,
  test.use = "MAST",
  logfc.threshold = 0.25,
  min.pct = 0.5
)

glut_markers <- glut_markers[glut_markers$p_val_adj < 0.05, ]
glut_markers$gene <- rownames(glut_markers)

# Rank by fold change
glut_sorted <- glut_markers[order(glut_markers$avg_log2FC, decreasing = TRUE), ]

top10_up   <- head(glut_sorted$gene, 10)
top10_down <- tail(glut_sorted$gene, 10)
highlight  <- c(top10_up, top10_down)


#7. Volcano plot ---------------------------------------------------------------

volc <- EnhancedVolcano(
  glut_markers,
  lab = glut_markers$gene,
  selectLab = highlight,
  x = "avg_log2FC",
  y = "p_val_adj",
  pCutoff = 0.05,
  FCcutoff = 0.5,
  drawConnectors = TRUE
)

ggsave("Glutamatergic_VolcanoPlot.png", volc, width = 10, height = 8)


#8. Feature plots -------------------------------------------------------------

feat_plot <- FeaturePlot(abc_hip, features = top10_up[1:4])
ggsave("Glut_top4_upregulated_featureplot.png",
       feat_plot, width = 10, height = 8, dpi = 600)


# Dot plot ---------------------------------------------------------------------

dot <- DotPlot(abc_hip, features = highlight, group.by = "subclass") + RotatedAxis()
ggsave("Glut_upregulated_dotplot.png",dot, width = 10, height = 8, dpi = 600)


#9 Heatmap for DE genes --------------------------------------------------------

DoHeatmap(abc_hip, features = highlight, group.by = "subclass") + NoLegend()




# Complete ATAC-seq processing and Integration pipeline ------------------------

#1 : Read data extracted from the h5ad file using python and creat counts
atac_counts <- ReadMtx(
  mtx = "atac_matrix.mtx",
  features = "atac_features.tsv",
  cells = "atac_barcodes.tsv",
  feature.column = 1
)


#2: Create genomic ranges ------------------------------------------------------
gr <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
# Convert to half-open intervals as bin are 0-500, 500-1000
gr_half_open <- gr
end(gr_half_open) <- end(gr_half_open) - 1
genome(gr_half_open) <- "mm10"
cat("Created", length(gr_half_open), "genomic bins\n\n")


#3: Create chromatin assay and seurat object----------------------------
chrom_assay <- suppressWarnings(
  CreateChromatinAssay(
    counts = atac_counts,
    ranges = gr_half_open,
    genome = "mm10",
    min.cells = 0,
    min.features = 0
  )
)

atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "ATAC"
)

#4 Adding metadata from the extract. ----------------------------
metadata <- read.csv("atac_metadata.csv", row.names = 1)
atac <- AddMetaData(atac, metadata)


#5 Create individual histograms for preprocess visualization -------------------
p1 <- ggplot(atac@meta.data, aes(x = nCount_ATAC)) +
  geom_histogram(bins = 100, fill = "blue", color = "black") +
  theme_minimal() + 
  ggtitle("nCount_ATAC") + 
  xlab("Total counts per cell")

p2 <- ggplot(atac@meta.data, aes(x = nFeature_ATAC)) +
  geom_histogram(bins = 100, fill = "orange", color = "black") +
  theme_minimal() + 
  ggtitle("nFeature_ATAC") + 
  xlab("Features per cell")

p3 <- ggplot(atac@meta.data, aes(x = tsse)) +
  geom_histogram(bins = 100, fill = "green", color = "black") +
  theme_minimal() + 
  ggtitle("TSS Enrichment") + 
  xlab("TSSe")

p4 <- ggplot(atac@meta.data, aes(x = frac_mito)) +
  geom_histogram(bins = 100, fill = "firebrick", color = "black") +
  theme_minimal() + 
  ggtitle("Mitochondrial Fraction") + 
  xlab("Fraction Mitochondrial")

# Combine plots
(p1 | p2) / (p3 | p4)
combined_plot <- (p1 | p2) / (p3 | p4)

# Save the combined plot
ggsave("ATAC_QC_Histograms.png", plot = combined_plot, width = 10, height = 8, dpi = 300)

#6: Add Gene annotations --------------------------------------------------------
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"

#7 Direct assignment to avoid Signac bug ---------------------------------------
atac[["ATAC"]]@annotation <- annotations
atac@misc$gene_annotations <- annotations

#8 : Standard atac pre processing. --------------------------------------------
atac <- RunTFIDF(atac)
atac <- FindTopFeatures(atac, min.cutoff = 10)
atac <- RunSVD(atac)
ElbowPlot(atac, reduction = "lsi")
atac <- RunUMAP(atac, reduction = "lsi", dims = 2:15) #chosen 15 after Elbow plot
atac <- FindNeighbors(atac, reduction = "lsi", dims = 2:15)
atac <- FindClusters(atac, resolution = 0.5)
DimPlot(atac, reduction = "umap", label = TRUE)
f1<- FeaturePlot(atac, "nCount_ATAC")
f2<-FeaturePlot(atac, "nFeature_ATAC")
f3<-FeaturePlot(atac, "tsse")
f1|f2|f3



#9. : Calculating manual gene activity matrix. ----------------------------
# Get peak information
peak_ranges <- granges(atac[["ATAC"]])
peak_counts <- GetAssayData(atac, assay = "ATAC", slot = "counts")

gene_promoters <- promoters(annotations, upstream = 2000, downstream = 0)

# Find overlaps
overlaps <- findOverlaps(peak_ranges, gene_promoters)
gene_names <- gene_promoters$gene_name[subjectHits(overlaps)]
peak_indices <- queryHits(overlaps)

# Remove NAs
valid_indices <- !is.na(gene_names)
gene_names <- gene_names[valid_indices]
peak_indices <- peak_indices[valid_indices]

# Create gene-peak mapping
gene_peak_list <- split(peak_indices, gene_names)
unique_genes <- names(gene_peak_list)

# Initialize gene activity matrix
gene_activity <- Matrix(0, 
                        nrow = length(unique_genes), 
                        ncol = ncol(atac),
                        dimnames = list(unique_genes, colnames(atac)),
                        sparse = TRUE)

# Sum counts for each gene
pb <- txtProgressBar(min = 0, max = length(unique_genes), style = 3)
for(i in seq_along(unique_genes)) {
  gene <- unique_genes[i]
  peak_ids <- gene_peak_list[[gene]]
  
  if(length(peak_ids) > 0) {
    gene_counts <- peak_counts[peak_ids, , drop = FALSE]
    if(nrow(gene_counts) > 0) {
      gene_activity[gene, ] <- Matrix::colSums(gene_counts)
    }
  }
  setTxtProgressBar(pb, i)
}
close(pb)

# Remove zero-count genes
nonzero_rows <- Matrix::rowSums(gene_activity) > 0
gene_activity <- gene_activity[nonzero_rows, ]

# Remove empty gene names
valid_rows <- !(rownames(gene_activity) == "" | is.na(rownames(gene_activity)))
gene_activity <- gene_activity[valid_rows, ]

# Add to Seurat object as GENE assay
cat("   Adding gene activity as GENE assay...\n")
atac[["GENE"]] <- CreateAssayObject(counts = gene_activity)

# Normalize gene activity
DefaultAssay(atac) <- "GENE"
atac <- NormalizeData(atac)
atac <- ScaleData(atac)

#10 Saving intermediate result for future use-----------------------------------
saveRDS(atac, "atac_processed_with_gene_activity.rds")


#------------------------------Starting the integration------------------------


#1. reading saved pre processed RNA seq and atac seq data
atac <- readRDS("atac_processed_with_gene_activity.rds")
rna <- readRDS("ABC_hip.rds")

#2.: Finding variable feature again and assigning default assay for anchor and label transfer
rna <- FindVariableFeatures(rna, nfeatures = 2000)
DefaultAssay(atac) <- "GENE"
DefaultAssay(rna) <- "RNA"

#3 : Identifying common genes between to types
common_genes <- intersect(rownames(atac[["GENE"]]), rownames(rna))
integration_genes <- intersect(common_genes, VariableFeatures(rna))
FeaturePlot(atac,features=c("Neurod1", "Slc17a7", "Gad2"))


#4: Creatig anchor 
anchors <- FindTransferAnchors(
  reference = rna,
  query = atac,
  reference.assay = "RNA",
  query.assay = "GENE",
  features = integration_genes,
  reduction = "cca",
  dims = 1:15
)

#5 Transferring anchor using class from rna-seq data
pred <- TransferData(
  anchorset = anchors,
  refdata = rna$class, 
  weight.reduction = atac[["lsi"]],
  dims = 2:15
)

#6 adding metadata back
atac <- AddMetaData(object = atac, metadata = pred)
DimPlot(atac, group.by="predicted.id", label=TRUE, repel=TRUE)
FeaturePlot(atac, features="prediction.score.max")

plot5 <- DimPlot(rna, group.by = 'class', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
plot6 <- DimPlot(atac, group.by = 'predicted.id', label = TRUE, repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
combined_integrated_plot <-plot5 + plot6
ggsave("combined_integrated_plot.png", plot = combined_integrated_plot, width = 10, height = 8, dpi = 300)

#7 Assigning ATAC as default assaY for identifying markes in atac seq
DefaultAssay(atac) <- "ATAC"
Idents(atac) <- atac$predicted.id

#8 Identifying markers against GABA from GLUt in atac after integration using linear regression
da_peaks <- FindMarkers(
  object = atac,
  ident.1 = "04 DG-IMN Glut",
  ident.2 = NULL,
  test.use = "LR",
  min.pct = 0.05
)



#9 Subsetting result and visualization
all_peaks <- da_peaks[da_peaks$p_val_adj < 0.05 & da_peaks$avg_log2FC > 0,]

#Open chromatin region in glut and gaba
open_glut <- rownames(da_peaks[da_peaks$avg_log2FC > 0 & da_peaks$p_val_adj < 0.05, ])
open_GABA <- rownames(da_peaks[da_peaks$avg_log2FC < 0 & da_peaks$p_val_adj < 0.05, ])

#Identifying closest region of open chromatin
closest_glut <- ClosestFeature(atac, open_glut)
closest_GABA <- ClosestFeature(atac, closed_glut)

#violin plot of first three open chromatin region
VlnPlot(object = atac, 
        features = open_glut[1:3],
        pt.size =0.1,
        idents = c("06 CTX-CGE GABA","04 DG-IMN Glut"))


plot7 <- VlnPlot(
  object = atac,
  features = rownames(all_peaks)[3],
  pt.size = 0.1,
  idents = c("06 CTX-CGE GABA","04 DG-IMN Glut")
)

plot8 <- FeaturePlot(
  object = atac,
  features = rownames(all_peaks)[3],
  pt.size = 0.1,
  max.cutoff = 'q95'
)
peaks_combined <- plot7 | plot8
ggsave("bdnf_open_region_GLUT.png", plot = peaks_combined, width = 10, height = 8, dpi = 300)


#10 Result summary
peak_summary <- tibble(
  Peak = c(open_glut, open_GABA),
  Regulation = c(rep("Open in DG-IMN Glut", length(open_glut)), "Open in CGE-CTX GABA "),
  Gene = c(closest_glut$gene_name, closest_GABA$gene_name),
  Distance = c(closest_glut$distance, closest_GABA$distance)
)
kable(peak_summary, caption = "Peak-to-Gene Summary for DG-IMN Glut Cells")

#11. plotting the reuslt
peak_to_plot <- open_glut
VlnPlot(atac, features=peak_to_plot, group.by="predicted.id")
FeaturePlot(atac,features=peak_to_plot ,group.by="predicted.id")




