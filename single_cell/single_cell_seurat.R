# Load libraries--------------------------------------------------------------------------------
library(Seurat)        
library(SeuratObject)
library(tidyverse)
library(ggplot2)
library(sctransform)
library(patchwork)     # arrange plots
library(DoubletFinder)
library(glmGamPoi)
library(future)
library(dplyr)

# Load data----------------------------------------------------------------------------------
raw.data <- Read10X(
  data.dir = "filtered_feature_bc_matrix/"
)

# Inspect structure and dimensions of raw matrix------------------------------------------------
str(raw.data)
dim(raw.data)

# Create Seurat object-----------------------------------------------------------------------
se_ob <- CreateSeuratObject(
  counts    = raw.data,
  project   = "PBMC_10k",
  min.cells = 3,     
  min.features = 200 
)

cat(sprintf("  Cells: %d | Genes: %d\n", ncol(se_ob), nrow(se_ob)))

# Quality Control (QC)-----------------------------------------------------------------------

# Percentage of mitochondrial gene expression
se_ob[["percent.mt"]] <- PercentageFeatureSet(se_ob, pattern = "^MT-")

# Percentage of ribosomal genes
se_ob[["percent.rb"]] <- PercentageFeatureSet(se_ob, pattern = "^RP[SL]")

# Log10 genes per UMI
se_ob[["log10GenesPerUMI"]] <- log10(se_ob$nFeature_RNA) / log10(se_ob$nCount_RNA)

View(se_ob@meta.data)

# QC Visualisations-----------------------------------------------------------------------------
# Distribution of QC metrics across cells
VlnPlot(
  se_ob,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
  ncol = 4
)

# Relationship between counts and detected genes
FeatureScatter(se_ob, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") 

# Relationship between counts and mitochondrial content
FeatureScatter(se_ob, feature1 = "nCount_RNA", feature2 = "percent.mt") 

# QC Filtering----------------------------------------------------------------------------------
cat(" Cells before filtering :", ncol(se_ob))

se_ob <- subset(se_ob, subset = 
                  nFeature_RNA > 200 & 
                  nFeature_RNA < 5000 &
                  nCount_RNA < 50000 &
                  percent.mt < 20 &
                  percent.rb < 40 &
                  log10GenesPerUMI > 0.80)

cat(" Cells after filtering :", ncol(se_ob))

# Doublet Detection (DoubletFinder)-------------------------------------------------------------
se_ob_tmp  <- NormalizeData(se_ob) |>
              FindVariableFeatures() |>
              ScaleData() |>
              RunPCA(npcs = 20)

# Parameter sweep to determine optimal pK value
sweep.res <- paramSweep(se_ob_tmp, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn     <- find.pK(sweep.stats)

# Select pK with highest BCmetric
optimal_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
optimal_pK

# Estimate expected number of doublets
nExp <- round(0.075 * ncol(se_ob_tmp))

# Run DoubletFinder
se_ob_tmp <- doubletFinder(
  se_ob_tmp, 
  PCs = 1:20, 
  pN = 0.25,
  pK = optimal_pK, 
  nExp = nExp, 
  reuse.pANN = NULL,
  sct = FALSE)

# Identify classification column
doublet_col <- grep("DF.classifications", colnames(se_ob_tmp@meta.data), value = TRUE)

table(se_ob_tmp@meta.data[[doublet_col]])

# Remove predicted doublets
se_ob <- se_ob[, se_ob_tmp@meta.data[[doublet_col]] == "Singlet"]

cat("Cells after doublet removal: ", ncol(se_ob))

# Normalization (SCTransform)-------------------------------------------------------------------
se_ob <- SCTransform(
  se_ob,
  vars.to.regress = c("percent.mt"),
  verbose         = FALSE
)

# Dimensionality Reduction (PCA)----------------------------------------------------------------
se_ob <- RunPCA(se_ob,
               features = VariableFeatures(object = se_ob),
               npcs     = 50,
               verbose  = FALSE)

# Elbow plot to determine number of informative PCs
ElbowPlot(se_ob, ndims = 50)

# Top contributing genes per PC
print(se_ob[["pca"]], dims = 1:10, nfeatures = 5)

# Heatmap of top genes for first 10 PCs
DimHeatmap(se_ob, dims = 1:10, cells = 500, 
           balanced = TRUE, fast = FALSE) & 
  theme(axis.text.y = element_text(size = 6))

# Clustering------------------------------------------------------------------------------------

# Construct nearest-neighbor graph
se_ob <- FindNeighbors(se_ob, dims = 1:30, verbose = FALSE)

# Test multiple clustering resolutions
se_ob <- FindClusters(se_ob, resolution = c(0.3, 0.5, 0.7, 1.0), verbose = FALSE)

# Visualize clustering at different resolutions
DimPlot(se_ob, group.by = "SCT_snn_res.0.3") 
DimPlot(se_ob, group.by = "SCT_snn_res.0.5")
DimPlot(se_ob, group.by = "SCT_snn_res.0.7")
DimPlot(se_ob, group.by = "SCT_snn_res.1")

# Select final clustering resolution
se_ob <- FindClusters(se_ob, resolution = 0.5, verbose = FALSE)

# Run UMAP using selected PCs-------------------------------------------------------------------
se_ob <- RunUMAP(se_ob, dims = 1:30, verbose = FALSE)

# Visualize clusters in 2D space
DimPlot(se_ob, reduction = "umap", label = TRUE)

# Disable parallelization to avoid memory issues------------------------------------------------
options(future.globals.maxSize = 8000 * 1024^2)
plan("sequential")

# Marker Gene Identification--------------------------------------------------------------------

# Identify cluster-specific marker genes
all_markers <- FindAllMarkers(
  se_ob,
  only.pos  = TRUE,  # positive markers only
  min.pct   = 0.25,  # expressed in ≥25% of cells in cluster
  logfc.threshold = 0.25,
  verbose   = FALSE
)

# Select top markers with strong fold change
sig_markers <- all_markers |>
  group_by(cluster) |>
  dplyr::filter(avg_log2FC > 1) |>        # keep only strong markers
  top_n(5, wt = avg_log2FC) |>            # then take top 5 per cluster
  ungroup()

# Summary of marker genes
cat(sprintf("  High-confidence markers (log2FC > 1): %d genes across %d clusters\n",
            nrow(sig_markers), length(unique(sig_markers$cluster))))

# Visualize marker expression
DotPlot(se_ob, features = unique(sig_markers$gene)) + RotatedAxis()

# Annotate cell types---------------------------------------------------------------------------

# Assign biological identities to clusters
new_labels <- c(
  "0"  = "Classical Monocytes",
  "1"  = "CD4+ Naive T cells",
  "2"  = "Non-classical Monocytes",
  "3"  = "Naive B cells",
  "4"  = "NK cells",
  "5"  = "MAIT cells",
  "6"  = "Memory B cells",
  "7"  = "CD8+ T cells",
  "8"  = "Dendritic cells",
  "9"  = "Plasma cells",
  "10" = "Platelets",
  "11" = "Gamma-delta T cells",
  "12" = "Uncertain",
  "13" = "CD8+ Effector T cells",
  "14" = "HSC"
)

# Store cluster identities
se_ob <- RenameIdents(se_ob, new_labels)
se_ob$cell_type <- Idents(se_ob)
View(se_ob@meta.data)

# UMAP with annotated cell types
DimPlot(se_ob, reduction = "umap", label = TRUE, 
        repel = TRUE, pt.size = 0.5) +
  ggtitle("PBMC 10k - Cell Type Annotation")

# Feature Visualization-------------------------------------------------------------------------
FeaturePlot(se_ob, 
            features = c("CD14", "CD3D", "CD79A", 
                         "NKG7", "FCGR3A", "IL7R"),
            ncol = 3)
