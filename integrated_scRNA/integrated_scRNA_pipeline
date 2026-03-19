# Load libraries-----------------------------------------------------------------
library(Seurat)
library(harmony)
library(Matrix)
library(future)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(glmGamPoi)
library(pheatmap)

# Set folder path-------------------------------------------------------------------------
folder_path <- "C:/Users/Admin/Desktop/scRNA"

# List all .rds files---------------------------------------------------------------------
files <- list.files(folder_path, pattern = "\\.rds$", full.names = TRUE)

# Load all files simultaneously-----------------------------------------------------------
data_list <- lapply(files, readRDS)

# View structure of loaded files----------------------------------------------------------
for (i in seq_along(data_list)) {
  cat("=== Structure of", names(data_list)[i], "===\n")
  str(data_list[[i]], max.level = 2)
  cat("\n")
}

# Extract exon----------------------------------------------------------------------------
exon_list <- lapply(data_list, function(x) x$exon)

# Define sample information---------------------------------------------------------------
sample_info <- data.frame(
  file = c("GSM4557327_555_1_cell.counts.matrices.rds",
           "GSM4557328_555_2_cell.counts.matrices.rds",
           "GSM4557329_556_cell.counts.matrices.rds",
           "GSM4557330_557_cell.counts.matrices.rds",
           "GSM4557331_558_cell.counts.matrices.rds",
           "GSM4557332_559_cell.counts.matrices.rds",
           "GSM4557333_561_cell.counts.matrices.rds",
           "GSM4557334_HIP002_cell.counts.matrices.rds",
           "GSM4557335_HIP015_cell.counts.matrices.rds",
           "GSM4557336_HIP023_cell.counts.matrices.rds",
           "GSM4557337_HIP043_cell.counts.matrices.rds",
           "GSM4557338_HIP044_cell.counts.matrices.rds",
           "GSM4557339_HIP045_cell.counts.matrices.rds"),
  sample = c("COVID_555_1", "COVID_555_2", "COVID_556",
             "COVID_557", "COVID_558", "COVID_559", "COVID_561",
             "Healthy_HIP002", "Healthy_HIP015", "Healthy_HIP023",
             "Healthy_HIP043", "Healthy_HIP044", "Healthy_HIP045"),
  condition = c(rep("COVID", 7), rep("Healthy", 6))
)

# Create Seurat objects for each sample---------------------------------------------------
seurat_list <- lapply(1:nrow(sample_info), function(i) {
  mat <- readRDS(sample_info$file[i])
  so <- CreateSeuratObject(
    counts = mat$exon,
    project = sample_info$sample[i],
    min.cells = 3,
    min.features = 200
  )
  so$condition <- sample_info$condition[i]  # COVID ya Healthy label
  so$sample <- sample_info$sample[i]
  so
})

# Name each Seurat object in the list by sample-------------------------------------------
names(seurat_list) <- sample_info$sample
seurat_list  

# QC metrics------------------------------------------------------------------------------
seurat_list <- lapply(seurat_list, function(so) {
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  so[["percent.rb"]] <- PercentageFeatureSet(so, pattern = "^RP[SL]")
  so[["log10GenesPerUMI"]] <- log10(so$nFeature_RNA) / log10(so$nCount_RNA)
  so
})

# Plot QC metrics-----------------------------------------------------------------
pdf("QC_all_samples.pdf", width = 14, height = 10)

for (sample in names(seurat_list)) {
  so <- seurat_list[[sample]]
  
  # Violin plot
  p1 <- VlnPlot(so,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3) + 
    plot_annotation(title = sample)
  
  # Scatter plots
  p2 <- FeatureScatter(so, "nCount_RNA", "nFeature_RNA") +
    geom_smooth(method = "lm", se = FALSE)
  p3 <- FeatureScatter(so, "nCount_RNA", "percent.mt")
  
  print(p1)
  print(p2 + p3 + plot_annotation(title = sample))
}

dev.off()

# Cells filtering ----------------------------------------------------------------------------
cat("Cells before filtering:\n")
print(sapply(seurat_list, ncol))

seurat_list <- lapply(seurat_list, function(so) {
  subset(so, subset =
           nFeature_RNA    > 200  &
           nFeature_RNA    < 2500 &
           nCount_RNA      < 15000 &
           percent.mt      < 20   &
           log10GenesPerUMI > 0.80)
})

cat("Cells after filtering:\n")
cell_counts <- sapply(seurat_list, ncol)
print(cell_counts)
cat("Total cells:", sum(cell_counts), "\n")

# Merge all filtered Seurat objects into one combined object----------------------------------
covid_merged <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project      = "COVID_study"
)

covid_merged

# Normalization using SCTransform-------------------------------------------------------------
covid_merged <- SCTransform(covid_merged, 
                            vars.to.regress = "percent.mt",
                            method = "glmGamPoi",
                            verbose = FALSE)

# Dimensionality Reduction (PCA)--------------------------------------------------------------
covid_merged <- RunPCA(covid_merged, npcs = 30, verbose = FALSE)

covid_merged

# Elbow Plot------------------------------------------------------------------------------------
ElbowPlot(covid_merged)

# Harmony Integration- correct batch effects across sample--------------------------------------
covid_merged <- IntegrateLayers(
  covid_merged,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  grouping.var = "sample",
  verbose = FALSE
)

# Clustering- build neighbor graph and identify clusters----------------------------------------
covid_merged <- FindNeighbors(covid_merged, reduction = "harmony", dims = 1:20)
covid_merged <- FindClusters(covid_merged, resolution = 0.5, verbose = FALSE)

# UMAP- compute 2D visualization--------------------------------------------------------------- 
covid_merged <- RunUMAP(covid_merged, reduction = "harmony", dims = 1:30)

# UMAP Visualization----------------------------------------------------------------------------
# Plot 1 - color cells by condition (COVID vs Healthy)
DimPlot(covid_merged, reduction = "umap", group.by = "condition", pt.size = 0.3)

# Plot 2 - color cells by cluster and label clusters
DimPlot(covid_merged, reduction = "umap", label = TRUE, pt.size = 0.3)

# Setup for parallel processing (optional, here using sequential for safety)-------------
plan("sequential")
options(future.globals.maxSize = 8000 * 1024^2)

# Prepare Seurat object for SCT-based marker detection-----------------------------------
covid_merged <- PrepSCTFindMarkers(covid_merged, verbose = FALSE)

# Find all cluster-specific markers------------------------------------------------------
all_markers <- FindAllMarkers(
  covid_merged,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = FALSE
)

write.csv(all_markers, "all_markers_per_cluster.csv", row.names = FALSE)

# Select top 3 markers per cluster for visualization-------------------------------------
markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC) %>%
  arrange(cluster) %>%
  select(cluster, gene, avg_log2FC, pct.1, pct.2) #%>%
  #print(n = 70)

View (markers)
write.csv(markers, "top_markers_per_cluster.csv", row.names = FALSE)

cat("Total markers found:", nrow(all_markers), "\n")
cat("Clusters:", length(unique(all_markers$cluster)), "\n")

# Check how many rows in top markers
nrow(markers)  # dekho kitni rows hain

# Print all selected markers in console
print(markers, n = nrow(markers))

# Define cell type labels for each cluster (manual annotation based on marker genes)------------
new_labels <- c(
  "0"  = "CD4+ Naive T cells",
  "1"  = "CD8+ T cells",
  "2"  = "NK cells",
  "3"  = "Classical Monocytes",
  "4"  = "Low Quality",
  "5"  = "Naive B cells",
  "6"  = "Activated Monocytes",
  "7"  = "Erythrocytes",
  "8"  = "Non-classical Monocytes",
  "9"  = "Plasma cells IgG",
  "10" = "Memory B cells",
  "11" = "Plasma cells IgM",
  "12" = "Platelets",
  "13" = "Proliferating cells",
  "14" = "Plasma cells IgG2",
  "15" = "Plasma cells IgA",
  "16" = "Neutrophils",
  "17" = "Dendritic cells",
  "18" = "Plasma cells mixed",
  "19" = "Plasmacytoid DC",
  "20" = "Basophils",
  "21" = "Stressed cells",
  "22" = "Progenitor cells"
)

# Rename cluster IDs to biologically meaningful cell type labels--------------------------------
Idents(covid_merged) <- "seurat_clusters"
covid_merged <- RenameIdents(covid_merged, new_labels)
covid_merged$cell_type <- Idents(covid_merged)

# Visualize annotated cell types on UMAP--------------------------------------------------------
DimPlot(covid_merged, reduction = "umap", 
        label = TRUE, repel = TRUE, pt.size = 0.3) +
  ggtitle("COVID vs Healthy - Cell Types")

# Set cell type as the active identity class----------------------------------------------------
Idents(covid_merged) <- "cell_type"

# Calculate proportion of each cell type within each condition (COVID vs Healthy)---------------
prop_table <- prop.table(table(covid_merged$cell_type, 
                               covid_merged$condition), 
                         margin = 2)

# Plot stacked barplot showing cell type composition per condition
as.data.frame(prop_table) %>%
  ggplot(aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity") +
  labs(x = "Condition", y = "Proportion", 
       fill = "Cell Type",
       title = "Cell Type Proportions: COVID vs Healthy") +
  theme_classic()

# Switch to RNA assay for differential expression (DEGs)----------------------------------------
DefaultAssay(covid_merged) <- "RNA"
covid_merged <- NormalizeData(covid_merged, verbose = FALSE)

# DEGs------------------------------------------------------------------------------------------
all_degs <- list()
cell_types <- unique(covid_merged$cell_type)
cell_types <- cell_types[!cell_types %in% c("Low Quality", "Erythrocytes")]

for(ct in cell_types) {
  tryCatch({
    degs <- FindMarkers(
      covid_merged,
      ident.1 = "COVID",
      ident.2 = "Healthy",
      group.by = "condition",
      subset.ident = ct,
      min.pct = 0.25,
      logfc.threshold = 0.25,
      verbose = FALSE
    )
    degs$cell_type <- ct
    degs$gene <- rownames(degs)
    all_degs[[ct]] <- degs
    cat("Done:", ct, "->", nrow(degs), "DEGs\n")
  }, error = function(e) {
    cat("Skipped:", ct, "- Error:", conditionMessage(e), "\n")
  })
}

# Combine all results into one dataframe
all_degs_df <- do.call(rbind, all_degs)
rownames(all_degs_df) <- NULL

# Filter significant DEGs
sig_degs <- all_degs_df %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5)

cat("Total significant DEGs:", nrow(sig_degs), "\n")

# Save
write.csv(all_degs_df, "all_DEGs_COVID_vs_Healthy.csv", row.names = FALSE)
write.csv(sig_degs, "significant_DEGs.csv", row.names = FALSE)

# Volacno Plots----------------------------------------------------------------------------------------------
volcano_plots <- lapply(unique(all_degs_df$cell_type), function(ct) {
  data <- all_degs_df %>%
    filter(cell_type == ct) %>%
    mutate(
      significance = case_when(
        p_val_adj < 0.05 & avg_log2FC > 0.5  ~ "Up in COVID",
        p_val_adj < 0.05 & avg_log2FC < -0.5 ~ "Up in Healthy",
        TRUE ~ "Not significant"
      ),
      label = ifelse(abs(avg_log2FC) > 3 & p_val_adj < 0.05, gene, "")
    )
  
ggplot(data, aes(x = avg_log2FC, y = -log10(p_val_adj + 1e-300),
                   color = significance)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 15) +
    scale_color_manual(values = c("Up in COVID"     = "red",
                                  "Up in Healthy"   = "blue",
                                  "Not significant" = "grey70")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    labs(title = ct,
         x = "log2 Fold Change",
         y = "-log10 adj p-value",
         color = "") +
    theme_classic() +
    theme(plot.title = element_text(size = 9, face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(size = 7))
})

# Save the plot
pdf("volcano_all_celltypes.pdf", width = 20, height = 25)
wrap_plots(volcano_plots, ncol = 4)
dev.off()

# Heatmap------------------------------------------------------------------------------------------
# Top 3 genes per cell type
top_genes <- sig_degs %>%
  filter(avg_log2FC > 0) %>%
  group_by(cell_type) %>%
  top_n(3, avg_log2FC) %>%
  pull(gene) %>%
  unique()

covid_merged$cell_type_condition <- paste(covid_merged$cell_type, 
                                          covid_merged$condition, 
                                          sep = "_")
Idents(covid_merged) <- "cell_type_condition"

# Average expression matrix
avg_exp <- AverageExpression(covid_merged, 
                             features = top_genes,
                             assays = "RNA",
                             return.seurat = FALSE)$RNA

# Plot heatmap-----------------------------------------------------------------------------------
p <- pheatmap(avg_exp,
              scale = "row",
              color = colorRampPalette(c("blue", "white", "red"))(100),
              angle_col = 45,
              fontsize_row = 7,
              fontsize_col = 7,
              main = "Top DEGs: COVID vs Healthy across Cell Types")

png("pheatmap.png", width = 15, height = 15, units = "in", res = 600)
grid::grid.draw(p$gtable)
dev.off()


