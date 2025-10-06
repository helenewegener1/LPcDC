setwd("~/Documents/projects/project_cDC/LPcDC/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)

# Load data
seurat_integrated <- saveRDS("04_integration/out/seurat_integrated_v5_RNA.rds")

# Settings
DefaultAssay(seurat_integrated) <- "RNA"

Reductions(seurat_integrated)

reduction <- "RNA_umap.cca"
cluster.name <- "RNA_cca_clusters_res.0.4"

# Final UMAP
DimPlot(seurat_integrated, reduction = reduction, group.by = cluster.name, label = TRUE) +
  labs(title = glue("UMAP - post integration"),
       subtitle = cluster.name)

# Cluster to subset
target_cluster <- c(1, 5, 10)

# Subset cluster
# subset_seurat <- subset(seurat_integrated, subset = RNA_cca_clusters_res.0.4 == target_cluster)
subset_seurat <- subset(seurat_integrated, subset = RNA_cca_clusters_res.0.4 %in% target_cluster)

# Seurat workflow 
subset_seurat <- NormalizeData(subset_seurat)
subset_seurat <- FindVariableFeatures(subset_seurat, selection.method = "vst", nfeatures = 2000)	
subset_seurat <- ScaleData(subset_seurat)
subset_seurat <- RunPCA(subset_seurat, features = VariableFeatures(object = subset_seurat))
ElbowPlot(subset_seurat)
subset_seurat <- FindNeighbors(subset_seurat, dims = 1:20)
subset_seurat <- FindClusters(subset_seurat, resolution = 0.2)

subset_seurat <- RunUMAP(subset_seurat, dims = 1:15)
DimPlot(subset_seurat, reduction = "umap", label = TRUE)	

# Find markers 
subset_markers <- FindAllMarkers(subset_seurat, only.pos = TRUE)

