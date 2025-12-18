setwd("~/Documents/projects/project_cDC/LPcDC/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)

# Load data
seurat_object <- readRDS("04_SILP_integration/out/SILP_seurat_integrated_v5_RNA.rds")

# Settings
Reductions(seurat_object)

reduction <- "RNA_umap.cca"
cluster.name <- "RNA_cca_clusters_res.1"

# Final UMAP
DimPlot(seurat_object, reduction = reduction, group.by = cluster.name, label = TRUE) +
  labs(title = glue("UMAP - post integration"),
       subtitle = cluster.name)

# Dataset to remove
seurat_object$orig.ident %>% unique()
sample_name <- "E_MTAB_9522"

# Subset cluster
seurat_object_clean <- subset(seurat_object, subset = orig.ident != sample_name)
seurat_object_clean$orig.ident %>% unique()

# Re-seurat workflow and integrate 
seurat_object_clean <- NormalizeData(seurat_object_clean)
seurat_object_clean <- FindVariableFeatures(seurat_object_clean)
seurat_object_clean <- ScaleData(seurat_object_clean)
# seurat_object_clean <- SCTransform(seurat_object_clean)
seurat_object_clean <- RunPCA(seurat_object_clean)
ElbowPlot(seurat_object_clean)

# Integrate
methods <- c("cca", "harmony")

res <- 1.2
n_dim <- 30

for (method in methods){
  
  method <- "cca"
  
  new.reduction <- glue("RNA_{method}_clean") %>% as.character()
  cluster.name <- glue("RNA_{method}_clean_clusters_res.{res}") %>% as.character()
  umap_reduction <- glue("RNA_umap.{method}_clean") %>% as.character()
  
  if (method == "cca"){
    
    seurat_object_clean <- IntegrateLayers(
      object = seurat_object_clean, 
      method = CCAIntegration,
      orig.reduction = "pca", 
      new.reduction = new.reduction,
      assay = "RNA", 
      verbose = FALSE
    )
    
  } else if (method == "harmony"){
    
    seurat_object_clean <- IntegrateLayers(
      object = seurat_object_clean, 
      method = HarmonyIntegration,
      orig.reduction = "pca", 
      new.reduction = new.reduction,
      assay = "RNA", 
      verbose = FALSE
    )
    
  }
  
  Reductions(seurat_object_clean)
  
  DefaultAssay(seurat_object_clean)
  seurat_object_clean <- FindNeighbors(seurat_object_clean, reduction = new.reduction, dims = 1:n_dim)
  seurat_object_clean <- FindClusters(seurat_object_clean, reduction = new.reduction, resolution = res, cluster.name = cluster.name)
  seurat_object_clean <- RunUMAP(seurat_object_clean, reduction = new.reduction, dims = 1:n_dim, reduction.name = umap_reduction)
  
  # Cluster plot
  DimPlot(seurat_object_clean, reduction = umap_reduction, group.by = cluster.name, label = TRUE) +
    labs(title = glue("UMAP - {new.reduction}"),
         subtitle = cluster.name)
  # ggsave(glue("07_SILP_subclustering/plot/UMAP_clean_{method}_res.{res}.pdf"), width = 8, height = 7)
  
  # Feature plots
  FeaturePlot(seurat_object_clean, reduction = umap_reduction, features = c("Rorc", "Prdm16"), order = TRUE, pt.size = 1) +
    labs(subtitle = glue("UMAP - {new.reduction}"))
  # ggsave(glue("07_SILP_subclustering/plot/UMAP_clean_{method}_Rorc_Prdm16.pdf"), width = 14, height = 7)
  
}
