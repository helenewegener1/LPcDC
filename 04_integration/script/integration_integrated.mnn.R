setwd("~/Documents/projects/project_cDC/LPcDC/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(glmGamPoi)
library(dplyr)
library(glue)
library(ggplot2)
# remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
# remotes::install_github('satijalab/azimuth', ref = 'master')
library(Azimuth)

options(future.globals.maxSize = 1e9)

seurat_merged <- readRDS("04_integration/out/seurat_merged.rds")

########################## Integration using FastMNNIntegration ###########################

# Split SCT layers
Layers(seurat_merged[["SCT"]])
seurat_merged_split <- seurat_merged
seurat_merged_split[["SCT"]] <- split(seurat_merged_split[["SCT"]], f = seurat_merged_split$orig.ident)
Layers(seurat_merged_split[["SCT"]])

seurat_mnn <- IntegrateLayers(
  # object = seurat_merged, # does not work 
  object = seurat_merged_split, # works 
  method = FastMNNIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.mnn", 
  assay = "SCT", 
  verbose = FALSE
)

reduction <- "integrated.mnn"

seurat_integrated <- JoinLayers(seurat_mnn)

####################### Run UMAP using CCA Integration embedding #######################

seurat_integrated <- FindNeighbors(seurat_integrated, reduction = reduction, dims = 1:20)
seurat_integrated <- RunUMAP(seurat_integrated, reduction = reduction, dims = 1:20)

# Visualize with UMAP stratified by dataset - post harmony integration 
DimPlot(seurat_integrated, reduction = "umap", group.by = "orig.ident") +
  labs(title = glue("UMAP - post {reduction} integration")) + 
  theme(legend.text = element_text(size = 8))

ggsave(glue("04_integration/plot/UMAP_post_{reduction}_integration_orig.ident.pdf"), 
       width = 8, 
       height = 7)

DimPlot(seurat_integrated, reduction = "umap", split.by = "orig.ident", ncol = 3)
ggsave(glue("04_integration/plot/UMAP_post_{reduction}_integration_orig.ident_split.pdf"), 
       width = 12, 
       height = 8)

# Visualize with UMAP stratified by seurat clusters - post harmony integration 

res_list <- c(0.3, 0.5, 0.7)

for (res in res_list){
  
  # res <- 0.3
  
  seurat_integrated <- FindClusters(seurat_integrated, resolution = res)
  
  DimPlot(seurat_integrated, reduction = "umap", group.by = glue("SCT_snn_res.{res}"), label = TRUE) +
    labs(title = glue("UMAP - post {reduction} integration"),
         subtitle = glue("SCT_snn_res.{res}"))
  
  ggsave(glue(glue("04_integration/plot/UMAP_post_{reduction}_integration_SCT_snn_res_{res}.pdf")), 
         width = 8, 
         height = 7)
  
  DimPlot(seurat_integrated, reduction = "umap", group.by = glue("SCT_snn_res.{res}"), split.by = "orig.ident", ncol = 3) +
    labs(title = glue("UMAP - post {reduction} integration"),
         subtitle = glue("SCT_snn_res.{res}"))
  
  ggsave(glue("04_integration/plot/UMAP_post_{reduction}_integration_SCT_snn_res_{res}_split.by_orig.ident.pdf"), 
         width = 12, 
         height = 8)
  
  
}

# Feature plots: UMAP stratified by continuous variable 

features <- c("nFeature_RNA", "percent.mt", "percent.ribo")

lapply(features, function(x) {
  
  FeaturePlot(seurat_integrated, reduction = "umap", features = x)
  ggsave(glue("04_integration/plot/UMAP_post_{reduction}_integration_{x}.pdf"), 
         width = 8, 
         height = 7)
  
})

################### Export list of integrated Seurat objects ################### 

saveRDS(seurat_integrated, glue("04_integration/out/seurat_{reduction}_integrated.rds"))















