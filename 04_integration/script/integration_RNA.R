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

# options(future.globals.maxSize = 1e9)

# seurat_merged <- readRDS("04_integration/out/seurat_merged_SCT.rds")
# 
# ###################### Merge into one Seurat object - RNA ######################
# 
# DefaultAssay(seurat_merged) <- "RNA"
# seurat_merged <- NormalizeData(seurat_merged)
# seurat_merged <- FindVariableFeatures(seurat_merged)
# seurat_merged <- ScaleData(seurat_merged)
# seurat_merged <- RunPCA(seurat_merged)
# DefaultAssay(seurat_merged)
# 
# saveRDS(seurat_merged, "04_integration/out/seurat_merged_RNA.rds")
seurat_merged <- readRDS("04_integration/out/seurat_merged_RNA.rds")

seurat_merged <- FindNeighbors(seurat_merged, assay = "RNA",  dims = 1:30)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.4)
seurat_merged <- RunUMAP(seurat_merged, assay = "RNA", reduction = "pca", dims = 1:30)

# Visualize with UMAP stratified by dataset - pre integration 
DimPlot(seurat_merged, reduction = "umap", group.by = "orig.ident") + 
  labs(title = "UMAP - RNA - pre integration") + 
  theme(legend.text = element_text(size = 8))

ggsave("04_integration/plot/RNA/UMAP_RNA_pre_integration_orig.ident.pdf", 
       width = 8, 
       height = 7)

####################### Integration using CCAIntegration #######################

seurat_integrated <- IntegrateLayers(
  object = seurat_merged, 
  method = CCAIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.cca",
  assay = "RNA", 
  verbose = FALSE
)

seurat_integrated <- IntegrateLayers(
  object = seurat_integrated, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.harmony",
  assay = "RNA", 
  verbose = FALSE
)

# seurat_integrated <- IntegrateLayers(
#   object = seurat_integrated,
#   method = RPCAIntegration,
#   orig.reduction = "pca",
#   new.reduction = "integrated.rpca",
#   assay = "RNA",
#   verbose = FALSE
# )

seurat_integrated <- IntegrateLayers(
  object = seurat_integrated,
  method = FastMNNIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.mnn",
  assay = "RNA",
  verbose = FALSE
)

# seurat_integrated <- IntegrateLayers(
#   object = seurat_integrated,
#   method = scVIIntegration,
#   new.reduction = "integrated.scvi",
#   conda_env = "../miniconda3/envs/scvi-env",
#   assay = "RNA",
#   verbose = FALSE
# )

seurat_integrated@reductions %>% names()

################### Export list of integrated Seurat objects ################### 

# saveRDS(seurat_integrated, "04_integration/out/seurat_integrated_RNA.rds")
seurat_integrated <- readRDS("04_integration/out/seurat_integrated_RNA.rds")

####################### Run UMAP using CCA Integration embedding #######################

reductions <- list(
  c("integrated.cca", "umap.caa", "cca_clusters"), 
  c("integrated.harmony", "umap.harmony", "harmony_clusters"), 
  c("integrated.mnn", "umap.mnn", "mnn_clusters")
)


for (red in reductions){
  
  reduction <- red[[1]]
  umap_reduction.name <- red[[2]]
  cluster.name <- red[[3]]

  seurat_integrated <- FindNeighbors(seurat_integrated, reduction = reduction, dims = 1:30)
  seurat_integrated <- FindClusters(seurat_integrated, resolution = 2, cluster.name = cluster.name)
  seurat_integrated <- RunUMAP(seurat_integrated, reduction = reduction, dims = 1:20, reduction.name = umap_reduction.name)
  seurat_integrated <- FindNeighbors(seurat_integrated, reduction = reduction, dims = 1:20)
  
  # Visualize with UMAP stratified by dataset - post harmony integration 
  DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = "orig.ident") +
    labs(title = glue("UMAP - RNA - post {reduction} integration")) + 
    theme(legend.text = element_text(size = 8))
  
  ggsave(glue("04_integration/plot/RNA/UMAP_RNA_post_{reduction}_integration_orig.ident.pdf"), 
         width = 8, 
         height = 7)
  
  DimPlot(seurat_integrated, reduction = umap_reduction.name, split.by = "orig.ident", ncol = 3) +
    labs(title = glue("UMAP - post RNA {reduction} integration")) + 
    theme(legend.text = element_text(size = 8))
  
  ggsave(glue("04_integration/plot/RNA/UMAP_RNA_post_{reduction}_integration_orig.ident_split.pdf"), 
         width = 12, 
         height = 8)
  
  # Visualize with UMAP stratified by seurat clusters - post harmony integration 
  
  res_list <- c(0.3, 0.5, 0.7)
  
  for (res in res_list){
    
    # res <- 0.3
    
    seurat_integrated <- FindClusters(seurat_integrated, resolution = res)
    
    DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = glue("RNA_snn_res.{res}"), label = TRUE) +
      labs(title = glue("UMAP - post RNA {reduction} integration"),
           subtitle = glue("RNA_snn_res.{res}"))
    
    ggsave(glue(glue("04_integration/plot/RNA/UMAP_RNA_post_{reduction}_integration_RNA_snn_res_{res}.pdf")), 
           width = 8, 
           height = 7)
    
    DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = glue("RNA_snn_res.{res}"), split.by = "orig.ident", ncol = 3) +
      labs(title = glue("UMAP - post RNA {reduction} integration"),
           subtitle = glue("RNA_snn_res.{res}"))
    
    ggsave(glue("04_integration/plot/RNA/UMAP_RNA_post_{reduction}_integration_RNA_snn_res_{res}_split.by_orig.ident.pdf"), 
           width = 12, 
           height = 8)
    
    
  }
  
  # Feature plots: UMAP stratified by continuous variable 
  
  features <- c("nFeature_RNA", "percent.mt", "percent.ribo")
  
  lapply(features, function(x) {
    
    FeaturePlot(seurat_integrated, reduction = umap_reduction.name, features = x)
    ggsave(glue("04_integration/plot/RNA/UMAP_RNA_post_{reduction}_integration_{x}.pdf"), 
           width = 8, 
           height = 7)
    
  })
  
  
}

######################### Save as h5ad file for python ######################### 

library(SeuratDisk)
library(rhdf5)

# list of reductions to benchmark
reductions <- c("integrated.cca", "integrated.mnn", "integrated.harmony")

for (red in reductions) {
  
  # red <- reductions[[1]]
  
  message("Saving: ", red)
  
  # create a copy of the object with this reduction as default
  obj_tmp <- seurat_integrated
  DefaultAssay(obj_tmp) <- "RNA"
  obj_tmp@reductions <- obj_tmp@reductions[red]  # keep only relevant reduction
  
  obj_tmp@meta.data[] <- lapply(obj_tmp@meta.data, function(x) {
    if (is.factor(x)) as.character(x) else x
  })
  
  # Caviet to save for the Convert 
  obj_tmp[["RNA3"]] <- as(object = obj_tmp[["RNA"]], Class = "Assay")
  DefaultAssay(obj_tmp) <- "RNA3"
  obj_tmp[["RNA"]] <- NULL
  obj_tmp <- RenameAssays(object = obj_tmp, RNA3 = 'RNA')
  
  filename <- glue("04_integration/out/mydata_RNA_{red}.h5Seurat")
  SaveH5Seurat(obj_tmp, filename = filename, overwrite = TRUE)
  Convert(filename, dest = "h5ad", overwrite = TRUE, verbose = FALSE)
  
}


# All in one file 
obj_tmp <- seurat_integrated 
obj_tmp[["RNA3"]] <- as(object = obj_tmp[["RNA"]], Class = "Assay")
DefaultAssay(obj_tmp) <- "RNA3"
obj_tmp[["RNA"]] <- NULL
obj_tmp <- RenameAssays(object = obj_tmp, RNA3 = 'RNA')

filename <- glue("04_integration/out/mydata_RNA.h5Seurat")
SaveH5Seurat(obj_tmp, filename = filename, overwrite = TRUE)
Convert(filename, dest = "h5ad", overwrite = TRUE, verbose = FALSE)













