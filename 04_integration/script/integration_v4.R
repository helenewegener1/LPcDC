setwd("~/Documents/projects/project_cDC/LPcDC/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(glmGamPoi)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)
library(harmony)
# remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)
# remotes::install_github('satijalab/azimuth', ref = 'master')
library(Azimuth)

# Load data
seurat_obj_list <- readRDS("03_QC/out/seurat_obj_finalQC_list.rds")

############################ RNA integration prep ############################

# Merge
seurat_merged <- merge(seurat_obj_list[[1]], y = seurat_obj_list[-1])

# Split SCT layer
seurat_merged[["SCT"]] <- split(seurat_merged[["SCT"]], f = seurat_merged$orig.ident)

Layers(seurat_merged[["SCT"]])

# Control all cells from each dataset are included
seurat_merged@meta.data$orig.ident %>% table(useNA = "always")

# Seurat workflow on merged data
seurat_merged <- SCTransform(seurat_merged, assay = "RNA", layer = "counts", verbose = FALSE, vst.flavor = "v2")
seurat_merged <- RunPCA(seurat_merged, assay = "SCT", verbose = FALSE) # HW: Yes - SCT assay for PCA

# ^^ Needs the part above to run harmony ^^

# This is just for plotting UMAP
ElbowPlot(seurat_merged)
seurat_merged <- FindNeighbors(seurat_merged, assay = "SCT",  dims = 1:20)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.4)
seurat_merged <- RunUMAP(seurat_merged, assay = "SCT", reduction = "pca", dims = 1:20)

# Save merged object
# saveRDS(seurat_merged, "04_integration/out/seurat_merged_SCT.rds")
# seurat_merged <- readRDS("04_integration/out/seurat_merged_SCT.rds")

DefaultAssay(seurat_merged)

# Visualize with UMAP stratified by dataset - pre integration
DimPlot(seurat_merged, reduction = "umap", group.by = "orig.ident") +
  labs(title = "UMAP - SCT - pre integration") +
  theme(legend.text = element_text(size = 8))

# ggsave("04_integration/plot/SCT/UMAP_SCT_pre_integration_orig.ident.pdf",
#        width = 8,
#        height = 7)

############################# SCT integration prep #############################
# 
# DefaultAssay(seurat_merged) <- "RNA"
# seurat_merged <- NormalizeData(seurat_merged, assay = "RNA")
# seurat_merged <- FindVariableFeatures(seurat_merged, assay = "RNA")
# seurat_merged <- ScaleData(seurat_merged, assay = "RNA")
# seurat_merged <- RunPCA(seurat_merged, assay = "RNA")
# 
# ElbowPlot(seurat_merged)
# seurat_merged <- FindNeighbors(seurat_merged, assay = "RNA",  dims = 1:20)
# seurat_merged <- FindClusters(seurat_merged, resolution = 0.4)
# seurat_merged <- RunUMAP(seurat_merged, assay = "RNA", reduction = "pca", dims = 1:20)
# 
# DimPlot(seurat_merged, reduction = "umap", group.by = "orig.ident") +
#   labs(title = "UMAP - RNA - pre integration") +
#   theme(legend.text = element_text(size = 8))
# 
# # ggsave("04_integration/plot/RNA/UMAP_RNA_pre_integration_orig.ident.pdf",
# #        width = 8,
# #        height = 7)
# 
# DefaultAssay(seurat_merged)

############################## Integration on SCT ##############################

# SCT assay: Uses sctransform for variance stabilization, which better handles technical noise and sequencing depth differences

seurat_integrated <- seurat_merged

DefaultAssay(seurat_integrated) <- "SCT"

Layers(seurat_integrated[["SCT"]])

# Integrate using harmony
seurat_integrated <- IntegrateLayers(
  object = seurat_integrated, 
  method = HarmonyIntegration, 
  orig.reduction = "pca", 
  new.reduction = "SCT_integrated.harmony", 
  # normalization.method = "SCT",
  assay = "SCT",
  verbose = FALSE
)

# Run Harmony on the merged object
seurat_integrated <- RunHarmony(
  seurat_integrated,
  group.by.vars = "orig.ident",  # metadata column indicating dataset
  dims = 1:20
)

seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "harmony", dims = 1:20)
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "harmony", dims = 1:20)

# Visualize with UMAP stratified by dataset - post harmony integration 
DimPlot(seurat_integrated, reduction = "umap", group.by = "orig.ident") +
  labs(title = "UMAP - post harmonny integration")
# ggsave("04_integration/plot/UMAP_post_harmony_integration_orig.ident.pdf")


############################## Integration on RNA ##############################

# DefaultAssay(seurat_integrated) <- "RNA"
# 
# Layers(seurat_integrated[["RNA"]])
# 
# seurat_integrated <- IntegrateLayers(
#   object = seurat_integrated, 
#   method = HarmonyIntegration,
#   orig.reduction = "pca", 
#   new.reduction = "RNA_integrated.harmony",
#   assay = "RNA", 
#   verbose = FALSE
# )
# 
# Reductions(seurat_integrated)

####################### Run UMAP using Harmony embedding #######################

reductions <- list(

  # c("RNA_integrated.cca", "RNA_umap.cca", "RNA_cca_clusters"),
  # c("RNA_integrated.harmony", "RNA_umap.harmony", "RNA_harmony_clusters"),
  # c("RNA_integrated.mnn", "RNA_umap.mnn", "RNA_mnn_clusters"),
  # c("RNA_integrated.rpca", "RNA_umap.rpca", "RNA_rpca_clusters"),
  # 
  c("SCT_integrated.harmony", "SCT_umap.harmony", "SCT_harmony_clusters")
  # c("SCT_integrated.mnn", "SCT_umap.mnn", "SCT_mnn_clusters"),
  # c("SCT_integrated.rpca", "SCT_umap.rpca", "SCT_rpca_clusters")
  
)

for (red in reductions){
  
  reduction <- red[[1]]
  umap_reduction.name <- red[[2]]
  cluster.name <- red[[3]]
  
  # Either RNA or SCT 
  assay <- str_split_i(reduction, "_", 1)
  
  # Set default assay 
  DefaultAssay(seurat_integrated) <- assay
  
  seurat_integrated <- FindNeighbors(seurat_integrated, reduction = reduction, dims = 1:20)
  seurat_integrated <- FindClusters(seurat_integrated, resolution = 2, cluster.name = cluster.name)
  seurat_integrated <- RunUMAP(seurat_integrated, reduction = reduction, dims = 1:20, reduction.name = umap_reduction.name)
  seurat_integrated <- FindNeighbors(seurat_integrated, reduction = reduction, dims = 1:20)
  
  # Visualize with UMAP stratified by dataset - post harmony integration 
  DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = "orig.ident") +
    labs(title = glue("UMAP - post {reduction}")) + 
    theme(legend.text = element_text(size = 8))
  
  ggsave(glue("04_integration/plot/v4/{assay}/UMAP_{reduction}_orig.ident.pdf"), 
         width = 8, 
         height = 7)
  
  DimPlot(seurat_integrated, reduction = umap_reduction.name, split.by = "orig.ident", ncol = 3) +
    labs(title = glue("UMAP - post {reduction}")) + 
    theme(legend.text = element_text(size = 8))
  
  ggsave(glue("04_integration/plot/v4/{assay}/UMAP_{reduction}_orig.ident_split.pdf"), 
         width = 12, 
         height = 8)
  
  # Visualize with UMAP stratified by seurat clusters - post harmony integration 
  
  res_list <- c(0.3, 0.5, 0.7)
  
  for (res in res_list){
    
    # res <- 0.3
    
    seurat_integrated <- FindClusters(seurat_integrated, resolution = res)
    
    DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = glue("{assay}_snn_res.{res}"), label = TRUE) +
      labs(title = glue("UMAP - post {reduction}"),
           subtitle = glue("{assay}_snn_res.{res}"))
    
    ggsave(glue(glue("04_integration/plot/v4/{assay}/UMAP_{reduction}_{assay}_snn_res_{res}.pdf")), 
           width = 8, 
           height = 7)
    
    DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = glue("{assay}_snn_res.{res}"), split.by = "orig.ident", ncol = 3) +
      labs(title = glue("UMAP - post {reduction}"),
           subtitle = glue("{assay}_snn_res.{res}"))
    
    ggsave(glue("04_integration/plot/v4/{assay}/UMAP_{reduction}_{assay}_snn_res_{res}_split.by_orig.ident.pdf"), 
           width = 12, 
           height = 8)
    
    
  }
  
  # Feature plots: UMAP stratified by continuous variable 
  
  features <- c("nFeature_RNA", "percent.mt", "percent.ribo")
  
  lapply(features, function(x) {
    
    FeaturePlot(seurat_integrated, reduction = umap_reduction.name, features = x)
    ggsave(glue("04_integration/plot/v4/{assay}/UMAP_{reduction}_{x}.pdf"), 
           width = 8, 
           height = 7)
    
  })
  
  
}

######################### Save as h5ad file for python ######################### 

Reductions(seurat_integrated)

saveRDS(seurat_integrated, "04_integration/out/seurat_integrated_sct_harmony.rds")
seurat_integrated <- readRDS("04_integration/out/seurat_integrated_sct_harmony.rds")


######################### Save as h5ad file for python ######################### 

library(SeuratDisk)
library(rhdf5)

# All in one file 
obj_tmp <- seurat_integrated 

DefaultAssay(obj_tmp)

Layers(obj_tmp[["SCT"]])

# IMPORTANT: Join layers before export (h5ad doesn't support split layers)
# obj_tmp[["SCT"]] <- JoinLayers(obj_tmp[["SCT"]])

obj_tmp[["SCT3"]] <- as(object = obj_tmp[["SCT"]], Class = "Assay")
DefaultAssay(obj_tmp) <- "SCT3"
obj_tmp[["SCT"]] <- NULL
obj_tmp <- RenameAssays(object = obj_tmp, SCT3 = 'SCT')

# # Check new names
# reductions <- Reductions(obj_tmp)[str_detect(Reductions(obj_tmp), "integrated")]
# 
# # 1. Extract the embeddings
# embeddings <- lapply(reductions, function(x) Embeddings(obj_tmp, reduction = x))
# names(embeddings) <- reductions
# 
# # 2. Add the clean matrix as a new reduction
# for(x in reductions) {
#   
#   new_name <- str_replace(x, "\\.", "_")
#   
#   new_key <- str_replace(x, "_integrated\\.", "") %>% 
#     str_to_upper() %>% 
#     paste0("_")
#   
#   obj_tmp[[new_name]] <- CreateDimReducObject(
#     embeddings = embeddings[[x]], 
#     key = new_key,
#     assay = DefaultAssay(obj_tmp),
#     global = TRUE
#   )
#   
#   # Remove the problematic reduction to keep things clean
#   obj_tmp[[x]] <- NULL
#   
# }
# 
# # Handle PCA 
# obj_tmp[["PCA"]] <- CreateDimReducObject(
#   embeddings = Embeddings(obj_tmp, reduction = "pca"), 
#   key = "PCA_",
#   assay = DefaultAssay(obj_tmp),
#   global = TRUE
# )
# 
# # Remove the problematic reduction to keep things clean
# obj_tmp[["pca"]] <- NULL
# 
# umap_reductions <- Reductions(obj_tmp)[str_detect(Reductions(obj_tmp), "umap")]
# for (x in umap_reductions) {
#   obj_tmp[[x]] <- NULL
# }

# Verify the new reduction name (optional)
Reductions(obj_tmp)

filename <- glue("04_integration/out/mydata_v4.h5Seurat")
SaveH5Seurat(obj_tmp, filename = filename, overwrite = TRUE)
Convert(filename, dest = "h5ad", overwrite = TRUE, verbose = FALSE)














