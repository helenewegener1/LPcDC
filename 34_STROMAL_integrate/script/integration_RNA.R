# setwd("~/Documents/projects/project_cDC/LPcDC/")

getwd()

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
library(clustree)

# Load data
seurat_obj <- readRDS("31_STROMAL_make_seurat/out/STROMAL_object.rds")

Reductions(seurat_obj)
seurat_obj@meta.data %>% head()

############################ Check batch effect ############################

DimPlot(seurat_obj, reduction = "umap", group.by = "Cluster", split.by = "orig.ident") + 
  labs(title = "UMAP - RNA - pre integration") +
  theme(legend.text = element_text(size = 8))

ggsave("34_STROMAL_integrate/plot/UMAP_RNA_pre_integration_annotation.pdf", width = 12, height = 7)

DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident") + 
  labs(title = "UMAP - RNA - pre integration") +
  theme(legend.text = element_text(size = 8))

ggsave("34_STROMAL_integrate/plot/UMAP_RNA_pre_integration_sample.pdf", width = 12, height = 7)

############################ RNA integration prep ############################

seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$orig.ident)

Layers(seurat_obj[["RNA"]]) # already split

seurat_merged <- seurat_obj

DefaultAssay(seurat_merged) <- "RNA"

# Seurat workflow
seurat_merged <- NormalizeData(seurat_merged)
seurat_merged <- FindVariableFeatures(seurat_merged)
seurat_merged <- ScaleData(seurat_merged)
seurat_merged <- RunPCA(seurat_merged)

ElbowPlot(seurat_merged)
seurat_merged <- FindNeighbors(seurat_merged,  dims = 1:20)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.4)
seurat_merged <- RunUMAP(seurat_merged, reduction = "pca", dims = 1:20)

# Save merged object
# saveRDS(seurat_merged, "14_inflammation_integration/out/seurat_merged_SCT.rds")
# seurat_merged <- readRDS("14_inflammation_integration/out/seurat_merged_SCT.rds")

DefaultAssay(seurat_merged)


############################## Integration on RNA ##############################

seurat_integrated <- seurat_merged

DefaultAssay(seurat_integrated) <- "RNA"

Layers(seurat_integrated[["RNA"]])

seurat_integrated <- IntegrateLayers(
  object = seurat_integrated, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "RNA_integrated.harmony",
  assay = "RNA", 
  verbose = FALSE
)

seurat_integrated <- IntegrateLayers(
  object = seurat_integrated, 
  method = CCAIntegration,
  orig.reduction = "pca", 
  new.reduction = "RNA_integrated.cca",
  assay = "RNA", 
  verbose = FALSE
)

################### Export list of integrated Seurat objects ################### 

Reductions(seurat_integrated)

####################### Run UMAP using Harmony embedding #######################

reductions <- list(
  
  c("RNA_integrated.cca", "RNA_umap.cca", "RNA_cca_clusters"),
  c("RNA_integrated.harmony", "RNA_umap.harmony", "RNA_harmony_clusters")
  # c("RNA_integrated.mnn", "RNA_umap.mnn", "RNA_mnn_clusters"),
  # c("RNA_integrated.rpca", "RNA_umap.rpca", "RNA_rpca_clusters")
  
)

# red <- c("RNA_integrated.cca", "RNA_umap.cca", "RNA_cca_clusters")

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
  
  ggsave(glue("34_STROMAL_integrate/plot/{assay}/UMAP_{reduction}_sample.pdf"), 
         width = 8, 
         height = 7)
  
  DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = "Cluster") +
    labs(title = glue("UMAP - post {reduction}")) + 
    theme(legend.text = element_text(size = 8))
  
  ggsave(glue("34_STROMAL_integrate/plot/{assay}/UMAP_{reduction}_Cluster.pdf"), 
         width = 8, 
         height = 7)
  
  # DimPlot(seurat_integrated, reduction = umap_reduction.name, split.by = "orig.ident", group.by = "Cluster") +
  #   labs(title = glue("UMAP - post {reduction}")) + 
  #   theme(legend.text = element_text(size = 8))
  # 
  # ggsave(glue("34_STROMAL_integrate/plot/{assay}/UMAP_{reduction}_sample.pdf"), 
  #        width = 12, 
  #        height = 8)
  
  # Visualize with UMAP stratified by seurat clusters - post harmony integration 
  
  res_list <- seq(0.1, 1.5, by = 0.1)
  
  for (res in res_list){
    
    # res <- 0.3
    
    seurat_integrated <- FindClusters(seurat_integrated, resolution = res, cluster.name = glue("{cluster.name}_res.{res}"))
    
    # WHAT TO NAME THEM?
    DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = glue("{cluster.name}_res.{res}"), label = TRUE) +
      labs(title = glue("UMAP - post {reduction}"),
           subtitle = glue("{cluster.name}_res.{res}"))
    
    ggsave(glue(glue("34_STROMAL_integrate/plot/{assay}/UMAP_{reduction}_{assay}_snn_res_{res}.pdf")), 
           width = 8, 
           height = 7)
    
    # DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = glue("{cluster.name}_res.{res}"), split.by = "orig.ident", ncol = 3) +
    #   labs(title = glue("UMAP - post {reduction}"),
    #        subtitle = glue("{cluster.name}_res.{res}"))
    # 
    # ggsave(glue("34_STROMAL_integrate/plot/{assay}/UMAP_{reduction}_{assay}_snn_res_{res}_split.by_sample.pdf"), 
    #        width = 12, 
    #        height = 8)
    
    
  }

  
  
  # Feature plots: UMAP stratified by continuous variable 
  
  # features <- c("nFeature_RNA", "percent.mt", "percent.ribo")
  # 
  # lapply(features, function(x) {
  #   
  #   FeaturePlot(seurat_integrated, reduction = umap_reduction.name, features = x)
  #   ggsave(glue("14_inflammation_integration/plot/{assay}/UMAP_{reduction}_{x}.pdf"), 
  #          width = 8, 
  #          height = 7)
  #   
  # })
  
  
}

# clustree

cluster.name <- "RNA_cca_clusters"
pdf(file = glue("34_STROMAL_integrate/plot/{assay}/clustree_{cluster.name}.pdf"), width = 12, height = 10)
clustree(seurat_integrated, assay = "RNA", return = "plot", prefix = glue("{cluster.name}_res."))
dev.off()

cluster.name <- "RNA_harmony_clusters"
pdf(file = glue("34_STROMAL_integrate/plot/{assay}/clustree_{cluster.name}.pdf"), width = 12, height = 10)
clustree(seurat_integrated, assay = "RNA", return = "plot", prefix = glue("{cluster.name}_res."))
dev.off()

# cluster.name <- "RNA_rpca_clusters"
# pdf(file = glue("14_inflammation_integration/plot/{assay}/clustree_{cluster.name}.pdf"), width = 12, height = 10)
# clustree(seurat_integrated, assay = "RNA", return = "plot", prefix = glue("{cluster.name}_res."))
# dev.off()
# 
# cluster.name <- "RNA_mnn_clusters"
# pdf(file = glue("14_inflammation_integration/plot/{assay}/clustree_{cluster.name}.pdf"), width = 12, height = 10)
# clustree(seurat_integrated, assay = "RNA", return = "plot", prefix = glue("{cluster.name}_res."))
# dev.off()

######################### Save as h5ad file for python ######################### 

Reductions(seurat_integrated)

saveRDS(seurat_integrated, "34_STROMAL_integrate/out/STROMAL_seurat_integrated_v5_RNA.rds")
seurat_integrated <- readRDS("34_STROMAL_integrate/out/STROMAL_seurat_integrated_v5_RNA.rds")

######################## Save as h5ad file for python #########################

# library(SeuratDisk)
# library(rhdf5)
# 
# # All in one file 
# obj_tmp <- seurat_integrated 
# 
# DefaultAssay(obj_tmp)
# 
# Layers(obj_tmp[["RNA"]])
# 
# # IMPORTANT: Join layers before export (h5ad doesn't support split layers)
# obj_tmp[["RNA"]] <- JoinLayers(obj_tmp[["RNA"]])
# 
# obj_tmp[["RNA3"]] <- as(object = obj_tmp[["RNA"]], Class = "Assay")
# DefaultAssay(obj_tmp) <- "RNA3"
# obj_tmp[["RNA"]] <- NULL
# obj_tmp <- RenameAssays(object = obj_tmp, RNA3 = 'RNA')
# 
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
# 
# # Verify the new reduction name (optional)
# Reductions(obj_tmp)
# 
# filename <- glue("34_STROMAL_integrate/out/mydata_RNA_v5.h5Seurat")
# SaveH5Seurat(obj_tmp, filename = filename, overwrite = TRUE)
# Convert(filename, dest = "h5ad", overwrite = TRUE, verbose = FALSE)














