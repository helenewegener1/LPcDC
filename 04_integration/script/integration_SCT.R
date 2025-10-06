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

# Check that layers are split
Layers(seurat_merged[["SCT"]])

# https://github.com/satijalab/seurat/issues/7542
seurat_merged <- SCTransform(seurat_merged, vst.flavor = "v2")
DefaultAssay(seurat_merged) <- "SCT"
seurat_merged <- RunPCA(seurat_merged, npcs = 30, verbose = FALSE)

# This is just for plotting UMAP
ElbowPlot(seurat_merged)
seurat_merged <- FindNeighbors(seurat_merged,  dims = 1:30)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.4)
seurat_merged <- RunUMAP(seurat_merged, reduction = "pca", dims = 1:30)

# Visualize with UMAP stratified by dataset - pre integration
DimPlot(seurat_merged, reduction = "umap", group.by = "orig.ident") +
  labs(title = "UMAP - SCT - pre integration") +
  theme(legend.text = element_text(size = 8))

ggsave("04_integration/plot/SCT/UMAP_SCT_pre_integration_orig.ident.pdf",
       width = 8,
       height = 7)


############################## Integration on SCT ##############################

# SCT assay: Uses sctransform for variance stabilization, which better handles technical noise and sequencing depth differences

# Initiate integrated seurat obejct
seurat_integrated <- seurat_merged

DefaultAssay(seurat_integrated) 

Layers(seurat_integrated[["SCT"]])

Reductions(seurat_integrated)

# Integrate using harmony
seurat_integrated <- IntegrateLayers(
  object = seurat_integrated, 
  method = HarmonyIntegration, 
  orig.reduction = "pca", 
  new.reduction = "SCT_integrated.harmony", 
  normalization.method = "SCT",
  assay = "SCT",
  group.by.vars = "orig.ident",
  verbose = FALSE
)

# Does not work well for SCT layer
seurat_integrated <- IntegrateLayers(
  object = seurat_integrated,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "SCT_integrated.cca",
  normalization.method = "SCT",
  # group.by.vars = "orig.ident",
  assay = "SCT",
  verbose = TRUE
)

options(future.globals.maxSize = 20 * 1024^3)  # 8 GB
seurat_integrated <- IntegrateLayers(
  object = seurat_integrated,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "SCT_integrated.rpca",
  normalization.method = "SCT",
  assay = "SCT",
  verbose = FALSE
)

# Integrate using FastMNNIntegration
# seurat_integrated <- IntegrateLayers(
#   object = seurat_integrated,
#   method = FastMNNIntegration,
#   orig.reduction = "pca",
#   new.reduction = "SCT_integrated.mnn",
#   normalization.method = "SCT",
#   assay = "SCT",
#   verbose = FALSE
# )


####################### Run UMAP using Harmony embedding #######################

Reductions(seurat_integrated)

reductions <- list(
  
  c("SCT_integrated.harmony", "SCT_umap.harmony", "SCT_harmony_clusters"),
  # c("SCT_integrated.mnn", "SCT_umap.mnn", "SCT_mnn_clusters"),
  c("SCT_integrated.rpca", "SCT_umap.rpca", "SCT_rpca_clusters"),
  c("SCT_integrated.cca", "SCT_umap.cca", "SCT_cca_clusters")
  
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
  
  ggsave(glue("04_integration/plot/{assay}/UMAP_{reduction}_orig.ident.pdf"), 
         width = 8, 
         height = 7)
  
  DimPlot(seurat_integrated, reduction = umap_reduction.name, split.by = "orig.ident", ncol = 3) +
    labs(title = glue("UMAP - post {reduction}")) + 
    theme(legend.text = element_text(size = 8))
  
  ggsave(glue("04_integration/plot/{assay}/UMAP_{reduction}_orig.ident_split.pdf"), 
         width = 12, 
         height = 8)
  
  # Visualize with UMAP stratified by seurat clusters - post harmony integration 
  
  res_list <- c(0.1, 0.3, 0.5)
  
  for (res in res_list){
    
    # res <- 0.3
    
    seurat_integrated <- FindClusters(seurat_integrated, resolution = res)
    
    DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = glue("{assay}_snn_res.{res}"), label = TRUE) +
      labs(title = glue("UMAP - post {reduction}"),
           subtitle = glue("{assay}_snn_res.{res}"))
    
    ggsave(glue(glue("04_integration/plot/{assay}/UMAP_{reduction}_{assay}_snn_res_{res}.pdf")), 
           width = 8, 
           height = 7)
    
    DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = glue("{assay}_snn_res.{res}"), split.by = "orig.ident", ncol = 3) +
      labs(title = glue("UMAP - post {reduction}"),
           subtitle = glue("{assay}_snn_res.{res}"))
    
    ggsave(glue("04_integration/plot/{assay}/UMAP_{reduction}_{assay}_snn_res_{res}_split.by_orig.ident.pdf"), 
           width = 12, 
           height = 8)
    
    
  }
  
  # Feature plots: UMAP stratified by continuous variable 
  # 
  # features <- c("nFeature_RNA", "percent.mt", "percent.ribo")
  # 
  # lapply(features, function(x) {
  #   
  #   FeaturePlot(seurat_integrated, reduction = umap_reduction.name, features = x)
  #   ggsave(glue("04_integration/plot/{assay}/UMAP_{reduction}_{x}.pdf"), 
  #          width = 8, 
  #          height = 7)
  #   
  # })
  
  
}

######################### Save as h5ad file for python ######################### 

Reductions(seurat_integrated)

saveRDS(seurat_integrated, "04_integration/out/seurat_integrated_v5_SCT.rds")
# seurat_integrated <- readRDS("04_integration/out/seurat_integrated_v5_SCT.rds")


######################### Save as h5ad file for python ######################### 

library(SeuratDisk)
library(rhdf5)

# All in one file 
obj_tmp <- seurat_integrated 

DefaultAssay(obj_tmp)

Layers(obj_tmp[["SCT"]])

# IMPORTANT: Join layers before export (h5ad doesn't support split layers)
# obj_tmp[["SCT"]] <- JoinLayers(obj_tmp[["SCT"]]) # already joined

obj_tmp[["SCT3"]] <- as(object = obj_tmp[["SCT"]], Class = "Assay")
DefaultAssay(obj_tmp) <- "SCT3"
obj_tmp[["SCT"]] <- NULL
obj_tmp <- RenameAssays(object = obj_tmp, SCT3 = 'SCT')

# Check new names
reductions <- Reductions(obj_tmp)[str_detect(Reductions(obj_tmp), "integrated")]

# 1. Extract the embeddings
embeddings <- lapply(reductions, function(x) Embeddings(obj_tmp, reduction = x))
names(embeddings) <- reductions

# 2. Add the clean matrix as a new reduction
for(x in reductions) {
  
  new_name <- str_replace(x, "\\.", "_")
  
  new_key <- str_replace(x, "_integrated\\.", "") %>% 
    str_to_upper() %>% 
    paste0("_")
  
  obj_tmp[[new_name]] <- CreateDimReducObject(
    embeddings = embeddings[[x]], 
    key = new_key,
    assay = DefaultAssay(obj_tmp),
    global = TRUE
  )
  
  # Remove the problematic reduction to keep things clean
  obj_tmp[[x]] <- NULL
  
}

# Handle PCA 
obj_tmp[["PCA"]] <- CreateDimReducObject(
  embeddings = Embeddings(obj_tmp, reduction = "pca"), 
  key = "PCA_",
  assay = DefaultAssay(obj_tmp),
  global = TRUE
)

# Remove the problematic reduction to keep things clean
obj_tmp[["pca"]] <- NULL

umap_reductions <- Reductions(obj_tmp)[str_detect(Reductions(obj_tmp), "umap")]
for (x in umap_reductions) {
  obj_tmp[[x]] <- NULL
}

# Verify the new reduction name (optional)
Reductions(obj_tmp)

filename <- glue("04_integration/out/mydata_SCT_v5.h5Seurat")
SaveH5Seurat(obj_tmp, filename = filename, overwrite = TRUE)
Convert(filename, dest = "h5ad", overwrite = TRUE, verbose = FALSE)














