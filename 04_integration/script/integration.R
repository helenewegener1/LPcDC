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

# # Load data
# seurat_obj_finalQC_list <- readRDS("03_QC/out/seurat_obj_finalQC_list.rds")
# 
# ############### SCTransform individual datasets after filtering ################
# 
# seurat_obj_list <- list()
# 
# for (sample_name in names(seurat_obj_finalQC_list)) {
# 
#   # Define seurat object of sample
#   seurat_obj_finalQC <- seurat_obj_finalQC_list[[sample_name]]
# 
#   # SCTransform gives a warning because the SCTransform assay from doublets were removed is overwritten
#   seurat_obj_tmp <- SCTransform(seurat_obj_finalQC, assay = "RNA", layer = "counts", verbose = FALSE)
# 
#   # Save in list
#   seurat_obj_list[[sample_name]] <- seurat_obj_tmp
# 
#   # Clean up
#   rm(seurat_obj_tmp)
# 
# }
# 
# ############################ RNA integration prep ############################
# 
# # Merge
# seurat_merged <- merge(seurat_obj_list[[1]], y = seurat_obj_list[-1])
# 
# # Control all cells from each dataset are included
# seurat_merged@meta.data$orig.ident %>% table(useNA = "always")
# 
# # Remove redundant layers
# Layers(seurat_merged[["RNA"]])
# 
# seurat_merged[["RNA"]]$`counts.Gene Expression.GSM9122899` <- NULL
# seurat_merged[["RNA"]]$`counts.Antibody Capture.GSM9122899` <- NULL
# 
# Layers(seurat_merged[["RNA"]])
# 
# # Seurat workflow on merged data
# seurat_merged <- SCTransform(seurat_merged, assay = "RNA", layer = "counts", verbose = FALSE)
# seurat_merged <- RunPCA(seurat_merged, assay = "SCT", verbose = FALSE) # HW: Yes - SCT assay for PCA
# 
# # ^^ Needs the part above to run harmony ^^
# 
# # This is just for plotting UMAP
# ElbowPlot(seurat_merged)
# seurat_merged <- FindNeighbors(seurat_merged, assay = "SCT",  dims = 1:30)
# seurat_merged <- FindClusters(seurat_merged, resolution = 0.4)
# seurat_merged <- RunUMAP(seurat_merged, assay = "SCT", reduction = "pca", dims = 1:30)
# 
# # Save merged object
# saveRDS(seurat_merged, "04_integration/out/seurat_merged_SCT.rds")
seurat_merged <- readRDS("04_integration/out/seurat_merged_SCT.rds")

DefaultAssay(seurat_merged)

# Visualize with UMAP stratified by dataset - pre integration 
DimPlot(seurat_merged, reduction = "umap", group.by = "orig.ident") + 
  labs(title = "UMAP - pre integration") + 
  theme(legend.text = element_text(size = 8))

ggsave("04_integration/plot/SCT/UMAP_SCT_pre_integration_orig.ident.pdf", 
       width = 8, 
       height = 7)

############################# SCT integration prep #############################

DefaultAssay(seurat_merged) <- "RNA"
seurat_merged <- NormalizeData(seurat_merged)
seurat_merged <- FindVariableFeatures(seurat_merged)
seurat_merged <- ScaleData(seurat_merged)
seurat_merged <- RunPCA(seurat_merged)
DefaultAssay(seurat_merged)

############################## Integration on SCT ##############################

# Initiate integrated seurat obejct
seurat_integrated <- seurat_merged

DefaultAssay(seurat_integrated) <- "SCT"

# Integrate using harmony
seurat_integrated <- IntegrateLayers(
  object = seurat_integrated, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "SCT_integrated.harmony",
  assay = "SCT", 
  verbose = FALSE
)

# Split SCT layers for FastMNN Integration
seurat_integrated[["SCT"]] <- split(seurat_integrated[["SCT"]], f = seurat_integrated$orig.ident)

# Integrate using FastMNNIntegration
seurat_integrated <- IntegrateLayers(
  object = seurat_integrated,
  method = FastMNNIntegration,
  orig.reduction = "pca",
  new.reduction = "SCT_integrated.mnn",
  assay = "SCT",
  verbose = FALSE
)

# Merge back the SCT Layers
seurat_integrated[["SCT"]] <- JoinLayers(seurat_integrated[["SCT"]])

############################## Integration on RNA ##############################

DefaultAssay(seurat_integrated) <- "RNA"

seurat_integrated <- IntegrateLayers(
  object = seurat_integrated, 
  method = CCAIntegration,
  orig.reduction = "pca", 
  new.reduction = "RNA_integrated.cca",
  assay = "RNA", 
  verbose = FALSE
)

seurat_integrated <- IntegrateLayers(
  object = seurat_integrated, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "RNA_integrated.harmony",
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
  new.reduction = "RNA_integrated.mnn",
  assay = "RNA",
  verbose = FALSE
)


################### Export list of integrated Seurat objects ################### 

seurat_integrated@reductions %>% 

# seurat_integrated[["SCT"]]

# saveRDS(seurat_integrated, "04_integration/out/seurat_integrated_SCT.rds")
# seurat_integrated <- readRDS("04_integration/out/seurat_integrated_SCT.rds")

saveRDS(seurat_integrated, "04_integration/out/seurat_integrated_all.rds")

####################### Run UMAP using Harmony embedding #######################

reductions <- list(
  c("RNA_integrated.cca", "RNA_umap.cca", "RNA_cca_clusters"),
  c("RNA_integrated.harmony", "RNA_umap.harmony", "RNA_harmony_clusters"),
  c("RNA_integrated.mnn", "RNA_umap.mnn", "RNA_mnn_clusters"),
  
  c("SCT_integrated.harmony", "SCT_umap.harmony", "SCT_harmony_clusters"),
  c("SCT_integrated.mnn", "SCT_umap.mnn", "SCT_mnn_clusters") 
)

for (red in reductions){
  
  reduction <- red[[1]]
  umap_reduction.name <- red[[2]]
  cluster.name <- red[[3]]
  
  # Either RNA or SCT 
  assay <- str_split_i(reduction, "_", 1)
  
  # Set default assay 
  DefaultAssay(seurat_integrated) <- assay
  
  seurat_integrated <- FindNeighbors(seurat_integrated, reduction = reduction, dims = 1:30)
  seurat_integrated <- FindClusters(seurat_integrated, resolution = 2, cluster.name = cluster.name)
  seurat_integrated <- RunUMAP(seurat_integrated, reduction = reduction, dims = 1:20, reduction.name = umap_reduction.name)
  seurat_integrated <- FindNeighbors(seurat_integrated, reduction = reduction, dims = 1:20)
  
  # Visualize with UMAP stratified by dataset - post harmony integration 
  DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = "orig.ident") +
    labs(title = glue("UMAP - {assay} - post {reduction} integration")) + 
    theme(legend.text = element_text(size = 8))
  
  ggsave(glue("04_integration/plot/{assay}/UMAP_{assay}_post_{reduction}_integration_orig.ident.pdf"), 
         width = 8, 
         height = 7)
  
  DimPlot(seurat_integrated, reduction = umap_reduction.name, split.by = "orig.ident", ncol = 3) +
    labs(title = glue("UMAP - post {assay} {reduction} integration")) + 
    theme(legend.text = element_text(size = 8))
  
  ggsave(glue("04_integration/plot/{assay}/UMAP_{assay}_post_{reduction}_integration_orig.ident_split.pdf"), 
         width = 12, 
         height = 8)
  
  # Visualize with UMAP stratified by seurat clusters - post harmony integration 
  
  res_list <- c(0.3, 0.5, 0.7)
  
  for (res in res_list){
    
    # res <- 0.3
    
    seurat_integrated <- FindClusters(seurat_integrated, resolution = res)
    
    DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = glue("{assay}_snn_res.{res}"), label = TRUE) +
      labs(title = glue("UMAP - post {assay} {reduction} integration"),
           subtitle = glue("{assay}_snn_res.{res}"))
    
    ggsave(glue(glue("04_integration/plot/{assay}/UMAP_{assay}_post_{reduction}_integration_{assay}_snn_res_{res}.pdf")), 
           width = 8, 
           height = 7)
    
    DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = glue("{assay}_snn_res.{res}"), split.by = "orig.ident", ncol = 3) +
      labs(title = glue("UMAP - post {assay} {reduction} integration"),
           subtitle = glue("{assay}_snn_res.{res}"))
    
    ggsave(glue("04_integration/plot/{assay}/UMAP_{assay}_post_{reduction}_integration_{assay}_snn_res_{res}_split.by_orig.ident.pdf"), 
           width = 12, 
           height = 8)
    
    
  }
  
  # Feature plots: UMAP stratified by continuous variable 
  
  features <- c("nFeature_RNA", "percent.mt", "percent.ribo")
  
  lapply(features, function(x) {
    
    FeaturePlot(seurat_integrated, reduction = umap_reduction.name, features = x)
    ggsave(glue("04_integration/plot/{assay}/UMAP_{assay}_post_{reduction}_integration_{x}.pdf"), 
           width = 8, 
           height = 7)
    
  })
  
  
}


######################### Save as h5ad file for python ######################### 

library(SeuratDisk)
library(rhdf5)

# All in one file 
obj_tmp <- seurat_integrated 

obj_tmp[["RNA3"]] <- as(object = obj_tmp[["RNA"]], Class = "Assay")
DefaultAssay(obj_tmp) <- "RNA3"
obj_tmp[["RNA"]] <- NULL
obj_tmp <- RenameAssays(object = obj_tmp, RNA3 = 'RNA')

obj_tmp[["SCT3"]] <- as(object = obj_tmp[["SCT"]], Class = "Assay")
DefaultAssay(obj_tmp) <- "SCT3"
obj_tmp[["SCT"]] <- NULL
obj_tmp <- RenameAssays(object = obj_tmp, SCT3 = 'SCT')

filename <- glue("04_integration/out/mydata_SCT.h5Seurat")
SaveH5Seurat(obj_tmp, filename = filename, overwrite = TRUE)
Convert(filename, dest = "h5ad", overwrite = TRUE, verbose = FALSE)













