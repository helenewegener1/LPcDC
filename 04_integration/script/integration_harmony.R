setwd("~/Documents/projects/project_cDC/LPcDC/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(glmGamPoi)
library(dplyr)
library(glue)
library(ggplot2)
library(harmony)

# Load data
seurat_obj_finalQC_list <- readRDS("03_QC/out/seurat_obj_finalQC_list.rds")

############### SCTransform individual datasets after filtering ################ 

seurat_obj_list <- list()

for (sample_name in names(seurat_obj_finalQC_list)) {
  
  # Define seurat object of sample 
  seurat_obj_finalQC <- seurat_obj_finalQC_list[[sample_name]]
  
  # SCTransform gives a warning because the SCTransform assay from doublets were removed is overwritten
  seurat_obj_tmp <- SCTransform(seurat_obj_finalQC, assay = "RNA", layer = "counts", verbose = FALSE)
  
  # Save in list 
  seurat_obj_list[[sample_name]] <- seurat_obj_tmp
  
  # Clean up
  rm(seurat_obj_tmp)
  
}

######################### Merge into one Seurat object #########################

# Merge
seurat_merged <- merge(seurat_obj_list[[1]], y = seurat_obj_list[-1])

# Control all cells from each dataset are included 
seurat_merged@meta.data$orig.ident %>% table(useNA = "always")

# Find variable features 
# features <- SelectIntegrationFeatures(object.list = seurat_obj_finalQC_list, nfeatures = 3000)
# VariableFeatures(seurat_merged) <- features
rm(seurat_obj_finalQC_list)

seurat_merged <- FindVariableFeatures(seurat_merged, assay = "SCT", nfeatures = 3000)

# Seurat workflow on merged data 
seurat_merged <- RunPCA(seurat_merged, assay = "SCT", verbose = FALSE) # HW: Yes - SCT assay for PCA 

# ^^ Needs the part above to run harmony ^^

# This is just for plotting UMAP
ElbowPlot(seurat_merged)
seurat_merged <- FindNeighbors(seurat_merged, dims = 1:20)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.4)
seurat_merged <- RunUMAP(seurat_merged, reduction = "pca", dims = 1:20)

# Visualize with UMAP stratified by dataset - pre integration 
DimPlot(seurat_merged, reduction = "umap", group.by = "orig.ident") + 
  labs(title = "UMAP - pre integration")
ggsave("04_integration/plot/UMAP_pre_integration_orig.ident.pdf")

########################## Integration using Harmony ###########################

# Run Harmony on the merged object
seurat_integrated <- RunHarmony(
  seurat_merged,
  group.by.vars = "orig.ident",  # metadata column indicating dataset
  dims = 1:20
)

####################### Run UMAP using Harmony embedding #######################

seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "harmony", dims = 1:20)
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "harmony", dims = 1:20)

# Visualize with UMAP stratified by dataset - post harmony integration 
DimPlot(seurat_integrated, reduction = "umap", group.by = "orig.ident") +
  labs(title = "UMAP - post harmonny integration")
ggsave("04_integration/plot/UMAP_post_harmony_integration_orig.ident.pdf")

# Visualize with UMAP stratified by seurat clusters - post harmony integration 

res_list <- c(0.3, 0.5, 0.7)

for (res in res_list){
  
  # res <- 0.2
  
  seurat_integrated <- FindClusters(seurat_integrated, resolution = res)
  DimPlot(seurat_integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
    labs(title = "UMAP - post harmony integration",
         subtitle = glue("res: {res}"))
  
  ggsave(glue("04_integration/plot/UMAP_post_harmony_integration_seurat_clusters_res_{res}.pdf"))
  
}
















