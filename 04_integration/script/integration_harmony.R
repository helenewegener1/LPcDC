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

# Save merged object
saveRDS(seurat_merged, "04_integration/out/seurat_merged.rds")

# Visualize with UMAP stratified by dataset - pre integration 
DimPlot(seurat_merged, reduction = "umap", group.by = "orig.ident") + 
  labs(title = "UMAP - pre integration") + 
  theme(legend.text = element_text(size = 8))

ggsave("04_integration/plot/UMAP_pre_integration_orig.ident.pdf", 
       width = 8, 
       height = 7)

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
  labs(title = "UMAP - post harmonny integration") + 
  theme(legend.text = element_text(size = 8))
ggsave("04_integration/plot/UMAP_post_harmony_integration_orig.ident.pdf", 
       width = 8, 
       height = 7)

DimPlot(seurat_integrated, reduction = "umap", split.by = "orig.ident", ncol = 3)

# Visualize with UMAP stratified by seurat clusters - post harmony integration 

res_list <- c(0.3, 0.5, 0.7)

for (res in res_list){
  
  # res <- 0.3
  
  seurat_integrated <- FindClusters(seurat_integrated, resolution = res)
  
  DimPlot(seurat_integrated, reduction = "umap", group.by = glue("SCT_snn_res.{res}"), label = TRUE) +
    labs(title = "UMAP - post harmony integration",
         subtitle = glue("SCT_snn_res.{res}"))
  
  ggsave(glue("04_integration/plot/UMAP_post_harmony_integration_SCT_snn_res_{res}.pdf"), 
         width = 8, 
         height = 7)
  
  DimPlot(seurat_integrated, reduction = "umap", group.by = glue("SCT_snn_res.{res}"), split.by = "orig.ident", ncol = 3) +
    labs(title = "UMAP - post harmony integration",
         subtitle = glue("SCT_snn_res.{res}"))
  
  ggsave(glue("04_integration/plot/UMAP_post_harmony_integration_SCT_snn_res_{res}_split.by_orig.ident.pdf"), 
         width = 12, 
         height = 8)
  
  
}

# Feature plots: UMAP stratified by continuous variable 

features <- c("nFeature_RNA", "percent.mt", "percent.ribo")

lapply(features, function(x) {
  
  FeaturePlot(seurat_integrated, reduction = "umap", features = x)
  ggsave(glue("04_integration/plot/UMAP_post_harmony_integration_{x}.pdf"), 
         width = 8, 
         height = 7)
  
})

################### Export list of integrated Seurat objects ################### 

saveRDS(seurat_integrated, "04_integration/out/seurat_integrated.rds")


 












