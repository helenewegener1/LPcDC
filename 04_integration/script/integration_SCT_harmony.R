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

######################### Merge into one Seurat object #########################

# Merge
seurat_obj_list <- seurat_obj_finalQC_list
seurat_merged <- merge(seurat_obj_list[[1]], y = seurat_obj_list[-1])
DefaultAssay(seurat_merged) <- "SCT"

# Control all cells from each dataset are included 
seurat_merged@meta.data$orig.ident %>% table(useNA = "always")

# Find variable features 
# features <- SelectIntegrationFeatures(object.list = seurat_obj_finalQC_list, nfeatures = 3000)
# VariableFeatures(seurat_merged) <- features
# rm(seurat_obj_finalQC_list)

# Seurat workflow on merged data


# features <- SelectIntegrationFeatures(seurat_obj_list)
# seurat_obj_list <- PrepSCTIntegration(seurat_obj_list, anchor.features = features)
# seurat_merged <- merge(seurat_obj_list[[1]], y = seurat_obj_list[-1])

seurat_merged <- SCTransform(seurat_merged, assay = "RNA", layer = "counts", verbose = FALSE)
# VariableFeatures(seurat_merged) <- SelectIntegrationFeatures(
#   object.list = SplitObject(seurat_merged, split.by = "orig.ident"),
#   nfeatures = 3000
# )

# Seurat workflow on merged data 
seurat_merged <- RunPCA(seurat_merged, assay = "SCT", verbose = FALSE, features = features) # HW: Yes - SCT assay for PCA 

# ^^ Needs the part above to run harmony ^^

# This is just for plotting UMAP
ElbowPlot(seurat_merged)
seurat_merged <- FindNeighbors(seurat_merged, dims = 1:30)
seurat_merged <- FindClusters(seurat_merged, resolution = 0.4)
seurat_merged <- RunUMAP(seurat_merged, reduction = "pca", dims = 1:30)

# Visualize with UMAP stratified by dataset - pre integration 
DimPlot(seurat_merged, reduction = "umap", group.by = "orig.ident") + 
  labs(title = "UMAP - pre integration")
ggsave("04_integration/plot/v4/SCT/UMAP_pre_integration_orig.ident.pdf")

########################## Integration using Harmony ###########################

DefaultAssay(seurat_merged)

# Run Harmony on the merged object
seurat_integrated <- RunHarmony(
  seurat_merged,
  group.by.vars = "orig.ident",  # metadata column indicating dataset
  dims = 1:30
)

####################### Run UMAP using Harmony embedding #######################

seurat_integrated <- FindNeighbors(seurat_integrated, reduction = "harmony", dims = 1:30)
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "harmony", dims = 1:30)

# Visualize with UMAP stratified by dataset - post harmony integration 
DimPlot(seurat_integrated, reduction = "umap", group.by = "orig.ident") +
  labs(title = "UMAP - post harmony integration")
ggsave("04_integration/plot/v4/SCT/UMAP_post_harmony_integration_orig.ident.pdf")

saveRDS(seurat_integrated, "04_integration/out/seurat_SCT_harmony.rds")


# Visualize with UMAP stratified by seurat clusters - post harmony integration 

res_list <- c(0.3, 0.5, 0.7)

for (res in res_list){
  
  # res <- 0.2
  
  seurat_integrated <- FindClusters(seurat_integrated, resolution = res)
  DimPlot(seurat_integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
    labs(title = "UMAP - post harmony integration",
         subtitle = glue("res: {res}"))
  
  ggsave(glue("04_integration/plot/v4/SCT/UMAP_post_harmony_integration_seurat_clusters_res_{res}.pdf"))
  
}

######################### Save as h5ad file for python ######################### 

library(SeuratDisk)
library(rhdf5)

# All in one file 
obj_tmp <- seurat_integrated 

DefaultAssay(obj_tmp)

Layers(obj_tmp[["SCT"]])

obj_tmp[["SCT3"]] <- as(object = obj_tmp[["SCT"]], Class = "Assay")
DefaultAssay(obj_tmp) <- "SCT3"
obj_tmp[["SCT"]] <- NULL
obj_tmp <- RenameAssays(object = obj_tmp, SCT3 = 'SCT')


# Extract the embeddings
embedding <- Embeddings(obj_tmp, reduction = "harmony")

obj_tmp[["SCT_harmony"]] <- CreateDimReducObject(
  embeddings = embedding, 
  key = "SCT_harmony",
  assay = DefaultAssay(obj_tmp),
  global = TRUE
)


# Remove the problematic reduction to keep things clean
obj_tmp[["harmony"]] <- NULL

# Verify the new reduction name (optional)
Reductions(obj_tmp)

filename <- glue("04_integration/out/mydata_v4_sct.h5Seurat")
SaveH5Seurat(obj_tmp, filename = filename, overwrite = TRUE)
Convert(filename, dest = "h5ad", overwrite = TRUE, verbose = FALSE)



