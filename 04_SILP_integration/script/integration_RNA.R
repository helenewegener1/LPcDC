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
library(clustree)

# Load data
seurat_obj_list <- readRDS("03_SILP_QC_filtering/out/SILP_seurat_obj_QC_filtered_list.rds")

############################ RNA integration prep ############################

# Merge
seurat_merged <- merge(seurat_obj_list[[1]], y = seurat_obj_list[-1])

Layers(seurat_merged[["RNA"]]) # already split

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
# saveRDS(seurat_merged, "04_integration/out/seurat_merged_SCT.rds")
# seurat_merged <- readRDS("04_integration/out/seurat_merged_SCT.rds")

DefaultAssay(seurat_merged)

# Visualize with UMAP stratified by dataset - pre integration
DimPlot(seurat_merged, reduction = "umap", group.by = "orig.ident", split.by = "orig.ident") +
  NoLegend() + 
  labs(title = "UMAP - RNA - pre integration") +
  theme(legend.text = element_text(size = 8))

ggsave("04_SILP_integration/plot/RNA/UMAP_RNA_pre_integration_orig.ident.pdf",
       width = 15,
       height = 5)

DefaultAssay(seurat_merged)


############################## Integration on RNA ##############################

seurat_integrated <- seurat_merged

DefaultAssay(seurat_integrated) <- "RNA"

Layers(seurat_integrated[["RNA"]])

# seurat_integrated <- seurat_merged
# 
# DefaultAssay(seurat_integrated) <- "SCT"
# 
# Layers(seurat_integrated[["SCT"]])
# 
# Reductions(seurat_integrated)

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

# options(future.globals.maxSize = 8000 * 1024^2)  # 8 GB
# seurat_integrated <- IntegrateLayers(
#   object = seurat_integrated,
#   method = RPCAIntegration,
#   orig.reduction = "pca",
#   new.reduction = "RNA_integrated.rpca",
#   assay = "RNA",
#   verbose = FALSE
# )
# 
# seurat_integrated <- IntegrateLayers(
#   object = seurat_integrated,
#   method = FastMNNIntegration,
#   orig.reduction = "pca",
#   new.reduction = "RNA_integrated.mnn",
#   assay = "RNA",
#   verbose = FALSE
# )

################### Export list of integrated Seurat objects ################### 

Reductions(seurat_integrated)

####################### Run UMAP using Harmony embedding #######################

reductions <- list(
  
  c("RNA_integrated.cca", "RNA_umap.cca", "RNA_cca_clusters"),
  c("RNA_integrated.harmony", "RNA_umap.harmony", "RNA_harmony_clusters")
  # c("RNA_integrated.mnn", "RNA_umap.mnn", "RNA_mnn_clusters"),
  # c("RNA_integrated.rpca", "RNA_umap.rpca", "RNA_rpca_clusters")
  
)

for (red in reductions){
  
  red <- c("RNA_integrated.cca", "RNA_umap.cca", "RNA_cca_clusters")
  
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
  
  ggsave(glue("04_SILP_integration/plot/{assay}/UMAP_{reduction}_orig.ident.pdf"), 
         width = 8, 
         height = 7)
  
  Idents(seurat_integrated) <- "orig.ident"
  DimPlot(seurat_integrated, reduction = umap_reduction.name, split.by = "orig.ident", ncol = 3) +
    labs(title = glue("UMAP - post {reduction}")) + 
    theme(legend.text = element_text(size = 8))
  
  ggsave(glue("04_SILP_integration/plot/{assay}/UMAP_{reduction}_orig.ident_split.pdf"), 
         width = 12, 
         height = 8)
  
  # Visualize with UMAP stratified by seurat clusters - post harmony integration 
  
  res_list <- seq(0.1, 1.5, by = 0.1)
  
  for (res in res_list){
    
    res <- 1
    
    seurat_integrated <- FindClusters(seurat_integrated, resolution = res, cluster.name = glue("{cluster.name}_res.{res}"))
    
    # WHAT TO NAME THEM?
    DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = glue("{cluster.name}_res.{res}"), label = TRUE) +
      labs(title = glue("UMAP - post {reduction}"),
           subtitle = glue("{cluster.name}_res.{res}"))
    
    ggsave(glue(glue("04_SILP_integration/plot/{assay}/UMAP_{reduction}_{assay}_snn_res_{res}.pdf")), 
           width = 8, 
           height = 7)
    
    DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = glue("{cluster.name}_res.{res}"), split.by = "orig.ident", ncol = 3) +
      labs(title = glue("UMAP - post {reduction}"),
           subtitle = glue("{cluster.name}_res.{res}"))
    
    ggsave(glue("04_SILP_integration/plot/{assay}/UMAP_{reduction}_{assay}_snn_res_{res}_split.by_orig.ident.pdf"), 
           width = 12, 
           height = 8)
    
    
  }

  
  
  # Feature plots: UMAP stratified by continuous variable 
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

# clustree

cluster.name <- "RNA_cca_clusters"
pdf(file = glue("04_SILP_integration/plot/{assay}/clustree_{cluster.name}.pdf"), width = 12, height = 12)
clustree(seurat_integrated, assay = "RNA", return = "plot", prefix = glue("{cluster.name}_res."))
dev.off()

cluster.name <- "RNA_harmony_clusters"
pdf(file = glue("04_SILP_integration/plot/{assay}/clustree_{cluster.name}.pdf"), width = 12, height = 12)
clustree(seurat_integrated, assay = "RNA", return = "plot", prefix = glue("{cluster.name}_res."))
dev.off()

# cluster.name <- "RNA_rpca_clusters"
# pdf(file = glue("04_SILP_integration/plot/{assay}/clustree_{cluster.name}.pdf"), width = 12, height = 12)
# clustree(seurat_integrated, assay = "RNA", return = "plot", prefix = glue("{cluster.name}_res."))
# dev.off()
# 
# cluster.name <- "RNA_mnn_clusters"
# pdf(file = glue("04_SILP_integration/plot/{assay}/clustree_{cluster.name}.pdf"), width = 12, height = 12)
# clustree(seurat_integrated, assay = "RNA", return = "plot", prefix = glue("{cluster.name}_res."))
# dev.off()

######################### Save as h5ad file for python ######################### 

Reductions(seurat_integrated)

saveRDS(seurat_integrated, "04_SILP_integration/out/SILP_seurat_integrated_v5_RNA.rds")
# seurat_integrated <- readRDS("04_SILP_integration/out/seurat_integrated_v5_RNA.rds")

######################## Save as h5ad file for python #########################

library(SeuratDisk)
library(rhdf5)

# All in one file 
obj_tmp <- seurat_integrated 

DefaultAssay(obj_tmp)

Layers(obj_tmp[["RNA"]])

# IMPORTANT: Join layers before export (h5ad doesn't support split layers)
obj_tmp[["RNA"]] <- JoinLayers(obj_tmp[["RNA"]])

obj_tmp[["RNA3"]] <- as(object = obj_tmp[["RNA"]], Class = "Assay")
DefaultAssay(obj_tmp) <- "RNA3"
obj_tmp[["RNA"]] <- NULL
obj_tmp <- RenameAssays(object = obj_tmp, RNA3 = 'RNA')

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

filename <- glue("04_SILP_integration/out/mydata_RNA_v5.h5Seurat")
SaveH5Seurat(obj_tmp, filename = filename, overwrite = TRUE)
Convert(filename, dest = "h5ad", overwrite = TRUE, verbose = FALSE)














