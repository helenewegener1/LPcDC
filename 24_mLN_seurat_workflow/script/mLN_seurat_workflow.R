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
seurat_obj <- readRDS("23_mLN_QC_filtering/out/mLN_seurat_obj_QC_filtered.rds")

umap_reduction.name <- "umap"

DimPlot(seurat_obj, reduction = umap_reduction.name) + 
  labs(title = "UMAP - RNA - pre integration") +
  theme(legend.text = element_text(size = 8))

ggsave("14_INF_integration/plot/UMAP_RNA.pdf", width = 12, height = 7)


# Either RNA or SCT 
assay <- "RNA"

# Set default assay 
DefaultAssay(seurat_obj) <- assay

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 2)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

res_list <- seq(0.1, 1.5, by = 0.1)

for (res in res_list){

  # res <- 0.3
  
  cluster.name <- glue("RNA_snn_res.{res}")
  seurat_obj <- FindClusters(seurat_obj, resolution = res, cluster.name = cluster.name)
  
  # WHAT TO NAME THEM?
  DimPlot(seurat_obj, reduction = umap_reduction.name, group.by = cluster.name, label = TRUE)
  
  ggsave(glue(glue("24_mLN_seurat_workflow/plot/{assay}/UMAP_{cluster.name}.pdf")), 
         width = 8, 
         height = 7)

}

# clustree

# cluster.name <- "RNA_cca_clusters"
# pdf(file = glue("14_inflammation_integration/plot/{assay}/clustree_{cluster.name}.pdf"), width = 12, height = 10)
# clustree(seurat_obj, assay = "RNA", return = "plot", prefix = glue("{cluster.name}_res."))
# dev.off()

cluster.name <- "RNA"
pdf(file = glue("24_mLN_seurat_workflow/plot/{assay}/clustree_{cluster.name}.pdf"), width = 12, height = 10)
clustree(seurat_obj, assay = "RNA", return = "plot", prefix = glue("RNA_snn_res."))
dev.off()

# cluster.name <- "RNA_rpca_clusters"
# pdf(file = glue("14_inflammation_integration/plot/{assay}/clustree_{cluster.name}.pdf"), width = 12, height = 10)
# clustree(seurat_obj, assay = "RNA", return = "plot", prefix = glue("{cluster.name}_res."))
# dev.off()
# 
# cluster.name <- "RNA_mnn_clusters"
# pdf(file = glue("14_inflammation_integration/plot/{assay}/clustree_{cluster.name}.pdf"), width = 12, height = 10)
# clustree(seurat_obj, assay = "RNA", return = "plot", prefix = glue("{cluster.name}_res."))
# dev.off()

######################### Save as h5ad file for python ######################### 

Reductions(seurat_obj)

saveRDS(seurat_obj, "24_mLN_seurat_workflow/out/mLN_seurat_obj_v5_RNA.rds")
# seurat_obj <- readRDS("14_inflammation_integration/out/seurat_obj_v5_RNA.rds")

