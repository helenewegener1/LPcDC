setwd("~/Documents/projects/project_cDC/LPcDC/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(stringr)
library(tibble)

# Load files 
sample_name <- "STROMAL_GSM5819064"

# Read data
counts <- ReadMtx(
  mtx = glue("00_data/{sample_name}/filtered_feature_bc_matrix/matrix.mtx.gz"),
  features = glue("00_data/{sample_name}/filtered_feature_bc_matrix/features.tsv.gz"),
  cells = glue("00_data/{sample_name}/filtered_feature_bc_matrix/barcodes.tsv.gz")
)

seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)

# # Quick workflow
# seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
# VlnPlot(seurat_obj, features = c("nFeature_RNA", "percent.mt"))
# 
# # Filter cells based on QC plots
# filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 5)
# seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)
# 
# # Seurat workflow
# seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
# seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
# seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
# DefaultAssay(seurat_obj)
# seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
# ElbowPlot(seurat_obj)
# seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20, verbose = FALSE)
# seurat_obj <- FindClusters(seurat_obj, resolution = 0.8, verbose = FALSE)
# seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:20, verbose = FALSE)
# 
# DimPlot(seurat_obj)

# RDS file
SeuratObject_mLN <- readRDS("00_data/STROMAL_GSM5819064/SeuratObject_mLN_day_0_10_24_56_300_nonENDO_SC.Rds")
SeuratObject_mLN <- UpdateSeuratObject(SeuratObject_mLN)

# Check QC
range(SeuratObject_mLN$percent.mito)
range(SeuratObject_mLN$nFeature_RNA)

# Different days 
SeuratObject_mLN$orig.ident %>% table()

head(SeuratObject_mLN@meta.data)

n_cells <- ncol(SeuratObject_mLN)

DimPlot(SeuratObject_mLN)
DimPlot(SeuratObject_mLN, group.by = "orig.ident")

# Make UMAP
SeuratObject_mLN <- NormalizeData(SeuratObject_mLN, verbose = FALSE)
SeuratObject_mLN <- FindVariableFeatures(SeuratObject_mLN, verbose = FALSE)
SeuratObject_mLN <- ScaleData(SeuratObject_mLN, verbose = FALSE)
# DefaultAssay(seurat_obj)
SeuratObject_mLN <- RunPCA(SeuratObject_mLN, verbose = FALSE)
ElbowPlot(SeuratObject_mLN)
SeuratObject_mLN <- FindNeighbors(SeuratObject_mLN, dims = 1:20, verbose = FALSE)
SeuratObject_mLN <- FindClusters(SeuratObject_mLN, resolution = 0.8, verbose = FALSE)
SeuratObject_mLN <- RunUMAP(SeuratObject_mLN, reduction = "pca", dims = 1:20, verbose = FALSE)

DimPlot(SeuratObject_mLN)
DimPlot(SeuratObject_mLN, group.by = "orig.ident")
DimPlot(SeuratObject_mLN, group.by = "Cluster")

Reductions(SeuratObject_mLN)

saveRDS(SeuratObject_mLN, "31_STROMAL_make_seurat/out/STROMAL_object.rds")

# QC is dones
# We need to integrate

########################################## Export of Seurat object ##########################################

saveRDS(seurat_obj, "31_STROMAL_make_seurat/out/STROMAL_seurat_obj.rds")
