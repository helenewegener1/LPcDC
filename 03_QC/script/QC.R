setwd("~/Documents/projects/project_cDC/LPcDC/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)

# Load data
seurat_obj_roughQC_list <- readRDS("02_roughQC/out/seurat_obj_roughQC_list.rds")

datasets <- names(seurat_obj_list)

dataset <- datasets[3]

# for (dataset in datasets){}
seurat_obj <- seurat_obj_list[[dataset]]

# Calculate QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rps|^Rpl")

# if large data set, subset cells for visualization purposes (it takes too long to plot all cells)
n_cells <- ncol(seurat_obj)

if (n_cells > 100000){
  # subset to 5k cells
  subset_cells <- sample(colnames(seurat_obj), 5000)  
  seurat_obj_subset <- seurat_obj[, subset_cells]
  
  # Plot QC metrics in violin plots
  VlnPlot(seurat_obj_subset, features = c("percent.mt", "percent.ribo"), ncol = 2)
  VlnPlot(seurat_obj_subset, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
  
} else {
  
  # Plot QC metrics in violin plots
  VlnPlot(seurat_obj, features = "percent.mt") + geom_hline(yintercept = 5, color = "black")
  VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
  
}

# Filter cells based on QC plots
seurat_obj_filtered <- subset(seurat_obj, subset = nFeature_RNA < 6000 & percent.mt < 5)

ncol(seurat_obj_filtered)

# Plot QC metrics in violin plots after filtering
VlnPlot(seurat_obj_filtered, features = c("percent.mt", "percent.ribo"), ncol = 2)
VlnPlot(seurat_obj_filtered, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)


# 3. Normalize
seurat_obj <- NormalizeData(seurat_obj)

# 4. Identify variable features
seurat_obj <- FindVariableFeatures(seurat_obj)

# 5. Optional: Scale and run PCA (helps for integration)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)


