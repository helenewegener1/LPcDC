# setwd("~/Documents/projects/project_cDC/LPcDC/")
getwd()

# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)
library(patchwork)

source("02_roughQC/script/functions.R")

# Load data
seurat_obj <- readRDS("11_inflammation_make_seurat/out/INF_seurat_obj.rds")

# Calculate QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rps|^Rpl")
seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Hb[ab]")

# if large data set, subset cells for visualization purposes (it takes too long to plot all cells)
n_cells_raw <- ncol(seurat_obj) 

# Plot QC metrics in violin plots
plot_qc(seurat_obj = seurat_obj, 
        sample_name = "INF", 
        n_cells = n_cells_raw, 
        version = "raw", 
        filtering = "")

ggsave("12_inflammation_roughQC/plot/raw_QC_plot_INF.pdf", 
       width = 9, 
       height = 8)

# Extra plots 
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 2000 & nFeature_RNA < 7500 & percent.mt < 5)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

n_cells_filtered <- ncol(seurat_obj_filtered)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, 
        sample_name = "INF", 
        n_cells = n_cells_filtered, 
        version = "filtered", 
        filtering = rlang::expr_text(filtering_expr))

ggsave("12_inflammation_roughQC/plot/filtered_QC_plot_INF.pdf", 
       width = 9, 
       height = 8)


saveRDS(seurat_obj_filtered, "12_inflammation_roughQC/out/INF_seurat_obj_roughQC.rds")
