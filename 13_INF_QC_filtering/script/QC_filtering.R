# setwd("~/Documents/projects/project_cDC/LPcDC/")
getwd()

# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)
library(DropletUtils)
library(glmGamPoi)

source("13_INF_QC_filtering/script/functions.R")

# Load data
sample_name <- "INF"
seurat_obj <- readRDS(glue("12_INF_QC/out/{sample_name}_seurat_obj_QC_metrics.rds") )

Idents(seurat_obj) <- "orig.ident"
seurat_obj <- subset(seurat_obj, subset = scDblFinder.class == "singlet")

# Calculate QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rps|^Rpl")
seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Hb[ab]")

# if large data set, subset cells for visualization purposes (it takes too long to plot all cells)
n_cells_raw <- ncol(seurat_obj) 

# Real raw n_cells: 2156566
# subset to 5k cells
# subset_cells <- sample(colnames(seurat_obj), 5000)  
# seurat_obj_subset <- seurat_obj[, subset_cells]

# Plot QC metrics in violin plots
plot_qc(seurat_obj = seurat_obj, 
        sample_name = sample_name, 
        n_cells = n_cells_raw, 
        version = "raw", 
        filtering = "")

# Extra plots 
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA > 1800 & nFeature_RNA < 6000 & percent.mt < 5)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

n_cells_filtered <- ncol(seurat_obj_filtered)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, 
        sample_name = sample_name, 
        n_cells = n_cells_filtered, 
        version = "filtered", 
        filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_finalQC <- seurat_obj_filtered

########################################## Export list of filtered Seurat objects ##########################################

saveRDS(seurat_obj_finalQC, "13_INF_QC_filtering/out/INF_seurat_obj_QC_filtered_list.rds")





