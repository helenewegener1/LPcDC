setwd("~/Documents/projects/project_cDC/LPcDC/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)
library(patchwork)

source("03_SILP_QC_filtering/script/functions.R")

# Load data
seurat_obj_list <- readRDS("02_SILP_QC/out/SILP_seurat_obj_QC_metrics.rds")

# Initialize filtered list
seurat_obj_QC_filtered_list <- list()

############################################ Data set 1. Caspar Ohnmacht ############################################

# Sample: GSM7789315 control SI-LP steady state
sample_name <- "GSM7789315"

seurat_obj <- seurat_obj_list[[sample_name]]
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
filtering_expr <- expr(nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 5)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

n_cells_filtered <- ncol(seurat_obj_filtered)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, 
        sample_name = sample_name, 
        n_cells = n_cells_filtered, 
        version = "filtered", 
        filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

############################################ Data set 2. Vuk Cerovic ############################################

# Sample: GSM8672515	Small intestinal dendritic cells CCR7gfp SI LP DCs
sample_name <- "GSM8672515"

seurat_obj <- seurat_obj_list[[sample_name]]
Idents(seurat_obj) <- "orig.ident"
seurat_obj <- subset(seurat_obj, subset = scDblFinder.class == "singlet")

# Calculate QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rps|^Rpl")

# if large data set, subset cells for visualization purposes (it takes too long to plot all cells)
n_cells_raw <- ncol(seurat_obj) # 1963

# Plot QC metrics in violin plots
plot_qc(seurat_obj = seurat_obj, 
        sample_name = sample_name, 
        n_cells = n_cells_raw, 
        version = "raw", 
        filtering = "")

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA < 6000 & nFeature_RNA > 400 & percent.mt < 4)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

n_cells_filtered <- ncol(seurat_obj_filtered) 

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, 
        sample_name = sample_name, 
        n_cells = n_cells_filtered, 
        version = "filtered", 
        filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Sample: GSM9122899	Small intestinal dendritic cells from WT mice
sample_name <- "GSM9122899"

seurat_obj <- seurat_obj_list[[sample_name]]
Idents(seurat_obj) <- "orig.ident"
seurat_obj <- subset(seurat_obj, subset = scDblFinder.class == "singlet")

# Calculate QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rps|^Rpl")

# if large data set, subset cells for visualization purposes (it takes too long to plot all cells)
n_cells_raw <- ncol(seurat_obj) # 11414

# Plot QC metrics in violin plots
plot_qc(seurat_obj = seurat_obj, 
        sample_name = sample_name, 
        n_cells = n_cells_raw, 
        version = "raw", 
        filtering = "")

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA < 6000 & nFeature_RNA > 400 & percent.mt < 2.5)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

n_cells_filtered <- ncol(seurat_obj_filtered) 

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, 
        sample_name = sample_name, 
        n_cells = n_cells_filtered, 
        version = "filtered", 
        filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

############################################ Data set 3. Vuk Cerovic ############################################

# THIS IS BULK DATA!

############################################ Data set 4. George Kollias ############################################

# GSE255350_myeloid_raw_counts.txt.gz
sample_name <- "GSE255350"

seurat_obj <- seurat_obj_list[[sample_name]]
Idents(seurat_obj) <- "orig.ident"
seurat_obj <- subset(seurat_obj, subset = scDblFinder.class == "singlet")

# Calculate QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rps|^Rpl")

# if large data set, subset cells for visualization purposes (it takes too long to plot all cells)
n_cells_raw <- ncol(seurat_obj) # 2228

# Plot QC metrics in violin plots
plot_qc(seurat_obj = seurat_obj, 
        sample_name = sample_name, 
        n_cells = n_cells_raw, 
        version = "raw", 
        filtering = "")

# Filter cells based on QC plots
filtering_expr <- expr(nFeature_RNA < 6000 & nFeature_RNA > 2000 & percent.mt < 3.5)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

n_cells_filtered <- ncol(seurat_obj_filtered) 

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, 
        sample_name = sample_name, 
        n_cells = n_cells_filtered, 
        version = "filtered", 
        filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

############################################ Data set 5. Fiona Powrie ############################################

# sample_name <- "CRAM1"
# 
# seurat_obj <- seurat_obj_list[[sample_name]]
# 
# # Calculate QC metrics
# seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
# seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rps|^Rpl")
# 
# # if large data set, subset cells for visualization purposes (it takes too long to plot all cells)
# n_cells_raw <- ncol(seurat_obj) # 224334
# 
# # # Real raw n_cells: 224334
# # # subset to 5k cells
# # subset_cells <- sample(colnames(seurat_obj), 5000)
# # seurat_obj_subset <- seurat_obj[, subset_cells]
# 
# # Plot QC metrics in violin plots
# plot_qc(seurat_obj = seurat_obj, 
#         sample_name = sample_name, 
#         n_cells = n_cells_raw, 
#         version = "pre_filtered", 
#         filtering = "")
# 
# # # Filter cells based on QC plots
# # filtering_expr <- expr(nFeature_RNA < 2500 & nFeature_RNA > 400 & percent.mt < 2.5)
# # seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)
# # 
# # n_cells_filtered <- ncol(seurat_obj_filtered) 
# # 
# # # Plot QC metrics in violin plots after filtering
# # plot_qc(seurat_obj = seurat_obj_filtered, 
# #         sample_name = sample_name, 
# #         n_cells = n_cells_filtered, 
# #         version = "filtered", 
# #         filtering = rlang::expr_text(filtering_expr))
# 
# # Save filtered seurat object
# seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj
# 
# # Clean up
# rm(seurat_obj, n_cells_raw)
# 
# ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# 
# sample_name <- "CRAM2"
# 
# seurat_obj <- seurat_obj_list[[sample_name]]
# 
# # Calculate QC metrics
# seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
# seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rps|^Rpl")
# 
# # if large data set, subset cells for visualization purposes (it takes too long to plot all cells)
# n_cells_raw <- ncol(seurat_obj) # 281591
# 
# # Real raw n_cells: 281591
# # subset to 5k cells
# # subset_cells <- sample(colnames(seurat_obj), 5000)
# # seurat_obj_subset <- seurat_obj[, subset_cells]
# 
# # Plot QC metrics in violin plots
# plot_qc(seurat_obj = seurat_obj, 
#         sample_name = sample_name, 
#         n_cells = n_cells_raw, 
#         version = "pre_filtered", 
#         filtering = "")
# 
# # # Filter cells based on QC plots
# # filtering_expr <- expr(nFeature_RNA < 2500 & nFeature_RNA > 400 & percent.mt < 2.5)
# # seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)
# # 
# # n_cells_filtered <- ncol(seurat_obj_filtered) 
# # 
# # # Plot QC metrics in violin plots after filtering
# # plot_qc(seurat_obj = seurat_obj_filtered, 
# #         sample_name = sample_name, 
# #         n_cells = n_cells_filtered, 
# #         version = "filtered", 
# #         filtering = rlang::expr_text(filtering_expr))
# 
# # Save filtered seurat object
# seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj
# 
# # Clean up
# rm(seurat_obj, n_cells_raw)

############################################ Data set 6. Rivera CA ############################################

# Sample: GSM5678427 CD103+CD11b+ dendritic cells from lamina propria
sample_name <- "GSM5678427"

seurat_obj <- seurat_obj_list[[sample_name]]
Idents(seurat_obj) <- "orig.ident"
seurat_obj <- subset(seurat_obj, subset = scDblFinder.class == "singlet")

# Calculate QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rps|^Rpl")
seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Hb[ab]")

# if large data set, subset cells for visualization purposes (it takes too long to plot all cells)
n_cells_raw <- ncol(seurat_obj) 

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
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 4000 & percent.mt < 10)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

n_cells_filtered <- ncol(seurat_obj_filtered)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, 
        sample_name = sample_name, 
        n_cells = n_cells_filtered, 
        version = "filtered", 
        filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

############################################ Data set 7. ############################################

# Sample: 

sample_name <- "E_MTAB_9522"

seurat_obj <- seurat_obj_list[[sample_name]]
Idents(seurat_obj) <- "orig.ident"
seurat_obj <- subset(seurat_obj, subset = scDblFinder.class == "singlet")

# Calculate QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rps|^Rpl")
seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Hb[ab]")

# if large data set, subset cells for visualization purposes (it takes too long to plot all cells)
n_cells_raw <- ncol(seurat_obj) 

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
filtering_expr <- expr(nFeature_RNA > 400 & nFeature_RNA < 6000 & percent.mt < 10)
seurat_obj_filtered <- subset(seurat_obj, subset = !!filtering_expr)

n_cells_filtered <- ncol(seurat_obj_filtered)

# Plot QC metrics in violin plots after filtering
plot_qc(seurat_obj = seurat_obj_filtered, 
        sample_name = sample_name, 
        n_cells = n_cells_filtered, 
        version = "filtered", 
        filtering = rlang::expr_text(filtering_expr))

# Save filtered seurat object
seurat_obj_QC_filtered_list[[sample_name]] <- seurat_obj_filtered

# Clean up
rm(seurat_obj, seurat_obj_filtered, n_cells_raw, n_cells_filtered)

########################################## Export list of filtered Seurat objects ##########################################

names(seurat_obj_QC_filtered_list)

saveRDS(seurat_obj_QC_filtered_list, "03_SILP_QC_filtering/out/SILP_seurat_obj_QC_filtered_list.rds")
