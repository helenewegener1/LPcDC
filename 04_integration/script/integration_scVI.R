setwd("~/Documents/projects/project_cDC/LPcDC/")

# scVI doesn’t fully support Seurat v5 objects yet.

# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)
# install.packages("reticulate")
library(reticulate)
# devtools::install_github("cellgeni/sceasy")
library(sceasy)

# reticulate::py_config()
# reticulate::py_install(c("scanpy", "scvi-tools"))
sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)

# Load data
seurat_obj_finalQC_list <- readRDS("03_QC/out/seurat_obj_finalQC_list.rds")
seurat_merged <- readRDS("04_integration/out/seurat_merged.rds")

# Integration using scVI (R wrapper)
# Following this tutorial: https://docs.scvi-tools.org/en/1.0.0/tutorials/notebooks/scvi_in_R.html - Integrating datasets with scVI

adata <- sceasy::convertFormat(seurat_merged, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)


adata <- sceasy::convertFormat(seurat_obj_finalQC_list$GSM7789315, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)

# ChatGPT workaround 

library(Seurat)

# Extract counts from the layer
counts_mat <- seurat_obj_finalQC_list$GSM7789315[["RNA"]]@layers[["counts"]]

# Recreate a classic Assay object
assay_v4 <- CreateAssayObject(counts = counts_mat)

# Replace the assay in a temporary Seurat object
seurat_v4 <- seurat_obj_finalQC_list$GSM7789315
seurat_v4[["RNA"]] <- assay_v4

# Now try conversion
adata <- sceasy::convertFormat(
  seurat_v4,
  from = "seurat",
  to = "anndata",
  main_layer = "counts",
  drop_single_values = FALSE
)

