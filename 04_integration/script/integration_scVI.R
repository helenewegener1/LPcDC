setwd("~/Documents/projects/project_cDC/LPcDC/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)

# Load data
seurat_obj_finalQC_list <- readRDS("03_QC/out/seurat_obj_finalQC_list.rds")

# Integration using scVI (R wrapper)
