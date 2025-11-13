# setwd("~/Documents/projects/project_cDC/LPcDC/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(stringr)
library(tibble)

# List data files for overview
data_files <- list.files("00_data/")
data_files <- data_files[startsWith(data_files, prefix = "G")]
data_files

############################################ Data set 1. Caspar Ohnmacht mLN ############################################

# Sample: GSM7789317 control mLN steady state

sample_name <- "GSM7789317"

# Read data
counts <- ReadMtx(
  mtx = glue("00_data/{sample_name}.matrix.mtx.gz"),
  features = glue("00_data/{sample_name}.features.tsv.gz"),
  cells = glue("00_data/{sample_name}.barcodes.tsv.gz")
)

# Create seurat object 
# This sample has 2.1 million cells 
# Many of the cells have a very low number of genes detected 
# So we do a rough filtering here
seurat_obj <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)

########################################## Export Seurat object ##########################################

saveRDS(seurat_obj, "21_mLN_make_seurat/out/mLN_seurat_obj.rds")
