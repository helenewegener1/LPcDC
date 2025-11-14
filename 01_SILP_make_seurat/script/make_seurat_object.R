setwd("~/Documents/projects/project_cDC/LPcDC/")

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

seurat_obj_list <- list()

############################################ Data set 1. Caspar Ohnmacht ############################################

# Sample: GSM7789315 control SI-LP steady state

sample_name <- "GSM7789315"

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
seurat_obj_list[[sample_name]] <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)

rm(counts, sample_name)

############################################ Data set 2. Vuk Cerovic ############################################

# Sample: GSM8672515	Small intestinal dendritic cells CCR7gfp SI LP DCs

sample_name <- "GSM8672515"

# Read data
counts <- ReadMtx(
  mtx = glue("00_data/{sample_name}.matrix.mtx.gz"),
  features = glue("00_data/{sample_name}.features.tsv.gz"),
  cells = glue("00_data/{sample_name}.barcodes.tsv.gz")
)

# Create seurat object 
seurat_obj_list[[sample_name]] <- CreateSeuratObject(counts = counts, project = sample_name)
rm(counts, sample_name)

# Sample: GSM9122899	Small intestinal dendritic cells from WT mice

sample_name <- "GSM9122899"

# Read data
counts <- Read10X_h5(glue("00_data/{sample_name}.feature.matrix.h5"))

# Create seurat object 
seurat_obj_list[[sample_name]] <- CreateSeuratObject(counts = counts, project = sample_name)
rm(counts, sample_name)

############################################ Data set 3. Vuk Cerovic ############################################

# THIS IS BULK DATA!

# Download all DC read counts
# Data should be filtered for the gut samples (exclude duct and node)
# gut:	Intestinal lamina propria	(LP)
# duct: Mesenteric lymphatic duct	(mL)
# node: Mesenteric lymph node	(mLN)

# sample_name <- "GSE160156"
# 
# 
# # Read data
# counts <- read.csv(gzfile(glue("00_data/{sample_name}.read_counts.csv.gz")), sep = '\t', row.names = 1)
# counts_mat <- as.matrix(counts)
# 
# # Create seurat object 
# seurat_obj_list[[sample_name]] <- CreateSeuratObject(counts = counts_mat, project = sample_name)
# rm(counts, sample_name)

############################################ Data set 4. George Kollias ############################################

# GSE255350_myeloid_raw_counts.txt.gz

sample_name <- "GSE255350"

# Read data
counts <- read.csv(gzfile(glue("00_data/{sample_name}_myeloid_raw_counts.txt.gz")), sep = '\t')
counts_mat <- as.matrix(counts)

# Create seurat object 
seurat_obj_list[[sample_name]] <- CreateSeuratObject(counts = counts_mat, project = sample_name)

# Metadata
# 2 conditions. TNFDARE= Inflamed, WT= Wild type
metadata_GSE255350 <- read.csv(gzfile(glue("00_data/{sample_name}_metadata.txt.gz")), sep = '\t')

# Check sample names of myeloid cells in meta data
metadata_GSE255350_myeloid <- metadata_GSE255350 %>% filter(cell_type == "myeloid")
myeloid_sample_names <- metadata_GSE255350_myeloid$Cell_id
myeloid_sample_names %>% length() # 8492
myeloid_sample_names %>% unique() %>% length() # 8021

# Check sample names in seurat object
seurat_obj_list[[sample_name]]@meta.data %>% nrow() # 8492
seurat_obj_list[[sample_name]]@meta.data %>% rownames() %>% unique() %>% length() # 8492

# In metadata, some ID's are used in both wt and dare (sample)
metadata_GSE255350_myeloid %>%
  group_by(Cell_id) %>%
  filter(n() > 1) %>% 
  arrange(Cell_id)

# Filter for wild-type (WT) sample (wt-12w)
metadata_GSE255350_myeloid_wt <- metadata_GSE255350_myeloid %>% 
  filter(sample == "wt-12w")

# Check sample names of myeloid, wt cells
myeloid_wt_sample_names <- metadata_GSE255350_myeloid_wt$Cell_id
myeloid_wt_sample_names %>% length() # 2228
myeloid_wt_sample_names %>% unique() %>% length() # 2228

# Investigate sample names in seurat object
seurat_obj_list[[sample_name]]@meta.data %>% rownames() %>% str_split_i('\\.', 2) %>% table() # 1_1: 2228,  1_2: 6264

# How numbers match:

# 1_1: 2228 --> myeloid, wt-12w cells
metadata_GSE255350_myeloid_wt %>% nrow() # 2228

# 1_2: 6264 --> myeloid, dare-12w cells
metadata_GSE255350_myeloid %>% filter(sample == "dare-12w") %>% nrow() # 6264

# Filter seurat object for the cell with "1_1" suffix
mask_cells_to_keep <- str_split_i(colnames(seurat_obj_list[[sample_name]]), "\\.", 2) == "1_1"
cells_to_keep <- colnames(seurat_obj_list[[sample_name]])[mask_cells_to_keep]
length(cells_to_keep)

# Subset the Seurat object
seurat_obj_list[[sample_name]] <- subset(seurat_obj_list[[sample_name]], cells = cells_to_keep)

# Rename cells in metadata to match metadata
metadata_GSE255350_myeloid_wt$Cell_id <- metadata_GSE255350_myeloid_wt$Cell_id %>% str_replace_all("-1", ".1_1")

# Check that cell names match 
table(metadata_GSE255350_myeloid_wt$Cell_id == colnames(seurat_obj_list[[sample_name]])) # TRUE

# Add metadata to seurat object 
metadata_GSE255350_myeloid_wt <- metadata_GSE255350_myeloid_wt %>% column_to_rownames("Cell_id")
seurat_obj_list[[sample_name]] <- AddMetaData(seurat_obj_list[[sample_name]], metadata = metadata_GSE255350_myeloid_wt)

seurat_obj_list[[sample_name]]@meta.data %>% head()

############################################ Data set 5. Fiona Powrie ############################################
# # 2025-10-15, email from Faidra: "The Fiona Powrie dataset will have to wait a bit since we might not be able to use it".
# # Myeloid. Contain multiple site such as Lamina Propria, MLN and others. 
# # Two conditions: Helicobacter infection = Inflammation, and WT
# 
# sample_name <- "CRAM1"
# 
# # Read data
# # counts <- ReadMtx(
# #   mtx = glue("00_data/{sample_name}/raw_feature_bc_matrix/matrix.mtx.gz"),
# #   features = glue("00_data/{sample_name}/raw_feature_bc_matrix/features.tsv.gz"),
# #   cells = glue("00_data/{sample_name}/raw_feature_bc_matrix/barcodes.tsv.gz")
# # )
# counts <- ReadMtx(
#   mtx = glue("00_data/{sample_name}_filtered/filtered_feature_bc_matrix/matrix.mtx.gz"),
#   features = glue("00_data/{sample_name}_filtered/filtered_feature_bc_matrix/features.tsv.gz"),
#   cells = glue("00_data/{sample_name}_filtered/filtered_feature_bc_matrix/barcodes.tsv.gz")
# )
# 
# seurat_obj_list[[sample_name]] <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
# 
# # Create seurat object 
# seurat_obj_list[[sample_name]] <- CreateSeuratObject(counts = counts, project = sample_name)
# rm(counts, sample_name)
# 
# sample_name <- "CRAM2"
# 
# # Read data
# # counts <- ReadMtx(
# #   mtx = glue("00_data/{sample_name}/raw_feature_bc_matrix/matrix.mtx.gz"),
# #   features = glue("00_data/{sample_name}/raw_feature_bc_matrix/features.tsv.gz"),
# #   cells = glue("00_data/{sample_name}/raw_feature_bc_matrix/barcodes.tsv.gz")
# # )
# counts <- ReadMtx(
#   mtx = glue("00_data/{sample_name}_filtered/filtered_feature_bc_matrix/matrix.mtx.gz"),
#   features = glue("00_data/{sample_name}_filtered/filtered_feature_bc_matrix/features.tsv.gz"),
#   cells = glue("00_data/{sample_name}_filtered/filtered_feature_bc_matrix/barcodes.tsv.gz")
# )
# 
# seurat_obj_list[[sample_name]] <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
# 
# # Create seurat object 
# seurat_obj_list[[sample_name]] <- CreateSeuratObject(counts = counts, project = sample_name)
# rm(counts, sample_name)

############################################ Data set 6. Rivera CA ############################################

# Sample: GSM5678427 CD103+CD11b+ dendritic cells from lamina propria

sample_name <- "GSM5678427"

# Read data
counts <- Read10X_h5(glue("00_data/{sample_name}_LP_raw_gene_bc_matrices_h5.h5"))

# Create seurat object 
seurat_obj_list[[sample_name]] <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)
rm(counts, sample_name)

############################################ Data set 7. ############################################

# Sample: 

sample_name <- "E_MTAB_9522"

# Read data
counts <- ReadMtx(
  mtx = glue("00_data/{sample_name}/filtered_feature_bc_matrix/matrix.mtx.gz"),
  features = glue("00_data/{sample_name}/filtered_feature_bc_matrix/features.tsv.gz"),
  cells = glue("00_data/{sample_name}/filtered_feature_bc_matrix/barcodes.tsv.gz")
)

seurat_obj_list[[sample_name]] <- CreateSeuratObject(counts = counts, project = sample_name, min.cells = 3, min.features = 200)

rm(counts, sample_name)

########################################## Export list of Seurat objects ##########################################

names(seurat_obj_list)

saveRDS(seurat_obj_list, "01_SILP_make_seurat/out/SILP_seurat_obj_list.rds")
