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

############################################ Data set 4. George Kollias ############################################

# GSE255350_myeloid_raw_counts.txt.gz
GSE255350
sample_name <- "GSE255350"

# Read data
counts <- read.csv(gzfile(glue("00_data/{sample_name}_myeloid_raw_counts.txt.gz")), sep = '\t')
counts_mat <- as.matrix(counts)

# Create seurat object 
seurat_obj <- CreateSeuratObject(counts = counts_mat, project = sample_name)

# Metadata
# 2 conditions. TNFDARE = Inflamed, WT = Wild type
metadata_GSE255350 <- read.csv(gzfile(glue("00_data/{sample_name}_metadata.txt.gz")), sep = '\t')

# Check sample names of myeloid cells in meta data
metadata_GSE255350_myeloid <- metadata_GSE255350 %>% filter(cell_type == "myeloid")
myeloid_sample_names <- metadata_GSE255350_myeloid$Cell_id
myeloid_sample_names %>% length() # 8492
myeloid_sample_names %>% unique() %>% length() # 8021

# Check sample names in seurat object
seurat_obj@meta.data %>% nrow() # 8492
seurat_obj@meta.data %>% rownames() %>% unique() %>% length() # 8492

# In metadata, some ID's are used in both wt and dare (sample)
metadata_GSE255350_myeloid %>%
  group_by(Cell_id) %>%
  filter(n() > 1) %>% 
  arrange(Cell_id)

# Check sample names of myeloid, wt cells
myeloid_wt_sample_names <- metadata_GSE255350_myeloid_wt$Cell_id
myeloid_wt_sample_names %>% length() # 2228
myeloid_wt_sample_names %>% unique() %>% length() # 2228

# Investigate sample names in seurat object
seurat_obj@meta.data %>% rownames() %>% str_split_i('\\.', 2) %>% table() # 1_1: 2228,  1_2: 6264

# How numbers match with metadata:

# 1_1: 2228 --> myeloid, wt-12w cells
metadata_GSE255350_myeloid %>% filter(sample == "wt-12w") %>% nrow() # 2228

# 1_2: 6264 --> myeloid, dare-12w cells
metadata_GSE255350_myeloid %>% filter(sample == "dare-12w") %>% nrow() # 6264

# Rename cells in metadata to match metadata
metadata_GSE255350_myeloid <- metadata_GSE255350_myeloid %>% 
  mutate(Cell_id = case_when(sample == "wt-12w" ~ Cell_id %>% str_split_i("-", i = 1) %>% paste0(".1_1"),
                             sample == "dare-12w" ~ Cell_id %>% str_split_i("-", i = 1) %>% paste0(".1_2"))
         )
                                 
# Check that cell names match 
table(metadata_GSE255350_myeloid$Cell_id == colnames(seurat_obj)) # TRUE

# Add metadata to seurat object 
seurat_obj <- AddMetaData(seurat_obj, metadata = metadata_GSE255350_myeloid)

seurat_obj@meta.data %>% head()


########################################## Export list of Seurat objects ##########################################

saveRDS(seurat_obj, "11_inflammation_make_seurat/out/INF_seurat_obj.rds")
