# setwd("~/Documents/projects/project_cDC/LPcDC/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)
library(DropletUtils)
library(scDblFinder)
library(glmGamPoi)
# devtools::install_github("constantAmateur/SoupX", ref='devel')
library(SoupX)
library(multtest)

# Load data
seurat_obj <- readRDS("21_mLN_make_seurat/out/mLN_seurat_obj.rds") 

# Get count matrix
counts <- seurat_obj@assays$RNA$counts 

# Run scDblFinder
sce <- scDblFinder(counts, dbr=0.1)

# Access doublets and make metadata
doublet_metadata <- data.frame(scDblFinder.class = sce$scDblFinder.class,
                               scDblFinder.score = sce$scDblFinder.score)

# Add doublet analysis to metadata
seurat_obj <- AddMetaData(seurat_obj, doublet_metadata)

# Number of singlet and doublet - Add to plot
result <- table(seurat_obj@meta.data$scDblFinder.class, useNA = "ifany")

# Seurat workflow so I can UMAP
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
#Doublet detection
seurat_obj <- RunPCA(seurat_obj)
ElbowPlot(seurat_obj) #to determine dimentions used for following steps in doublet detection. Adjust dims. 
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
DimPlot(seurat_obj)

# Plot
DimPlot(seurat_obj, reduction = 'umap', group.by = "scDblFinder.class") + 
  labs(title = "scDblFinder", subtitle = glue("N doublets: {result[[2]]}, N singlets: {result[[1]]}"))
ggsave("22_mLN_QC/plot/scDblFinder.pdf", width = 7, height = 6)

# FeaturePlot(seurat_obj, reduction = 'umap', features = "MKI67") + 
#   labs(title = "MKI67", subtitle = "MKI67 is a proliferation marker")
# ggsave(glue("02_SILP_QC/plot/scDblFinder/MKI67_{sample_name}.pdf"), width = 7, height = 6)

# We do not do this here, but maybe later if it makes sense: Subset the Seurat object to include only "Singlet" cells
# seurat_obj_finalQC <- seurat_obj[, seurat_obj@meta.data[[DF.classification]] == "Singlet"]
# seurat_obj_finalQC_list[[sample_name]] <- seurat_obj_finalQC

seurat_obj_QC <- seurat_obj


#################### Export list of Seurat objects with QC metrices in metadata #################### 

saveRDS(seurat_obj_QC, "22_mLN_QC/out/mLN_seurat_obj_QC_metrics.rds")
