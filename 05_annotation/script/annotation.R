setwd("~/Documents/projects/project_cDC/LPcDC/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)
library(patchwork)
library(stringr)

source("05_annotation/script/functions.R")

# Load data
seurat_integrated <- readRDS("04_integration/out/seurat_integrated.rds")

# Cluster plot
res <- 0.5
dp <- DimPlot(seurat_integrated, reduction = "umap", group.by = glue("SCT_snn_res.{res}"), label = TRUE) +
  labs(title = "UMAP - post harmony integration",
       subtitle = glue("SCT_snn_res.{res}"))

########################### CD45+ -> All immune cells ########################## 
grep("CD45", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # Can't find
grep("Ptprc", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # chatGPT alternative to 

# Define marker cell information 
marker_genes <- c("Ptprc")
description <- "CD45+ -> All immune cells"

# Make plots 
explore_annotation_plot(seurat_obj = seurat_integrated, dimplot = dp, marker_genes = marker_genes, description = description)

# Clean up
rm(marker_genes, description)

##################### MHCII -> all antigen presenting cells ##################### 
grep("H2-A", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE)

# Define marker cell information 
marker_genes <- c("H2-Ab1", "H2-Aa")
description <- "MHCII -> all antigen presenting cells"

# Make plots 
explore_annotation_plot(seurat_obj = seurat_integrated, dimplot = dp, marker_genes = marker_genes, description = description)

# Clean up
rm(marker_genes, description)

################################# XCR1 -> cDC1 ################################# 
grep("XCR1", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE)

# Define marker cell information 
marker_genes <- c("Xcr1")
description <- "XCR1 -> cDC1"

# Make plots 
explore_annotation_plot(seurat_obj = seurat_integrated, dimplot = dp, marker_genes = marker_genes, description = description)

# Clean up
rm(marker_genes, description)

#####################  SIRPa, CD11b -> Macrophages and cDC2 ##################### 
grep("SIRPa", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE)

grep("CD11b", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # Can't find
grep("Itgam", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # chatGPT alternative to 

# Define marker cell information 
marker_genes <- c("Sirpa", "Itgam")
description <- "SIRPa, CD11b -> Macrophages and cDC2"

# Make plots 
explore_annotation_plot(seurat_obj = seurat_integrated, dimplot = dp, marker_genes = marker_genes, description = description)

# Clean up
rm(marker_genes, description)

########################## CD11c -> Macrophages + DCs ########################## 
grep("CD11c", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # Can't find
grep("Itgax", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # chatGPT alternative to 

# Define marker cell information 
marker_genes <- c("Itgax")
description <- "CD11c -> Macrophages + DCs"

# Make plots 
explore_annotation_plot(seurat_obj = seurat_integrated, dimplot = dp, marker_genes = marker_genes, description = description)

# Clean up
rm(marker_genes, description)

########################## CD64, F4/80 -> Macrophages ########################## 
grep("CD64", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # Can't find
grep("Fcgr1", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # chatGPT alternative 

grep("F4/80", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # Can't find
grep("Adgre1", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # chatGPT alternative 

# Define marker cell information 
marker_genes <- c("Fcgr1", "Adgre1")
description <- "CD64, F4/80 -> Macrophages"

# Make plots 
explore_annotation_plot(seurat_obj = seurat_integrated, dimplot = dp, marker_genes = marker_genes, description = description)

# Clean up
rm(marker_genes, description)

############################# CD19, B220 -> B cells ############################ 
grep("CD19", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) 

grep("B220", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # Can't find
grep("Ptprc", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # chatGPT alternative

# Define marker cell information 
marker_genes <- c("Cd19", "Ptprc")
description <- "CD19, B220 -> B cells"

# Make plots 
explore_annotation_plot(seurat_obj = seurat_integrated, dimplot = dp, marker_genes = marker_genes, description = description)

# Clean up
rm(marker_genes, description)

############################# CD3, TCRb -> T cells ############################# 
grep("CD3", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) 

grep("TCRb", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # Can't find

# Define marker cell information 
marker_genes <- c("Cd3e", "Cd3d", "Cd3g")
description <- "CD3, TCRb -> T cells"

# Make plots 
explore_annotation_plot(seurat_obj = seurat_integrated, dimplot = dp, marker_genes = marker_genes, description = description)

# Clean up
rm(marker_genes, description)

############################# Ly6G -> Neutrophils ############################# 
grep("Ly6G", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) 

# Define marker cell information 
marker_genes <- c("Ly6g")
description <- "Ly6G -> Neutrophils"

# Make plots 
explore_annotation_plot(seurat_obj = seurat_integrated, dimplot = dp, marker_genes = marker_genes, description = description)

# Clean up
rm(marker_genes, description)

############################### NK1.1 -> NK cells ############################## 
grep("NK1", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # Can't find
grep("Klrb1c", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # chatGPT alternative 

# Define marker cell information 
marker_genes <- c("Klrb1c")
description <- "NK1.1 -> NK cells"

# Make plots 
explore_annotation_plot(seurat_obj = seurat_integrated, dimplot = dp, marker_genes = marker_genes, description = description)

# Clean up
rm(marker_genes, description)

############################ SiglecF -> Eosinophils ############################ 
grep("SiglecF", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE)

# Define marker cell information 
marker_genes <- c("Siglecf")
description <- "SiglecF -> Eosinophils"

# Make plots 
explore_annotation_plot(seurat_obj = seurat_integrated, dimplot = dp, marker_genes = marker_genes, description = description)

# Clean up
rm(marker_genes, description)

#################################################################################
# Markers for subsets in the intesine specifically: cDC2: CD103+ CD11b+ or CD103- CD11b+, cDC1: CD103+ CD11b-
grep("Itgae", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # CD103
grep("Itgam", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # CD11b

# Define marker cell information 
marker_genes <- c("Itgae", "Itgam")
description <- "Markers for subsets in the intesine specifically: cDC2: CD103+ CD11b+ or CD103- CD11b+, cDC1: CD103+ CD11b-"

# Make plots 
explore_annotation_plot(seurat_obj = seurat_integrated, dimplot = dp, marker_genes = marker_genes, description = description)

# Clean up
rm(marker_genes, description)




