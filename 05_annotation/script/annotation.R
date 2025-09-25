setwd("~/Documents/projects/project_cDC/LPcDC/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)
library(patchwork)

# Load data
seurat_integrated <- readRDS("04_integration/out/seurat_integrated.rds")

# CD45+ -> All immune cells
grep("CD45", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # Can't find
grep("Ptprc", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # chatGPT alternative to 

FeaturePlot(seurat_integrated, features = "Ptprc") + plot_annotation(subtitle = "CD45+ -> All immune cells")
ggsave("05_annotation/plot/UMAP_Ptprc.pdf", width = 7, height = 6)

# MHCII -> all antigen presenting cells
grep("H2-A", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE)

FeaturePlot(seurat_integrated, features = c("H2-Ab1", "H2-Aa")) + plot_annotation(subtitle = "MHCII -> all antigen presenting cells")
ggsave("05_annotation/plot/UMAP_H2.Ab1_H2.Aa.pdf", width = 12, height = 6)

# XCR1 -> cDC1
grep("XCR1", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE)

FeaturePlot(seurat_integrated, features = "Xcr1") + plot_annotation(subtitle = "XCR1 -> cDC1")
ggsave("05_annotation/plot/UMAP_Xcr1.pdf", width = 7, height = 6)

# SIRPa, CD11b -> Macrophages and cDC2
grep("SIRPa", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE)

grep("CD11b", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # Can't find
grep("Itgam", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # chatGPT alternative to 

FeaturePlot(seurat_integrated, features = c("Sirpa", "Itgam")) + plot_annotation(subtitle = "SIRPa, CD11b -> Macrophages and cDC2")
ggsave("05_annotation/plot/UMAP_Sirpa_Itgam.pdf", width = 12, height = 6)

# CD11c -> Macrophages + DCs
grep("CD11c", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # Can't find
grep("Itgax", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # chatGPT alternative to 

FeaturePlot(seurat_integrated, features = "Itgax") + plot_annotation(subtitle = "CD11c -> Macrophages + DCs")
ggsave("05_annotation/plot/UMAP_Itgax.pdf", width = 7, height = 6)

# CD64, F4/80 -> Macrophages
grep("CD64", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # Can't find
grep("Fcgr1", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # chatGPT alternative 

grep("F4/80", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # Can't find
grep("Adgre1", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # chatGPT alternative 

FeaturePlot(seurat_integrated, features = c("Fcgr1", "Adgre1")) + plot_annotation(subtitle = "CD64, F4/80 -> Macrophages")
ggsave("05_annotation/plot/UMAP_Fcgr1_Adgre1.pdf", width = 12, height = 6)

# CD19, B220 -> B cells
grep("CD19", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) 

grep("B220", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # Can't find
grep("Ptprc", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # chatGPT alternative 

FeaturePlot(seurat_integrated, features = c("Cd19", "Ptprc")) + plot_annotation(subtitle = "CD19, B220 -> B cells")
ggsave("05_annotation/plot/UMAP_Cd19_Ptprc.pdf", width = 12, height = 6)

# CD3, TCRb -> T cells
grep("CD3", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) 

grep("TCRb", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # Can't find

FeaturePlot(seurat_integrated, features = c("Cd3e", "Cd3d", "Cd3g")) + plot_annotation(subtitle = "CD3, TCRb -> T cells")
ggsave("05_annotation/plot/UMAP_Cd3e_Cd3d_Cd3g.pdf", width = 12, height = 12)

# Ly6G -> Neutrophils
grep("Ly6G", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) 

FeaturePlot(seurat_integrated, features = "Ly6g") + plot_annotation(subtitle = "Ly6G -> Neutrophils")
ggsave("05_annotation/plot/UMAP_Ly6g.pdf", width = 7, height = 6)

# NK1.1 -> NK cells
grep("NK1", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # Can't find
grep("Klrb1c", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # chatGPT alternative 

FeaturePlot(seurat_integrated, features = "Klrb1c") + plot_annotation(subtitle = "NK1.1 -> NK cells")
ggsave("05_annotation/plot/UMAP_Klrb1c.pdf", width = 7, height = 6)

# SiglecF -> Eosinophils
grep("SiglecF", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE)

FeaturePlot(seurat_integrated, features = "Siglecf") + plot_annotation(subtitle = "SiglecF -> Eosinophils")
ggsave("05_annotation/plot/UMAP_Siglecf.pdf", width = 7, height = 6)

# Markers for subsets in the intesine specifically: cDC2: CD103+ CD11b+ or CD103- CD11b+, cDC1: CD103+ CD11b-
grep("Itgae", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # CD103
grep("Itgam", rownames(seurat_integrated), value = TRUE, ignore.case = TRUE) # CD11b

FeaturePlot(seurat_integrated, features = c("Itgae", "Itgam")) + plot_annotation(subtitle = "Markers for subsets in the intesine specifically: cDC2: CD103+ CD11b+ or CD103- CD11b+, cDC1: CD103+ CD11b-")
ggsave("05_annotation/plot/UMAP_Itgae_Itgam.pdf", width = 12, height = 6)





