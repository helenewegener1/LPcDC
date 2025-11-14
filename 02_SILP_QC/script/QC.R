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
seurat_obj_list <- readRDS("01_SILP_make_seurat/out/SILP_seurat_obj_list.rds") # cellranger filtered

# Clean up Sample GSM9122899
seurat_obj_list$GSM9122899@assays$RNA$counts <- seurat_obj_list$GSM9122899@assays$RNA$`counts.Gene Expression`
seurat_obj_list$GSM9122899@assays$RNA$`counts.Gene Expression` <- NULL
seurat_obj_list$GSM9122899@assays$RNA$`counts.Antibody Capture` <- NULL
Layers(seurat_obj_list$GSM9122899)


# # Investigate need for removal of empty droplets 
# # https://bioconductor.org/packages/release/bioc/vignettes/DropletUtils/inst/doc/DropletUtils.html
# for (sample_name in names(seurat_obj_list)){
#   
#   ######################## CHECK EMPTY DROPLETS IN RAW #########################
#   
#   # # Load raw 
#   seurat_obj_raw <- seurat_obj_list[[sample_name]]
#   # 
#   # # Transform to SingleCellExperiment object
#   sc_exp_raw <- as.SingleCellExperiment(seurat_obj_raw)
#   # 
#   # # Run emptyDrops
#   # set.seed(100)
#   # e.out <- emptyDrops(counts(sc_exp_raw))
#   # 
#   # # See ?emptyDrops for an explanation of why there are NA values.
#   # is.cell <- table(e.out$FDR <= 0.01, useNA = "ifany")
#   # 
#   # # Computing barcode ranks
#   br.out <- barcodeRanks(counts(sc_exp_raw))
# 
#   # Making Barcode Rank Plot.
#   pdf(glue("02_QC/plot/emptyDrops/emptyDrops_{sample_name}_raw.pdf"), width = 8, height = 6)
#   
#   plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total", 
#        main=glue("{sample_name} raw (prefilter) Barcode Rank Plot")) 
#        # sub = glue("Cells: {is.cell[['TRUE']]}, Not cells: {is.cell[['FALSE']]}, NAs: {is.cell[[3]]}"))
#   o <- order(br.out$rank)
#   lines(br.out$rank[o], br.out$fitted[o], col="red")
#   
#   abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
#   abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
#   legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
#          legend=c("knee", "inflection"))
#   
#   dev.off()
#   
#   # Counts not raw enough to SoupX
#   
# }

################################ DoubletFinder on cellranger filtered ################################ 

# Initialize final QC list 
# seurat_obj_DoubletFinder <- list()
# 
# for (sample_name in names(seurat_obj_list)){
# 
#   # Define sample
#   seurat_obj <- seurat_obj_list[[sample_name]]
#   
#   # Seurat workflow
#   seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
#   seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
#   seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
#   #Doublet detection
#   seurat_obj <- RunPCA(seurat_obj)
#   ElbowPlot(seurat_obj) #to determine dimentions used for following steps in doublet detection. Adjust dims. 
#   seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
#   seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
#   seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
#   DimPlot(seurat_obj)
#   
#   #define expected number of doublets (10x Genomics) based on the number of cells
#   # Get the number of cells in the Seurat object
#   num_cells <- ncol(seurat_obj)
#   print(num_cells)
#   
#   #identify the Pk parameter for each sample, using no ground truth strategy
#   ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
#   sweep.res.list <- paramSweep(seurat_obj, PCs = 1:20, sct = TRUE)
#   sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
#   bcmvn<- find.pK(sweep.stats)
#   #visualize plot to find Pk parameter (highest peak)
#   ggplot(bcmvn, aes(pK, BCmetric, group = 1)) +
#     geom_point()+
#     geom_line()
#   
#   #store max Pk as Pk variable 
#   pK <- bcmvn %>% 
#     filter(BCmetric == max(BCmetric)) %>%
#     select(pK)
#   pK <- as.numeric(as.character(pK[[1]]))
#   
#   #homotypic doublet proportions
#   annotations <- seurat_obj@meta.data$seurat_clusters
#   homotypic.prop <- modelHomotypic(annotations)
#   nExp_poi <- round(0.061*nrow(seurat_obj@meta.data))
#   nExp.poi.adj <- round(nExp_poi*(1-homotypic.prop))
#   
#   #Run Doublet finder
#   seurat_obj <- doubletFinder(seurat_obj,
#                               PCs = 1:20,
#                               pN = 0.25,
#                               pK = pK,
#                               nExp = nExp.poi.adj,
#                               # reuse.pANN = FALSE, 
#                               sct = TRUE)
#   
#   
#   #view metadata to see single vs. double in DF.classicication
#   # View(seurat_obj@meta.data)
#   #get the name of the coloumn, copy DF.classification name  
#   # names(seurat_obj@meta.data)
#   
#   # visualize
#   DF.classification <- colnames(seurat_obj@meta.data)[colnames(seurat_obj@meta.data) %>% str_starts('DF')]
#   
#   # Number of singlet and doublet - Add to plot
#   result <- table(seurat_obj@meta.data[DF.classification], useNA = "ifany")
#   
#   # Plot
#   DimPlot(seurat_obj, reduction = 'umap', group.by = DF.classification) + 
#     labs(title = "DoubletFinder", subtitle = glue("N doublets: {result[[1]]}, N singlets: {result[[2]]}"))
#   ggsave(glue("02_QC/plot/DoubletFinder/DoubletFinder_{sample_name}.pdf"), width = 7, height = 6)
#   
#   FeaturePlot(seurat_obj, reduction = 'umap', features = "Mki67") + 
#     labs(title = "MKI67", subtitle = "MKI67 is a proliferation marker")
#   ggsave(glue("02_QC/plot/DoubletFinder/Mki67_{sample_name}.pdf"), width = 7, height = 6)
#   
#   
#   # We do not do this here, but maybe later if it makes sense: Subset the Seurat object to include only "Singlet" cells
#   # seurat_obj_finalQC <- seurat_obj[, seurat_obj@meta.data[[DF.classification]] == "Singlet"]
#   # seurat_obj_finalQC_list[[sample_name]] <- seurat_obj_finalQC
#   
#   seurat_obj_DoubletFinder[[sample_name]] <- seurat_obj
#   
# }



# Initialize final QC list 
seurat_obj_QC <- list()

for (sample_name in names(seurat_obj_list)){
  
  # Define sample
  seurat_obj <- seurat_obj_list[[sample_name]]
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
  ggsave(glue("02_SILP_QC/plot/scDblFinder/scDblFinder_{sample_name}.pdf"), width = 7, height = 6)
  
  # FeaturePlot(seurat_obj, reduction = 'umap', features = "MKI67") + 
  #   labs(title = "MKI67", subtitle = "MKI67 is a proliferation marker")
  # ggsave(glue("02_SILP_QC/plot/scDblFinder/MKI67_{sample_name}.pdf"), width = 7, height = 6)
  
  # We do not do this here, but maybe later if it makes sense: Subset the Seurat object to include only "Singlet" cells
  # seurat_obj_finalQC <- seurat_obj[, seurat_obj@meta.data[[DF.classification]] == "Singlet"]
  # seurat_obj_finalQC_list[[sample_name]] <- seurat_obj_finalQC
  
  seurat_obj_QC[[sample_name]] <- seurat_obj
  
}

#################### Export list of Seurat objects with QC metrices in metadata #################### 

saveRDS(seurat_obj_QC, "02_SILP_QC/out/SILP_seurat_obj_QC_metrics.rds")
