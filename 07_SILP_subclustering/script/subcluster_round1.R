# setwd("~/Documents/projects/project_cDC/LPcDC/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)

# Load data
seurat_object <- readRDS("04_SILP_integration/out/SILP_seurat_integrated_v5_RNA.rds")

# Settings
Reductions(seurat_object)

reduction <- "RNA_umap.cca"
cluster.name <- "RNA_cca_clusters_res.1"

# Final UMAP
DimPlot(seurat_object, reduction = reduction, group.by = cluster.name, label = TRUE) +
  labs(title = glue("UMAP - post integration"),
       subtitle = cluster.name)

# Clusters to remove 
target_cluster <- c(9, 24, 21, 22, 19) # CCA integration, res = 1, decided in meeting with Faidra

# Create a new logical column
seurat_object$keep_bool <- ifelse(seurat_object$RNA_cca_clusters_res.1 %in% target_cluster, FALSE, TRUE)
table(seurat_object$keep_bool)

# Subset cluster
seurat_object_clean <- subset(seurat_object, subset = keep_bool == TRUE)

# Check layers (already split)
# Layers(seurat_object_clean[["RNA"]])
# seurat_object_clean[["RNA"]] <- split(seurat_object_clean[["RNA"]], f = seurat_object_clean$orig.ident)

# Re-seurat workflow and integrate 
seurat_object_clean <- NormalizeData(seurat_object_clean)
seurat_object_clean <- FindVariableFeatures(seurat_object_clean)
seurat_object_clean <- ScaleData(seurat_object_clean)
# seurat_object_clean <- SCTransform(seurat_object_clean)
seurat_object_clean <- RunPCA(seurat_object_clean)
ElbowPlot(seurat_object_clean)

# Integrate
methods <- c("cca", "harmony")

# res <- 3
n_dim <- 20
# method <- "cca"

for (method in methods){
    
    # method <- "harmony"
    
    new.reduction <- glue("RNA_{method}_clean") %>% as.character()
    umap_reduction <- glue("RNA_umap.{method}_clean") %>% as.character()
    
    if (method == "cca"){
      
      seurat_object_clean <- IntegrateLayers(
        object = seurat_object_clean, 
        method = CCAIntegration,
        orig.reduction = "pca", 
        new.reduction = new.reduction,
        assay = "RNA", 
        verbose = FALSE
      )
      
    } else if (method == "harmony"){
      
      seurat_object_clean <- IntegrateLayers(
        object = seurat_object_clean, 
        method = HarmonyIntegration,
        orig.reduction = "pca", 
        new.reduction = new.reduction,
        assay = "RNA", 
        verbose = FALSE
      )
      
    }
    
    # Reductions(seurat_object_clean)
    
    for (res in c(1, 1.1, 1.2, 1.3, 1.4, 1.5)){
      
      cluster.name <- glue("RNA_{method}_clean_clusters_res.{res}") %>% as.character()
      
      DefaultAssay(seurat_object_clean)
      # ElbowPlot(seurat_object_clean, ndims = 30)
      seurat_object_clean <- FindNeighbors(seurat_object_clean, reduction = new.reduction, dims = 1:n_dim)
      seurat_object_clean <- FindClusters(seurat_object_clean, reduction = new.reduction, resolution = res, cluster.name = cluster.name)
      seurat_object_clean <- RunUMAP(seurat_object_clean, reduction = new.reduction, dims = 1:n_dim, reduction.name = umap_reduction)
      
      # Cluster plot
      DimPlot(seurat_object_clean, reduction = umap_reduction, group.by = cluster.name, label = TRUE, pt.size = 1.5) +
        labs(title = glue("UMAP - {new.reduction}"),
             subtitle = cluster.name)
      ggsave(glue("07_SILP_subclustering/plot/UMAP_clean_{method}_res.{res}.pdf"), width = 8, height = 7)
      
      # # Cluster plot by orig.ident
      # DimPlot(seurat_object_clean, reduction = umap_reduction, group.by = cluster.name, label = TRUE, split.by = "orig.ident") +
      #   NoLegend() + 
      #   labs(title = glue("UMAP - {new.reduction}"),
      #        subtitle = cluster.name)
      # ggsave(glue("07_SILP_subclustering/plot/UMAP_clean_{method}_res.{res}_split.pdf"), width = 24, height = 7)
      # 
    }
    
      # Feature plots
      FeaturePlot(seurat_object_clean, reduction = umap_reduction, features = c("Rorc", "Prdm16"), order = TRUE, pt.size = 1) +
        labs(subtitle = glue("UMAP - {new.reduction}"))
      ggsave(glue("07_SILP_subclustering/plot/UMAP_clean_{method}_Rorc_Prdm16.pdf"), width = 14, height = 7)
      
    
}

# Export object 
saveRDS(seurat_object_clean, "07_SILP_subclustering/out/seurat_object_clean.rds")

# seurat_object_clean <- readRDS("07_SILP_subclustering/out/seurat_object_clean.rds")

##################################### DEGs #####################################

res <- c(1, 1.1, 1.2, 1.3, 1.4, 1.5)
method <- c("harmony", "cca")

# Create all combinations
combos <- expand.grid(method = method, res = res)

# Generate cluster names
cluster.names <- glue("RNA_{combos$method}_clean_clusters_res.{combos$res}") %>% as.character()

for (cluster.name in cluster.names){
  
  # Join layers before DE analysis 
  seurat_object_clean[["RNA"]] <- JoinLayers(seurat_object_clean[["RNA"]])
  Layers(seurat_object_clean[["RNA"]])
  
  # Find all markers 
  all_markers <- FindAllMarkers(seurat_object_clean, 
                                group.by = cluster.name)
  
  # # Ensure the gene names are in a column (if they were in row names)
  df_all_markers <- all_markers #%>%
  # rownames_to_column("gene")
  
  # Process markers
  top_markers_list <- df_all_markers %>%
    filter(p_val_adj < 0.05) %>% # Filter for significant markers
    group_by(cluster) %>%
    # Sort by adjusted p-value (most significant) and then avg_log2FC (highest expression)
    arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 100) %>%
    ungroup() %>%
    # The split() function creates the named list needed for openxlsx
    split(., .$cluster)
  
  # Export xlsx file 
  out_file <- glue("07_SILP_subclustering/out/SILP_clean_Top_100_DE_Markers_{cluster.name}.xlsx")
  
  # Use openxlsx::write.xlsx, which takes the named list and writes
  # each element as a separate sheet (sheet name = list name, i.e., Cluster ID)
  openxlsx::write.xlsx(
    x = top_markers_list,
    file = out_file,
    overwrite = TRUE # Overwrite the file if it already exists
  )
  
}




