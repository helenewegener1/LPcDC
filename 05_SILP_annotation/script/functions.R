
explore_annotation_plot <- function(seurat_obj, dimplot, marker_genes, description = "", subfolder) {

  # Color UMAP by marker genes
  fp <- FeaturePlot(seurat_obj, features = marker_genes, reduction = umap_reduction.name) + plot_annotation(subtitle = description)
  
  # Collapse list of marker genes to name (a single string)
  marker_genes_name <- paste(marker_genes, collapse = '_')
  
  # Save UMAP
  n_marker_genes <- length(marker_genes)
  
  if (n_marker_genes == 1){
    
    ggsave(glue("05_SILP_annotation/plot/{subfolder}/UMAP_{marker_genes_name}.pdf"), plot = fp, width = 7, height = 6)
    
  } else if (n_marker_genes == 2){
    
    ggsave(glue("05_SILP_annotation/plot/{subfolder}/UMAP_{marker_genes_name}.pdf"), plot = fp, width = 12, height = 6)
    
  } else if (n_marker_genes == 3){
    
    ggsave(glue("05_SILP_annotation/plot/{subfolder}/UMAP_{marker_genes_name}.pdf"), plot = fp, width = 12, height = 12)
    
  } else {
    
    ggsave(glue("05_SILP_annotation/plot/{subfolder}/UMAP_{marker_genes_name}.pdf"), plot = fp)
    
  }
  
  
  # Wrap fp with UMAP colored by clusters and save 
  if (n_marker_genes == 1){
    
    wp <- dimplot + fp + plot_annotation(subtitle = description)
    ggsave(glue("05_SILP_annotation/plot/{subfolder}/UMAP_clusters_{marker_genes_name}.pdf"), plot = wp, width = 12, height = 6)
    
  } else if (n_marker_genes == 2){
    
    wp <- (dimplot + ggplot()) / fp + plot_annotation(subtitle = description)
    ggsave(glue("05_SILP_annotation/plot/{subfolder}/UMAP_clusters_{marker_genes_name}.pdf"), plot = wp, width = 14, height = 12)
    
  } else if (n_marker_genes == 3){
    
    wp <- dimplot + fp + plot_annotation(subtitle = description)
    ggsave(glue("05_SILP_annotation/plot/{subfolder}/UMAP_clusters_{marker_genes_name}.pdf"), plot = wp, width = 14, height = 7)
    
  } else {
    
    wp <- dimplot + fp + plot_annotation(subtitle = description)
    ggsave(glue("05_SILP_annotation/plot/{subfolder}/UMAP_clusters_{marker_genes_name}.pdf"), plot = wp) 
    
  }
  
  
}

# Basic function to convert human to mouse gene names
# https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
convertHumanGeneList <- function(x){
  
  # https://github.com/Huber-group-EMBL/biomaRt/issues/61
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  # print(head(humanx))
  return(humanx)
}

