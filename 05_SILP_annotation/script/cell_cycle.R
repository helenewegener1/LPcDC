# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(glue)
library(ggplot2)
library(patchwork)
library(stringr)
library(biomaRt)

source("05_SILP_annotation/script/functions.R")

# assay <- "SCT"
assay <- "RNA"

int_method <- "harmony"
# int_method <- "cca"

# Load data
seurat_integrated <- readRDS(glue("04_SILP_integration/out/SILP_seurat_integrated_v5_{assay}.rds"))

DefaultAssay(seurat_integrated) <- assay

# N cells 
n_cells <- ncol(seurat_integrated)

Reductions(seurat_integrated)

# Define reductions for integration method
res <- 1
red <- c(glue("RNA_integrated.{int_method}"), glue("RNA_umap.{int_method}"), glue("RNA_{int_method}_clusters_res.{res}"))

reduction <- red[[1]]
umap_reduction.name <- red[[2]]
cluster.name <- red[[3]]

# Get cell cycle genes (human)
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes

# Convert to mouse: https://github.com/satijalab/seurat/issues/2493
m.s.genes <- convertHumanGeneList(s.genes)
m.g2m.genes <- convertHumanGeneList(g2m.genes)

# Look at mouse genes 
length(s.genes) - length(m.s.genes) # 2 genes are not translated
length(g2m.genes) - length(m.g2m.genes) # 1 gene is not translated 

m.s.genes %in% rownames(seurat_integrated) %>% table()
m.g2m.genes %in% rownames(seurat_integrated) %>% table()

# Calculate mouse cell cycling score (mouse)
seurat_integrated <- CellCycleScoring(seurat_integrated,
                                      s.features = m.s.genes,
                                      g2m.features = m.g2m.genes)


# Plot Cell cycle
features <- c("S.Score", "G2M.Score")
for (feature in features){
  # feature <- "S.Score"
  FeaturePlot(seurat_integrated, 
              reduction = umap_reduction.name,
              features = feature) + 
    labs(
      title = feature, 
      caption = glue("N cells: {n_cells}")
    )
  
  ggsave(glue("05_SILP_annotation/plot/cell_cycle/UMAP_normal_{feature}.pdf"), width = 8, height = 7)
  
}

DimPlot(seurat_integrated, 
        reduction = umap_reduction.name,
        group.by = "Phase") + 
  labs(
    title = "Cell cycle phase", 
    caption = glue("N cells: {n_cells}")
  )
ggsave("05_SILP_annotation/plot/cell_cycle/UMAP_normal_Phase.pdf", width = 8, height = 7)

DimPlot(seurat_integrated, reduction = umap_reduction.name, group.by = cluster.name, label = TRUE) +
  labs(title = glue("UMAP - {reduction}"),
       subtitle = glue("{cluster.name}_res.{res}"))
ggsave("05_SILP_annotation/plot/cell_cycle/UMAP_normal.pdf", width = 8, height = 7)

######################## Regressing out cell cycle genes #######################

seurat_obj <- seurat_integrated
rm(seurat_integrated)

n_dim <- 20
res <- 1
reduction <- "RNA_integrated.harmony_cell_cycle_out"
umap_reduction <- glue("UMAP_{reduction}")
cluster.name <- glue("res_{res}_{reduction}")

DefaultAssay(seurat_obj)
# RNA already split
# seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$orig.ident)

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c(m.s.genes, m.g2m.genes))
# seurat_obj <- SCTransform(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
ElbowPlot(seurat_obj)

# # Get standard deviations of each PC
# pca_sd <- seurat_obj[["pca"]]@stdev
# 
# # Calculate variance explained per PC
# pca_var <- pca_sd^2
# pca_var_percent <- 100 * pca_var / sum(pca_var)
# 
# # Create a data frame for plotting
# pc_df <- data.frame(
#   PC = 1:length(pca_var_percent),
#   VariancePercent = pca_var_percent
# )
# barplot(pc_df$VariancePercent,
#         names.arg = pc_df$PC,
#         xlab = "PC",
#         ylab = "Variance Explained (%)",
#         main = "Variance Explained by PCs",
#         col = "steelblue")

seurat_integrated <- IntegrateLayers(
  object = seurat_obj, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = reduction,
  assay = "RNA", 
  verbose = FALSE
)

Reductions(seurat_integrated)

DefaultAssay(seurat_integrated)

seurat_integrated <- FindNeighbors(seurat_integrated, reduction = reduction, dims = 1:n_dim)
seurat_integrated <- FindClusters(seurat_integrated, reduction = reduction, resolution = res, cluster.name = cluster.name)
seurat_integrated <- RunUMAP(seurat_integrated, reduction = reduction, dims = 1:n_dim, reduction.name = umap_reduction)

# Cluster plot
DimPlot(seurat_integrated, reduction = umap_reduction, group.by = cluster.name, label = TRUE) +
  labs(title = glue("UMAP - {reduction}"),
       subtitle = glue("{cluster.name}_res.{res}"))
ggsave("05_SILP_annotation/plot/cell_cycle/UMAP_regressed_out.pdf", width = 8, height = 7)

# Plot Cell cycle
features <- c("S.Score", "G2M.Score")
for (feature in features){
  # feature <- "S.Score"
  FeaturePlot(seurat_integrated, 
              reduction = umap_reduction,
              features = feature) + 
    labs(
      title = feature, 
      caption = glue("N cells: {n_cells}")
    )
  
  ggsave(glue("05_SILP_annotation/plot/cell_cycle/UMAP_regressed_out_{feature}.pdf"), width = 8, height = 7)
  
}

DimPlot(seurat_integrated,
        reduction = umap_reduction,
        group.by = "Phase") +
  labs(
    title = "Cell cycle phase",
    caption = glue("N cells: {n_cells}")
  )
ggsave("05_SILP_annotation/plot/cell_cycle/UMAP_regressed_out_Phase.pdf", width = 8, height = 7)

saveRDS(seurat_integrated, "05_SILP_annotation/out/seurat_integrated_cell_cycle_genes_regressed_out.rds")

rm(seurat_integrated)

############################ Remove cell cycle genes ###########################

# Reload data
seurat_obj <- readRDS(glue("04_SILP_integration/out/SILP_seurat_integrated_v5_{assay}.rds"))

DefaultAssay(seurat_obj)

n_dim <- 20
res <- 1
reduction <- "RNA_integrated.harmony_cell_cycle_rm"
umap_reduction <- glue("UMAP_{reduction}")
cluster.name <- glue("res_{res}_{reduction}")

# Remove cell cycle genes
genes_to_remove <- c(cc.genes$s.genes, cc.genes$g2m.genes)

seurat_obj[["RNA"]] <- subset(
  seurat_obj[["RNA"]],
  features = setdiff(rownames(seurat_obj), genes_to_remove)
)

# seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$orig.ident)
# seurat_obj

seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
# seurat_obj <- SCTransform(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
ElbowPlot(seurat_obj)

seurat_integrated <- IntegrateLayers(
  object = seurat_obj, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = reduction,
  assay = "RNA", 
  verbose = FALSE
)

Reductions(seurat_integrated)

DefaultAssay(seurat_integrated)
seurat_integrated <- FindNeighbors(seurat_integrated, reduction = reduction, dims = 1:n_dim)
seurat_integrated <- FindClusters(seurat_integrated, reduction = reduction, resolution = res, cluster.name = cluster.name)
seurat_integrated <- RunUMAP(seurat_integrated, reduction = reduction, dims = 1:n_dim, reduction.name = umap_reduction)

# Cluster plot
# res <- 0.8
DimPlot(seurat_integrated, reduction = umap_reduction, group.by = cluster.name, label = TRUE) +
  labs(title = glue("UMAP - {reduction}"),
       subtitle = glue("{cluster.name}_res.{res}"))

ggsave("05_SILP_annotation/plot/cell_cycle/UMAP_removed.pdf", width = 8, height = 7)

saveRDS(seurat_integrated, "05_SILP_annotation/out/seurat_integrated_cell_cycle_genes_removed.rds")













