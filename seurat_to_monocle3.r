library(SeuratWrappers)
library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)

seurat <- readRDS("/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Brain_Organoids_Cocaine_Control_Only_15_07_2021/merged/merged_Cerebro.rds")


project.dir = gsub("//", "/",seurat@misc[["filepaths"]][["project.dir"]])
save.dir = gsub("//", "/",seurat@misc[["filepaths"]][["save.dir"]])
QC.dir = gsub("//", "/",seurat@misc[["filepaths"]][["QC.dir"]])
Analysis.dir = gsub("//", "/",seurat@misc[["filepaths"]][["Analysis.dir"]])
Clustering.dir = gsub("//", "/",seurat@misc[["filepaths"]][["Clustering.dir"]])
UMAPS.dir = gsub("//", "/",seurat@misc[["filepaths"]][["UMAPS.dir"]])

seed.use = seurat@misc[["parameters"]][["seed.use"]]


nPC <- 40
monocle.dir <- file.path(save.dir,"monocle")
dir.create(monocle.dir)

cds_from_seurat <- as.cell_data_set(seurat, assay = "spliced")

## Calculate size factors using built-in function in monocle3
cds_from_seurat <- estimate_size_factors(cds_from_seurat)

## Add gene names into CDS
cds_from_seurat@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(seurat[["spliced"]])

## Process Using Monocle3

cds_from_seurat <- cluster_cells(cds = cds_from_seurat, reduction_method = "UMAP")
cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = TRUE)

cds_from_seurat <- order_cells(cds_from_seurat)

## Can plot Pseudotime

p1 <- plot_cells(
  cds = cds_from_seurat,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)  + facet_wrap(~orig.ident, nrow = 2)

ggsave(filename = paste(monocle.dir,"Trajectory UMAP.svg", sep = "/"), plot = p1, width = 6, height = 5, units = "in", dpi = 600)
ggsave(filename = paste(monocle.dir,"Trajectory UMAP.png", sep = "/"), plot = p1, width = 6, height = 5, units = "in", dpi = 600)



## Subset the cells and do branch expression/pseudotime analysis

cds_from_seurat_pr_test_res = graph_test(cds_from_seurat, neighbor_graph="knn", cores=6)

pr_deg_ids = row.names(subset(cds_from_seurat_pr_test_res, q_value < 0.05))

gene_list_1 <- pr_deg_ids[1:16]

### Plotting individual genes on the map

plot_cells(cds_from_seurat, 
           genes = gene_list_1,
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE) #+ facet_wrap(~orig.ident, nrow = 2)

### Now we can collect the trajectory-variable genes into modules for co-reg elements

cds_from_seurat <- preprocess_cds(cds_from_seurat, num_dim = 100)
gene_module_df = find_gene_modules(cds_from_seurat[pr_deg_ids,], resolution=1e-7)

cell_group_df = tibble::tibble(cell=row.names(colData(cds_from_seurat)), cell_group=cds_from_seurat@colData@listData[["orig.ident"]])
agg_mat = aggregate_gene_expression(cds_from_seurat, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
p1 <- pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

ggsave(filename = paste(monocle.dir,"GeneModulesHM.svg", sep = "/"), plot = p1, width = 8, height = 12, units = "in", dpi = 600)
ggsave(filename = paste(monocle.dir,"GeneModulesHM.png", sep = "/"), plot = p1, width = 8, height = 12, units = "in", dpi = 600)

### Now we can plot these modules

plot_cells(cds_from_seurat,
           genes=gene_module_df %>% filter(module %in% c(1)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE) + facet_wrap(~orig.ident, nrow = 2)


# ---- Finding genes that vary over Pseudotime ----

cds_pr_test_res <- graph_test(cds_from_seurat, neighbor_graph="principal_graph", cores=6)
pr_deg_ids <- row.names(subset(cds_pr_test_res, q_value < 0.05))

gene_list_1 <- pr_deg_ids[1:16]

gene_module_df = find_gene_modules(cds_from_seurat[pr_deg_ids,], resolution=1e-7)

cell_group_df = tibble::tibble(cell=row.names(colData(cds_from_seurat)), cell_group=cds_from_seurat@colData@listData[["orig.ident"]])
agg_mat = aggregate_gene_expression(cds_from_seurat, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
p1 <- pheatmap::pheatmap(agg_mat,
                         scale="column", clustering_method="ward.D2")

p1 <- plot_genes_in_pseudotime(cds_from_seurat[rowData(cds_from_seurat)$gene_short_name %in% c("HES1", "HTR2C", "SOX2","HES4","HES5","HES6", "SOX11")],
                         min_expr=0.5) 

ggsave(filename = paste(monocle.dir,"GenesInPseudo.svg", sep = "/"), plot = p1, width = 8, height = 12, units = "in", dpi = 600)
ggsave(filename = paste(monocle.dir,"GenesInPseudo.png", sep = "/"), plot = p1, width = 8, height = 12, units = "in", dpi = 600)

+ facet_wrap(~orig.ident, nrow = 2)
