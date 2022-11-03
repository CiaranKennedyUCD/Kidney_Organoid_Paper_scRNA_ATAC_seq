seurat <- FindSubCluster(seurat,
                         cluster = "Stroma 2/3_0",
               algorithm = 3,
               resolution = 0.5,
               graph.name = "peaks_snn")


DimPlot(object = seurat, label = TRUE, repel = TRUE, group.by = "sub.cluster")
DimPlot(object = seurat, label = TRUE, repel = TRUE)
DimPlot(object = seurat, label = TRUE, repel = TRUE, group.by = "orig.ident")

DimPlot(object = subset(seurat, idents = "Stroma 2/3_0"), label = TRUE, repel = TRUE, group.by = "sub.cluster")
DimPlot(object = subset(seurat, idents = "Stroma 2/3_0"), label = TRUE, repel = TRUE, group.by = "orig.ident")

seurat <- RunUMAP(seurat, dims = 2:50, reduction = "harmony", reduction.name = "umap_harmony")

seurat <- RenameIdents(seurat, "Stroma 2/3_0_1" = "Stroma 3", "Stroma 2/3_0_2" = "Stroma 3", "Stroma 2/3_0_3" = "Stroma 2", "Stroma 2/3_0_0" = "Stroma 3")
Idents(seurat) <- seurat$sub.cluster

seurat$clusters <- Idents(seurat)

levels <- c("Stroma 1", "Stroma 2", "Stroma 3", "Stroma 4", "Podocyte", "Mixed Distal Endothelial", "Kidney Progenitor 1", "Kidney Progenitor 2", "Muscle Progenitor", "Glial Progenitor", "Neural Progenitor 1", "Neural Progenitor 2", "Neural Progenitor 3",)
