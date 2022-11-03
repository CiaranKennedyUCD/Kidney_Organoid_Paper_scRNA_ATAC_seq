library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(ggplot2)

combined <- readRDS("/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Kidney Organoids ATAC 2_14_03_21/combined/combined_ChromVar.rds")

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

combined <- AddMotifs(
  object = combined,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

combined <- RunChromVAR(
  object = combined,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

saveRDS(combined, "/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Kidney Organoids ATAC 2_14_03_21/combined/combined_ChromVar.rds")

DefaultAssay(combined) <- "chromvar"

#Originals From Figure
SMAD5 <- FeaturePlot(combined, "MA1557.1", min.cutoff = 0.1, max.cutoff = 3)
SMAD3 <- FeaturePlot(combined, "MA0795.1", min.cutoff = 0, max.cutoff = 2.5)
FOSL1 <- FeaturePlot(combined, "MA0477.2", min.cutoff = 0, max.cutoff = 3)
ZNF384 <- FeaturePlot(combined, "MA1125.1", min.cutoff = 0, max.cutoff = 3)
JUND <- FeaturePlot(combined, "MA0491.2", min.cutoff = 0, max.cutoff = 3)
FOSL2JUN <- FeaturePlot(combined, "MA1130.1", min.cutoff = 0, max.cutoff = 2.5)

#Top 6 Enriched Stroma2_Control vs Stroma3_TGFB
HOXB13 <- FeaturePlot(combined, "MA0901.2", min.cutoff = 0, max.cutoff = 3)
ZNF384 <- FeaturePlot(combined, "MA1125.1", min.cutoff = 0, max.cutoff = 3)
ZBTB18 <- FeaturePlot(combined, "MA0698.1", min.cutoff = 0, max.cutoff = 8)
IRF1 <- FeaturePlot(combined, "MA0050.2", min.cutoff = 0, max.cutoff = 2.5)
MEF2A <- FeaturePlot(combined, "MA0052.4", min.cutoff = 0, max.cutoff = 3)
SPI1 <- FeaturePlot(combined, "MA0080.5", min.cutoff = 0, max.cutoff = 2.5)

#6 Relevant Significant Motifs Stroma2_Control vs Stroma3_TGFB
ZEB1 <- FeaturePlot(combined, "MA0103.3", min.cutoff = 0, max.cutoff = 4)
TWIST1 <- FeaturePlot(combined, "MA1123.2", min.cutoff = 0, max.cutoff = 9)
STAT3 <- FeaturePlot(combined, "MA0144.2", min.cutoff = 0, max.cutoff = 3)
SNAI2 <- FeaturePlot(combined, "MA0745.2", min.cutoff = 0, max.cutoff = 4)
SNAI1 <- FeaturePlot(combined, "MA1558.1", min.cutoff = 0, max.cutoff = 4.5)
SMAD2SMAD3SMAD4 <- FeaturePlot(combined, "MA0513.1", min.cutoff = 0, max.cutoff = 2.5)

ls <- objects()[c(6,7,8,9,10,12,22,23,24,25,26,27,28,29,31,32,33)]

for(i in 1:length(ls)){
  ggsave(filename = paste(save.dir,ls[i],"_Motif_UMAP.png", sep = ""), plot = get(ls[i]), width = 9, height = 8, units = "in", dpi = 600)
  ggsave(filename = paste(save.dir,ls[i],"_Motif_UMAP.svg", sep = ""), plot = get(ls[i]), width = 9, height = 8, units = "in", dpi = 600)
}


FOSL2 <- FeaturePlot(combined, "MA0478.1", min.cutoff = 0, max.cutoff = 2)

ggsave(filename = paste(save.dir,"EGR1_Motif_UMAP.png", sep = ""), plot = EGR1, width = 9, height = 8, units = "in", dpi = 600)
ggsave(filename = paste(save.dir,"EGR1_Motif_UMAP.svg", sep = ""), plot = EGR1, width = 9, height = 8, units = "in", dpi = 600)


DefaultAssay(combined) <- "peaks"

combined <- Footprint(
  object = combined,
  motif.name = "SMAD5",
  genome = BSgenome.Hsapiens.UCSC.hg38
)


ggsave(filename = paste(save.dir,"Footprints_UMAP_2.png", sep = ""), plot = p1 + plot_layout(ncol=1), width = 9, height = 20, units = "in", dpi = 600)
ggsave(filename = paste(save.dir,"Footprints_UMAP_2.svg", sep = ""), plot = p1 + plot_layout(ncol=1), width = 9, height = 20, units = "in", dpi = 600)

Idents(combined) <- combined@meta.data[["celltype.stim"]]

p1 <- PlotFootprint(combined, idents = c("Stroma 2_Control", "Stroma 3_ControlTGFB", "Stroma 3_GSKTGFB"), c("SMAD3", "SMAD5", "FOSL2", "JUNB", "STAT3", "IRF1", "EGR1"), label = FALSE)
