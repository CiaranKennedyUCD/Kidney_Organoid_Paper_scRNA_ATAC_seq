library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(cicero)
library(SeuratWrappers)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(DoMultiBarHeatmap)
set.seed = 13357

#Enable parallelisation / multicore processing

plan("multiprocess", workers = 4)
options(future.globals.maxSize = 30000 * 1024^2) # for 30 Gb RAM

#Set Up directories for the analysis

project.name = "Kidney Organoids ATAC 2"
date = format(Sys.time(), "%d_%m_%Y")
date = gsub("2021","21",date)
sample.name = "combined"
project.dir = (file.path("/media/sf_G_DRIVE/Shared drives/Crean Group Data Share", paste(project.name,"_",date, sep="" )))
dir.create(project.dir, showWarnings = FALSE)
save.dir = (file.path(project.dir, sample.name))
dir.create(save.dir, showWarnings = FALSE)
dir.create(file.path(save.dir,"QC"), showWarnings = FALSE)
QC.dir = file.path(save.dir,"QC")
dir.create(file.path(save.dir,"Analysis"), showWarnings = FALSE)
Analysis.dir = file.path(save.dir,"Analysis")
dir.create(file.path(save.dir,"Clustering"), showWarnings = FALSE)
Clustering.dir = file.path(save.dir,"Clustering")
dir.create(file.path(Clustering.dir,"UMAPS"), showWarnings = FALSE)
UMAPS.dir = file.path(Clustering.dir,"UMAPS")

#Set up some nice custom color pallets


custom_colors <- list()

colors_dutch <- c(
  '#FFC312','#C4E538','#12CBC4','#FDA7DF','#ED4C67',
  '#F79F1F','#A3CB38','#1289A7','#D980FA','#B53471',
  '#EE5A24','#009432','#0652DD','#9980FA','#833471',
  '#EA2027','#006266','#1B1464','#5758BB','#6F1E51'
)

colors_spanish <- c(
  '#40407a','#706fd3','#f7f1e3','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62'
)

custom_colors$discrete <- c(colors_dutch, colors_spanish)

custom_colors$cell_cycle <- setNames(
  c('#45aaf2', '#f1c40f', '#e74c3c', '#7f8c8d'),
  c('G1',      'S',       'G2M',     '-')
)

#Extract gene annotation

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"


#Load in data and create seurat object

counts_control <- Read10X_h5(filename = "/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Kidney_Organoid_ATAC/CellRanger_Outs/Control/outs/filtered_peak_bc_matrix.h5")
metadata_control <- read.csv(
  file = "/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Kidney_Organoid_ATAC/CellRanger_Outs/Control/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

counts_controlTGFB <- Read10X_h5(filename = "/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Kidney_Organoid_ATAC/CellRanger_Outs/TGFB/outs/filtered_peak_bc_matrix.h5")
metadata_controlTGFB <- read.csv(
  file = "/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Kidney_Organoid_ATAC/CellRanger_Outs/TGFB/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

counts_GSK <- Read10X_h5(filename = "/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Kidney_Organoid_ATAC/CellRanger_Outs/GSK/outs/filtered_peak_bc_matrix.h5")
metadata_GSK <- read.csv(
  file = "/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Kidney_Organoid_ATAC/CellRanger_Outs/GSK/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

counts_GSKTGFB <- Read10X_h5(filename = "/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Kidney_Organoid_ATAC/CellRanger_Outs/GSKTGFB/outs/filtered_peak_bc_matrix.h5")
metadata_GSKTGFB <- read.csv(
  file = "/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Kidney_Organoid_ATAC/CellRanger_Outs/GSKTGFB/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)
##
chrom_assay_control <- CreateChromatinAssay(
  counts = counts_control,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = "/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Kidney_Organoid_ATAC/CellRanger_Outs/Control/outs/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200,
  annotation = annotation
)

chrom_assay_controlTGFB <- CreateChromatinAssay(
  counts = counts_controlTGFB,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = "/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Kidney_Organoid_ATAC/CellRanger_Outs/TGFB/outs/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200,
  annotation = annotation
)

chrom_assay_GSK <- CreateChromatinAssay(
  counts = counts_GSK,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = "/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Kidney_Organoid_ATAC/CellRanger_Outs/GSK/outs/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200,
  annotation = annotation
)

chrom_assay_GSKTGFB <- CreateChromatinAssay(
  counts = counts_GSKTGFB,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = "/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Kidney_Organoid_ATAC/CellRanger_Outs/GSKTGFB/outs/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200,
  annotation = annotation
  
)
##
control <- CreateSeuratObject(
  counts = chrom_assay_control,
  assay = "peaks",
  meta.data = metadata_control
)

controlTGFB <- CreateSeuratObject(
  counts = chrom_assay_controlTGFB,
  assay = "peaks",
  meta.data = metadata_controlTGFB
)

GSK <- CreateSeuratObject(
  counts = chrom_assay_GSK,
  assay = "peaks",
  meta.data = metadata_GSK
)

GSKTGFB <- CreateSeuratObject(
  counts = chrom_assay_GSKTGFB,
  assay = "peaks",
  meta.data = metadata_GSKTGFB
)
##

#Add Information to identify samples

control$dataset <- "ATAC"
controlTGFB$dataset <- "ATAC"
GSK$dataset <- "ATAC"
GSKTGFB$dataset <- "ATAC"

control@meta.data[["orig.ident"]] = "Control"
controlTGFB@meta.data[["orig.ident"]] = "ControlTGFB"
GSK@meta.data[["orig.ident"]] = "GSK"
GSKTGFB@meta.data[["orig.ident"]] = "GSKTGFB"

control@meta.data[["TGFB"]] = "NoTGFB"
controlTGFB@meta.data[["TGFB"]] = "TGFB"
GSK@meta.data[["TGFB"]] = "NoTGFB"
GSKTGFB@meta.data[["TGFB"]] = "TGFB"

control@meta.data[["GSK"]] = "NoGSK"
controlTGFB@meta.data[["GSK"]] = "NoGSK"
GSK@meta.data[["GSK"]] = "GSK"
GSKTGFB@meta.data[["GSK"]] = "GSK"

saveRDS(control, file = file.path(save.dir, "control_signac.rds"))
saveRDS(controlTGFB, file = file.path(save.dir, "controlTGFB_signac.rds"))
saveRDS(GSK, file = file.path(save.dir, "GSK_signac.rds"))
saveRDS(GSKTGFB, file = file.path(save.dir, "GSKTGFB_signac.rds"))

#Merge all datasets, adding a cell ID to make sure Cell names are unique

combined <- merge(x = control, 
                  y= c(controlTGFB, GSK, GSKTGFB),
                  add.cell.ids = c("control", "controlTGFB", "GSK", "GSKTGFB")
)

saveRDS(combined, file = file.path(save.dir, "combined_signac.rds"))

#Perform basic QC

# compute nucleosome signal score per cell
combined <- NucleosomeSignal(object = combined)

# compute TSS enrichment score per cell
combined <- TSSEnrichment(object = combined, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
combined$pct_reads_in_peaks <- combined$peak_region_fragments / combined$passed_filters * 100
combined$blacklist_ratio <- combined$blacklist_region_fragments / combined$peak_region_fragments

#plot TSS enrichment
combined$high.tss <- ifelse(combined$TSS.enrichment > 2, 'High', 'Low')
tss_plot <- TSSPlot(combined, group.by = 'high.tss') + theme_bw() + NoLegend()
tss_plot_sep <- TSSPlot(combined, group.by = 'orig.ident') + theme_bw() + NoLegend() 

tss_plots <- tss_plot + tss_plot_sep

ggsave(filename = paste(QC.dir,"TSS.svg", sep = ""), plot = tss_plots, width = 16, height = 8, units = "in", dpi = 600)
ggsave(filename = paste(QC.dir,"TSS.png", sep = ""), plot = tss_plots, width = 16, height = 8, units = "in", dpi = 600)

#Plot fragment size
combined$nucleosome_group <- ifelse(combined$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

fragment_plot <- FragmentHistogram(object = combined, group.by = 'nucleosome_group') + theme_bw()
fragment_plot_sep <- FragmentHistogram(object = combined, group.by = 'orig.ident') + theme_bw()

fragment_plots <- fragment_plot + fragment_plot_sep

ggsave(filename = paste(QC.dir,"fragments.svg", sep = ""), plot = fragment_plots, width = 16, height = 8, units = "in", dpi = 600)
ggsave(filename = paste(QC.dir,"fragments.png", sep = ""), plot = fragment_plots, width = 16, height = 8, units = "in", dpi = 600)

#Violin plots of various QC stats
combined_plots <- VlnPlot(
  object = combined,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
) + theme_bw()



p1 <- ggplot(FetchData(combined, vars = c("pct_reads_in_peaks", "orig.ident")), aes(x = orig.ident, y = pct_reads_in_peaks, fill = orig.ident)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(combined$orig.ident))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Percent Reads in Peaks', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p2 <- ggplot(FetchData(combined, vars = c("peak_region_fragments", "orig.ident")), aes(x = orig.ident, y = peak_region_fragments, fill = orig.ident)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(combined$orig.ident))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Peak Region Fragments', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p3 <- ggplot(FetchData(combined, vars = c("TSS.enrichment", "orig.ident")), aes(x = orig.ident, y = TSS.enrichment, fill = orig.ident)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(combined$orig.ident))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'TSS Enrichment', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p4 <- ggplot(FetchData(combined, vars = c("blacklist_ratio", "orig.ident")), aes(x = orig.ident, y = blacklist_ratio, fill = orig.ident)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(combined$orig.ident))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Blacklist Ratio', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p5 <- ggplot(FetchData(combined, vars = c("nucleosome_signal", "orig.ident")), aes(x = orig.ident, y = nucleosome_signal, fill = orig.ident)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(combined$orig.ident))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Nucleasome Signal', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

ggsave(
  file.path(QC.dir,"qc_histograms.png"),
  p1 + p3 + p5 +
    p2 + p4 + plot_layout(ncol = 3),
  height = 10, width = 15
)
ggsave(
  file.path(QC.dir,"qc_istograms.svg"),
  p1 + p3 + p5 +
    p2 + p4 + plot_layout(ncol = 3),
  height = 10, width = 15
)

saveRDS(combined, file = file.path(save.dir, "combined_signac_QC.rds"))

#Subset based on these stats
combined <- subset(
  x = combined,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

#Basic Processing, normalization and linear dimensional reduction
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = "lsi", seed.use = set.seed)

depth.cor <- DepthCor(combined, n = 20) + theme_bw()

ggsave(filename = paste(QC.dir,"depthcor.svg", sep = ""), plot = depth.cor, width = 16, height = 8, units = "in", dpi = 600)
ggsave(filename = paste(QC.dir,"depthcor.png", sep = ""), plot = depth.cor, width = 16, height = 8, units = "in", dpi = 600)

#Visualise to make sure everything went to plan 
DimPlot(combined, 
        group.by = "orig.ident",
        pt.size = 0.1)

CoveragePlot(
  object = combined,
  group.by = "predicted.id",
  region = "chr16-68735292-68837541",
  features = "CDH1"
)

saveRDS(combined, file = file.path(save.dir, "combined_signac_processed.rds"))


#Create a gene activity matrix
gene.activities <- GeneActivity(combined)
save(gene.activities, file = file.path(save.dir, "Gene_Activities.RData"))

# add the gene activity matrix to the Seurat object as a new assay and normalize it
combined[['RNA']] <- CreateAssayObject(counts = gene.activities)


saveRDS(combined, file = file.path(save.dir, "combined_signac_Gene_Matrix.rds"))

combined <- NormalizeData(
  object = combined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(combined$nCount_RNA)
)

saveRDS(combined, file = file.path(save.dir, "combined_signac_Gene_Matrix_Norm.rds"))

# Load the pre-processed scRNA-seq data for combined Kidney Organoids
combined_rna <- readRDS("/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Kidney_Organoids_2019_22_10_20/merged/Labelling/Merged Kidney Organoids 2019_cerebro_labelled.rds")
combined_rna$dataset <- "RNA"

DefaultAssay(combined) <- "RNA"

#Compute transfer anchors and transfer labels using harmony integration
transfer.anchors <- FindTransferAnchors(
  reference = combined_rna,
  query = combined,
  reduction = 'cca',
  normalization.method = "SCT"
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = combined_rna$seurat_clusters,
  weight.reduction = combined[['lsi']],
  dims = 2:50
)

combined<- AddMetaData(object = combined, metadata = predicted.labels)

hist(combined$prediction.score.max)
abline(v = 0.5, col = "red")

p1 <- DimPlot(combined, 
        group.by = "predicted.id",
        pt.size = 0.2,
        cols = custom_colors$discrete)

p2 <- DimPlot(combined, 
        group.by = "orig.ident",
        pt.size = 0.2,
        cols = custom_colors$discrete)

p1+p2

DefaultAssay(combined) <- "peaks"

saveRDS(combined, file = file.path(save.dir, "combined_signac_Labelled.rds"))


#Plotting 

CoveragePlot(
  object = combined,
  group.by = "predicted.id",
  region = "chr16-68735292-68835541",
  features = "CDH1",
)

#Find Co-accesable links

combined.cds <- as.cell_data_set(x = combined)
combined.cicero <- make_cicero_cds(combined.cds, reduced_coordinates = reducedDims(combined.cds)$UMAP)

# get the chromosome sizes from the Seurat object
genome <- seqlengths(annotation)

# use chromosome 1 to save some time
# omit this step to run on the whole genome
#genome <- genome[1]

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns <- run_cicero(combined.cicero, genomic_coords = genome.df, sample_num = 100)

#Convert pairwise co-assability links to co-accessability networks
ccans <- generate_ccans(conns)

#Convert connections to links and add to the seurat object

links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(combined) <- linked@assays[["peaks"]]@links

saveRDS(combined, file = file.path(save.dir, "combined_signac_Labelled_Links.rds"))

#Sample integration using Harmony

library(harmony)

hm.integrated <- RunHarmony(
  object = combined,
  group.by.vars = 'orig.ident',
  reduction = 'lsi',
  assay.use = 'peaks',
  project.dim = FALSE
)

# re-compute the UMAP using corrected LSI embeddings
hm.integrated <- RunUMAP(hm.integrated, dims = 2:30, reduction = 'harmony', seed.use = set.seed)

p1 <- DimPlot(hm.integrated, group.by = 'orig.ident', pt.size = 0.2, cols = custom_colors$discrete) + labs(title = NULL)
p2 <- DimPlot(hm.integrated, group.by = 'predicted.id', pt.size = 0.2, cols = custom_colors$discrete, label = TRUE, repel = TRUE) + labs(title = NULL) 
p3 <- DimPlot(hm.integrated, pt.size = 0.2, cols = custom_colors$discrete, label = TRUE, repel = TRUE) + labs(title = NULL) 

harmony_integrated <- wrap_plots(p1 + p2) + plot_annotation(title = "Harmony Integrated", theme = theme(plot.title = element_text(hjust = 0.5)))

p4 <- DimPlot(combined, group.by = 'orig.ident', pt.size = 0.2, cols = custom_colors$discrete)+ labs(title = NULL)
p5 <- DimPlot(combined, group.by = 'predicted.id', pt.size = 0.2, cols = custom_colors$discrete, label = TRUE, repel = TRUE) + labs(title = NULL)
p6 <- DimPlot(combined, pt.size = 0.2, cols = custom_colors$discrete, label = TRUE, repel = TRUE) + labs(title = NULL)

non_integrated <- wrap_plots(p3 + p4) + plot_annotation(title = "Non-Integrated", theme = theme(plot.title = element_text(hjust = 0.5)))



harmony_clust <- wrap_plots(p3 + p6) + plot_annotation(title = "Harmony Clustered", theme = theme(plot.title = element_text(hjust = 0.5)))

ggsave(filename = "/media/sf_G_DRIVE/My drive/harmony_clust_0.3.svg", plot = harmony_clust, width = 16, height = 8, units = "in", dpi = 600)
ggsave(filename = "/media/sf_G_DRIVE/My drive/harmony_clust_0.3.png", plot = harmony_clust, width = 16, height = 8, units = "in", dpi = 600)

ggsave(filename = paste(save.dir,"/non_integrated.svg", sep = ""), plot = non_integrated, width = 16, height = 8, units = "in", dpi = 600)
ggsave(filename = paste(save.dir,"/non_integrated.png", sep = ""), plot = non_integrated, width = 16, height = 8, units = "in", dpi = 600)

ggsave(filename = paste(save.dir,"/harmony_integrated.svg", sep = ""), plot = harmony_integrated, width = 16, height = 8, units = "in", dpi = 600)
ggsave(filename = paste(save.dir,"/harmony_integrated.png", sep = ""), plot = harmony_integrated, width = 16, height = 8, units = "in", dpi = 600)

saveRDS(hm.integrated, file = file.path(save.dir, "combined_signac_Labelled_Links_Harmony.rds"))

combined <= hm.integrated

# table <- table(combined$predicted.id, combined@meta.data$orig.ident)
# write.table(table, file = paste(save.dir,"table.tsv"), sep = ",")

#Clustering

seurat <- FindNeighbors(
  object = seurat,
  reduction = 'harmony',
  dims = 2:50
)

seurat <- FindClusters(
  object = seurat,
  algorithm = 3,
  resolution = 0.3,
  verbose = FALSE,
  random.seed = set.seed
)

DimPlot(object = hm.integrated, label = TRUE)
DimPlot(object = hm.integrated, label = TRUE, repel = TRUE, group.by = "predicted.id")

##Find Motifs and add to Seurat object

#Download motif positions from JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = "Homo sapiens", all_versions = FALSE)
)

# add motif information
combined <- AddMotifs(
  object = combined,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

saveRDS(combined, file = file.path(save.dir, "combined_signac_Labelled_Links_Harmony_Motifs.rds"))

#Find TF foothprints

combined <- Footprint(
  object = combined,
  motif.name = c("SIX2"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

p2 <- PlotFootprint(combined, features = "SIX2", group.by = "predicted.id")
p3 <- PlotFootprint(combined, features = "SIX2")
p2 + plot_layout(ncol = 1)

## ---- Store Variables in Seurat Object, include = False--------------
combined@misc$experiment <- list(
  experiment_name = project.name,
  organism = 'hg',
  date_of_analysis = Sys.Date(),
  genome = "hg38"
)

combined@misc$filepaths <- list(
  project.dir = project.dir,
  save.dir = save.dir,
  QC.dir = QC.dir,
  Analysis.dir = Analysis.dir,
  Clustering.dir = Clustering.dir,
  UMAPS.dir = UMAPS.dir
  
)

combined@misc$parameters <- list(
  gene_nomenclature = 'gene_name',
  seed.use = set.seed
)


sink(file.path(save.dir,"sessionInfo.txt"))
sessionInfo()
sink()

