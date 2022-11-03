## ----include = FALSE-------------------------------------------------

library(dplyr)
library(Seurat)
library(devEMF)
library(future)
library(scDblFinder)
library(SeuratWrappers)
library(velocyto.R)
library(ggplot2)
library(patchwork)


## ---- Setup a number of variables and file directories ----
#Set number of cores to use and object size limit 

# plan("multisession", workers = 4) // Parrallisation seems to cause issues with a few of these steps so disable it for this script
options(future.globals.maxSize= 89128960000)

#Set Seed for visualisations
seed.use = 824

#List your samples
samples.list = c("control", "cocaine")

#Set up variables for the cell cycle assignments
s.genes <- cc.genes$s.genes 
g2m.genes <- cc.genes$g2m.genes

#Setup Sample name and output folders

project.name = "Brain_Organoids_Cocaine_Control_Only"
sample.name = "merged"
date = format(Sys.time(), "%d_%m_%Y")
date = gsub("2020","20",date)
project.dir = (file.path("/media/sf_G_DRIVE/Shared drives/Crean Group Data Share", paste(project.name,"_",date, sep="")))
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




## ----Setup Colour Pallet ----------------------------

###Setup a pallet of custom colours (credit to https://romanhaa.github.io/projects/scsplicedseq_workflow/ for the code)

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


## ---- loading in data and converting to seurat objects, as the data was pre-processed for RNA velocity analysis using Velocyto, cell names were edited to 
## ----- adhere to the seurat naming convention --------------------------

control.v <- ReadVelocity(file = "/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Brain_Organoid_Cellranger_Outputs_RNAseq/Brain/ControlBrainCount.loom")

control.v[["spliced"]]@Dimnames[[2]] <- gsub("ControlBrainCount:", "", control.v[["spliced"]]@Dimnames[[2]])
control.v[["spliced"]]@Dimnames[[2]] <- gsub("x", "-1", control.v[["spliced"]]@Dimnames[[2]])

control.v[["unspliced"]]@Dimnames[[2]] <- gsub("ControlBrainCount:", "", control.v[["unspliced"]]@Dimnames[[2]])
control.v[["unspliced"]]@Dimnames[[2]] <- gsub("x", "-1", control.v[["unspliced"]]@Dimnames[[2]])

control.v[["ambiguous"]]@Dimnames[[2]] <- gsub("ControlBrainCount:", "", control.v[["ambiguous"]]@Dimnames[[2]])
control.v[["ambiguous"]]@Dimnames[[2]] <- gsub("x", "-1", control.v[["ambiguous"]]@Dimnames[[2]])

control <- as.Seurat(control.v)

cocaine.v <- ReadVelocity(file = "/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Brain_Organoid_Cellranger_Outputs_RNAseq/Brain/CokeBrainCount.loom")

cocaine.v[["spliced"]]@Dimnames[[2]] <- gsub("CokeBrainCount:", "", cocaine.v[["spliced"]]@Dimnames[[2]])
cocaine.v[["spliced"]]@Dimnames[[2]] <- gsub("x", "-1", cocaine.v[["spliced"]]@Dimnames[[2]])

cocaine.v[["unspliced"]]@Dimnames[[2]] <- gsub("CokeBrainCount:", "", cocaine.v[["unspliced"]]@Dimnames[[2]])
cocaine.v[["unspliced"]]@Dimnames[[2]] <- gsub("x", "-1", cocaine.v[["unspliced"]]@Dimnames[[2]])

cocaine.v[["ambiguous"]]@Dimnames[[2]] <- gsub("CokeBrainCount:", "", cocaine.v[["ambiguous"]]@Dimnames[[2]])
cocaine.v[["ambiguous"]]@Dimnames[[2]] <- gsub("x", "-1", cocaine.v[["ambiguous"]]@Dimnames[[2]])

cocaine <- as.Seurat(cocaine.v)

##---- Assign an original identity variable for each sample for later identification ----

control@meta.data[["orig.ident"]] = "Control"
cocaine@meta.data[["orig.ident"]] = "Cocaine"

## ---- Initial merge of the processed objects for plotting some QC stats ----

sample.list <- c(control, cocaine)

seurat <- merge(x = control, y = cocaine, add.cell.ids = samples.list, project = "Cerebral Organoids Cocaine Control Only")


## ---- Assign Mitochondrial DNA %, include = FALSE--------------------

#Add colum to object matadata to store Mitochondrial percentage of each cell
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")

#View Mitochondrial % from first 5 cells
head(seurat@meta.data, 5)



## ---- Setup QC thresholds, include = FALSE---------------------------
median_nCount <- median(seurat$nCount_spliced)

sd_nCount <- sd(seurat$nCount_spliced)

median_nFeature <- median(seurat$nFeature_spliced)

sd_nFeature <- sd(seurat$nFeature_spliced)


median_percent_MT <- median(seurat@meta.data[["percent.mt"]])

sd_percent_MT <- sd(seurat@meta.data[["percent.mt"]])


thresholds_nCount <- c(0, median_nCount + 2*sd_nCount)

thresholds_nFeature <- c(0, median_nFeature + 2*sd_nFeature)

thresholds_percent_MT <- c(0, median_percent_MT + 2*sd_percent_MT)



## ---- Plot QC, results = FALSE, message = FALSE, warning= FALSE, error = FALSE, echo = FALSE, fig.cap = "QC metrics plotted for the organoids scRNAseq run. Black verticle line represents the median value across all samples. Red lines represent the filtering threshold (5 x Median Absolute Deviation)."----
  p1 <- ggplot(FetchData(seurat, vars = c("nCount_spliced", "orig.ident")), aes(x = orig.ident, y = nCount_spliced, fill = orig.ident)) +
    geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
    geom_hline(yintercept = median_nCount, color = 'black') +
    geom_hline(yintercept = thresholds_nCount, color = 'red') +
    theme_bw() +
    scale_fill_manual(values = custom_colors$discrete) +
    scale_x_discrete(limits = rev(levels(seurat$orig.ident))) +
    scale_y_continuous(labels = scales::comma) +
    labs(title = 'Number of transcripts', subtitle = 'linear scale') +
    theme(
      axis.title = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = 'none'
    ) +
    coord_flip()

p2 <- ggplot(FetchData(seurat, vars = c("nCount_spliced", "orig.ident")), aes(x = orig.ident, y = nCount_spliced, fill = orig.ident)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_hline(yintercept = median_nCount, color = 'black') +
  geom_hline(yintercept = thresholds_nCount, color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(seurat$orig.ident))) +
  scale_y_log10(labels = scales::comma) +
  labs(title = 'Number of transcripts', subtitle = 'log-scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()
p3 <- ggplot(FetchData(seurat, vars = c("nFeature_spliced", "orig.ident")), aes(x = orig.ident, y = nFeature_spliced, fill = orig.ident)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_hline(yintercept = median_nFeature, color = 'black') +
  geom_hline(yintercept = thresholds_nFeature, color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(seurat$orig.ident))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()
p4 <- ggplot(FetchData(seurat, vars = c("nFeature_spliced", "orig.ident")), aes(x = orig.ident, y = nFeature_spliced, fill = orig.ident)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_hline(yintercept = median_nFeature, color = 'black') +
  geom_hline(yintercept = thresholds_nFeature, color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(seurat$orig.ident))) +
  scale_y_log10(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'log-scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()
p5 <- ggplot(FetchData(seurat, vars = c("percent.mt", "orig.ident")), aes(x = orig.ident, y = percent.mt, fill = orig.ident)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  geom_hline(yintercept = median_percent_MT, color = 'black') +
  geom_hline(yintercept = thresholds_percent_MT, color = 'red') +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(seurat$orig.ident))) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = 'Percent MT transcripts', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()
ggsave(
  file.path(QC.dir,"qc_histogram_nCount_nFeature_percentMT_thresholds.png"),
  p1 + p3 + p5 +
  p2 + p4 + plot_layout(ncol = 3),
  height = 10, width = 15
)
ggsave(
  file.path(QC.dir,"qc_histogram_nCount_nFeature_percentMT_thresholds.svg"),
  p1 + p3 + p5 +
  p2 + p4 + plot_layout(ncol = 3),
  height = 10, width = 15
)

p1 + p3 + p5 +
p2 + p4


## ---- Plotting some initial QC figures, eval = FALSE, echo = FALSE----

VlnPlot(seurat, features = c("nFeature_spliced", "nCount_spliced", "percent.mt"), group.by = "orig.ident", assay = "spliced")
plot1 <- FeatureScatter(seurat, feature1 = "nCount_spliced", feature2 = "percent.mt", group.by = "orig.ident")
plot1
plot2 <- FeatureScatter(seurat, feature1 = "nCount_spliced", feature2 = "nFeature_spliced", group.by = "orig.ident") +
  geom_hline(yintercept = thresholds_nFeature, color = 'red') +
  geom_vline(xintercept = thresholds_nCount, color = 'red')
plot2

#PlotQC data
png(file = file.path(QC.dir,"Feature,Count,MT.png"), width=8, height = 4, units= "in", res= 600)
VlnPlot(seurat, features = c("nFeature_spliced", "nCount_spliced", "percent.mt"), group.by = "orig.ident", assay = "spliced")
dev.off()
svg(file = file.path(QC.dir,"Feature,Count,MT.svg"))
VlnPlot(seurat, features = c("nFeature_spliced", "nCount_spliced", "percent.mt"), group.by = "orig.ident", assay = "spliced")
dev.off()


#Create scatter plots of QC data
png(file = file.path(QC.dir,"MTvsCount.png"), width=8, height = 4, units= "in", res= 300)
plot1 <- FeatureScatter(seurat, feature1 = "nCount_spliced", feature2 = "percent.mt", group.by = "orig.ident")
plot1
dev.off()
svg(file = file.path(QC.dir,"MTvsCount.svg"))
plot1 <- FeatureScatter(seurat, feature1 = "nCount_spliced", feature2 = "percent.mt", group.by = "orig.ident")
plot1
dev.off()

#Label the top10 most DR genes
png(file = file.path(QC.dir,"FeaturevsCount.png"), width=8, height = 4, units= "in", res= 300)
plot2 <- FeatureScatter(seurat, feature1 = "nCount_spliced", feature2 = "nFeature_spliced", group.by = "orig.ident")
plot2
dev.off()
svg(file = file.path(QC.dir,"FeaturevsCount.svg"))
plot2 <- FeatureScatter(seurat, feature1 = "nCount_spliced", feature2 = "nFeature_spliced", group.by = "orig.ident")
plot2
dev.off()

png(file = file.path(QC.dir,"FeatureScatterCombined.png"), width=16, height = 4, units= "in", res= 300)
CombinePlots(plots = list(plot1,plot2))
dev.off()
svg(file = file.path(QC.dir,"FeatureScatterCombined.svg"), width = 16, height = 4)
CombinePlots(plots = list(plot1,plot2))
dev.off()


## ---- Doublet Detection, results = FALSE, message = FALSE, warning= FALSE, error = FALSE, fig.cap = "scDblFinder was used to automatically detect doublets. Doublet rate detected is around the expected 1% per thousand cells for 10X chromium."----
dblt <- as.Seurat(scDblFinder(as.SingleCellExperiment(seurat)))
seurat@meta.data[["scDblFinder.class"]] <- dblt@meta.data[["scDblFinder.class"]]

doublets <- table(dblt@meta.data[["scDblFinder.class"]])

p1 <- ggplot(FetchData(seurat, vars = c("scDblFinder.class", "nCount_spliced")), aes(x = scDblFinder.class, y = nCount_spliced, fill = scDblFinder.class)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(seurat$scDblFinder.class))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of transcripts', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p2 <- ggplot(FetchData(seurat, vars = c("scDblFinder.class", "nCount_spliced")), aes(x = scDblFinder.class, y = nCount_spliced, fill = scDblFinder.class)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(seurat$scDblFinder.class))) +
  scale_y_log10(labels = scales::comma) +
  labs(title = 'Number of transcripts', subtitle = 'log-scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p3 <- ggplot(FetchData(seurat, vars = c("scDblFinder.class", "nFeature_spliced")), aes(x = scDblFinder.class, y = nFeature_spliced, fill = scDblFinder.class)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(seurat$scDblFinder.class))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'linear scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

p4 <- ggplot(FetchData(seurat, vars = c("scDblFinder.class", "nFeature_spliced")), aes(x = scDblFinder.class, y = nFeature_spliced, fill = scDblFinder.class)) +
  geom_violin(draw_quantiles = c(0.5), scale = 'area', trim = FALSE) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(seurat$scDblFinder.class))) +
  scale_y_log10(labels = scales::comma) +
  labs(title = 'Number of expressed genes', subtitle = 'log-scale') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) +
  coord_flip()

ggsave(
  file.path(QC.dir,'qc_nCount_nFeature_by_multiplet_class.png'),
  p1 + p3 +
  p2 + p4 + plot_layout(ncol = 2),
  height = 7, width = 10
)

ggsave(
  file.path(QC.dir,'qc_nCount_nFeature_by_multiplet_class.svg'),
  p1 + p3 +
  p2 + p4 + plot_layout(ncol = 2),
  height = 7, width = 10
)

p1 + p3 + 
p2 + p4 



## ---- QC and Process each sample separetly -----------------------------------

##  ----- Control -----

    ## ---- Set up some QC stats and thresholds ----

control[["percent.mt"]] <- PercentageFeatureSet(control, pattern = "^MT-")

median_nCount <- median(control$nCount_spliced)

sd_nCount <- sd(control$nCount_spliced)

median_nFeature <- median(control$nFeature_spliced)

sd_nFeature <- sd(control$nFeature_spliced)


median_percent_MT <- median(control@meta.data[["percent.mt"]])

sd_percent_MT <- sd(control@meta.data[["percent.mt"]])


thresholds_nCount <- c(0, median_nCount + 2*sd_nCount)

thresholds_nFeature <- c(0, median_nFeature + 2*sd_nFeature)

thresholds_percent_MT <- c(0, median_percent_MT + 2*sd_percent_MT)

    ## ---- Use scDblFinder to label each cell as doublet or singlet ----

dblt <- as.Seurat(scDblFinder(as.SingleCellExperiment(control)))
control@meta.data[["scDblFinder.class"]] <- dblt@meta.data[["scDblFinder.class"]]

    ## ---- Filter out any cells outside of the QC cutoffs and remove doublets ----

control <- subset(control, subset = nCount_spliced < thresholds_nCount[2] & nFeature_spliced < thresholds_nFeature[2] & percent.mt < thresholds_percent_MT[2] & percent.mt < 25 & scDblFinder.class == "singlet")

    ## ---- Normalize to allow cellcyclescoring to be run

control <- NormalizeData(control)

control <- CellCycleScoring(control, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
control$CC.Difference <- control$S.Score - control$G2M.Score

    ## ---- SCTransform the data regressing out the difference in cell cycle score and the percent of mitochondrial reads in each cell ----

vars.to.regress = c("CC.Difference", "percent.mt")
control <- SCTransform(control, assay = "spliced", verbose = FALSE, vars.to.regress = vars.to.regress)

##  ----- Cocaine -----

    ## ---- Set up some QC stats and thresholds ----

cocaine[["percent.mt"]] <- PercentageFeatureSet(cocaine, pattern = "^MT-")

median_nCount <- median(cocaine$nCount_spliced)

sd_nCount <- sd(cocaine$nCount_spliced)

median_nFeature <- median(cocaine$nFeature_spliced)

sd_nFeature <- sd(cocaine$nFeature_spliced)


median_percent_MT <- median(cocaine@meta.data[["percent.mt"]])

sd_percent_MT <- sd(cocaine@meta.data[["percent.mt"]])


thresholds_nCount <- c(0, median_nCount + 2*sd_nCount)

thresholds_nFeature <- c(0, median_nFeature + 2*sd_nFeature)

thresholds_percent_MT <- c(0, median_percent_MT + 2*sd_percent_MT)

    ## ---- Use scDblFinder to label each cell as doublet or singlet ----

dblt <- as.Seurat(scDblFinder(as.SingleCellExperiment(cocaine)))
cocaine@meta.data[["scDblFinder.class"]] <- dblt@meta.data[["scDblFinder.class"]]

    ## ---- Filter out any cells outside of the QC cutoffs and remove doublets  ----

cocaine <- subset(cocaine, subset = nCount_spliced < thresholds_nCount[2] & nFeature_spliced < thresholds_nFeature[2] & percent.mt < thresholds_percent_MT[2] & percent.mt < 25 & scDblFinder.class == "singlet")

    ## ---- Normalize to allow cellcyclescoring to be run

cocaine <- NormalizeData(cocaine)

cocaine <- CellCycleScoring(cocaine, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
cocaine$CC.Difference <- cocaine$S.Score - cocaine$G2M.Score


    ## ---- SCTransform the data regressing out the difference in cell cycle score and the percent of mitochondrial reads in each cell ----

vars.to.regress = c("CC.Difference", "percent.mt")
cocaine <- SCTransform(cocaine, assay = "spliced", verbose = FALSE, vars.to.regress = vars.to.regress)

## ---- Merge the processed objects together ----

seurat <- merge(x = control, y = cocaine, add.cell.ids = samples.list, project = "Cerebral Organoids Cocaine Control Only")

## ---- Must manually set the variable features to the variable genes returned by SCTransform ----

obj.features <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 2000)
VariableFeatures(seurat[["SCT"]]) <- obj.features

## ---- Run PCA -------------------------------------------
PCA_npcs = 100
seurat <- RunPCA(seurat, npcs = PCA_npcs, )


## ---- Visualise Dimensions ----

p1 = VizDimLoadings(seurat, dims = 1:2, reduction = "pca")
p2 = DimPlot(seurat, reduction = "pca", group.by = "orig.ident") + theme(legend.position = 'none')
p3 = DimHeatmap(seurat, dims = 1:5, cells = 500, balanced = TRUE, fast = FALSE)

p3 | (p1 / p2)

ggsave(
  file = file.path(Analysis.dir,"DimVis.png"),
 p3 | (p1 / p2),
  height = 12, width = 15
)
ggsave(
  file = file.path(Analysis.dir,"DimVis.svg"),
 p3 | (p1 / p2),
  height = 12, width = 15
)


## ---- Elbowplot, to visualise significant PCA dimensions for downstream use ----

p1 = ElbowPlot(seurat, ndims = 100)

p1


ggsave(
  file = file.path(Analysis.dir,"Elbowplot.png"), 
  p1,
  height = 7, width = 10
       )

ggsave(
  file = file.path(Analysis.dir,"Elbowplot.svg"), 
  p1,
  height = 7, width = 10
       )



## ---- Store Variables in Seurat Object --------------

seurat@misc$experiment <- list(
  experiment_name = project.name,
  organism = 'hg',
  date_of_analysis = Sys.Date(),
  sample = "Merged Cerebral Organoids"
)

seurat@misc$filepaths <- list(
  project.dir = project.dir,
  save.dir = save.dir,
  QC.dir = QC.dir,
  Analysis.dir = Analysis.dir,
  Clustering.dir = Clustering.dir,
  UMAPS.dir = UMAPS.dir
  
)

seurat@misc$parameters <- list(
  gene_nomenclature = 'gene_name',
  discard_genes_expressed_in_fewer_cells_than = "3",
  keep_mitochondrial_genes = TRUE,
  PCA_Dims = PCA_npcs,
  number_PCs = "N/A",
  cluster_resolution = "N/A",
  seed.use = seed.use
)



## ---- Save File -------------------------------------
saveRDS(seurat, file= file.path(save.dir, paste(sample.name, ".rds", sep="")))

## ---- Save session info ----
sink(file.path(save.dir,"sessionInfoInitial.txt"))
     sessionInfo()
     sink()
