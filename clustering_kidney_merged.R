## ---- Load Libraries --------------------------------
library(plyr)
library(Seurat)
library(future)
library(clustree)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)

## ---- Set parallelisation ---------------------------
plan("multisession", workers = 8)
options(future.globals.maxSize = 30000 * 1024^2)


## ---- Load file -------------------------------------
RDS.file.location <- args[1]
seurat <- readRDS(RDS.file.location)


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
  '#40407a','#706fd3','#1e272e','#34ace0','#33d9b2',
  '#2c2c54','#474787','#aaa69d','#227093','#218c74',
  '#ff5252','#ff793f','#d1ccc0','#ffb142','#ffda79',
  '#b33939','#cd6133','#84817a','#cc8e35','#ccae62'
)

custom_colors$discrete <- c(colors_dutch, colors_spanish)

custom_colors$cell_cycle <- setNames(
  c('#45aaf2', '#f1c40f', '#e74c3c', '#7f8c8d'),
  c('G1',      'S',       'G2M',     '-')
)


## ---- Extrace File paths and variables --------------
project.dir = gsub("//", "/",seurat@misc[["filepaths"]][["project.dir"]])
save.dir = gsub("//", "/",seurat@misc[["filepaths"]][["save.dir"]])
QC.dir = gsub("//", "/",seurat@misc[["filepaths"]][["QC.dir"]])
Analysis.dir = gsub("//", "/",seurat@misc[["filepaths"]][["Analysis.dir"]])
Clustering.dir = gsub("//", "/",seurat@misc[["filepaths"]][["Clustering.dir"]])
UMAPS.dir = gsub("//", "/",seurat@misc[["filepaths"]][["UMAPS.dir"]])

seed.use = seurat@misc[["parameters"]][["seed.use"]]


## ---- Set Dims Number -------------------------------
dims.number = 1:args[2]
seurat@misc[["parameters"]][["number_PCs"]] <- dims.number


## ---- Find Neigbhours --------------------------------
seurat <- FindNeighbors(seurat, dims = dims.number)


## ---- Calculate and save clustering and umaps at multiple resolutions ----

res_list <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2.0)

plots <- list()

seurat <- RunUMAP(seurat, dims = dims.number, seed.use = seed.use)

for (i in seq_along(res_list)) 
  {
res.dir <- file.path(UMAPS.dir,res_list[i],"/")
dir.create(res.dir, showWarnings = FALSE)

seurat <- FindClusters(seurat, resolution = res_list[i] , random.seed = seed.use, reduction = "lsi")

#UMAP

p = DimPlot(seurat, reduction = "umap", cols = custom_colors$discrete)

plots[[i]] = p + ggtitle(paste("Resolution = ", res_list[i], sep = "")) +
          theme(plot.title = element_text(hjust = 0.5))

ggsave(filename = paste(res.dir,"UMAP.svg", sep = ""), plot = p, width = 12, height = 8, units = "in", dpi = 600)
ggsave(filename = paste(res.dir,"UMAP.png", sep = ""), plot = p, width = 12, height = 8)
} 




## ---- Stitch together all of the UMAPS into a figure ----
wrap_plots(plots, ncol = 4) + plot_annotation(tag_levels = 'A')

ggsave(file.path(UMAPS.dir,"UMAPS_Combined.png"), scale = 2, dpi = 400, width = 20, height = 15, units = "in",
wrap_plots(plots, ncol = 4) + plot_annotation(tag_levels = 'A')
)
ggsave(file.path(UMAPS.dir,"UMAPS_Combined.svg"), scale = 2, dpi = 400, width = 20, height = 15, units = "in",
wrap_plots(plots, ncol = 4) + plot_annotation(tag_levels = 'A')
)


## ---- Create and save the clustree trees-----------------------------
  p1 <- clustree(seurat)   
 
 p1
 
  ggsave(file.path(UMAPS.dir,"clustree.png"), 
         p1,
         width = 10, height = 24 
         )
  ggsave(file.path(UMAPS.dir,"clustree.svg"), 
         p1,
         width = 10, height = 24
         )



## --------- Plot a few useful figures ----------
p1 <- ggplot(seurat@meta.data, aes(x = orig.ident, fill = Phase)) + 
    geom_bar(position = "dodge", ) +
    scale_fill_manual(values = custom_colors$discrete) +
    scale_x_discrete(limits = levels(seurat$orig.ident)) +
    scale_y_continuous() +
    labs(title = 'Proportion of cells in each cell state', subtitle = 'By Sample') +
    theme(
      axis.title = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
    ) 
    
p2 <- DimPlot(seurat, reduction = "umap", group.by = "orig.ident") + 
  theme(
      axis.title = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
    ) 
p3 <- DimPlot(seurat, reduction = "umap", group.by = "TGFB") +
  theme(
      axis.title = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
    ) 
p4 <- DimPlot(seurat, reduction = "umap", group.by = "GSK") +
  theme(
      axis.title = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
    )

p1 + p2 + p3 + p4

ggsave(file.path(UMAPS.dir,"breakdown.png"), width = 15, height = 10, units = "in", dpi = 600,
   wrap_plots(p1 + p2 + p3 + p4,  widths = 10, 
    heights = 10 ) +
    plot_annotation(tag_levels = 'A')
   
)

ggsave(file.path(UMAPS.dir,"breakdown.svg"), width = 15, height = 10, units = "in",
    wrap_plots(p1 + p2 + p3 + p4,  widths = 10, 
    heights = 10) +
    plot_annotation(tag_levels = 'A')

)



## -------- Run Phate Analysis -----
seurat <- RunPHATE(object = seurat, assay = "RNA")

p1 <- DimPlot(seurat, reduction = "phate")

ggsave(file.path(UMAPS.dir,"phate.svg"), width = 20, height = 10, units = "in",
       p1
)
ggsave(file.path(UMAPS.dir,"phate.png"), width = 20, height = 10, units = "in",
       p1
)

## ------ Sae file--------------
saveRDS(seurat, RDS.file.location)

