library(Seurat)
library(patchwork)
library(ggplot2)

##Adapted from https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )+
    ylab(feature) +
    theme(legend.position = "none",
          plot.title= element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          axis.text.y = element_text(size = rel(1), vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(0,0,0,0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle=45,hjust=1, vjust = 0.9), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x +
                            scale_y_continuous(breaks = c(y)) +
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1) +
    plot_annotation(theme = theme(plot.margin = margin(0.5,0.5,0.5,0.5,"cm")))
  
  return(p)
}

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


features<- c("PITX2", "AKR1C2", "LUM", "EDNRA", "RBM47", "PAX2", "SOX10", "CDH19", "PLP1", "RIPPLY3", 
             "ACTA1", "MYBPH", "HES5", "OLIG3", "STMN2", "BARHL2", "LHX2", "STMN4", "CELF3", 
             "CDK1", "MKI67", "PBK"
             )

vln <- StackedVlnPlot(obj = seurat, features = features)

ggsave(filename = file.path(save.dir,"UMAP_Prediction.svg"), plot = vln, width = 10, height = 10, units = "in", dpi = 600)
ggsave(filename = file.path(save.dir,"UMAP_Prediction.png"), plot = vln, width = 10, height = 10)


ggsave(filename = file.path(save.dir,"UMAP_Labelled.svg"), plot = umap, width = 10, height = 10, units = "in", dpi = 600)
ggsave(filename = file.path(save.dir,"UMAP_Labelled.png"), plot = umap, width = 10, height = 10)

heatmap.genes <- c("PAX2", "PAX8", "EMX2", "KRT19", "KRT18",
                   "SIX1", "PITX2", "PDGFC", "EYA1", "SIX2",
                   "RSPO3", "IGFBP7", "PCSK1N", "RSPO3", "SALL1",
                   "NPHS2", "PODXL", "MAFB", "WT1", "SOST",
                   "TYMS", "NUSAP1", "CLSPN", "CENPF",
                   "UBE2C", "CENPF", "TOP2A", "CCNB2", "PTTG1",
                   "BST2", "NR2F1", "NR2F2", "ATP2B1", "IGDCC3",
                   "PLP1", "TFAP2B", "S100B", "PMEL", "EDNRB",
                   "STMN2", "ELAVL4", "TAGLN3", "PPP1R17", "ELAVL2",
                   "COL9A1", "MGP", "CNMD", "COL9A3", "FIBIN",
                   "MEOX2", "CXCL14", "LUM", "KIFAP3", "PITX1",
                   "TSHZ1", "RND3", "AMOT", "EFNA5", "TCF7L2",
                   "COL1A1", "DNM3OS", "COL3A1", "COL1A2", "KCNQ1OT1",
                   "IGF1", "COL21A1", "LUM", "CXCL14", "TWIST1"
)

heatmap <- DoHeatmap(subset(seurat, downsample = 100), features = heatmap.genes, size = 3)

ggsave(filename = paste("G:/Shared drives/Crean Group Data Share/Kidney_Organoids_Hydrogel_09_10_20/merged/Labelling/","Heatmap.svg", sep = ""), plot = heatmap, width = 12, height = 10, units = "in", dpi = 600)
ggsave(filename = paste("G:/Shared drives/Crean Group Data Share/Kidney_Organoids_Hydrogel_09_10_20/merged/Labelling/","Heatmap.png", sep = ""), plot = heatmap, width = 12, height = 10)

##Reordering levels for labels
# levels(seurat) <- c("Nephron Progenitor I", "Nephron Progenitor II", "Nephron Progenitor III", "Podocyte", "Cell Cycle I", "Mesangial", "Melanocyte", "Neural Progenitor", "Stroma", "Stroma Progenitor I", "Stroma Progenitor II", "Stroma Progenitor III", "Stroma Progenitor IV")  
