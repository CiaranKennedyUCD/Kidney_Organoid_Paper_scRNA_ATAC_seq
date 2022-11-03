library(ggtree)
library(Seurat)
library(dplyr)
library(ggplot2)

save.dir = seurat@misc[["filepaths"]][["save.dir"]]

seurat <- BuildClusterTree(seurat, assay = "SCT")

tree <- seurat@tools[["BuildClusterTree"]]


# ____Custom colours_____

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

colors_german <- c(
  "#ef5777", "#fd9644", "#fed330", "#26de81", "#2bcbba",
  "#eb3b5a", "#fa8231", "#f7b731", "#20bf6b", "#0fb9b1",
  "#45aaf2", "#4b7bec", "#a55eea", "#d1d8e0", "#778ca3",
  "#2d98da", "#3867d6", "#8854d0", "#a5b1c2", "#4b6584"
)
custom_colors$discrete <- c(colors_dutch, colors_spanish, colors_german)

custom_colors$cell_cycle <- setNames(
  c('#45aaf2', '#f1c40f', '#e74c3c', '#7f8c8d'),
  c('G1',      'S',       'G2M',     '-')
)


d = data.frame(label = tree[["tip.label"]], cluster = levels(seurat@meta.data[["clusters"]]))
tree <- full_join(tree, d, by='label')
tree$tip.label <- paste0("Cluster ", tree$tip.label)

p1 <- ggtree::ggtree(tree) +
  #scale_y_reverse() +
  ggtree::geom_tree() +
  ggtree::theme_tree() +
  ggtree::geom_tiplab(offset = 10, geom = "text", size = 5) + xlim(0, 300) +
  ggtree::geom_tippoint(color = custom_colors$discrete[1:length(tree@phylo[["tip.label"]])], shape = "square", size = 8)+
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(0,2.5,0,0), 'cm')) +
  theme(legend.position = "none") 




ggsave(file.path(save.dir,"h-clustree_text.png"), 
       p1,
       width = 5, height = 6 
)
ggsave(file.path(save.dir,"h-clustree_text.svg"),
       p1,
       width = 5, height = 6
)
