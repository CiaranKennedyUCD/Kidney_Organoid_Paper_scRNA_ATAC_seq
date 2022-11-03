p1 <- DotPlot(combined3, assay = "chromvar", features = c("MA0795.1", "MA1557.1", "MA0478.1", "MA0490.2", "MA0144.2", "MA0050.2", "MA0162.4")) + scale_y_discrete(limits = rev(c("Stroma 2_Control", "Stroma 3_ControlTGFB", "Stroma 3_GSKTGFB")))+ theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(filename = paste(save.dir,"Plots_For_Paper/","dotplot.png", sep = ""), plot = p1, width = 8, height = 6, units = "in")
ggsave(filename = paste(save.dir,"Plots_For_Paper/","dotplot.svg", sep = ""), plot = p1, width = 8, height = 6, units = "in", dpi = 600, limitsize = FALSE)


DefaultAssay(combined) <- "chromvar"
p2 <- FeaturePlot(combined, features = "MA0478.1", min.cutoff = 0.1, max.cutoff = 2)

ggsave(filename = paste(save.dir,"Plots_For_Paper/","FOSL2.png", sep = ""), plot = p2, width = 9, height = 8, units = "in")
ggsave(filename = paste(save.dir,"Plots_For_Paper/","FOSL2.svg", sep = ""), plot = p2, width = 9, height = 8, units = "in", dpi = 600, limitsize = FALSE)
