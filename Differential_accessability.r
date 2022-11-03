library(openxlsx)
library(BSgenome.Hsapiens.UCSC.hg38)
library(enrichR)

combined <- readRDS("/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Kidney Organoids ATAC 2_14_03_21/combined/combined_Relabelled.rds")

project.dir = gsub("//", "/",combined@misc[["filepaths"]][["project.dir"]])
save.dir = gsub("//", "/",combined@misc[["filepaths"]][["save.dir"]])
QC.dir = gsub("//", "/",combined@misc[["filepaths"]][["QC.dir"]])
Analysis.dir = gsub("//", "/",combined@misc[["filepaths"]][["Analysis.dir"]])
Clustering.dir = gsub("//", "/",combined@misc[["filepaths"]][["Clustering.dir"]])
UMAPS.dir = gsub("//", "/",combined@misc[["filepaths"]][["UMAPS.dir"]])
de.dir <- paste(save.dir, "DE", sep = "")

Idents(combined) <- combined$orig.ident
DefaultAssay(combined) <- 'peaks'

sample.list <- unique(combined@meta.data[["orig.ident"]])
samples.pairwise <- combn(sample.list, 2, simplify = FALSE)

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = "Homo sapiens", all_versions = FALSE)
)

combined <- AddMotifs(
  object = combined,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

dbs <- c("GO_Biological_Process_2021","GO_Cellular_Component_2021","GO_Molecular_Function_2021",
         "Elsevier_Pathway_Collection", "MSigDB_Hallmark_2020", "Allen_Brain_Atlas_10x_scRNA_2021","Descartes_Cell_Types_and_Tissue_2021",
         "KEGG_2021_Human", "WikiPathway_2021_Human", "MGI_Mammalian_Phenotype_Level_4_2021","CellMarker_Augmented_2021",
         "PanglaoDB_Augmented_2021","Azimuth_Cell_Types_2021", "BioPlanet_2019","ClinVar_2019", "TRRUST_Transcription_Factors_2019",
         "TF_Perturbations_Followed_by_Expression", "DSigDB" )


for (j in 1:length(samples.pairwise)){
  sample1 = as.data.frame((samples.pairwise))[[j]][[1]]
  sample2 = as.data.frame((samples.pairwise))[[j]][[2]]

  dir.create(file.path(de.dir,paste(sample1,"vs",sample2, sep = " ")))
  samp.dir = file.path(de.dir,paste(sample1,"vs",sample2, sep = " "))

da_peaks <- FindMarkers(
  object = combined,
  ident.1 = sample1,
  ident.2 = sample2,
  min.pct = 0.05,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

# head(da_peaks)

open <- rownames(da_peaks[da_peaks$p_val_adj <= 0.05, ])

da_peaks <- tibble::rownames_to_column(da_peaks, var = "Peaks")

closest_gene <- ClosestFeature(combined, regions = open)

enrich <- enrichr(closest_gene[[2]], dbs)

merged <- merge(x = da_peaks, y = closest_gene, by.x = "Peaks", by.y = "query_region", all.x = TRUE)

list <- c(list(merged), enrich)
names(list)[1] <- "DE Results"

write.xlsx(setNames(as.list(lapply(list, data.frame)), strtrim(names(list), 30)), file = file.path(samp.dir, paste(sample1, sample2, "DA_Peaks_Closest_Gene_Enrichr.xlsx", sep = "_")), overwrite = TRUE)

sig_da_peaks <- subset(da_peaks, subset = p_val_adj <= 0.05)

combined <- RegionStats(combined, BSgenome.Hsapiens.UCSC.hg38, assay = "peaks")

enriched.motifs <- FindMotifs(
  object = combined,
  features = open
)

mp <- MotifPlot(
  object = combined,
  motifs = head(rownames(enriched.motifs)),
  assay = 'peaks',
  ncol = 2
) 


write.xlsx(enriched.motifs, file.path(samp.dir, paste(sample1, sample2, "DA_Peaks_Motif_Enrichment.xlsx", sep = "_")), overwrite = TRUE)

ggsave(filename = paste(samp.dir,paste(sample1, sample2, "Enriched_Motifs.svg", sep = "_"), sep = "/"), plot = mp, width = 8, height = 5, units = "in", dpi = 600)
ggsave(filename = paste(samp.dir,paste(sample1, sample2, "Enriched_Motifs.png", sep = "_"), sep = "/"), plot = mp, width = 8, height = 5, units = "in", dpi = 600)
}

combined$celltype.stim <- paste(combined$clusters, combined$orig.ident, sep = "_")

Idents(combined) <- "celltype.stim"

for (j in 1:length(samples.pairwise)){
  sample1 = "Kidney Progenitor 1_Control"
  sample2 = "Kidney Progenitor 2_ControlTGFB"
  
  dir.create(file.path(de.dir,paste(sample1,"vs",sample2, sep = " ")))
  samp.dir = file.path(de.dir,paste(sample1,"vs",sample2, sep = " "))
  
  da_peaks <- FindMarkers(
    object = combined,
    ident.1 = sample1,
    ident.2 = sample2,
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'peak_region_fragments'
  )
  
  # head(da_peaks)
  
  open <- rownames(da_peaks[da_peaks$p_val_adj <= 0.05, ])
  
  da_peaks <- tibble::rownames_to_column(da_peaks, var = "Peaks")
  
  closest_gene <- ClosestFeature(combined, regions = open)
  
  enrich <- enrichr(closest_gene[[2]], dbs)
  
  merged <- merge(x = da_peaks, y = closest_gene, by.x = "Peaks", by.y = "query_region", all.x = TRUE)
  
  list <- c(list(merged), enrich)
  names(list)[1] <- "DE Results"
  
  write.xlsx(setNames(as.list(lapply(list, data.frame)), strtrim(names(list), 30)), file = file.path(samp.dir, paste(sample1, sample2, "DA_Peaks_Closest_Gene_Enrichr.xlsx", sep = "_")), overwrite = TRUE)
  
  sig_da_peaks <- subset(da_peaks, subset = p_val_adj <= 0.05)
  
  combined <- RegionStats(combined, BSgenome.Hsapiens.UCSC.hg38, assay = "peaks")
  
  enriched.motifs <- FindMotifs(
    object = combined,
    features = open
  )
  
  mp <- MotifPlot(
    object = combined,
    motifs = head(rownames(enriched.motifs)),
    assay = 'peaks',
    ncol = 2
  ) 
  
  
  write.xlsx(enriched.motifs, file.path(samp.dir, paste(sample1, sample2, "DA_Peaks_Motif_Enrichment.xlsx", sep = "_")), overwrite = TRUE)
  
  ggsave(filename = paste(samp.dir,paste(sample1, sample2, "Enriched_Motifs.svg", sep = "_"), sep = "/"), plot = mp, width = 8, height = 5, units = "in", dpi = 600)
  ggsave(filename = paste(samp.dir,paste(sample1, sample2, "Enriched_Motifs.png", sep = "_"), sep = "/"), plot = mp, width = 8, height = 5, units = "in", dpi = 600)
}
