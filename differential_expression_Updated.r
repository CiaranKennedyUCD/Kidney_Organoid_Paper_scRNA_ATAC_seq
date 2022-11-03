#Import the needed libraries, MAST and DESeq2 are only needed if you want to use those tests, can be installed from bioconductor

library(Seurat)
library(MAST)
library(DESeq2)
library(openxlsx)
library(enrichR)

#Open the seurat file if not already open

seurat <- readRDS("your/file/path")

# Set what folder you want to save everything in, need to create folder first if it doesn't exist yet

DE.dir = "/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Brain_Organoids_Cocaine_Control_Only_15_07_2021/merged/Differential Expression"

#Set the default assay to "RNA"" (or "spliced" if using Velocyto pre-processed data) and make sure its normalized. Change this to "RNA" if the seurat object was made directly from 10x data

DefaultAssay(seurat) <- "spliced"

#Counts need to be normalized for a lot of the tests
seurat <- NormalizeData(seurat, verbose = FALSE)

#Extract number of clusters from the assigned identities

num.clusts = nlevels(Idents(seurat))

#Create a list of Differential Expression tools you want to use
de.list <- c("wilcox", "bimod", "t", "MAST")

#Create a list of samples, when creating a Seurat object I always save the sample as orig.ident (Original identity) so I can keep track of samples that I merge or integrate so
#I can pull a list of the different identities out of it. you could also save them as seurat@meta.data[["sample"]] if that works easier for you.

sample.list <- unique(seurat@meta.data[["orig.ident"]])

#Creates a list of each pairwise combination within the list of samples
samples.pairwise <- combn(sample.list, 2, simplify = FALSE)

dbs <- c("GO_Biological_Process_2021","GO_Cellular_Component_2021","GO_Molecular_Function_2021",
         "Elsevier_Pathway_Collection", "MSigDB_Hallmark_2020", "Allen_Brain_Atlas_10x_scRNA_2021","Descartes_Cell_Types_and_Tissue_2021",
         "KEGG_2021_Human", "WikiPathway_2021_Human", "MGI_Mammalian_Phenotype_Level_4_2021","CellMarker_Augmented_2021",
         "PanglaoDB_Augmented_2021","Azimuth_Cell_Types_2021", "BioPlanet_2019","ClinVar_2019", "TRRUST_Transcription_Factors_2019",
         "TF_Perturbations_Followed_by_Expression", "DSigDB" )

seurat$labels <- Idents(seurat)

#Creates a folder within your differential expression folder for each DE tool and then outputs a list of all marker genes for every cluster using each of the tools
for (x in seq_along(de.list)){
  DETool = de.list[[x]]
  dir.create(file.path(DE.dir,DETool))
  detool.dir = file.path(DE.dir,DETool)

  all.markers <- FindAllMarkers(seurat, min.pct = 0.1, logfc.threshold = 0.25, test.use = de.list[[x]], assay = "spliced", slot = "data", only.pos = TRUE)

  write.xlsx(tibble::rownames_to_column(all.markers, var = "Genes"), file = file.path(detool.dir,paste("All_Markers_", DETool, ".xlsx", sep = "")), overwrite = TRUE)
  
  Idents(seurat) <- seurat$orig.ident
  sample_DE <- FindMarkers(seurat, ident.1 = "Control", ident.2 = "Cocaine", test.use = de.list[[x]], assay = "spliced", min.pct = 0.1, logfc.threshold = 0.25, pseudocount.use = 1, verbose = TRUE)
  Idents(seurat) <- seurat$labels
  
  sample_DE2 <-tibble::rownames_to_column(sample_DE, var = "Genes")
  sig_de <- subset(sample_DE, subset = p_val_adj <= 0.05)
  enrich <- enrichr(rownames(sig_de), dbs)
  
  list <- c(list(sample_DE2), enrich)
  names(list)[1] <- "DE Results"
  
  write.xlsx(setNames(as.list(lapply(list, data.frame)), strtrim(names(list), 30)), file = file.path(detool.dir,paste("Control_Cocaine_DE_", DETool, "_enrichr.xlsx", sep = "")), overwrite = TRUE)
  
  }  

#Labelling each cell group in the form "cluster Name"_"Sample". This separates out the clusters by sample and labels the cells with both current label and sample. Might need to assign the correct identity first or change it to read from the slot you want to use

seurat$celltype <- Idents(seurat)
seurat$celltype.stim <- paste(Idents(seurat), seurat$orig.ident, sep = "_")

Idents(seurat) <- "celltype.stim"

#Run these three nested for-loops together. There are definitly more elegant ways to do this, but this is what I know how to do.

#First loop will go through each of the differential expression tools you want to use and set it from the list from earlier.
for (x in seq_along(de.list)){
DETool = de.list[[x]]
dir.create(file.path(DE.dir,DETool))
detool.dir = file.path(DE.dir,DETool)

#Second loop will iterate all of the pairwise combinations of samples as set them as sample 1 and sample 2 and create a folder for each comparison within each DE tools folder
for (j in 1:length(samples.pairwise)){
sample1 = as.data.frame((samples.pairwise))[[j]][[1]]
sample2 = as.data.frame((samples.pairwise))[[j]][[2]]

dir.create(file.path(detool.dir,paste(sample1,"vs",sample2, sep = " ")))
samp.dir = file.path(detool.dir,paste(sample1,"vs",sample2, sep = " "))

#Third loop will go through each cluster and do the differential expression within the cluster between sample 1 and sample 2
for(i in 1:num.clusts){ 
  
  #The if statement can be used if there's a sample that's missing one or more of the cell types, basically just tells the script to skip 
  #any comparisons that would use the missing group of cells or it would just break. Comment it and one of the brackets at the end if its not needed.
  
  #if (paste(levels(seurat$celltype)[[i]], "_", sample1, sep = "") != "Nephron Progenitor_Control" & paste(levels(seurat$celltype)[[i]], "_", sample2, sep = "") != "Nephron Progenitor_Control") &paste(levels(seurat$celltype)[[i]], "_", sample1, sep = "") != "Nephron Progenitor_Control" & paste(levels(seurat$celltype)[[i]], "_", sample2, sep = "") != "Nephron Progenitor_Control"){
    
  #Just prints a message to the console telling you what the script is doing now
  print(paste("Calculating cluster",levels(seurat$celltype)[[i]],"markers now in", sample1, "vs", sample2, "with", DETool))
    
  #The differential expression test
  clust <- try(FindMarkers(seurat, ident.1 = paste(levels(seurat$celltype)[[i]], "_", sample1, sep = ""), ident.2 = paste(levels(seurat$celltype)[[i]], "_", sample2, sep = ""), test.use = de.list[[x]], assay = "spliced", min.pct = 0.1, logfc.threshold = 0.25, pseudocount.use = 1, verbose = TRUE))
  
  #Move on to next item if a cluster isnt found
  if(isTRUE(class(clust)=="try-error")) { next }
  
  # Perform enrichR testing on the results
  
  clust2 <-tibble::rownames_to_column(clust, var = "Genes")
 
  sig_de <- subset(clust, subset = p_val_adj <= 0.05)
  enrich <- enrichr(rownames(sig_de), dbs)
  
  list <- c(list(clust2), enrich)
  names(list)[1] <- "DE Results"
  
  # Save all the results in a file 

  list <- c(list(clust2), enrich)
  names(list)[1] <- "DE Results"
  
  ct <- gsub("/", " ", paste(levels(seurat$celltype)[[i]]))
  
  write.xlsx(setNames(as.list(lapply(list, data.frame)), strtrim(names(list), 30)), file = file.path(samp.dir,paste(ct, "markers", sample1, "vs", sample2, DETool, "enrichR.xlsx", sep = "_")), overwrite = TRUE)
 
}
}
}
#}
