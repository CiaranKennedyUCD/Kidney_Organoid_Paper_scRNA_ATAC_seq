library(Signac)
library(monocle3)
library(cicero)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(SeuratWrappers)
library(future)

# This script splits a Signac object by sample and calculates cis-co-accessability 
# networks through Cicero and adds them to the individual objects

args = commandArgs(trailingOnly=TRUE)

save.dir <- "/scratch/10309667/SeuratOutputs/Kidney_Organoids"
combined <- readRDS("/scratch/10309667/SeuratOutputs/Kidney_Organoids/combined_Relabelled.rds")

plan("multiprocess", workers = 4)
options(future.globals.maxSize = 30000 * 1024^2) # for 30 Gb RAM

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"

samples <- c("Control", "ControlTGFB", "GSK", "GSKTGFB")
combined_split <- SplitObject(combined, split.by = "orig.ident")

i = args[1]  

sample = samples[[i]]

split <- combined_split[[i]]

combined.cds <- as.cell_data_set(x = split, assay = "peaks")
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

saveRDS(split, file = file.path(save.dir, paste(sample,"cicero_conns.rds", sep = "_")))

#Convert pairwise co-assability links to co-accessability networks
ccans <- generate_ccans(conns)

#Convert connections to links and add to the seurat object

links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(split) <- links

saveRDS(split, file = file.path(save.dir, paste(sample,"Signac_Links.rds", sep = "_")))

