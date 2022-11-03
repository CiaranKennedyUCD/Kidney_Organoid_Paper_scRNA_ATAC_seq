library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biomaRt)
library(Seurat)
library(ggplot2)
library(ggpubr) 
library(patchwork)

# ---- Setup ----


#Import The data

signac <- readRDS("/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Kidney Organoids ATAC 2_14_03_21/combined/combined_Relabelled.rds")

DefaultAssay(signac) <- "peaks"

## Set up annotations

genes <- genes(EnsDb.Hsapiens.v86)
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
bm <- getBM(attributes = c("external_gene_name",'ensembl_gene_id', "start_position", "end_position"), values=names(genes),filters ='ensembl_gene_id', mart = mart)
names(genes) <- bm$external_gene_name[match(genes$gene_id,bm$ensembl_gene_id)]

seqlevelsStyle(genes) <- "UCSC"

## Get save directories stored in the Seurat object

project.dir = gsub("//", "/",signac@misc[["filepaths"]][["project.dir"]])
save.dir = gsub("//", "/",signac@misc[["filepaths"]][["save.dir"]])
QC.dir = gsub("//", "/",signac@misc[["filepaths"]][["QC.dir"]])
Analysis.dir = gsub("//", "/",signac@misc[["filepaths"]][["Analysis.dir"]])
Clustering.dir = gsub("//", "/",signac@misc[["filepaths"]][["Clustering.dir"]])
UMAPS.dir = gsub("//", "/",signac@misc[["filepaths"]][["UMAPS.dir"]])

seed.use = signac@misc[["parameters"]][["seed.use"]]


## Set up discrete Custom Colour Pallet
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

# ---- Gene list ----

##Set list of features to use, and a name for the modules

module.name <- "Matrisome"

features <- c("ESX1", "CACNA2D2", "TEX35", "CACNA2D3", "ISL2", "BTBD11", "SCO2", "RNF220", "PTF1A", "HCN4", "CRABP1", 
                "DLL1", "DLL4", "DPYSL5", "SPOCK3", "DPYSL4", "SPOCK1", "PKNOX2", "SAMD5", "JAG1", "CHST7", "SAP30BP", 
                "CHST8", "NEBL", "NR2F2", "PBX3", "NR2F1", "TEX14", "C12ORF56", "AIF1L", "CHST2", "FAT4", "SFXN3", "STRN", 
                "LRRC10B", "CDCA4", "NCAPG2", "COL13A1", "NR2E1", "PDGFB", "GATA6", "ATP5C1", "GATA5", "TEX26", "GATA4", 
                "GATA3", "HTR6", "GATA2", "HTR7", "T", "EPB41L2", "FAM110C", "LYPD1", "NKX3-1", "MTERF1", "MNX1", "FAM19A5", 
                "PRKACA", "NKX3-2", "DUSP5", "JUN", "DUSP4", "DUSP2", "SPHK1", "ADAM33", "WNT9A", "WNT9B", "ADAM30", "DCBLD1", 
                "FXYD1", "CENPE", "CRLF1", "MMP15", "SLC7A8", "MMP16", "BHLHE40", "BHLHE41", "RAPGEF3", "ADGRL3", "MIXL1", 
                "ACHE", "KCNC2", "KCNC4", "ZC3HAV1", "HOXC13", "C12ORF43", "HOXC12", "MED16", "SHB", "HOXC11", "SHE", 
                "HOXC10", "PANX2", "MED15", "SHH", "MECOM", "DLGAP2", "AP4S1", "DLGAP4", "ZNF467", "SLC18A3", "EGLN1", 
                "PPARA", "EGLN3", "DGAT2", "DGAT1", "IFT140", "WNT7A", "CNPY1", "ADAM12", "WNT7B", "ADAM11", "UNCX", 
                "MYBBP1A", "LARS2", "PTRH2", "LGR6", "LGR5", "CRB3", "KCNA1", "KCNA2", "WRB", "PRKAG2", "KCNA6", "FADS3", 
                "FADS6", "KIAA1147", "NEFL", "NEFM", "JAK3", "TP53I11", "FADS2", "ZSWIM4", "MAP3K2", "RASL10B", "CTNNBIP1", 
                "FBXO17", "PCDH8", "RASL10A", "TBR1", "ATP11A", "WT1", "MSH5", "PKP4", "LY6H", "PKP1", "KCNG1", "KCNG3", 
                "ALDH1L1", "CRCP", "NPAS1", "NPAS2", "NPAS3", "NPAS4", "FRZB", "ADORA1", "ZNF423", "KCND3", "FST", "FBXO32", 
                "SMO", "HPSE", "CERS1", "TRIM71", "FBXO25", "LIN7A", "ELFN2", "ARHGAP4", "CORO1C", "ELFN1", "GOLIM4", "MMP27", 
                "G3BP1", "ZNF408", "TRIM67", "SYBU", "KCNJ6", "KCNJ9", "COL26A1", "SPACA1", "HTR1B", "C10ORF11", "NRG2", 
                "HTR1A", "NRG3", "BICC1", "SEZ6L", "SP5", "ITGA8", "SP9", "DCTPP1", "SP8", "TMBIM6", "KCNK1", "KCNK2", 
                "KCNK4", "CD244", "DLX1", "SNAP25", "TRIM50", "DLX2", "LCN9", "DLX3", "SHC3", "SHC4", "DLX4", "KCNK9", "DLX5", 
                "ONECUT3", "ONECUT1", "ONECUT2", "MLLT3", "ADAMTS15", "BOC", "PPP6R2", "CCDC85C", "ADAMTS19", "TRIM47", 
                "B4GALNT2", "B4GALNT1", "PRLHR", "KCNH3", "NTRK2", "TBX2", "TBX1", "NTRK3", "SRCIN1", "KCNH6", "NTRK1", 
                "UNC5B", "GRID1", "UNC5C", "KCNK12", "KCNK13", "UNC5A", "MYO10", "CABLES1", "HCRTR1", "KCNK15", "KIAA1549L", 
                "SYT16", "HDGFRP3", "UNC5D", "TBX4", "TSPAN15", "MYO1C", "MCM3AP", "HAND1", "SYT12", "HAND2", "NAGS", 
                "KCTD12", "EIF3D", "ITGB1", "ITGB4", "NRGN", "FOXRED2", "ZNRF2", "EFR3B", "ZNRF4", "BAHD1", "KCTD20", 
                "IGF2BP2", "IGF2BP1", "PCDH7", "IGF2BP3", "UBQLN4", "PCDH1", "KCNJ4", "CCDC79", "ACTN3", "FBXO2", "MOGS", 
                "ISLR2", "MYO7B", "SMARCA2", "MED13L", "APLN", "SOAT2", "CRISPLD2", "FLYWCH1", "NINJ2", "BCORL1", "CERK", 
                "CPNE7", "CCDC102A", "ADGRA2", "BSG", "E2F3", "KDR", "EGR2", "EGR3", "RTBDN", "NTNG1", "EGR4", "FAM78A", 
                "MYO5B", "IGF2", "TP73", "AIMP1", "SDC1", "GRAP", "TTBK1", "SLC45A3", "NPFFR2", "HHIP", "LPAR3", "LPAR5", 
                "IKZF3", "DMBX1", "MCIDAS", "CRTAC1", "APOF", "LGALS7", "APOBEC4", "GABBR2", "GABBR1", "GREB1", "GAB2", 
                "ZHX1", "ZHX2", "VEGFA", "TFAP4", "KITLG", "KCNS3", "PABPC4", "FAM193B", "FJX1", "PXT1", "DDR1", "NRN1", 
                "GRIK4", "GRIK3", "GRIK1", "IGSF21", "KALRN", "IGSF23", "PPP3CA", "C8G", "SEC14L1", "NFIL3", "PCK2", 
                "SLC32A1", "GAD2", "GAD1", "KAZALD1", "MYO9B", "KIN", "GAK", "SH3RF3", "KCNQ1", "KCNQ2", "SGSM3", "PNPLA7", 
                "NRP1", "NRP2", "SYNM", "PRICKLE1", "ENPP1", "MC4R", "RAC2", "FILIP1L", "MFSD2A", "GIT1", "TFAP2B", "TFAP2C", 
                "TFAP2D", "TFAP2E", "CACHD1", "HS3ST6", "HS3ST5", "GRIN1", "HS3ST2", "SHQ1", "GPRIN2", "CCPG1", "SIDT1", 
                "PRR15L", "EMID1", "YPEL4", "HHIPL1", "TFAP2A", "CRTC3", "CELSR1", "ODF3B", "SFSWAP", "PRDM6", "NEXN", "OTP", 
                "ASAP3", "CELSR3", "RELB", "GRIP1", "RELN", "PRR14L", "ASAP1", "PEBP4", "TBX21", "TBC1D14", "TBX20", "SLIT2", 
                "FGF20", "MSX1", "MSX2", "SLIT3", "SYT6", "POU3F1", "BHLHA9", "SYT2", "AUTS2", "POU3F2", "POU3F3", "GNAL", 
                "SDK2", "FGF16", "IFT22", "FGF19", "GNAS", "MICU2", "FGF12", "FBN1", "FBN2", "PKDCC", "PCDH11X", "HTRA4", 
                "CACNA1A", "TACO1", "CACNA1B", "PLAU", "BHLHE22", "CACNA1E", "HTRA1", "BHLHE23", "ETS1", "CACNA1G", "HOXB13", 
                "CACNA1I", "CACNA1H", "CNN2", "PTRF", "TGM7", "IL6R", "TTC16", "MKL1", "OLIG1", "OLIG2", "SHROOM1", "OLIG3", 
                "GGN", "SHANK1", "ADGRD1", "CYP26A1", "ACTA1", "TMEM215", "ACOX1", "ALPL", "FAM72A", "PLEC", "NAT8L", 
                "SHANK3", "DNAH2", "AGPAT9", "DNAH8", "TNFAIP2", "AGPAT1", "PHTF2", "TTC28", "LMNB2", "RAI1", "CA2", "CA4", 
                "GRK5", "SNORD37", "CA7", "C9ORF172", "EYA4", "ZIC5", "TBX18", "FAM155B", "SNTB1", "CPED1", "TBX15", "GREM1", 
                "ERCC4", "ACKR1", "GFOD1", "ADGRB1", "ADGRB2", "DENND6A", "ERCC8", "PTX4", "TRIB2", "TCHH", "RND1", "GRM1", 
                "GRM8", "GRM7", "ZIC1", "ZIC2", "FAM73A", "DIP2A", "HOMER2", "AFAP1", "EVA1A", "SCRIB", "HNRNPM", "CD5", 
                "POLR3A", "CCS", "ASTN2", "TES", "ADCY1", "ADCY5", "RNF8", "ADCY9", "ASCL2", "ADCY8", "ASCL1", "ADCY6", 
                "FBXL21", "MTA1", "FBXL20", "MXI1", "GNG7", "FLRT2", "PXDC1", "MTA3", "ST3GAL1", "PRSS3", "PROK2", "VPS18", 
                "PRKCH", "PLK5", "ICMT", "PHC3", "MCAM", "PRKCE", "PRKCB", "RMDN2", "FOSL2", "COL2A1", "BCAN", "SNAI3", 
                "FREM2", "FARSA", "XKR7", "ZNF280D", "CDK5RAP2", "ALX4", "ALX3", "PROM1", "GABRA1", "RPS6KL1", "GABRA5", 
                "NFATC2", "ASH1L", "NFATC1", "WDR33", "PHOX2B", "PHOX2A", "NCOR2", "GAS1", "CLCF1", "HLA-DQB1", "RARA", 
                "GAS6", "MGLL", "GAS8", "GRAMD1B", "MESP1", "ARPC1B", "PRDM16", "PRDM14", "AEBP1", "PRDM13", "PRDM12", "GJC2", 
                "RASSF10", "PLLP", "JARID2", "TEAD3", "ECEL1", "CDKN2C", "AHSG", "CDKN2A", "COL23A1", "RRP1B", "NFIA", 
                "WDR13", "ELP6", "DHRS3", "FARP1", "FAM212B", "PLEKHO1", "NFIC", "SMOC2", "NDUFA4L2", "TMEM229A", "CDKN1C", 
                "RARG", "GRP", "SETD7", "SETD8", "FAM212A", "PPP1R8", "GSC", "OBFC1", "RRP12", "DEPTOR", "CSGALNACT1", 
                "PLXNC1", "SKIDA1", "DIO3", "HIC1", "LRFN5", "LTB4R2", "ZNF319", "FKBP8", "TMEM178A", "PLXNB2", "CNIH3", 
                "NFIX", "PLXND1", "CSF1", "FAM129B", "FGD2", "LAMP5", "GOLGA7B", "RAB20", "ABLIM2", "CLCC1", "TPM1", 
                "SLC39A8", "MB21D1", "TMEM132E", "C1QL1", "TMEM132D", "OSBPL6", "TMEM132A", "TRPC3", "FAM217B", "NSG1", 
                "C1QL4", "SCRT2", "C1QL3", "SCRT1", "C1QL2", "SYNJ2", "LHX4", "LHX3", "LHX2", "PITPNM3", "LHX1", "LHX8", 
                "LHX6", "PGR", "LHX5", "TOX", "LHX9", "BRSK2", "RAB3A", "INSIG2", "LARP1B", "FGF3", "TM7SF2", "FGF4", "FGF5", 
                "RIMS2", "FGF8", "RIMS4", "FGF9", "SALL3", "SALL4", "SRGAP3", "PCDHAC2", "PLXNA4", "SECISBP2L", "APOA4", 
                "KBTBD12", "BCL6", "APOA1", "BCL2", "MCM2", "TNRC6C", "COL14A1", "FAM149B1", "FAM114A1", "LBH", "PSTPIP2", 
                "RRP7A", "DENND3", "P2RY1", "EMILIN3", "PHACTR2", "EMILIN2", "TNPO2", "PRKG2", "CCDC30", "PRKG1", "COL27A1", 
                "CR2", "UGT1A1", "CEP131", "FBXL18", "FBXL19", "LEMD1", "CTDSPL", "MAB21L1", "DNAJA4", "POLR1E", "TLL1", 
                "CEBPA", "PLEKHH3", "WIPF1", "WDR86", "SYNE1", "TULP4", "PIK3C2G", "TULP1", "COCH", "ADAMTS5", "THBD", 
                "ADAMTS3", "HECTD2", "ADAMTS1", "ADAMTS2", "HECTD4", "ADAMTS7", "ADAMTS8", "SRRM3", "OTOP1", "P3H1", "FRAS1", 
                "CYP2U1", "ZNF839", "IL23A", "MTHFD2L", "CTH", "BCOR", "NBL1", "AMER2", "STXBP1", "STXBP5", "LMF1", "AQP5", 
                "AQP3", "HS6ST2", "HS6ST3", "LMF2", "STRIP2", "NIPAL4", "STXBP6", "RPRML", "JAZF1", "DRGX", "LONP1", "HLA-A", 
                "CBFA2T3", "TRERF1", "WBSCR17", "KLHDC8A", "SIM2", "VGLL2", "SIM1", "CRHR1", "MTRR", "CASZ1", "GFI1", 
                "SLC35D3", "TWIST2", "CRHR2", "ZNHIT3", "WBSCR28", "DZIP1L", "TNRC18", "GSC2", "FEZF2", "CCR10", "EDN1", 
                "STRADB", "HEG1", "LHFPL3", "SNAP91", "MTSS1", "CLDN11", "P4HA2", "VWA2", "RGS7BP", "HMX1", "MEGF8", "HMX3", 
                "HMX2", "RTN4R", "TENM4", "CITED1", "FAM43A", "CITED4", "MEOX2", "CHRAC1", "DUOX2", "SH3PXD2A", "MLEC", 
                "FAM43B", "RAB6B", "KIAA1804", "PSD", "MPP3", "VWC2", "PROX1", "SORBS1", "ZDHHC14", "GP1BB", "PPHLN1", 
                "GFRA2", "VASH2", "GFRA1", "VAX1", "THOC6", "GFRA4", "VAX2", "PARP8", "LIMA1", "CYP2W1", "NBR1", "SLC22A23", 
                "HIP1", "PTGER3", "PTGER4", "GDPD5", "PODNL1", "CLDN23", "EGFLAM", "LMO2", "CPSF1", "WNT5A", "CPSF4", 
                "KLHDC1", "SHISA6", "ALG12", "LIFR", "SHISA3", "GRHL2", "SHISA8", "MYCL", "SHISA2", "HRASLS5", "CHRNA5", 
                "ATP1A3", "NOL4", "NKD1", "EFNA5", "HOXA13", "HOXA11", "NT5E", "COMP", "HOXA10", "EFNB2", "SCUBE2", "NLRP5", 
                "EMC3", "BMPER", "NLRP6", "NLRP4", "GAS2L2", "WNT3A", "HELT", "MPST", "HMHA1", "EFNA2", "HGS", "PRPF8", 
                "EFNA3", "LOX", "TLX3", "TLX2", "CDH11", "COL6A3", "FAM163A", "TLX1", "THOP1", "MATN4", "FGFR3", "FGFR2", 
                "WNT2B", "LPL", "SPPL2B", "CDH20", "CDH22", "SVEP1", "CBX8", "SVOP", "CBX4", "CYHR1", "MYRIP", "PPP2R5D", 
                "SEPT5", "DAP", "RASA3", "MYH9", "FAM84A", "SOWAHA", "SOWAHB", "UCN", "MYH8", "MEGF11", "NOP2", "NRXN2", 
                "CMAHP", "ALOX12", "HAPLN4", "FAM83H", "LRRTM1", "SCFD2", "WNT11", "HRH3", "THPO", "CASP7", "TOP1MT", 
                "LONRF3", "KDM6B", "IL15RA", "KDM2A", "TCF7L2", "TCF7L1", "TMEM30B", "PDE10A", "PPP2R2C", "ALDH1A3", "IGDCC3", 
                "ALDH1A2", "PFDN4", "HLF", "COLEC12", "LTK", "SIX1", "SLC4A1", "SIX2", "SLC4A4", "ATP5G2", "OASL", "MAPK9", 
                "HLX", "SIX6", "CHSY3", "SIX3", "NRCAM", "COLGALT2", "MAPK6", "HES7", "RBFOX1", "RBFOX3", "CADM2", "CADM3", 
                "RFX4", "TSPEAR", "C20ORF24", "GRIN2C", "HPCAL4", "PSMB9", "GRIN3B", "OTX2-AS1", "ENTHD1", "NGEF", "LEF1", 
                "STAU2", "RETSAT", "FZD10", "SLC5A1", "CLDN3", "PRAP1", "PREX1", "MICALL2", "HES2", "HES3", "PPARGC1B", 
                "HES6", "HES5", "ABCA2", "ABCA3", "PROSER1", "SRMS", "EFCC1", "FIBCD1", "ABCA4", "KLF10", "VPS37B", "FLI1", 
                "AR", "LATS2", "ABTB2", "SETX", "PRKAR1B", "KCNMA1", "LYPD6B", "SNORA43", "DHH", "C5ORF28", "SLC30A10", 
                "TRIOBP", "AAAS", "LIPG", "GRIN2B", "LRAT", "ZBTB7C", "RHBDL3", "LRRC26", "LYN", "LRRC24", "RIPK4", "LRRC29", 
                "SAMD4A", "LETM1", "GPRC5C", "TTLL6", "GPRC5B", "UBAP2", "HRK", "GPR27", "WFIKKN1", "PTH1R", "LITAF", "LRRC7", 
                "LRSAM1", "EXTL1", "RIPK2", "ATXN1", "MEF2D", "GSX2", "ARG2", "GSX1", "PICK1", "ELOVL4", "SPSB4", "SPSB1", 
                "BDH2", "TMEM56", "DSCAML1", "ZNF213", "NCAN", "TDRP", "JPH1", "FSTL3", "SLC9C1", "RCSD1", "FSTL4", "FOXQ1", 
                "JPH3", "ZFYVE28", "HEY1", "NKX1-2", "SRSF10", "CSMD1", "SLC16A2", "SRD5A2", "PRRX1", "ZBTB16", "GPR50", 
                "FOXP2", "BYSL", "SLC9A3", "OLFM2", "SLC9A5", "TFCP2", "CXCL12", "MAPKAPK5", "MAPKAPK3", "TUBGCP4", "SQSTM1", 
                "PTPN5", "SEC24D", "CCNT1", "HR", "ZBTB24", "FOXO6", "SEZ6", "CXCL14", "GGA1", "SLC9A2", "GBX1", "BCL2L11", 
                "GBX2", "HAS2", "HAS1", "FLNA", "FLNB", "FLNC", "SLC16A9", "PITX1", "NPTX1", "PITX3", "PITX2", "AICDA", 
                "C11ORF87", "SLC25A51", "FAM229A", "FOXN4", "MAF", "MLXIPL", "IL1RAPL1", "TROAP", "BIRC6", "CDC42EP1", 
                "LRRC3C", "ST6GALNAC6", "PVALB", "C11ORF96", "ZNF512B", "GPR88", "THADA", "TACSTD2", "RORA", "GLIS1", 
                "CYP19A1", "GLIS3", "STK10", "CTIF", "GAREM", "FBRSL1", "ATCAY", "EDIL3", "SLC38A6", "KL", "DGCR2", "EMX1", 
                "EMX2", "BDNF", "RAD23B", "FOXL2", "RAD23A", "FOXL1", "UST", "SNPH", "MDGA2", "SHOX2", "CRYBG3", "CNTNAP1", 
                "L3MBTL3", "LRRK1", "RASGRF1", "CELF4", "FOXK1", "RAP1GAP", "DST", "DHODH", "NELL2", "NELL1", "FLT1", 
                "PPP1R13B", "MDGA1", "FLT3", "PPM1M", "COL15A1", "SHCBP1L", "ZBTB46", "CXXC4", "PPM1L", "BARX1", "BARX2", 
                "PPM1E", "PPM1H", "CECR2", "CECR6", "SLITRK2", "PTPRZ1", "EPHB1", "LRRC8C", "SLC12A3", "EPHA6", "CCNJL", 
                "EPHA5", "EPHA8", "EPHA7", "LRRN4", "PLEKHA2", "SFMBT2", "PISD", "VWA5B2", "DCLK2", "MAPK11", "SFRP5", 
                "MAPK12", "SFRP2", "SPATA13", "KIF20A", "LRRN2", "EPHA4", "LRP2", "LRP1", "SATB2", "FOXG1", "NTN5", "TIFA", 
                "CCNF", "PPP1R3B", "NTN1", "SFRP1", "NTN4", "NCLN", "DHX37", "RAB11FIP3", "CALN1", "N4BP1", "PAQR9", "ARSB", 
                "XRCC5", "ST8SIA1", "FOXF2", "ST8SIA2", "FOXF1", "MIR9-3", "COMMD3", "HNF1B", "PPP1R1B", "WNK2", "WNK4", 
                "INTS7", "RGL3", "VSX2", "CCNO", "LPGAT1", "SPON1", "FOXE3", "DOC2B", "EPHA10", "FOXE1", "PCDH10", "PCDH18", 
                "PCDH17", "RASSF2", "PGBD5", "RCC1", "ARSI", "TSPAN5", "HIVEP2", "KIF21A", "KIF21B", "CEP170B", "SCARF2", 
                "FOXD4", "FOXD3", "WARS", "FOXD2", "FOXD1", "SORCS3", "RGMA", "NR5A1", "TG", "ID3", "MKX", "ID4", "SFI1", 
                "CD24", "STK40", "PTPRU", "BARHL2", "DOT1L", "FOXC2", "PTPRT", "FOXC1", "AATK", "PTPRM", "IKBKAP", "SLC7A14", 
                "GPAA1", "PTPRJ", "BMI1", "KIF17", "SOCS1", "PENK", "KIF1A", "STK39", "ST8SIA6", "CCR9", "SOCS7", "PAG1", 
                "FOXB2", "FOXB1", "RFTN1", "GCGR", "MIAT", "MN1", "ANO1", "PER2", "DPY19L3", "NLRP11", "KIF26A", "KIF26B", 
                "PER1", "FAM171B", "EPPK1", "TMEM151B", "POPDC3", "FOXA2", "TMEM151A", "FOXA1", "DIRAS1", "FHOD1", "DIRAS2", 
                "SLC43A1", "CRYBA2", "CAMKK1", "CAMK2N1", "NCS1", "EP300", "HEATR3", "EOMES", "PDGFRA", "KIRREL3", "NPNT", 
                "OSR2", "OSR1", "PDAP1", "CYP24A1", "COL5A1", "WBP1L", "SUPT6H", "CD83", "SLC1A2", "IHH", "IBA57", "ATOH8", 
                "PROKR2", "CCND2", "PPME1", "DUOXA1", "DUOXA2", "SCN7A", "GPC1", "GPC3", "HSF4", "GPC6", "BANF2", "NXPH1", 
                "CAPN2", "ARHGEF40", "PAK6", "ATOH1", "FDPS", "ADPRHL2", "SLC30A3", "ARHGEF17", "DDX54", "EVX1", "EVX2", 
                "INHBB", "PROC", "APPBP2", "MIDN", "ARHGEF26", "ARHGEF28", "JADE2", "CHAT", "JADE1", "KLHL14", "PPP1R9B", 
                "RASGRP2", "RGS6", "CALHM1", "CALCR", "CAMTA1", "MSC", "A4GALT", "SBF1", "ARID2", "VARS2", "ZNF385A", 
                "HS3ST3B1", "ZFHX4", "GALNT6", "GALNT3", "SLC52A3", "HMGA2", "SLC52A2", "SSTR1", "LAMB2", "BMP7", "BMP6", 
                "BMP4", "VCAN", "BMP3", "BMP2", "BMP1", "NDUFS7", "IRF5", "SGK223", "IRF2", "BRMS1L", "IRF8", "SCN4B", 
                "QRICH2", "HS3ST3A1", "MAST1", "YWHAG", "RTN1", "DDX3X", "NPR1", "NPR2", "NPR3", "LAMC1", "CHD3", "FUT4", 
                "AFF3", "WNT6", "PDLIM2", "CNPPD1", "UBE2Q2", "BCL7A", "WNT1", "PRPH", "WNT3", "GALNT9", "CCDC177", "LMTK3", 
                "PCOLCE", "ATP2B3", "FNDC3B", "THAP1", "TBC1D1", "PTHLH", "TBC1D4", "NPHS2", "VGF", "ABHD17A", "FAM131B", 
                "TAF8", "INA", "TAF5", "DDX27", "LAMA5", "DDX23", "PAPL", "DBX1", "ACVR1B", "CLCN2", "SOBP", "EFS", "EVA1C", 
                "MAT2A", "SCN5A", "MGC50722", "SCARB1", "CHGB", "UNC13D", "NKX2-8", "SOX14", "SOX13", "SOX11", "DKK2", 
                "NISCH", "CPEB2", "CPEB1", "ASB8", "TNFSF9", "DUS3L", "CDK14", "C9ORF50", "TRRAP", "HOXB1", "COL12A1", "RAX", 
                "TPGS1", "DCAF5", "BACH2", "ASAH1", "HOXA9", "MYB", "CHST12", "HENMT1", "GPR150", "NKX2-2", "HOXA3", "NKX2-1", 
                "HOXA2", "CYP1B1", "CHST15", "ANKS1B", "NKX2-6", "ANKS1A", "HOXA7", "NKX2-5", "HOXA6", "NKX2-4", "HOXA5", 
                "NKX2-3", "FZD2", "FZD1", "ZNF143", "WNT10B", "FZD4", "WNT10A", "TIGD5", "COL25A1", "DNAH17", "FZD8", "GPN3", 
                "BACE2", "RIPPLY2", "FHDC1", "RIPPLY3", "BIN1", "HOXB5", "HOXB4", "HOXB3", "PRR18", "HOXB2", "HOXB9", 
                "DNAH10", "RLTPR", "HOXB8", "PRR15", "HOXB7", "HOXB6", "TRIO", "HOXD3", "TRIL", "TMEM181", "SOX21", "RXFP3", 
                "SLC6A4", "SLC6A5", "GPR176", "SOX18", "GALR3", "GPR173", "GALR2", "SOX17", "GALR1", "HOXC6", "HOXC5", 
                "HOXC4", "IER5", "PRIMA1", "HOXC9", "HOXC8", "ANXA2", "TESC", "PAX5", "TET3", "PAX7", "PAX6", "PRSS12", 
                "PAX1", "TRMT1", "TTC9B", "PAX3", "PAX2", "TMEM179", "PAX9", "PAX8", "ITPKA", "ITPKB", "ANKRD63", "HOXD9", 
                "HOXD8", "RET", "ADM", "RASAL3", "TMEM163", "TRIM9", "SYCE3", "SNTG1", "CHL1", "RNF150", "EN1", "GPR1", "EN2", 
                "WWP2", "TTC7A", "SLC2A10", "SRP68", "AES", "COL4A1", "CLIP1", "NMNAT2", "COL4A2", "NOTUM", "ZNF579", "PLCB2", 
                "SF3A1", "PHLPP2", "HOXD13", "PRSS23", "HOXD12", "GLI1", "HOXD11", "PRSS21", "MCTS2P", "HOXD10", "LFNG", 
                "RSPO1", "RIOK3", "RSPO2", "BRD8", "AOC2", "CD160", "F12", "SEMA4D", "MITF", "TEX2", "DISP2", "SYPL2", 
                "LZTR1", "SEMA5B", "CHRM3", "DGKG", "SERTM1", "RNA45S5", "SEMA5A", "PHLDB1", "UHRF1BP1L", "RNF38", "CTNND2", 
                "DGKA", "METRNL", "C1ORF115", "JAKMIP1", "CDH8", "DACH1", "CDH4", "CDH3", "CDH2", "GRB14", "GRB10", "MISP", 
                "NKX6-2", "NKX6-1", "ESRRB", "SEMA6D", "ANKRD34B", "PNISR", "SEMA6A", "SEMA6B", "IGFBP5", "ATP8B1", "KREMEN2", 
                "MAFB", "MAFA", "ADRA2C", "ABHD15", "ADRA2B", "DUSP26", "ADRA2A", "CACNB2", "CACNB3", "CACNB4", "TRAF4", 
                "ZNF536", "ANKRD34A", "EIF4G3", "CNGB1", "DGKI", "ANKRD33B", "MYOM2", "DPP10", "NOTCH1", "C3ORF70", "IRS1", 
                "CLSTN1", "CLSTN2", "CHRD", "AK5", "IRS4", "RHOBTB1", "PHF21B", "BAG6", "CEP350", "EFHD1", "ZNF521", "PTPN18", 
                "CBLN1", "KCNIP1", "DNMT3A", "DNMT3B", "TXNRD1", "CDC42BPB", "PTPN14", "ZZEF1", "CDK8", "RTN4RL1", "RTN4RL2", 
                "DMRTA2", "ESRP1", "DMRTA1", "NES", "ZNF513", "VSTM2B", "VSTM2A", "SERINC5", "SOS1", "BNC1", "ITPR3", 
                "ADRA1D", "ADRA1A", "MROH1", "HMBOX1", "MGAT5B", "RREB1", "ZNF503", "AMH", "IRX2", "SKOR1", "IRX3", "LIMCH1", 
                "IRX1", "IRX6", "VWF", "RPH3AL", "IRX4", "IRX5", "WRAP73", "NFASC", "MGAT4B", "INSM1", "INSM2", "SLC24A2", 
                "SMARCD3", "SLC24A4", "SEMA3F", "CHADL", "STK32A", "SOST", "CACNG2", "STK32C", "INSL3", "UGGT2", "CACNG4", 
                "DCAF15", "EBF1", "EBF2", "RUNX1", "MRPS18B", "EBF3", "EBF4", "PUS7L", "AFG3L2", "NIN", "RRAGD", "TBKBP1", 
                "TMCC3", "APP", "VIPR2", "WBP2", "NDNF", "ANTXR2", "ITLN1", "NHLH2", "ADAD2", "BCL11B", "SLC11A2", "BCL11A", 
                "SIAH2", "GADD45G", "PAPPA", "AVPR1A", "OVOL2", "GDF7", "TUBG2", "OVOL1", "NUAK1", "GDF6", "DOPEY2", "PARP12", 
                "OCLN", "ZFPM1", "NCEH1", "SPATS2L", "TSHZ1", "TTLL12", "ZFPM2", "NECAB2", "CDS1", "GPM6B", "OXTR", "HGH1", 
                "TSHZ3", "CTDSPL2", "TMTC1", "SLC6A20", "SOX3", "TNS1", "SOX4", "SOX1", "CYP27B1", "SOX2", "DPF3", "FFAR4", 
                "SOX9", "SOX7", "SOX8", "SOX5", "TNS2", "TSC22D3", "TFEB", "MAPT", "RDH5", "CAMK4", "FAM181B", "MIR5100", 
                "ACVRL1", "CLIC6", "OPRD1", "ADRB1", "TWF1", "VLDLR", "NEUROD2", "NKIRAS2", "PRKCDBP", "MAP6", "DISC1", 
                "ASIC2", "MAP4K3", "NOG", "DMRT2", "PRMT2", "DMRT3", "DMRT1", "AP3D1", "LBX2", "LBX1", "B3GNT8", "B3GNT6", 
                "STAC2", "LLGL1", "NEUROG1", "NEUROG2", "LLGL2", "NEUROG3", "CDK5R2", "RYR2", "FAM69C", "CTNNB1", "PRIM1", 
                "CILP2", "KIF13A", "CXCR4", "NPY", "CTNNA3", "TGFB1", "TGFB2", "GDF15", "CDX2", "B3GAT2", "ARID3C", 
                "ARHGAP27", "POU4F2", "POU4F3", "GDF10", "LIN28A", "ZNFX1", "RGS9BP", "TGFBI", "ABCG4", "GALNT12", "GALNT14", 
                "TNXB", "GALNT16", "RGS17", "NRK", "CASKIN1", "MYH10", "S1PR5", "PDE8A", "GHSR", "RGS20", "NSF", "OR10K2", 
                "STIL", "LRRC75B", "USP44", "EPHX3", "EPHX4", "NAALADL1", "TDRD6", "SMAD6", "ESR2", "NR4A2", "ARHGAP10", 
                "GDNF", "RADIL", "DOCK9", "ESPN", "LTBP1", "CSF2RA", "SEL1L3", "LTBP2", "SLC2A7", "NR4A1", "ARHGAP22", 
                "ARHGAP20", "PODXL", "ZCWPW1", "TMEM108", "SIRPA", "HSPA5", "HSPA2", "SENP1", "TGFBR2", "TGFBR3", "FES", 
                "MEIS1", "FEV", "PDE3B", "PDE3A", "DOCK1", "NCAM1", "AGAP2", "NCALD", "OTX1", "AOX1", "PDE4C", "PDE4A", 
                "NCAM2", "OTX2", "ZC3H7A", "CCBE1", "CAMK1D", "NEGR1", "KCTD1", "ARID1B", "DLK1", "SGPP2", "DLK2", "KCTD8", 
                "SGPP1", "LMX1B", "LMX1A", "HKDC1", "RGS10", "SSBP2", "B4GALT7", "HCN2", "KRT83", "CAMK2D", "CAMK2B", "ACAN", 
                "IGFBPL1", "RMND5B", "LRIG1", "SCYL1", "KLRG2")

# ---- Promoters ----

## Find location of promoters of these features and convert to correct format

pro <- GenomicFeatures::promoters(genes, downstream=0, upstream = 2000)

pro.sub <- pro[(names(pro) %in% features)]

## Find overlaps between these promoter peaks and the detected peaks

pro.ovrlap <- findOverlaps(query = signac, subject = pro.sub)

promoter.peaks <- queryHits(pro.ovrlap)

pro.ovrlap.peaks <- signac@assays[["peaks"]]@ranges[promoter.peaks]

## Convert from Granges object to a names list of strings to input to AddChromatinModule

pro.string <- GRangesToString(pro.ovrlap.peaks, sep = c("-", "-"))

pro.string <- list(paste(module.name, ".Promoters", sep = "") = pro.string)

## Add module score for these features to the seurat object

signac <- AddChromatinModule(signac, features = pro.string, genome = BSgenome.Hsapiens.UCSC.hg38, assay = "peaks", verbose = TRUE)

# ---- Enhancers ----
 
## Read in EnhancerAtlas data, from interactions between Enhancers and Genes for the chosen cell type

enhancers <- read.table("/media/sf_G_DRIVE/Shared drives/Crean Group Data Share/Kidney Organoids ATAC 2_14_03_21/H9_EP.txt", header = TRUE)

## Subset the interactions based on the involved genes and the genes of interest

enh.sub <- subset(enhancers, Gene %in% features)

enh.ranges <- StringToGRanges(as.list(enh.sub[1]), sep = c(":", "-"))

## FInd overlap between these enhancer areas and the detected peaks

enh.ovrlap <- findOverlaps(query = signac, subject = enh.ranges)

enh.peaks <- queryHits(enh.ovrlap)

enh.ovrlap.peaks <- signac@assays[["peaks"]]@ranges[enh.peaks]

## Convert from Granges object to list of strings

enh.string <- GRangesToString(enh.ovrlap.peaks, sep = c("-", "-"))

enh.string <- list(paste(module.name, ".Promoters", sep = "") = enh.string)

## Add the module to the seurat object

signac <- AddChromatinModule(signac, features = enh.string, genome = BSgenome.Hsapiens.UCSC.hg38, assay = "peaks", verbose = TRUE)

VlnPlot(signac, "Matrisome.Promoters", split.by = "orig.ident", idents = c("Stroma 2", "Stroma 3"))

FeaturePlot(signac, "Matrisome.Promoters",  split.by = "orig.ident")

VlnPlot(signac, "Matrisome.Enhancers", split.by = "orig.ident", idents = c("Stroma 2", "Stroma 3")) 

FeaturePlot(signac, "Matrisome.Enhancers",  split.by = "orig.ident")

# ---- Plotting ----

signac.sub <- subset(signac, idents = c("Stroma 2", "Stroma 3"))

dodge <- position_dodge(width = 0.4)

m.pro <- ggplot(FetchData(signac.sub, vars = c("Matrisome.Promoters", "orig.ident", "clusters")), aes(x = clusters, y = Matrisome.Promoters, fill = orig.ident))+
  scale_fill_viridis_d( option = "D")+
  ggbeeswarm::geom_quasirandom(shape = 21, size=1.5, dodge.width = .75, color="black",alpha=1,show.legend = F)+
  geom_violin(alpha=0.6, position = position_dodge(width = .75),size=1, color=NA) +
  geom_boxplot(notch = TRUE, position = position_dodge(width = .75), outlier.size = -1, color="black",lwd=1, alpha = 1,show.legend = F, width=.15)+
  # geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
  theme_minimal()+
  rremove("legend.title")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.ticks = element_line(size=1,color="black"),
    axis.text = element_text(color="black"),
    axis.ticks.length=unit(0.2,"cm"),
    #legend.position = c(0.95, 0.85),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank())+
  # font("xylab",size=15)+  
   font("xy",size=0)+ 
   font("xy.text", size = 12) +  
  # font("legend.text",size = 15) +
  guides(fill = guide_legend(override.aes = list(alpha = 1,color="black"))) +
  #guides(fill = "none") +
  ggtitle("Matrisome Promoter Score")

m.enh <- ggplot(FetchData(signac.sub, vars = c("Matrisome.Enhancers", "orig.ident", "clusters")), aes(x = clusters, y = Matrisome.Enhancers, fill = orig.ident))+
  scale_fill_viridis_d( option = "D")+
  ggbeeswarm::geom_quasirandom(shape = 21, size=1.5, dodge.width = .75, color="black",alpha=1,show.legend = F)+
  geom_violin(alpha=0.6, position = position_dodge(width = .75),size=1, color=NA) +
  geom_boxplot(notch = TRUE, position = position_dodge(width = .75), outlier.size = -1, color="black",lwd=1, alpha = 1,show.legend = F, width=.15)+
  # geom_point( shape = 21,size=2, position = position_jitterdodge(), color="black",alpha=1)+
  theme_minimal()+
  rremove("legend.title")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.ticks = element_line(size=1,color="black"),
        axis.text = element_text(color="black"),
        axis.ticks.length=unit(0.2,"cm"),
        #legend.position = c(0.95, 0.85),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
   #font("xylab",size=12)+  
   font("xy",size=0)+ 
   font("xy.text", size = 12) +  
  # font("legend.text",size = 15) +
  guides(fill = guide_legend(override.aes = list(alpha = 1,color="black"))) +
  #guides(fill = "none") +
  ggtitle("Matrisome Enhancer Score")

wrap_plots(m.pro / m.enh)

ggsave(paste(save.dir, module.name, "_Scores.png", sep = ""), width = 15, height = 10, units = "in",
       wrap_plots(m.pro / m.enh, guides = "collect"))
ggsave(paste(save.dir, module.name, "_Scores.svg", sep = ""), width = 15, height = 10, units = "in",
       wrap_plots(m.pro / m.enh, guides = "collect"))


# ---- Other plotting ----

ggplot(FetchData(subset(signac, idents = c("Stroma2", "Stroma 3") vars = c("Matrisome.Promoters", "orig.ident", "clusters")), aes(x = clusters, y = Matrisome.Promoters, fill = orig.ident)) +
  geom_violin(scale = 'area', trim = FALSE) +
  geom_boxplot(width=.1, outlier.shape = NA) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(signac$orig.ident))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Matrisome Promoter Score') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) 

ggplot(FetchData(signac, vars = c("Matrisome.Enhancers", "orig.ident", "")), aes(x = orig.ident, y = Matrisome.Enhancers, fill = orig.ident)) +
  geom_violin(scale = 'area', trim = FALSE) +
  geom_boxplot(width=.1, outlier.shape = NA) +
  theme_bw() +
  scale_fill_manual(values = custom_colors$discrete) +
  scale_x_discrete(limits = rev(levels(signac$orig.ident))) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = 'Matrisome Enhancer Score') +
  theme(
    axis.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = 'none'
  ) 
