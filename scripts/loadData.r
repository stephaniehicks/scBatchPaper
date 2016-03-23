#### Load libraries 
library(Biobase)
library(GenomicRanges)
library(vcd)
library(stringr)
library(reshape2)
library(tidyr)
library(plyr)
library(dplyr)
library(irlba)


#### Load Deng et al. (2014) RPKM and phenotypic data
library(scRNASeqMouseDengMonoAllelic)
data("embryoMouseRPKM")

eset <- exprs(embryoMouseRPKM)
pd <- pData(embryoMouseRPKM)

# calculate PDG
pd$PDG <- 1 - (apply(eset, 2, function(x){ length(which(x == 0)) }) / nrow(eset))
pdDeng <- pd
eDeng <- log2(eset + 1) # log transform RPKMs

# group by time
zz <- levels(pdDeng$time)
pdDeng$bin <- ifelse(pdDeng$time %in% zz[1:2], "Bin A", 
                     ifelse(pdDeng$time %in% zz[3:4], "Bin B", 
                            ifelse(pdDeng$time %in% zz[5:7], "Bin C", "Bin D")))

# calculate SVD
dat <- sweep(eDeng[, pdDeng$bin == "Bin A"], 1, rowMeans(eDeng[, pdDeng$bin == "Bin A"]), FUN = "-")
sDeng.A <- svd(dat)
dat <- sweep(eDeng[, pdDeng$bin == "Bin B"], 1, rowMeans(eDeng[, pdDeng$bin == "Bin B"]), FUN = "-")
sDeng.B <- svd(dat)
dat <- sweep(eDeng[, pdDeng$bin == "Bin C"], 1, rowMeans(eDeng[, pdDeng$bin == "Bin C"]), FUN = "-")
sDeng.C <- svd(dat)
dat <- sweep(eDeng[, pdDeng$bin == "Bin D"], 1, rowMeans(eDeng[, pdDeng$bin == "Bin D"]), FUN = "-")
sDeng.D <- svd(dat)

print("[Data loaded]: Deng et al. (2014)")



#### Load Jaitin et al. (2014) barcode data and phenotypic data
library(scRNASeqMouseJaitinSpleen)
data("spleenMouseMARSSeq")

eset <- exprs(spleenMouseMARSSeq)
pd <- pData(spleenMouseMARSSeq)

# calculate PDG
pd$PDG <- 1 - (apply(eset, 2, function(x){ length(which(x == 0)) }) / nrow(eset))
pd$batch <- factor(pd$batch)
levels(pd$batch) <- paste("Batch", 1:44)

# remove cells with a library size of 0 and outlier cells
keepMe <- (colSums(eset) != 0) & (pd$PDG < 0.15)
pd <- pd[keepMe, ]
eset <- eset[, keepMe]

# normalize by total number of barcoded molecules and multiple by 1e4 (see drop-seq paper)
eset = sweep(eset, 2, colSums(eset)/1e4, FUN = "/")

pdJaitin <- pd
eJaitin <- log2(eset + 1) # log transform barcoded molecules

# Compute first 3 PCs (approximate) using irlba pkg
dat <- sweep(eJaitin, 1, rowMeans(eJaitin), FUN = "-")
sJaitin <- irlba(dat, nv = 3)

print("[Data loaded]: Jaitin et al. (2014)")



#### Load Kumar et al. (2014) TPM and phenotypic data
library(scRNASeqMouseKumarPSC)
data("pscKumarMouseTPM")

pd <- pData(pscKumarMouseTPM)
eset <- exprs(pscKumarMouseTPM)

# calculate PDG
pd$PDG <- 1 - (apply(eset, 2, function(x){ length(which(x == 0)) }) / nrow(eset))
pd$batch <- factor(paste(pd$instrument, pd$runID, pd$fcID, pd$fcLane, sep="_"))
pd$batch <- factor(pd$batch, levels(pd$batch)[c(3,1,4,2)])
levels(pd$batch) <- paste("Batch", 1:4)
pd$culture <- factor(str_sub(pd$characteristics_ch1.1, start = 21))
pd$group <- paste(str_sub(pd$culture, start = 21), str_sub(pd$characteristics_ch1, start=12, end = 31), sep=" ")
pd$mouse <- pd$source_name_ch1

# group by mouse and condition
keepDgcr8 <- (pd$mouse != "v6.5 mouse embryonic stem cells")
keepv6.5serum <- (pd$mouse == "v6.5 mouse embryonic stem cells" & pd$culture == "serum+LIF")
keepv6.52i <- (pd$mouse == "v6.5 mouse embryonic stem cells" & pd$culture != "serum+LIF")
pd$bin <- ifelse(keepDgcr8, "Group A", ifelse(keepv6.5serum, "Group B", "Group C"))

pdKumar <- pd
eKumar <- eset # already on log(TPM + 1) scale

# calculate SVD 
dat <- sweep(eKumar[, pdKumar$bin == "Group A"], 1, rowMeans(eKumar[, pdKumar$bin == "Group A"]), FUN = "-")
sKumar.A <- svd(dat)
dat <- sweep(eKumar[, pdKumar$bin == "Group B"], 1, rowMeans(eKumar[, pdKumar$bin == "Group B"]), FUN = "-")
sKumar.B <- svd(dat)
dat <- sweep(eKumar[, pdKumar$bin == "Group C"], 1, rowMeans(eKumar[, pdKumar$bin == "Group C"]), FUN = "-")
sKumar.C <- svd(dat)

print("[Data loaded]: Kumar et al. (2014)")



#### Load Patel et al. (2014) TPM and phenotypic data
library(scRNASeqHumanPatelGlioblastoma)

data("glioHumanCount")
eset <- assay(glioHumanCount)
pd <- colData(glioHumanCount)
PDG <- 1 - (apply(eset, 2, function(x){ length(which(x == 0)) }) / nrow(eset))

data("glioHumanTPM")
eset <- exprs(glioHumanTPM)
pd <- pData(glioHumanTPM)
pd$PDG <- PDG[match(rownames(pd), names(PDG))]

ePatel <- eset[, pd$sampleType == "SC"]
pdPatel <- pd[pd$sampleType == "SC", ]

# calculate SVD
sPatel <- svd(ePatel) # already row standardized 

print("[Data loaded]: Patel et al. (2014)")



#### Load Shalek et al. (2014) TPM and phenotypic data
# ```{r loadData-Shalek2014}
library(scRNASeqMouseShalekDendritic)
data("dendriticMouseTPM")

eset <- exprs(dendriticMouseTPM)
pd <- pData(dendriticMouseTPM)
keepIDs <- (grepl(pattern = "BMDC \\(([0-9]+h)", pd$source_name_ch1) |  
                grepl(pattern = "(Unstimulation)", pd$source_name_ch1)) & 
    (!grepl(pattern = "IFN-B", pd$source_name_ch1) & 
         !grepl(pattern = "StimulationReplicate Experiment)", pd$source_name_ch1))

pd <- pd[keepIDs, ]
eset <- eset[, keepIDs]

pd$source_name_ch1 <- factor(pd$source_name_ch1)
pd$stim <- ifelse(!grepl(pattern = "(Unstimulation)", pd$source_name_ch1), 
                  str_sub(unlist(lapply(str_split(pd$source_name_ch1, " ", n = 4), 
                                        function(x) x[4])), start = 1, end = -2),
                  str_sub(unlist(lapply(str_split(pd$source_name_ch1, " ", n = 4), 
                                        function(x) x[2])), start = 2, end = -2))
pd$time <- ifelse(grepl(pattern = "(Unstimulation)", pd$source_name_ch1), NA,
                  str_sub(unlist(lapply(str_split(pd$source_name_ch1, " ", n = 2), 
                                        function(x) x[2])), start = 2, end = -18))
pd$cond <- ifelse(grepl(pattern = "(Unstimulation)",
                        pd$source_name_ch1), NA,
                  str_sub(unlist(lapply(str_split(pd$source_name_ch1, " ", n = 2), 
                                        function(x) x[2])), start = 5, end = -13))

# calculate PDG
pd$PDG <- 1 - (apply(eset, 2, function(x){ length(which(x == 0)) }) / nrow(eset))

# Subset for only the LPS experimental condition 
keepIDs <- grepl(pattern = "^LPS_([0-9]+h)_S", pd$title) 
eShalek <- log2(eset[, keepIDs] + 1) # log transform FPKMs
pdShalek <- pd[keepIDs, ]

pdShalek$batch <- factor(paste(pdShalek$runID, pdShalek$fcLane, sep = "_"))
levels(pdShalek$batch) <- paste0("Batch ", 1:4) 

# calculate SVD 
dat <- sweep(eShalek, 1, rowMeans(eShalek), FUN = "-")
sShalek <- svd(dat)

print("[Data loaded]: Shalek et al. (2014)")


#### Load Trapnell et al. (2014) FPKM and phenotypic data
library(scRNASeqHumanTrapnellMyoblast)
data("myoblastHumanFPKM")

pd <- pData(myoblastHumanFPKM)
eset <- exprs(myoblastHumanFPKM)

# remove control cells, debris cells, and wells with more than 1 cell
pd <- pd[pd$sampleType == "SC" & pd$control == "control well: FALSE" & 
             pd$debris == "debris: FALSE" & 
             pd$numcells == "cells in well: 1", ]
eset <- eset[, match(rownames(pd), colnames(eset))]

pd$batch <- factor(paste(pd$runID, pd$fcLane, sep = "_"))
levels(pd$batch) <- paste0("Batch ", c(2, 1, 3:4))
pd$batch <- factor(pd$batch, levels = levels(pd$batch)[c(2,1,3,4)])
pd$hour <- factor(pd$hour, levels = levels(pd$hour)[c(2,3,4,1)])
levels(pd$hour) <- (c("0h", "24h", "48h", "72h"))


# calculate PDG
pd$PDG <- 1 - (apply(eset, 2, function(x){ length(which(x == 0)) }) / nrow(eset))
pdTrapnell <- pd
eTrapnell <- log2(eset + 1) # log transform FPKMs

# calculate SVD 
dat <- sweep(eTrapnell, 1, rowMeans(eTrapnell), FUN = "-")
sTrapnell <- svd(dat)

print("[Data loaded]: Trapnell et al. (2014)")


#### Load Treutlein et al. (2014) FPKM and phenotypic data
library(scRNASeqMouseTreutleinLineage)
data("lungMouseFPKM")

pd <- pData(lungMouseFPKM)
eset <- exprs(lungMouseFPKM)

# remove bulk and no cell
pd <- pd[pd$sampleType == "SC", ]
eset <- eset[, match(rownames(pd), colnames(eset))]

pd$batch <- factor(paste(pd$runID, pd$fcLane, sep = "_"))
levels(pd$batch) <- c("Batch 4", "Batch 5", "Batch 6", 
                          "Batch 1", "Batch 2", "Batch 7", "Batch 3")
pd$batch <- factor(pd$batch, levels = levels(pd$batch)[c(4,5,7,1:3,6)])

# calculate PDG
pd$PDG <- 1 - (apply(eset, 2, function(x){ length(which(x == 0)) }) / nrow(eset))
pdTreutlein <- pd
eTreutlein <- log2(eset + 1) # log transform FPKMs
levels(pdTreutlein$day) <- c(levels(pdTreutlein$day)[1:3], "Adult")

# calculate SVD
dat <- sweep(eTreutlein, 1, rowMeans(eTreutlein), FUN = "-")
sTreutlein <- svd(dat)

print("[Data loaded]: Treutlein et al. (2014)")



#### Load Bose et al. (2015) UMI and phenotypic data
library(scRNASeqHumanBosePrinting)
data("printHumanUMI_PS041")

eset <- exprs(printHumanUMI_PS041)
pd <- pData(printHumanUMI_PS041)

# calculate PDG
pd$PDG <- 1 - (apply(eset, 2, function(x){ length(which(x == 0)) }) / nrow(eset))
pd$textFile <- laply(str_split(pd$title, "_"), function(x){ x[1] })
pd$source_name_ch1 <- factor(pd$source_name_ch1)
pd$batch <- factor(pd$description)
levels(pd$batch) <- paste("Batch", 1:5)
pd$group <- pd$source_name_ch1
levels(pd$group)[1] <- c("mix U87 and WI-38")

eset = sweep(eset, 2, colSums(eset)/1e4, FUN = "/")

pdBose <- pd
eBose <- log2(eset + 1) # log transform FPKMs

# calculate SVD
dat <- sweep(eBose, 1, rowMeans(eBose), FUN = "-")
sBose <- svd(dat)

print("[Data loaded]: Bose et al. (2015)")



#### Load Burns et al. (2015) TPM and phenotypic data
library(scRNASeqMouseBurnsInnerEar)
data("innerEarMouseTPM")

eset <- exprs(innerEarMouseTPM)
pd <- pData(innerEarMouseTPM)

# calculate PDG
pd$PDG <- 1 - (apply(eset, 2, function(x){ length(which(x == 0)) }) / nrow(eset))
pd$tissue <- factor(pd$source_name_ch1)
levels(pd$tissue) <- c("cochlear", "utricular")
pd$fluidicChip <- factor(str_sub(pd$characteristics_ch1.1, start = 11))
pd$batch <- factor(paste(pd$instrument, pd$runID, pd$fcID, pd$fcLane, sep="_"))
levels(pd$batch) <- paste("Batch", 1:length(levels(pd$batch)))
pd$characteristics_ch1.4 <- str_sub(pd$characteristics_ch1.4, start = 18)

# remove bulk cells and cells purified with FACS
pd$bulk <- (pd$characteristics_ch1.4 %in% 
                c("Cochlear epithelium bulk population", "Downsample Outlier", 
                  "Utricular epithelium bulk population", "SINGuLAR Outlier", 
                  "Cochlear SC bulk population"))
pdBurns <- pd[!pd$bulk, ]
eBurns <- log2(eset[, !pd$bulk] + 1) # log transform FPKMs

# remove outlier due to FACS (not cochlear)
eBurns <- eBurns[, pdBurns$characteristics_ch1.4 != "FACs HC"]
pdBurns <- pdBurns[pdBurns$characteristics_ch1.4 != "FACs HC", ]

pdBurns$celltypespecific = tmp <- pdBurns$characteristics_ch1.4
pdBurns$celltype <- factor(ifelse(grepl("NSC", tmp), "NSC", 
                                  ifelse(grepl("HC", tmp), "HC", 
                                         ifelse(grepl("TEC", tmp), "TEC", "SC"))))
pdBurns$group = pdBurns$bin <- factor(paste(pdBurns$tissue, pdBurns$celltype, sep="_"))
pdBurns$groupspecific <- factor(paste(pdBurns$tissue, pdBurns$celltypespecific, sep="_"))
levels(pdBurns$bin) <- paste("Group", LETTERS[1:6])

# remove FACS purified samples
eBurns <- eBurns[, pdBurns$characteristics_ch1 != "purification: FACS"]
pdBurns <- pdBurns[pdBurns$characteristics_ch1 != "purification: FACS", ]


# calculate SVD 
dat <- sweep(eBurns[, pdBurns$tissue == "cochlear"], 1, rowMeans(eBurns[, pdBurns$tissue == "cochlear"]), FUN = "-")
sBurns.A <- svd(dat)
dat <- sweep(eBurns[, pdBurns$group == "utricular_HC"], 1, rowMeans(eBurns[, pdBurns$group == "utricular_HC"]), FUN = "-")
sBurns.B <- svd(dat)
dat <- sweep(eBurns[, pdBurns$group == "utricular_SC"], 1, rowMeans(eBurns[, pdBurns$group == "utricular_SC"]), FUN = "-")
sBurns.C <- svd(dat)
dat <- sweep(eBurns[, pdBurns$group == "utricular_TEC"], 1, rowMeans(eBurns[, pdBurns$group == "utricular_TEC"]), FUN = "-")
sBurns.D <- svd(dat)

print("[Data loaded]: Burns et al. (2015)")


#### Load Guo et al. (2015) FPKM and phenotypic data
library(scRNASeqHumanGuoGermCells)
data("germCellsHumanFPKM")

eset <- exprs(germCellsHumanFPKM)
pd <- pData(germCellsHumanFPKM)

# calculate PDG
pd$PDG <- 1 - (apply(eset, 2, function(x){ length(which(x == 0)) }) / nrow(eset))
pd$source_name_ch1 <- factor(pd$source_name_ch1)
pd$embryo <- unlist(lapply(str_split(as.character(pd$title),"_",n=5), function(x) x[4]))
pd$sex <- factor(str_sub(pd$characteristics_ch1.1, start = 9))
pd$week <- factor(unlist(lapply(str_split(as.character(pd$title),"_",n=5), function(x) str_sub(x[3],end = -2))), 
                  levels = c(4, 7, 8, 10, 11, 17, 19))
pd$weekGroup <- factor(ifelse(pd$week == "4", "4", ifelse(pd$week %in% c("7", "8"), "7-8", 
                                                          ifelse(pd$week %in% c("10", "11"), "10-11", "17-19"))), 
                       levels = c("4", "7-8", "10-11", "17-19"))
pd$group <- factor(paste(pd$weekGroup, pd$sex, pd$embryo, sep = "_"))
pd$group <- factor(pd$group, levels(pd$group)[c(9:15,1:8)])
pd$batch <- factor(paste(pd$instrument, pd$runID, pd$fcID, pd$fcLane, sep="_"))

# remove somatic cells
pd$celltype <- ifelse(pd$source_name_ch1 == "Primordial Germ Cells", "GermCell", "SomaticCell")
keepMe <- pd$celltype == "GermCell" & pd$weekGroup != "17-19"

pdGuo <- pd[keepMe, ]
eGuo <- log2(eset[,keepMe] + 1) # log transform FPKMs

pdGuo$batch <- factor(pdGuo$batch)
pdGuo$weekGroup <- factor(pdGuo$weekGroup)
pdGuo$week <- factor(pdGuo$week) 

# calculate SVD 
dat <- sweep(eGuo[, pdGuo$weekGroup == "4"], 1, rowMeans(eGuo[, pdGuo$weekGroup == "4"]), FUN = "-")
sGuo.A <- svd(dat)
dat <- sweep(eGuo[, pdGuo$weekGroup == "7-8"], 1, rowMeans(eGuo[, pdGuo$weekGroup == "7-8"]), FUN = "-")
sGuo.B <- svd(dat)
dat <- sweep(eGuo[, pdGuo$weekGroup == "10-11"], 1, rowMeans(eGuo[, pdGuo$weekGroup == "10-11"]), FUN = "-")
sGuo.C <- svd(dat)

print("[Data loaded]: Guo et al. (2015)")



#### Load Kowalczyk et al. (2015) TPM and phenotypic data
library(scRNASeqMouseKowalczykAging)
data("ageStrainC57BL6MouseTPM")

eset <- exprs(ageStrainC57BL6MouseTPM)
pd <- pData(ageStrainC57BL6MouseTPM)

# calculate PDG
pd$PDG <- 1 - (apply(eset, 2, function(x){ length(which(x == 0)) }) / nrow(eset))

pd$age <- str_sub(pd$characteristics_ch1.1, start = 6, end = 15)
pd$celltype <- str_sub(pd$characteristics_ch1.2, start = 12)
pd$celltype <- ifelse(pd$celltype == "short term hematopoietic stem cell",
                      "shortTermHSC", 
                      ifelse(pd$celltype == "long term hematopoietic stem cell", 
                             "longTermHSC", "multipotentProgenitor"))
pd$batch <- factor(paste(pd$instrument, pd$runID, pd$fcID, pd$fcLane, sep="_"))
pd$replicate <- ifelse(!grepl("replicate", pd$title), "rep1", "rep2")
pd$group <- factor(paste(pd$celltype, pd$age, pd$replicate, sep="_"))

# filter for short-term HSC
keepMe <- pd$celltype == "shortTermHSC" 

pdKowalczyk <- pd[keepMe, ]
eKowalczyk <- eset[, keepMe] # already on log scale

# calculate SVD
dat <- sweep(eKowalczyk, 1, rowMeans(eKowalczyk), FUN = "-")
sKowalczyk <- svd(dat)

print("[Data loaded]: Kowalczyk et al. (2015)")



#### Load Leng et al. (2015) TPM and phenotypic data
library(scRNASeqHumanLengOscillatoryGenes)
data("oscillatoryGenesTPM")

eset <- exprs(oscillatoryGenesTPM)
pd <- pData(oscillatoryGenesTPM)

# calculate PDG
pd$PDG <- 1 - (apply(eset, 2, function(x){ length(which(x == 0)) }) / nrow(eset))
pd$celltype <- ifelse(pd$source_name_ch1 == "single H1 hESC", "H1-hESC", 
                ifelse(pd$source_name_ch1 == 
                "single H1-Fucci cell sorted from G1 phase of the cell cycle only",
                "H1-Fucci-G1phase", ifelse(pd$source_name_ch1 == 
                        "single H1-Fucci cell sorted from G2/M phase of the cell cycle only", 
                        "H1-Fucci-G2/Mphase", "H1-Fucci-Sphase")))
pd$nPassages <- factor(str_sub(pd$characteristics_ch1.1, start = 11))
pd$sorted <- factor(str_sub(pd$characteristics_ch1.2, start = 12))
pd$batch <- factor(paste(pd$runID, paste0("L", pd$fcLane), sep="_"))

# remove two outlier cells
keepMe <- pd$PDG > 0.40
pdLeng <- pd[keepMe, ]
eLeng <- log2(eset[, keepMe] + 1) # log transform FPKMs

# calculate SVD
dat <- sweep(eLeng, 1, rowMeans(eLeng), FUN = "-")
sLeng <- svd(dat)

print("[Data loaded]: Leng et al. (2015)")



#### Load Macosko et al. (2015) UMI and phenotypic data
load("/net/irizarryfs01/srv/export/irizarryfs01/share_root/shicks/dataPackages/scRNASeqMouseMacoskoRetina/data/retinaMouseUMI.rda")

pdMacosko <- pData(retinaMouseUMI)
rm(retinaMouseUMI)

# eset <- exprs(retinaMouseUMI)
# # normalize by total number of UMI molecules and multiple by 1e4 (see drop-seq paper)
# tmp = colSums(eset)/1e4
# eset = sweep(eset, 2, tmp, FUN = "/")
# eMacosko <- log2(eset + 1) # log transform normalized UMIs
# save(eMacosko, file = "/net/irizarryfs01/srv/export/irizarryfs01/share_root/shicks/dataPackages/scRNASeqMouseMacoskoRetina/data/eMacosko.rda")
load("/net/irizarryfs01/srv/export/irizarryfs01/share_root/shicks/dataPackages/scRNASeqMouseMacoskoRetina/data/eMacosko.rda")

# calculate PDG
pdMacosko$PDG <- 1 - (apply(eMacosko, 2, function(x){ length(which(x == 0)) }) / nrow(eMacosko))
pdMacosko$retina <- factor(str_sub(pdMacosko$title, start = 11))

# Compute first 3 PCs (approximate) using irlba pkg
tmp = rowMeans(eMacosko)
dat <- sweep(eMacosko, 1, tmp, FUN = "-")
sMacosko <- irlba(dat, nv = 3)

print("[Data loaded]: Macosko et al. (2015)")



#### Load Satija et al. (2015) UMI and phenotypic data
library(scRNASeqDaniSajitaSeurat)
data("spatialZebrafishUMI")

eset <- exprs(spatialZebrafishUMI)
pd <- pData(spatialZebrafishUMI)

# calculate PDG
pd$PDG <- 1 - (apply(eset, 2, function(x){ length(which(x == 0)) }) / nrow(eset))
pd$group <- laply(str_split(pd$title, "_"), function(x){ x[1] })
pd$plate <- laply(str_split(pd$title, "_"), function(x){ x[2] })
pd$wellID <- laply(str_split(pd$title, "_"), function(x){ x[3] })
pd$batch <- factor(str_sub(pd$characteristics_ch1, start = 21))

# normalize by total number of barcoded molecules and multiple by 1e4 (see drop-seq paper)
eset = sweep(eset, 2, colSums(eset)/1e4, FUN = "/")

eSajita <- log2(eset + 1) # log transform normalized UMIs
pdSajita <- pd

# Compute first 3 PCs (approximate) using irlba pkg
dat <- sweep(eSajita, 1, rowMeans(eSajita), FUN = "-")
sSajita <- irlba(dat, nv = 3)

print("[Data loaded]: Satija et al. (2015)")



#### Load Zeisel et al. (2015) UMI and phenotypic data
library(scRNASeqMouseZeiselCortex)
data("cortexMouseUMI")

eset <- exprs(cortexMouseUMI)
pd <- pData(cortexMouseUMI)

# calculate PDG
pd$PDG <- 1 - (apply(eset, 2, function(x){ length(which(x == 0)) }) / nrow(eset))
pd$fluC1Run <- laply(str_split(pd$description, "_"), function(x){ x[1] })
pd$wellID <- laply(str_split(pd$description, "_"), function(x){ x[2] })
pd$brainRegion <- factor(ifelse(pd$source_name_ch1 == "ca1hippocampus", "hippocampus", "cortex"))
pd$gender <- factor(str_sub(pd$characteristics_ch1, start = 6))
pd$age <- factor(str_sub(pd$characteristics_ch1.1, start = 6))
pd$batch <- factor(paste(pd$runID, pd$fcLane, sep="_"))
pd$group <- factor(paste(pd$brainRegion, pd$gender, pd$age, sep="_"))

# normalize by total number of barcoded molecules and multiple by 1e4 (see drop-seq paper)
eset = sweep(eset, 2, colSums(eset)/1e4, FUN = "/")

eZeisel <- log2(eset + 1) # log transform normalized UMIs
pdZeisel <- pd

# Compute first 3 PCs (approximate) using irlba pkg
dat <- sweep(eZeisel[, pdZeisel$brainRegion == "cortex"], 1, 
             rowMeans(eZeisel[, pdZeisel$brainRegion == "cortex"]), FUN = "-")
sZeisel.A <- irlba(dat, nv = 3)
dat <- sweep(eZeisel[, pdZeisel$brainRegion == "hippocampus"], 1, 
             rowMeans(eZeisel[, pdZeisel$brainRegion == "hippocampus"]), FUN = "-")
sZeisel.B <- irlba(dat, nv = 3)

print("[Data loaded]: Zeisel et al. (2015)")

