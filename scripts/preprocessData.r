
#### Load Zheng et al. (2017) UMI and phenotypic data
library(TENxGenomics)
data_path <- "/net/irizarryfs01/srv/export/irizarryfs01/share_root/shicks/TENxGenomics/1M_neurons"
HDF5_20K_file <- "1M_neurons_neuron20k.h5"

# function to calculate cell-level summary statistics
cell.summary <- function(x, is.mito) {
  ucidx <- unique(x$cidx)
  x$cidx <- match(x$cidx, ucidx)
  
  # Taking sum of only mitochondrial genes.
  mito.val <- x$value
  mito.val[!is.mito[x$ridx]] <- 0L
  
  data.frame(
    cidx = ucidx,
    n = tabulate(x$cidx, length(ucidx)),
    sum = vapply(split(x$value, x$cidx), sum, numeric(1), USE.NAMES=FALSE),
    mito = vapply(split(mito.val, x$cidx), sum, numeric(1), USE.NAMES=FALSE),
    CDR = vapply(split(x$value>0, x$cidx), sum, numeric(1), USE.NAMES=FALSE),
    CDRlt1 = vapply(split(x$value>1, x$cidx), sum, numeric(1), USE.NAMES=FALSE)
  )
}

# load HDF5 data into R
tenx.se <- tenxSummarizedExperiment(file.path(data_path,HDF5_20K_file))

# add gene-level annotation 
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
chr.loc <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=rownames(tenx.se), 
                  keytype="GENEID", column="CDSCHROM")
anno <- data.frame(SYMBOL=rowData(tenx.se)$Symbol, CHR=chr.loc)
rowData(tenx.se) <- cbind(rowData(tenx.se), anno)
head(anno)

# calculate cell-level summary statistics
library(BiocParallel)
register(MulticoreParam(progressbar=TRUE, workers=1))
tenx <- assay(tenx.se, withDimnames=FALSE)
result <- tenxiterate(tenx, cell.summary, is.mito=(anno$CHR=="chrM" & !is.na(anno$CHR)))
cols <- do.call("rbind", result)
colnames(cols) <- c("Cell", "Ngenes", "Libsize", "Mito", "CDR", "CDRlt1")
cols$CDR <- cols$CDR / nrow(tenx.se)
cols$CDRlt1 <- cols$CDRlt1 / nrow(tenx.se)

# use scater to create index of discard cells
library(scater)
low.libsize <- isOutlier(cols$Libsize, log=TRUE, nmad=3, batch=tenx.se$Library, type="lower")
low.ngenes <- isOutlier(cols$Ngenes, log=TRUE, nmad=3, batch=tenx.se$Library, type="lower")
discard <- low.libsize | low.ngenes 
data.frame(LowLib=sum(low.libsize), LowGenes=sum(low.ngenes), Lost=sum(discard))

# create DelayedMatrix object and reload 10X data excluding discard cells
tenxmat <- TENxMatrix(file.path(data_path,HDF5_20K_file))[,!discard]

# write count table to HDF5Matrix with column chunk
library(HDF5Array)
options(DelayedArray.block.size=2e8)
mat.out <- writeHDF5Array(tenxmat, file=file.path(data_path, "qc_counts.h5"), name="neurons", chunk_dim=c(1e4, 1))

# create new SummarizedExperiment object excluding discard cells
se.out <- SummarizedExperiment(list(counts=mat.out), rowData=rowData(tenx.se), 
                               colData=cbind(colData(tenx.se)[!discard,], cols[!discard,]))

# pick row/column access
ncells <- ncol(mat.out)
ngenes <- nrow(mat.out)
cache.size <- ncells*sqrt(ngenes)
chunk.nrow <- ceiling(cache.size/ncells)
chunk.ncol <- ceiling(cache.size/ngenes)
c(chunk.nrow, chunk.ncol)

# calculate log2-expression values and save to new HDF5Matrix
se.out$sfs <- colData(se.out)$Libsize/1e6
log.norm <- log2(t(t(mat.out) / t(t(se.out$sfs))) + 1)
expr.mat <- writeHDF5Array(log.norm, file=file.path(data_path, "norm_exprs.h5"), name="neurons", 
                           chunk=c(chunk.nrow, chunk.ncol))

# store normalized expression into SummarizedExperiment
assay(se.out, "exprs", withDimnames=FALSE) <- expr.mat

# save object for serialization
saveRDS(se.out, file=file.path(data_path, "qc_mat.rds"))

### FOR PCA
# row standardize
exprs.mat <- (expr.mat - rowMeans(expr.mat))
dat <- t(realize(exprs.mat))

# Compute first 3 PCs (approximate) using irlba pkg
sTENx <- irlba(as.matrix(t(dat)), nv = 3)
saveRDS(sTENx, file=file.path(data_path, "pca_results.rds"))


### FOR PCA (with colMeans Removed)
# first remove colMeans on log scale and then row standardize
expr.mat.corrected <- t(t(expr.mat) - t(t(colMeans(expr.mat))))

# row standardize
exprs.mat <- (expr.mat.corrected - rowMeans(expr.mat.corrected))
dat <- t(realize(exprs.mat))

# Compute first 3 PCs (approximate) using irlba pkg
sTENx <- irlba(as.matrix(t(dat)), nv = 3)
saveRDS(sTENx, file=file.path(data_path, "pca_results_colMeansRemoved.rds"))



