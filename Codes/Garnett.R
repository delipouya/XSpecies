
library(monocle)
library(DelayedMatrixStats)
library(garnett)

## import data as a Cell Data Set class 

# load in the data
# NOTE: the 'system.file' file name is only necessary to read in
# included package data
#
mat <- Matrix::readMM(system.file("extdata", "exprs_sparse.mtx", package = "garnett"))
fdata <- read.table(system.file("extdata", "fdata.txt", package = "garnett"))
pdata <- read.table(system.file("extdata", "pdata.txt", package = "garnett"),
                    sep="\t")
row.names(mat) <- row.names(fdata)
colnames(mat) <- row.names(pdata)

# create a new CDS object
pd <- new("AnnotatedDataFrame", data = pdata)
fd <- new("AnnotatedDataFrame", data = fdata)
pbmc_cds <- newCellDataSet(as(mat, "dgCMatrix"),
                           phenoData = pd,
                           featureData = fd)

# generate size factors for normalization later
pbmc_cds <- estimateSizeFactors(pbmc_cds)