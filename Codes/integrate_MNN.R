source('Codes/Functions.R')
Initialize()
library(scran)
library('SeuratWrappers')
library(batchelor)

## running on a screen test_2
# merged_samples <- readRDS('Results/preproc_rats/merged/merged_rat_samples.rds')

Samples <- readRDS('Results/preproc_rats/merged/rat_samples_list.rds')
Samples_sce <- lapply(Samples, as.SingleCellExperiment)

genes_list <- lapply(Samples, rownames)
features <- Reduce(intersect,genes_list)
Samples_sce <- lapply(Samples_sce, function(x) x[features,])

# chosen.hvgs <- VariableFeatures(merged_samples)

f.out <- fastMNN(Samples_sce$rat_DA_01, Samples_sce$rat_DA_02,
                 Samples_sce$rat_Lew_01 , Samples_sce$ rat_Lew_02,
                 assay.type='logcounts')



