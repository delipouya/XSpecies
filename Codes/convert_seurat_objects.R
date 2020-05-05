source('Codes/Functions.R')
Initialize()
library(Seurat)
library(reticulate)
use_condaenv("scanpy")

library(loomR)
library("hdf5r")
library(scater)
library(Seurat)
library(cowplot)
library(reticulate)


names(Samples)

Samples <- readRDS('Results/preproc_rats/merged/rat_samples_list.rds')
sample = Samples[[1]]

sample_h5ad <- Convert(from = sample, to = "anndata", filename = "~/rat_DA_01.h5ad")
sample.loom <- as.loom(sample, filename = "~/rat_DA_01.loom", verbose = TRUE)

