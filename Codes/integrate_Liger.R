
# https://github.com/MacoskoLab/liger
# tutorial on each function meaning:
# https://macoskolab.github.io/liger/walkthrough_pbmc.html
## seurat wrapper:
# https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/liger.html
# http://htmlpreview.github.io/?https://github.com/MacoskoLab/liger/blob/master/vignettes/Integrating_multi_scRNA_data.html
source('Codes/Functions.R')
Initialize()
library(liger)
library('SeuratWrappers')

merged_samples <- readRDS('~/Desktop/merged_rat_samples.rds')
merged_samples_liger <- seuratToLiger(merged_samples,
                                      combined.seurat = T,
                                      names = "use-projects",
                                      meta.var = 'sample_type',
                                      assays.use = 'SCT',
                                      raw.assay = "RNA",
                                      remove.missing = T,
                                      renormalize = T,
                                      use.seurat.genes = T,
                                      num.hvg.info = NULL,
                                      use.idents = T,
                                      use.tsne = F,
                                      cca.to.H = F)
merged_samples_liger <- scaleNotCenter(merged_samples_liger)
k.suggest <- suggestK(merged_samples_liger, num.cores = 4)
merged_samples_liger <- optimizeALS(merged_samples_liger, k=10, thresh = 5e-5, nrep = 3)
merged_samples_liger <- runTSNE(merged_samples_liger, use.raw = T)
p1 <- plotByDatasetAndCluster(merged_samples_liger, return.plots = T)
merged_samples_liger <- quantileAlignSNF(merged_samples_liger, resolution = 0.4, small.clust.thresh = 20)
merged_samples_liger <- runTSNE(merged_samples_liger)
p_a <- plotByDatasetAndCluster(merged_samples_liger, return.plots = T,clusters='cell_types') 


merged_samples_liger <- runUMAP(merged_samples_liger)
p_u <- plotByDatasetAndCluster(merged_samples_liger, return.plots = T, clusters='cell_types') 








