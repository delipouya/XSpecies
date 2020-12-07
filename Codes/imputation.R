source('Codes/Functions.R')
Initialize()
library(SAVER)

merged_samples <- readRDS('Results/preproc_rats/merged/merged_rat_samples_2.rds')
merged_samples_mat <- GetAssayData(merged_samples)
merged_samples_imputaed <- saver(merged_samples_mat, ncores = detectCores()-2)


merged_samples_imputed <- readRDS('Results/preproc_rats/merged/merged_samples_imputed.rds')
test <- as.matrix(merged_samples_imputed$estimate)
test <- CreateSeuratObject(merged_samples_imputed$estimate)
