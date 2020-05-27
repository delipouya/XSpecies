source('Codes/Functions.R')
library('fastICA')
Initialize()

dir = 'Results/preproc_rats'
seur_genes_df <- read.delim('Data/rat_DA_M_10WK_003/features.tsv.gz', header = F)

merged_samples <- readRDS(paste0(dir,'/merged/merged_rat_samples.rds'))
merged_samples$strain = ifelse(merged_samples$sample_type =='rat_DA_M_10WK_003' | merged_samples$sample_type=='rat_DA_01', 
                               'rat_DA', 'rat_Lew')
merged_samples$merged_clusters <- as.character(Idents(merged_samples))



### split the gene expression matrix based on the cell-type
cell_type_names = names(table(merged_samples$cell_type))
merged_samples_split <- sapply(1:length(cell_type_names), 
                               function(i){merged_samples[,merged_samples$cell_type == cell_type_names[i]]}, simplify = F)
names(merged_samples_split) = cell_type_names



var_genes <- lapply(merged_samples_split, 
                    function(x){
                      x = FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
                      return(VariableFeatures(x))})

#### checking if the variable genes are present in the gene expression matrix
sapply(1:length(merged_samples_split), function(i) sum(rownames(merged_samples_split[[i]]) %in% var_genes[[i]]))

### only including the highly variable genes
merged_samples_split_hvg = sapply(1:length(merged_samples_split), function(i){
  mat = merged_samples_split[[i]]
  mat_hvg = var_genes[[i]]
  return(GetAssayData(mat[mat_hvg,]))
})

names(merged_samples_split_hvg) = cell_type_names
lapply(merged_samples_split_hvg, head)


ICA_results = sapply(1:length(merged_samples_split), function(i){
  mat <- merged_samples_split_hvg[[i]]
  mat_t <- t(mat)
  return(fastICA(mat_t, 50))
}, simplify = F)

names(ICA_results) = cell_type_names




