source('Codes/Functions.R')
source('Codes/convert_human_to_ortholog_functions.R')
library(cluster, quietly = TRUE)

Initialize()

### importing old analyzed seurat object > identity of the clusters are known
INPUT_NAME = 'rat_DA'
seur <- readRDS('Results/preproc_rats/rat_DA_mito_30_lib_1500.rds')
cluster_names = paste0('cluster_', as.character(seur$SCT_snn_res.0.6))
seur_data <- GetAssayData(seur[['SCT']])
cluster_names_types = unique(cluster_names)
mapper = read.delim('~/XSpecies/Data/rat_DA/features.tsv.gz',header=F)
reduction <- "pca"
dims <- 1:30
dataset = seur

# silhouette metric
dist.matrix <- dist(x = Embeddings(seur, 'pca')[1:PC_NUMBER])
clusters <- cluster_names
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
dataset$sil <- sil[, 3]

dist.matrix <- dist(x = Embeddings(object = dataset[[reduction]])[, dims])
clusters <- dataset$
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
dataset$sil <- sil[, 3]
length(dataset$sil)
dim(dataset)
pdf('silhouette_plot.pdf', height = 16, width=10)
plot(sil, main='DA-01 clusters', col=1:length(cluster_names_types))
dev.off()



library(RColorBrewer)
n <- 15
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

