source('~/HumanLiver/Code/PCanalysisUtils.R')
source('Codes/Functions.R')
Initialize()

mapper = read.delim('~/XSpecies/Data/rat_DA/features.tsv.gz',header=F)

convert_features <- function(sample){
  sample_data <- GetAssayData(sample, assay = 'SCT')
  rownames(sample_data) <- mapper$V2[match(rownames(sample_data), mapper$V1)]
  return(CreateSeuratObject(counts = sample_data, assay = 'SCT'))
}



get_varimax_rotated <- function(gene_exp_matrix, loading_matrix){
  
  ## gene_exp_matrix: gene expression matrix. rows named as genes and columns named as UMIs.
  ##                  attention: Data SHOULD be SCALED.
  ##                  Potentially gained from Seurat GetAssayData() function
  ## loading_matrix: the PCA loading matrix with rows named as gene names and 
  ##                 columns named as PC numbers. Potentially gained from Seurat Loadings() function
  
  ######### Varimax rotation
  initial_data <- t(gene_exp_matrix[rownames(gene_exp_matrix) %in% rownames(loading_matrix),])
  
  ## apply varimax rotation on the loadings
  varimax_res <- varimax(loading_matrix)
  rotatedLoadings <- varimax_res$loadings
  ## calculating the PC scores matrix
  invLoadings     <- t(pracma::pinv(rotatedLoadings))
  scores          <- scale(initial_data) %*% invLoadings ## this second scaling is not necessary
  
  ## compacting the rotated loading and score matrices in a list
  rotated_data <- list(rotLoadings=rotatedLoadings, rotScores = scores)
  return(rotated_data)
}


# rat_LEW_M09_WK_009_3pr_v3 , rat_DA_M09_WK_008_3pr_v3
sample_name <- 'rat_LEW_M09_WK_009_3pr_v3'
data_file <- 'Results/rat_LEW_M09_WK_009_3pr_v3/clusters/seur_clustered_rat_LEW_M09_WK_009_3pr_v3_mito_40_lib_1500.rds'
# 'Results/rat_DA_M09_WK_008_3pr_v3/clusters/seur_clustered_rat_DA_M09_WK_008_3pr_v3_mito_40_lib_1500.rds'
sample_name <- 'merged_sample'
merged_data_file <- 'Results/preproc_rats/merged/merged_rat_samples_newAdded.rds'


merged_samples <- readRDS(merged_data_file)
sample_type_df <- data.frame(sample_type=rep(x = sample_name, times= ncol(seur)))
merged_samples$strain <- 'rat_Lew'
merged_samples$merged_clusters <- as.character(merged_samples$SCT_snn_res.2)
head(sample_type_df)


########

merged_samples <- RunPCA(merged_samples, features=rownames(merged_samples))
### Check the number of PCs to use
standardDev <- Stdev(object = merged_samples, reduction = 'pca')
elbow_plot(standardDev, title = 'merged rat samples')
ElbowPlot(merged_samples, reduction = "pca")

### needed data for the varimax PCA 
loading_matrix = Loadings(merged_samples, 'pca')
gene_exp_matrix = GetAssayData(merged_samples)

### Apply the varimax rotation
rot_data <-  get_varimax_rotated(gene_exp_matrix, loading_matrix)
rotatedLoadings <- rot_data$rotLoadings
scores <- data.frame(rot_data$rotScores)
colnames(scores) = paste0('PC_', 1:ncol(scores))


### which PCs to evaluate in the rotated data
rot_standardDev <- apply(rot_data$rotScores, 2, sd)
elbow_plot(rot_standardDev, title = 'merged rat samples')
ndims = 50
rot_percVar = (rot_standardDev^2 * 100)/sum(rot_standardDev^2)
names(rot_percVar) <- paste0(1:(length(rot_percVar)-1))
perc_variance_threshold = 2
PCs_to_check <- as.numeric(names(rot_percVar[rot_percVar>perc_variance_threshold  ]))
# 1  2  3  4  5  6  7  8  9 11 12 17 22 30


##### Visualizing the PCA plots before the rotation #####
## two PC scatter plots

top_pc = 10

pca_embedding_df <- data.frame(Embeddings(merged_samples,reduction = 'pca'))
pca_embedding_df$clusters = merged_samples$merged_clusters

plot_dir = 'plots/newSamples/'
dir.create(plot_dir)

pdf(paste0(plot_dir, sample_name, '_PCA_plots.pdf'))
for(i in 1:top_pc){
  df = data.frame(PC_1=pca_embedding_df$PC_1,
                  emb_val=pca_embedding_df[,i],
                  cluster=pca_embedding_df$cluster)
  
  p=ggplot(df, aes(x=PC_1, y=emb_val, color=cluster))+geom_point()+
    theme_classic()+ylab(paste0('PC_',i))
  print(p)
  
}
dev.off()



embedd_df_rotated <- data.frame(scores)
pdf(paste0(plot_dir, sample_name, '_VarimaxPCA_plots.pdf'))
for(i in PCs_to_check){ #1:top_pc
  pc_num = i
  rot_df <- data.frame(PC_1=embedd_df_rotated$PC_1,
                       emb_val=embedd_df_rotated[,pc_num],
                       cluster=merged_samples$merged_clusters)
  p=ggplot(rot_df, aes(x=PC_1, y=emb_val, color=cluster))+geom_point()+
    theme_classic()+ylab(paste0('PC_',i))
  print(p)
  }
dev.off()                      




#### checking the expression of specific genes in the data 
sample_name_Lew <- 'rat_LEW_M09_WK_009_3pr_v3'
sample_name_DA <- 'rat_DA_M09_WK_008_3pr_v3'
data_file_Lew <- 'Results/rat_LEW_M09_WK_009_3pr_v3/clusters/seur_clustered_rat_LEW_M09_WK_009_3pr_v3_mito_40_lib_1500.rds'
data_file_DA <- 'Results/rat_DA_M09_WK_008_3pr_v3/clusters/seur_clustered_rat_DA_M09_WK_008_3pr_v3_mito_40_lib_1500.rds'

## pre-processed files with updated QC 
data_file_Lew <- 'objects/rat_LEW_M09_WK_009_3pr_v3/2.seur_dimRed_rat_LEW_M09_WK_009_3pr_v3_mito_50_lib_1000.rds'
data_file_DA <-'objects/rat_DA_M09_WK_008_3pr_v3/2.seur_dimRed_rat_DA_M09_WK_008_3pr_v3_mito_50_lib_1000.rds'

sample_name <- sample_name_Lew
sample_name <- sample_name_DA

data_file <- data_file_Lew
if(sample_name==sample_name_DA) data_file = data_file_DA

sample <- readRDS(data_file)
sample_geneConvert <- convert_features(sample)


gene_list <- c('Ptprc', 'Cd68', 'Cd163', 'Marco', 'Lyz2')


pdf('plots/newSamples/checkImmuneMarkerExpressio_updatedQC.pdf')
for( i in 1:length(gene_list)){
  gene_name = gene_list[i]
  gene_expression <- GetAssayData(sample_geneConvert, assay='SCT')[gene_name,]
  
  df_tsne <- data.frame(tSNE_1=Embeddings(sample, 'tsne')[,1], 
                        tSNE_2=Embeddings(sample, 'tsne')[,2], 
                        gene_expression= gene_expression)
  #head(df_tsne)
  p=ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=gene_expression))+geom_point()+theme_classic()+
    scale_color_viridis(direction = -1)+
    labs(title = paste0('Marker: ', gene_name),subtitle = paste0('Sample: ', sample_name))+
    theme(plot.title = element_text(hjust = 0.5, color = "black", face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "blue"))
  print(p)
}

dev.off()






