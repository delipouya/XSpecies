source('~/HumanLiver/Code/PCanalysisUtils.R')
source('Codes/Functions.R')
Initialize()


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


## import data 
mapper = read.delim('~/XSpecies/Data/rat_DA/features.tsv.gz',header=F)
#merged_samples <- readRDS('Results/preproc_rats/merged/merged_rat_samples_all_features.rds')
merged_samples <- readRDS('Results/preproc_rats/merged/merged_rat_samples_2.rds')

########
merged_samples <- your_scRNAseq_data_object
sample_type_df <- data.frame(str_split_fixed(colnames(merged_samples), pattern = '_', n = 6)  )
sample_type_df$samples_type <- paste0(sample_type_df$X1,'_', sample_type_df$X2, '_', sample_type_df$X3)

sample_type_df$samples_type_2 <- ifelse(sample_type_df$samples_type == 'rat_DA_M', 'rat_DA_M_10WK_003', sample_type_df$samples_type)
table(sample_type_df$samples_type_2)
merged_samples$sample_type <- sample_type_df$samples_type_2

merged_samples$RNA_snn_res.0.05 <- sCVdata_list$RNA_snn_res.0.05@Clusters
########
merged_samples$strain = ifelse(merged_samples$sample_type =='rat_DA_M_10WK_003' | merged_samples$sample_type=='rat_DA_01', 
                               'rat_DA', 'rat_Lew')
merged_samples$merged_clusters <- as.character(Idents(merged_samples))
merged_samples$merged_clusters <- as.character(merged_samples$RNA_snn_res.0.05)

saveRDS(merged_samples, 'Results/preproc_rats/merged/merged_rat_samples_all_features.rds')

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
rot_data <- readRDS('Results/preproc_rats/merged/rotated_Rat_PCs.rds')
rotatedLoadings <- rot_data$rotLoadings
scores <- data.frame(rot_data$rotScores)
colnames(scores) = paste0('PC_', 1:ncol(scores))
#saveRDS(rot_data, 'Results/preproc_rats/merged/rot_data.rds')


### which PCs to evaluate in the rotated data
rot_standardDev <- apply(rot_data$rotScores, 2, sd)
elbow_plot(rot_standardDev, title = 'merged rat samples')
ndims = 50
rot_percVar = (rot_standardDev^2 * 100)/sum(rot_standardDev^2)
names(rot_percVar) <- paste0(1:(length(rot_percVar)-1))
PCs_to_check <- as.numeric(names(rot_percVar[rot_percVar>1.5]))
# 1  2  3  4  5  6  7  8  9 10 13 14 15 16 17 24 33


###### Visualizations ###### 

##### Visualizing the PCA plots before the rotation #####
## two PC scatter plots
top_pc = 10
pca_embedding_df <- getDimRed_df(merged_samples,
                                 top_pc = 1:top_pc,
                                 reduction = 'pca')

### visualizing labeled based on cell-type annotation
pdf('plots/general_cell_type_annot_final.pdf')
for(i in 1:top_pc){
  df = data.frame(PC_1=pca_embedding_df$PC_1,
                  emb_val=pca_embedding_df[,i],
                  cell_type=pca_embedding_df$cell_type,
                  sample=pca_embedding_df$sample)
  
  p=ggplot(df, aes(x=PC_1, y=emb_val, color=cell_type))+geom_point()+
    theme_classic()+scale_color_manual(values = colorPalatte)+ylab(paste0('PC_',i))
  print(p)
  ggplot(df, aes(x=PC_1, y=emb_val, color=sample))+geom_point()+
    theme_classic()+scale_color_manual(values = colorPalatte)+ylab(paste0('PC_',i))
  
}
dev.off()


####### Visualizing the PCA plots after the rotation ###### 
## Some genes expressions to check
Pparg_expression <- GetAssayData(merged_samples)[mapper$V1[mapper$V2=='Pparg'],]
Ppara_expression <- GetAssayData(merged_samples)[mapper$V1[mapper$V2=='Ppara'],]
Rxra_expression <- GetAssayData(merged_samples)[mapper$V1[mapper$V2=='Rxra'],]
Hnf1b_expression <- GetAssayData(merged_samples)[mapper$V1[mapper$V2=='Hnf1b'],]
Cyp2c12_expression <- GetAssayData(merged_samples)[mapper$V1[mapper$V2=='Cyp2c12'],]
Marco_expression <- GetAssayData(merged_samples)[mapper$V1[mapper$V2=='Marco'],]

gene_name = 'AC134224.3'
gene_expression <- GetAssayData(merged_samples)[mapper$V1[mapper$V2==gene_name],]



sCVdata_list <- readRDS('~/XSpecies/objects/merged_subclusters/cluster_4_subclust.rds')
table(sCVdata_list$res.0.02@Clusters)
lapply(sCVdata_list$res.0.02@DEvsRest, head)
subclust <- data.frame(sCVdata_list$res.0.02@Clusters)
colnames(subclust) = 'cluster'
subclust$UMI <- rownames(subclust)
all_UMIs <- data.frame(all_umi=colnames(merged_samples))

all_UMIs_2 <- merge(all_UMIs, subclust, by.x='all_umi' ,by.y= 'UMI',all.x=T,sort=F)


all_UMIs_2$cluster <- ifelse(is.na(all_UMIs_2$cluster),'other',all_UMIs_2$cluster )
all_UMIs_2 <- all_UMIs_2[match(colnames(merged_samples), all_UMIs_2$all_umi),]
sum(all_UMIs_2$all_umi != colnames(merged_samples))

embedd_df_rotated <- data.frame(scores)
pdf('plots/verimax_rotated_PCs_detailedAnnot_2.pdf')
for(i in PCs_to_check){ #1:top_pc
  pc_num = i
  rot_df <- data.frame(PC_1=embedd_df_rotated$PC_1,
                       emb_val=embedd_df_rotated[,pc_num],
                       cell_type=merged_samples$cell_type,
                       cluster=merged_samples$merged_clusters,
                       sample_type=merged_samples$sample_type,
                       strain=merged_samples$strain
                       #mito_perc=merged_samples$mito_perc,
                       #libSize=merged_samples$nCount_RNA,
                       #numDetectedGenes= merged_samples$nFeature_RNA,
                       #Marco = Marco_expression,
                       #gene_exp = gene_expression,
                       #subclusters = all_UMIs_2$cluster
                       #subclusters = ifelse(merged_samples$merged_clusters %in% c('4', '8'), 
                      #                      merged_samples$merged_clusters, 'other')
                      )
  
  
  ggplot(rot_df, aes(x=PC_1, y=emb_val, color=gene_expression))+geom_point(alpha=0.5)+
    theme_classic()+ylab(paste0('PC_',i))+scale_color_viridis(direction = -1)+ggtitle(gene_name)
  
  p1=ggplot(rot_df, aes(x=PC_1, y=emb_val, color=cell_type))+geom_point()+
    theme_classic()+scale_color_manual(values = colorPalatte)+ylab(paste0('PC_',i))
  p2=ggplot(rot_df, aes(x=PC_1, y=emb_val, color=sample_type))+geom_point()+
    theme_classic()+scale_color_manual(values = colorPalatte)+ylab(paste0('PC_',i))
  p3=ggplot(rot_df, aes(x=PC_1, y=emb_val, color=strain))+geom_point()+
    theme_classic()+scale_color_manual(values = colorPalatte)+ylab(paste0('PC_',i))
  p4=ggplot(rot_df, aes(x=PC_1, y=emb_val, color=cluster))+geom_point()+
    theme_classic()+scale_color_manual(values = colorPalatte)+ylab(paste0('PC_',i))
  
  ggplot(rot_df, aes(x=PC_1, y=emb_val, color=subclusters))+geom_point(alpha=0.8)+
    theme_classic()+scale_color_manual(values = colorPalatte)+ylab(paste0('PC_',i))
  
  ggplot(rot_df, aes(x=subclusters, y=emb_val, fill=subclusters))+
    geom_violin()+theme_classic()+scale_fill_manual(values = colorPalatte)+
    coord_flip()+ggtitle(paste0('PC_', pc_num))
  
  
  
  p5=ggplot(rot_df, aes(x=cell_type, y=emb_val, fill=cell_type))+
    geom_violin()+theme_classic()+scale_fill_manual(values = colorPalatte)+
    coord_flip()+ggtitle(paste0('PC_', pc_num))
  p6=ggplot(rot_df, aes(x=cluster, y=emb_val, fill=cluster))+
    geom_violin()+theme_classic()+scale_fill_manual(values = colorPalatte)+
    coord_flip()+ggtitle(paste0('PC_', pc_num))
  
  #ggplot(rot_df, aes(x=PC_1, y=emb_val, color=mito_perc))+geom_point(alpha=0.5)+
  #  theme_classic()+ylab(paste0('PC_',i))+scale_color_viridis(direction = -1)
  
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
}
dev.off()



pdf('plots/find_PC_for_cluster_4.pdf')
for(i in PCs_to_check){ #1:top_pc
  pc_num = i
  rot_df <- data.frame(PC_1=embedd_df_rotated$PC_1,
                       emb_val=embedd_df_rotated[,pc_num],
                       cell_type=merged_samples$cell_type,
                       cluster=merged_samples$merged_clusters,
                       sample_type=merged_samples$sample_type,
                       strain=merged_samples$strain,
                       subclusters = ifelse(merged_samples$merged_clusters %in% c('4'), 
                                             merged_samples$merged_clusters, 'other')
  )
  
  p1=ggplot(rot_df, aes(x=PC_1, y=emb_val, color=subclusters))+geom_point(alpha=0.8)+
    theme_classic()+scale_color_manual(values = colorPalatte)+ylab(paste0('PC_',i))
  
  p2=ggplot(rot_df, aes(x=subclusters, y=emb_val, fill=subclusters))+
    geom_violin()+theme_classic()+scale_fill_manual(values = colorPalatte)+
    coord_flip()+ggtitle(paste0('PC_', pc_num))
  
  print(p1)
  print(p2)
}
dev.off()





## Evaluating the 5th PC in detail
top_pc = 10
pdf('plots/verimax_rotated_PC5.pdf')
for(i in 1:top_pc){
  pc_num = i
  rot_df <- data.frame(PC_5=embedd_df_rotated$PC_5,
                       emb_val=embedd_df_rotated[,pc_num],
                       cell_type=merged_samples$cell_type,
                       cluster=merged_samples$merged_clusters,
                       sample_type=merged_samples$sample_type,
                       strain=merged_samples$strain)
  
  p1=ggplot(rot_df, aes(x=PC_5, y=emb_val, color=cell_type))+geom_point()+
    theme_classic()+scale_color_manual(values = colorPalatte)+ylab(paste0('PC_',i))
  p2=ggplot(rot_df, aes(x=PC_5, y=emb_val, color=sample_type))+geom_point()+
    theme_classic()+scale_color_manual(values = colorPalatte)+ylab(paste0('PC_',i))
  p3=ggplot(rot_df, aes(x=PC_5, y=emb_val, color=strain))+geom_point()+
    theme_classic()+scale_color_manual(values = colorPalatte)+ylab(paste0('PC_',i))
  p4=ggplot(rot_df, aes(x=PC_5, y=emb_val, color=cluster))+geom_point()+
    theme_classic()+scale_color_manual(values = colorPalatte)+ylab(paste0('PC_',i))
  
  print(p1)
  print(p2)
  print(p3)
  print(p4)
}
dev.off()


Lyz_expression <- GetAssayData(merged_samples)[mapper$V1[mapper$V2=='Lyz'],]
gene_name = 'S100a6' #'S100a11'
gene_expression <- GetAssayData(merged_samples)[mapper$V1[mapper$V2==gene_name],]


#### Visualizing the UMAP plot ####
umap_df <- data.frame(Embeddings(merged_samples, reduction = 'umap_h'),
                      cell_type=merged_samples$cell_type,
                      cluster=merged_samples$merged_clusters,
                      sample_type=merged_samples$sample_type,
                      strain=merged_samples$strain,
                      mito_perc=merged_samples$mito_perc,
                      libSize=merged_samples$nCount_RNA,
                      numDetectedGenes= merged_samples$nFeature_RNA)
                      #Marco_EXP = Marco_expression)
ggplot(umap_df, aes(x=umap_h_1, y=umap_h_2, color=gene_expression))+geom_point(alpha=0.5)+
  theme_classic()+ylab(paste0('PC_',i))+scale_color_viridis(direction = -1)+ggtitle(gene_name)

ggplot(umap_df, aes(x=umap_h_1, y=umap_h_2, color=strain))+geom_point()+
  theme_classic()+scale_color_manual(values = colorPalatte)
ggplot(umap_df, aes(x=umap_h_1, y=umap_h_2, color=mito_perc))+geom_point(alpha=0.5)+
  theme_classic()+ylab(paste0('PC_',i))+scale_color_viridis(direction = -1)

cluster_to_check <- 'cluster_7'
umap_df$a_cluster <- umap_df$cluster == cluster_to_check
ggplot(umap_df, aes(x=umap_h_1, y=umap_h_2, color=a_cluster))+geom_point()+
  theme_classic()+ggtitle(cluster_to_check)+xlab('UMAP_1')+ylab('UMAP_2')

umap_df$subclusters = all_UMIs_2$cluster
umap_df$subclusters = ifelse(merged_samples$merged_clusters %in% c('4', '8'), 
                             merged_samples$merged_clusters, 'other')
ggplot(umap_df, aes(x=umap_h_1, y=umap_h_2, color=subclusters))+geom_point()+
  theme_classic()+xlab('UMAP_1')+ylab('UMAP_2')


###### Saving the results ###### 

## save the rotated loadings for pathway analysis 
loadings_dir <- 'Results/ranked_files/varimax_rotated_loadings_2/rot_loadings/'
dir.create(loadings_dir)
rotatedLoadings <- varimax_res$loadings

for(i in PCs_to_check){
  rot_loadings_df <- data.frame(genes=rownames(rotatedLoadings), ranks = rotatedLoadings[,i])
  rot_loadings_df= merge(rot_loadings_df, mapper, by.x='genes',by.y='V1')
  rot_loadings_df <- rot_loadings_df[order(rot_loadings_df$ranks, decreasing = T),]
  rot_loadings_df <- data.frame(genes=rot_loadings_df$V2, ranks=rot_loadings_df$ranks)
  file_name = colnames(loadings_total)[i]
  write.table(rot_loadings_df,
              paste0(loadings_dir, file_name, '.rnk'),
              sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
  
}

#### saving results
embedd_df_rotated <- data.frame(scores, cell_type= merged_samples$cell_type)
colnames(embedd_df_rotated)[1:ncol(scores)] <- paste0('PC_', 1:ncol(scores))
#saveRDS(embedd_df_rotated, 'objects/PCAembedd_df_rotated.rds')



########## Wilcoxon test
embedd_df_rotated <- readRDS('objects/PCAembedd_df_rotated.rds')
embedd_df_rotated$clusters = merged_samples$merged_clusters
embedd_df_rotated$subcluster = all_UMIs_2$cluster 

assined_cell_types = list()
assined_clusters = list()
assined_pc = list()
p.val.list = list()


cell_types <- names(table(embedd_df_rotated$cell_type))
clusters = names(table(embedd_df_rotated$clusters))

p.val.thr = 0.01/(length(clusters)*length(PCs_to_check)) # Bonferroni-corrected P-value threshold
neg_pc <- c(8, 4, 1, 13, 6, 17)
pos_pc <- c(9, 3, 2 , 16, 10, 14, 7)


a_cell_type = 'DC_cell'
a_pc_num = 16
a_cluster = '0'

for(a_cluster in clusters ){ # a_cell_type in cell_types
  print(a_cluster)
  
  for( a_pc_num in PCs_to_check){
    print(a_pc_num)
    
    group_1 <- embedd_df_rotated[embedd_df_rotated$clusters==a_cluster, a_pc_num]
    group_2 <- embedd_df_rotated[embedd_df_rotated$clusters!=a_cluster, a_pc_num]
    
    data = rbind(data.frame(value=group_1, group= 'group_1'), 
                 data.frame(value=group_2, group= 'group_2'))
    ggplot(data, aes(value, color=group))+geom_density()+
      theme_classic()+ggtitle(paste0('cluster: ', a_cluster, ' PC:', a_pc_num))
    
    # greater means x(group 1) is shifted to the right of y(group 2)
    # less means x(group 1) is shifted to the left of y(group 2)
    
    alternative_hyp = ifelse(a_pc_num %in% neg_pc, 'less', 
           ifelse(a_pc_num %in% pos_pc, 'greater', 'none'))
    
    alternative_hyp
    
    if(alternative_hyp == 'none') next
    wilcox_res <- wilcox.test(x=group_1, 
                              y=group_2, 
                              alternative = alternative_hyp)
    wilcox_res$p.value
    if(wilcox_res$p.value < p.val.thr){
      #assined_cell_types[[length(assined_cell_types)+1]] = a_cell_type
      assined_clusters[[length(assined_clusters)+1]] = a_cluster
      assined_pc[[length(assined_pc)+1]] = a_pc_num
      p.val.list[[length(p.val.list)+1]] = wilcox_res$p.value
    }
    
  }
  print('------------------')
}

res_df <- data.frame(cluster=unlist(assined_clusters), 
           PC=unlist(assined_pc), 
           p.val=unlist(p.val.list))

split.data.frame(res_df, res_df$PC)
split.data.frame(res_df, res_df$cluster)






a_pc_num = 13
a_cluster = '2'
group_1 <- embedd_df_rotated[embedd_df_rotated$subcluster==a_cluster, a_pc_num]
group_2 <- embedd_df_rotated[embedd_df_rotated$subcluster!=a_cluster, a_pc_num]

data = rbind(data.frame(value=group_1, group= 'group_1'), 
             data.frame(value=group_2, group= 'group_2'))
ggplot(data, aes(value, color=group))+geom_density()+
  theme_classic()+ggtitle(paste0('cluster: ', a_cluster, ' PC:', a_pc_num))

# greater means x(group 1) is shifted to the right of y(group 2)
# less means x(group 1) is shifted to the left of y(group 2)

alternative_hyp = ifelse(a_pc_num %in% neg_pc, 'less', 
                         ifelse(a_pc_num %in% pos_pc, 'greater', 'none'))

if(alternative_hyp == 'none') next
wilcox_res <- wilcox.test(x=group_1, 
                          y=group_2, 
                          alternative = alternative_hyp)
wilcox_res$p.value
if(wilcox_res$p.value < p.val.thr){
  #assined_cell_types[[length(assined_cell_types)+1]] = a_cell_type
  assined_clusters[[length(assined_clusters)+1]] = a_cluster
  assined_pc[[length(assined_pc)+1]] = a_pc_num
  p.val.list[[length(p.val.list)+1]] = wilcox_res$p.value
}
