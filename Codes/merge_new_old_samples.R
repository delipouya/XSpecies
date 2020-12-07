source('Codes/Functions.R')
source('~/HumanLiver/Code/PCanalysisUtils.R')
Initialize()
library(plyr)

##### Functions #####
clean_name <- function(x){
  x=x[-length(x)] 
  if(x[1]=='merged') x=x[-1]
  return(paste(x, collapse = '_'))
}

mapper = read.delim('~/XSpecies/Data/rat_DA/features.tsv.gz',header=F)

convert_features <- function(sample){
  sample_data <- GetAssayData(sample, assay = 'RNA')
  rownames(sample_data) <- mapper$V2[match(rownames(sample_data), mapper$V1)]
  return(CreateSeuratObject(sample_data))
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

# ToDo:
#  remove the corrupted file: 'Results/preproc_rats/merged/merged_rat_samples_newAdded.rds'


##### Load the data #####
sample_name_Lew <- 'rat_LEW_M09_WK_009_3pr_v3'
sample_name_DA <- 'rat_DA_M09_WK_008_3pr_v3'
data_file_Lew <- 'Results/rat_LEW_M09_WK_009_3pr_v3/clusters/seur_clustered_rat_LEW_M09_WK_009_3pr_v3_mito_40_lib_1500.rds'
data_file_DA <- 'Results/rat_DA_M09_WK_008_3pr_v3/clusters/seur_clustered_rat_DA_M09_WK_008_3pr_v3_mito_40_lib_1500.rds'

sample_Lew <- convert_features(readRDS(data_file_Lew))
sample_DA <- convert_features(readRDS(data_file_DA))
merged_samples <- readRDS('Results/preproc_rats/merged/merged_rat_samples_2.rds')


###### checking the number of expressed gene in each sample after QC ####
a_sample <- readRDS(data_file_DA)
hist(a_sample$nFeature_RNA, xlab = 'number of detected genes', main='DA sample (after QC)')
summary(a_sample$nFeature_RNA)
dim(a_sample)


# Lewis sample number of detected genes(after QC) stats:
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  89     936    1983    2287    3797    7415 

# DA sample number of detected genes(after QC) stats:
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  73     887    1357    1609    2162    5273 


# DA sample number of detected genes(after QC, lib size = 1000) stats:
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#. 70     634    1195    1474    1993    5273 

# Lew sample number of detected genes(after QC, lib size = 1000) stats:
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 89     618    1471    1964    3291    7415 


### merge the count files
### normalize and scale the count files
### run PCA on the merged data set

### merging the samples
all_merged_samples <- merge(merged_samples, y = c(sample_DA, sample_Lew), 
                        add.cell.ids = c('merged', 'rat_DA_M09_WK_008', 'rat_LEW_M09_WK_009'), 
                        project = "rat_data", 
                        merge.data = TRUE)


rm(merged_samples);rm(sample_DA); rm(sample_Lew)
gc()
all_merged_samples_seur <- CreateSeuratObject(GetAssayData(all_merged_samples))

#### add the sample names ##### 
all_sample_names <- str_split(string = colnames(all_merged_samples_seur), pattern = '_')
all_sample_names <- unlist(lapply(all_sample_names, clean_name))
all_merged_samples_seur$sample_name <- all_sample_names
table(all_merged_samples_seur@meta.data$sample_name)

##### normalize and scale the data  #####
all_merged_samples_seur <- NormalizeData(all_merged_samples_seur)
all_merged_samples_seur <- ScaleData(all_merged_samples_seur)
all_merged_samples_seur <- FindVariableFeatures(all_merged_samples_seur)

##### Run PCA and Harmony on the data  #####
all_merged_samples_seur <- RunPCA(all_merged_samples_seur) # features=rownames(all_merged_samples_seur) >> ram issues
#all_merged_samples_seur <- readRDS('Results/preproc_rats/merged/all_merged_samples_seur.rds')
all_merged_samples_seur <- RunHarmony(all_merged_samples_seur, "sample_name",assay.use="RNA")
## check the number of PCs to use
ElbowPlot(all_merged_samples_seur)
PC_NUMBER = 17


#################   Cluster the merged samples   #####################
dir = 'Results/preproc_rats'
seur <- all_merged_samples_seur
rm(all_merged_samples_seur)
gc()

max_seurat_resolution <- 2.4 ## change this to higher values
FDRthresh <- 0.05 # FDR threshold for statistical tests
min_num_DE <- 5
seurat_resolution <- 0 # Starting resolution is this plus the jump value below.
seurat_resolution_jump <- 0.1

seur <- FindNeighbors(seur,reduction="harmony",dims=1:PC_NUMBER,verbose=F)

sCVdata_list <- list()
DE_bw_clust <- TRUE

while(DE_bw_clust) {
  if (seurat_resolution >= max_seurat_resolution) { break }
  seurat_resolution <- seurat_resolution + seurat_resolution_jump 
  # ^ iteratively incrementing resolution parameter 
  
  seur <- FindClusters(seur,resolution=seurat_resolution,verbose=F)
  
  message(" ")
  message("------------------------------------------------------")
  message(paste0("--------  res.",seurat_resolution," with ",
                 length(levels(Idents(seur)))," clusters --------"))
  message("------------------------------------------------------")
  
  if (length(levels(Idents(seur))) <= 1) { 
    message("Only one cluster found, skipping analysis.")
    next 
  } 
  # ^ Only one cluster was found, need to bump up the resolution!
  
  if (length(sCVdata_list) >= 1) {
    temp_cl <- length(levels(Clusters(sCVdata_list[[length(sCVdata_list)]])))
    if (temp_cl == length(levels(Idents(seur)))) { 
      temp_cli <- length(levels(interaction(
        Clusters(sCVdata_list[[length(sCVdata_list)]]),
        Idents(seur),
        drop=T
      )))
      if (temp_cli == length(levels(Idents(seur)))) { 
        message("Clusters unchanged from previous, skipping analysis.")
        next 
      }
    }
  }
  
  curr_sCVdata <- CalcSCV(
    inD=seur,
    assayType="RNA",
    assaySlot="counts",
    cl=Idents(seur), 
    # ^ your most recent clustering results get stored in the Seurat "ident" slot
    exponent=NA, 
    # ^ going to use the corrected counts from SCTransform
    pseudocount=NA,
    DRthresh=0.1,
    DRforClust="harmony", #pca
    calcSil=T,
    calcDEvsRest=T,
    calcDEcombn=T
  )
  
  DE_bw_NN <- sapply(DEneighb(curr_sCVdata,FDRthresh),nrow)
  # ^ counts # of DE genes between neighbouring clusters at your selected FDR threshold
  message(paste("Number of DE genes between nearest neighbours:",min(DE_bw_NN)))
  
  if (min(DE_bw_NN) < min_num_DE) { DE_bw_clust <- FALSE }
  # ^ If no DE genes between nearest neighbours, don't loop again.
  
  sCVdata_list[[paste0("res.",seurat_resolution)]] <- curr_sCVdata
}



saveRDS(sCVdata_list, '~/newSamples_merged_rats_sCVdata_list.rds') 
sCVdata_list <- readRDS('~/newSamples_merged_rats_sCVdata_list.rds')
rm(seur); gc()

######################################


######### Run Varimax on the merged samples #########

#all_merged_samples_seur <- readRDS('Results/preproc_rats/merged/all_merged_samples_seur.rds')
saveRDS(all_merged_samples_seur, 'Results/preproc_rats/merged/newSamples_merged_samples_seur.rds')
all_merged_samples_seur <- RunPCA(all_merged_samples_seur, features=rownames(all_merged_samples_seur))

### needed data for the varimax PCA 
loading_matrix = Loadings(all_merged_samples_seur, 'pca')
gene_exp_matrix = GetAssayData(all_merged_samples_seur, assay = 'RNA')
dim(gene_exp_matrix)
dim(loading_matrix)

### Apply the varimax rotation
rot_data <-  get_varimax_rotated(gene_exp_matrix, loading_matrix)
#saveRDS(rot_data, 'Results/preproc_rats/merged/varimax_results_newSamples_merged_rats.rds')
rot_data <- readRDS('Results/preproc_rats/merged/varimax_results_All_merged_rats.rds')

rotatedLoadings <- rot_data$rotLoadings
scores <- data.frame(rot_data$rotScores)
colnames(scores) = paste0('PC_', 1:ncol(scores))
dim(scores)
dim(rotatedLoadings)

################## Visualize the rotated PCs and evaluate mapping with clusters ####################

all_merged_samples_seur$clusters <- paste0('cluster_', as.character(sCVdata_list$res.0.1@Clusters))

### which PCs to evaluate in the rotated data
rot_standardDev <- apply(rot_data$rotScores, 2, sd)
elbow_plot(rot_standardDev, title = 'merged rat samples')
ndims = 50
rot_percVar = (rot_standardDev^2 * 100)/sum(rot_standardDev^2)
names(rot_percVar) <- paste0(1:(length(rot_percVar)-1))
perc_variance_threshold = 1
PCs_to_check <- as.numeric(names(rot_percVar[rot_percVar>perc_variance_threshold  ]))


##### Visualizing the PCA plots before the rotation #####
## two PC scatter plots

top_pc = 10

pca_embedding_df <- data.frame(Embeddings(all_merged_samples_seur,reduction = 'pca'))
pca_embedding_df$clusters = all_merged_samples_seur$clusters
dim(pca_embedding_df)


plot_dir = 'plots/newSamples/newSamplesMerged/'
dir.create(plot_dir)

pdf(paste0(plot_dir, 'PCA_plots.pdf'))
for(i in 1:top_pc){
  df = data.frame(PC_1=pca_embedding_df$PC_1,
                  emb_val=pca_embedding_df[,i],
                  cluster=pca_embedding_df$cluster,
                  Library_size=all_merged_samples_seur$nCount_RNA,
                  num_expressed_genes=all_merged_samples_seur$nFeature_RNA)
  
  p1=ggplot(df, aes(x=PC_1, y=emb_val, color=cluster))+geom_point()+
    theme_classic()+ylab(paste0('PC_',i))+scale_color_manual(values=colorPalatte)
  p2=ggplot(df, aes(x=PC_1, y=emb_val, color=Library_size))+geom_point()+
    theme_classic()+ylab(paste0('PC_',i))
  p3=ggplot(df, aes(x=PC_1, y=emb_val, color=num_expressed_genes))+geom_point()+
    theme_classic()+ylab(paste0('PC_',i))
  print(p1);print(p2);print(p3) 
  
}
dev.off()



embedd_df_rotated <- data.frame(scores)
dim(embedd_df_rotated)
all_merged_samples_seur$strain = ifelse(all_merged_samples_seur$sample_name %in% 
                                          c('rat_DA_01', 'rat_DA_M_10WK_003', 'rat_DA_M09_WK_008'), 'rat_DA', 'rat_LEW')
all_merged_samples_seur$data_label = ifelse(all_merged_samples_seur$sample_name %in% 
                                              c('rat_LEW_M09_WK_009', 'rat_DA_M09_WK_008'), 'new', 'old')

dim(all_merged_samples_seur)
pdf(paste0(plot_dir, 'VarimaxPCA_plots.pdf'))
for(i in PCs_to_check){ #1:top_pc PCs_to_check
  pc_num = i
  rot_df <- data.frame(PC_1=embedd_df_rotated$PC_1,
                       emb_val=embedd_df_rotated[,pc_num],
                       cluster=all_merged_samples_seur$clusters,
                       Library_size=all_merged_samples_seur$nCount_RNA,
                       num_expressed_genes=all_merged_samples_seur$nFeature_RNA,
                       strain=all_merged_samples_seur$strain,
                       label=all_merged_samples_seur$data_label,
                       sample_name = all_merged_samples_seur$sample_name)
  
  p1=ggplot(rot_df, aes(x=PC_1, y=emb_val, color=cluster))+geom_point()+
    theme_classic()+ylab(paste0('PC_',i))+scale_color_manual(values=colorPalatte)
  p2=ggplot(rot_df, aes(x=PC_1, y=emb_val, color=strain))+geom_point()+
    theme_classic()+ylab(paste0('PC_',i))
  p3=ggplot(rot_df, aes(x=PC_1, y=emb_val, color=sample_name))+geom_point()+
    theme_classic()+ylab(paste0('PC_',i))
  p4=ggplot(rot_df, aes(x=PC_1, y=emb_val, color=label))+geom_point()+
    theme_classic()+ylab(paste0('PC_',i))
  print(p1);print(p2);print(p3)#;print(p4)
}
dev.off()                      



##### checking the sample composition of each cluster
all_merged_samples_seur$lowResCluster <- paste0('cluster_', as.character(sCVdata_list$res.0.1@Clusters))
clusters_pred = sCVdata_list$res.0.8@Clusters

df <- data.frame(sample_type = all_merged_samples_seur$sample_name, 
                 cluster = as.character(clusters_pred))
rownames(df) = NULL
head(df)

#### based on sample-type
counts <- ddply(df, .(df$sample_type, df$cluster), nrow)
names(counts) <- c("sample_type", "cluster", "Freq")
counts$cluster= factor(counts$cluster, levels = as.character(0:(length(unique(df$cluster))-1)) ) 
ggplot(data=counts, aes(x=cluster, y=Freq, fill=sample_type)) +
  geom_bar(stat="identity",color='black')+theme_classic()+scale_fill_brewer(palette = "Set3")+
  ylab('Counts')

counts_split <- split( counts , f = counts$cluster )
counts_split_norm <- lapply(counts_split, function(x) {x$Freq=x$Freq/sum(x$Freq);x})
counts_norm <- do.call(rbind,counts_split_norm )
ggplot(data=counts_norm, aes(x=cluster, y=Freq, fill=sample_type)) +
  geom_bar(stat="identity",color='black')+theme_classic()+scale_fill_brewer(palette = "Set3")

## after clustering:
df_umap$cluster = as.character(clusters_pred)
head(df_umap)
#scale_colour_manual(values=colorPalatte)



###### UMAP ###### 

ElbowPlot(all_merged_samples_seur, 'pca')
PC_NUMBER = 23
all_merged_samples_seur <- RunUMAP(all_merged_samples_seur,dims=1:PC_NUMBER, reduction="harmony")
df_umap <- data.frame(UMAP_1=getEmb(all_merged_samples_seur, 'umap')[,1], 
                      UMAP_2=getEmb(all_merged_samples_seur, 'umap')[,2], 
                      library_size= all_merged_samples_seur$nCount_RNA, 
                      n_expressed=all_merged_samples_seur$nFeature_RNA,
                      clusters=all_merged_samples_seur$clusters, 
                      lowResCluster=all_merged_samples_seur$lowResCluster,
                      sample_name=all_merged_samples_seur$sample_name)

pdf(paste0(plot_dir, 'umap_plots.pdf'))
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=library_size))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=n_expressed))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=sample_name))+geom_point()+theme_classic()+scale_color_manual(values = colorPalatte)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2))+geom_point(aes(color=lowResCluster),alpha=0.7,size=2)+theme_classic()+
  scale_color_manual(values = colorPalatte)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=clusters))+geom_point()+theme_classic()+scale_color_manual(values = colorPalatte)
dev.off()


###### tSNE ###### 

all_merged_samples_seur <- RunTSNE(all_merged_samples_seur,dims=1:PC_NUMBER, reduction="harmony")
df_tsne <- data.frame(tSNE_1=getEmb(all_merged_samples_seur, 'tsne')[,1], 
                      tSNE_2=getEmb(all_merged_samples_seur, 'tsne')[,2], 
                      library_size= all_merged_samples_seur$nCount_RNA, 
                      n_expressed=all_merged_samples_seur$nFeature_RNA,
                      clusters=all_merged_samples_seur$clusters, 
                      lowResCluster=all_merged_samples_seur$lowResCluster,
                      sample_name=all_merged_samples_seur$sample_name)

pdf(paste0(plot_dir, 'tsne_plots.pdf'))
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=library_size))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=n_expressed))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=sample_name))+geom_point()+theme_classic()+scale_color_manual(values = colorPalatte)
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2))+geom_point(aes(color=lowResCluster),alpha=0.7,size=2)+theme_classic()+
  scale_color_manual(values = colorPalatte)
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=clusters))+geom_point()+theme_classic()+scale_color_manual(values = colorPalatte)
dev.off()


#saveRDS(all_merged_samples_seur, 'Results/preproc_rats/merged/all_merged_samples_seur.rds')





