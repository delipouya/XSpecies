source('Codes/Functions.R')
Initialize()
library(plyr)


#####
### importing data
dir = 'Results/preproc_rats'
Samples <- lapply(list.files(dir,pattern = '*.rds',full.names = T), readRDS)
names(Samples) <- c('rat_DA_M_10WK_003', 'rat_DA_01', 'rat_Lew_01', 'rat_Lew_02', 'rat_DA_02')
lapply(Samples, function(x) dim(x@assays$RNA))
num_cells <- lapply(Samples, function(x) ncol(x@assays$SCT))
genes <- unlist(lapply(Samples, function(x) rownames(x@assays$SCT)))

cluster_cell_type_df <- get_manual_labels()


### adding the final clusters to metadata
for(i in 1:length(Samples)){
         sample_name = names(Samples)[i]
         sample = Samples[[i]]
         
         if(sample_name=='rat_DA_02') 
           Idents(Samples[[i]]) <- paste0('cluster_',as.character(sample$SCT_snn_res.1))
         
         if(sample_name=='rat_Lew_01') 
           Idents(Samples[[i]]) <- paste0('cluster_',as.character(sample$SCT_snn_res.0.8))
         
         if(sample_name %in% c('rat_DA_01', 'rat_Lew_02', 'rat_DA_M_10WK_003')) 
           Idents(Samples[[i]]) <- paste0('cluster_',Idents(sample))
         
         Samples[[i]]$final_cluster <-Idents(Samples[[i]])
}

### adding predicted cell-type to meta-data
for( i in 1:length(Samples)){
  sample_name = names(Samples)[i]
  a_sample = Samples[[sample_name]]
  cluster_list = data.frame(umi=colnames(a_sample), cluster=a_sample$final_cluster, index=1:ncol(a_sample))
  map_df = cluster_cell_type_df[[sample_name]]
  merged = merge(cluster_list, map_df, by.x='cluster', by.y='cluster', all.x=T, sort=F)
  merged <- merged[order(merged$index, decreasing = F),]
  Samples[[sample_name]]$cell_type <- merged$cell_type
}


lapply(Samples, function(x) head(x$cell_type))
#saveRDS(Samples, paste0(dir,'/merged/rat_samples_list.rds'))
  
### merging the samples
merged_samples <- merge(Samples[[1]], y = c(Samples[[2]], Samples[[3]], Samples[[4]], Samples[[5]]), 
                        add.cell.ids = names(Samples), project = "rat_data", 
                        merge.data = TRUE)

### check-ups on merge
sum(data.frame(num_cells)) == ncol(merged_samples@assays$SCT)
length(genes[!duplicated(genes)]) == nrow(merged_samples@assays$SCT)

### finding variables genes and scaling the data 
merged_samples <- FindVariableFeatures(merged_samples)
merged_samples <- ScaleData(merged_samples)
sample_type <- unlist(lapply(str_split(colnames(merged_samples),'_'), 
                             function(x) paste(x[-length(x)],collapse = '_')))
merged_samples@meta.data$sample_type = sample_type
#saveRDS(merged_samples, paste0(dir,'/merged/merged_rat_samples.rds'))






#######################################
#### load merged rat samples:

merged_samples <- readRDS(paste0(dir,'/merged/merged_rat_samples.rds'))

## PCA
merged_samples <- RunPCA(merged_samples,verbose=T)
plot(100 * merged_samples@reductions$pca@stdev^2 / merged_samples@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")

PC_NUMBER = 18



### Harmony
merged_samples <- RunHarmony(merged_samples, "sample_type",assay.use="SCT")
head(Embeddings(merged_samples, reduction = "harmony"))

### approach 2
harmonized_pcs <- HarmonyMatrix(
  data_mat  = Embeddings(merged_samples, reduction = "pca"),
  meta_data = merged_samples@meta.data,
  vars_use  = "sample_type",
  do_pca    = FALSE
)


## UMAP
# running UMAP seperately: uwot, umapr, M3C
merged_samples <- RunUMAP(merged_samples,dims=1:PC_NUMBER, reduction="pca",perplexity=30)
merged_samples <- RunUMAP(merged_samples,dims=1:PC_NUMBER, reduction = "harmony",perplexity=30)


umap_emb <- data.frame(Embeddings(merged_samples, 'umap'))
umap_emb$sample_type = merged_samples$sample_type
umap_emb$cell_type = merged_samples$cell_type

ggplot(umap_emb, aes(x=UMAP_1, y=UMAP_2))+
  geom_point(aes(color=sample_type),alpha=0.7,size=2)+theme_classic()

ggplot(umap_emb, aes(x=UMAP_1, y=UMAP_2))+
  geom_point(aes(shape=sample_type, color=cell_type),alpha=0.7,size=2)+
  theme_classic()+scale_colour_brewer(palette = "Paired")




## tSNE
merged_samples <- RunTSNE(merged_samples,dims=1:PC_NUMBER,reduction="pca",perplexity=30)
merged_samples <- RunTSNE(merged_samples,dims=1:PC_NUMBER,reduction="harmony",perplexity=30)

tsne_emb <- data.frame(Embeddings(merged_samples, 'tsne'))
tsne_emb$sample_type = merged_samples$sample_type
tsne_emb$cell_type = merged_samples$cell_type

ggplot(tsne_emb, aes(x=tSNE_1, y=tSNE_2))+
  geom_point(aes(color=sample_type),alpha=0.7,size=2)+theme_classic()

ggplot(tsne_emb, aes(x=tSNE_1, y=tSNE_2))+
  geom_point(aes(shape=sample_type, color=cell_type),alpha=0.7,size=2)+
  theme_classic()+scale_colour_brewer(palette = "Paired")



########### clustering the merged dataset
res = 0.4
merged_samples <- FindNeighbors(merged_samples,reduction="harmony",dims=1:PC_NUMBER,verbose=T)
merged_samples <- FindClusters(merged_samples,resolution=res,verbose=T)

## this clustering has been done on the integrated data using the harmony's method 
sCVdata_list <- readRDS('~/sCVdata_list_merged_rat.rds') 
clusters_pred = sCVdata_list$res.0.15@Clusters

df <- data.frame(sample_type = merged_samples$sample_type, 
                 cell_type = merged_samples$cell_type,
                 cluster = as.character(clusters_pred))
rownames(df) = NULL

# merged_rat_clusters.rds, sCVdata_list_merged_rat.rds

#### based on sample-type

counts <- ddply(df, .(df$sample_type, df$cluster), nrow)
names(counts) <- c("sample_type", "cluster", "Freq")
counts$cluster= factor(counts$cluster, levels = as.character(0:(length(unique(df$cluster))-1)) ) 
ggplot(data=counts, aes(x=cluster, y=Freq, fill=sample_type)) +
  geom_bar(stat="identity",color='black')+theme_classic()+scale_fill_brewer(palette = "Set3")

counts_split <- split( counts , f = counts$cluster )
counts_split_norm <- lapply(counts_split, function(x) {x$Freq=x$Freq/sum(x$Freq);x})
counts_norm <- do.call(rbind,counts_split_norm )
ggplot(data=counts_norm, aes(x=cluster, y=Freq, fill=sample_type)) +
  geom_bar(stat="identity",color='black')+theme_classic()+scale_fill_brewer(palette = "Set3")

## after clustering:
umap_emb$cluster = as.character(clusters_pred)
ggplot(umap_emb, aes(x=UMAP_1, y=UMAP_2))+
  geom_point(aes(color=cluster),alpha=0.7,size=2)+theme_classic()+
  scale_colour_manual(values=colorPalatte)



###### based on cell-type

counts <- ddply(df, .(df$cell_type, df$cluster), nrow)
names(counts) <- c("cell_type", "cluster", "Freq")
counts$cluster= factor(counts$cluster, levels = as.character(0:(length(unique(df$cluster))-1)) ) 
ggplot(data=counts, aes(x=cluster, y=Freq, fill=cell_type)) +
  geom_bar(stat="identity",color='black')+theme_classic()+scale_fill_brewer(palette = "Set3")

counts_split <- split( counts , f = counts$cluster )
counts_split_norm <- lapply(counts_split, function(x) {x$Freq=x$Freq/sum(x$Freq);x})
counts_norm <- do.call(rbind,counts_split_norm )
ggplot(data=counts_norm, aes(x=cluster, y=Freq, fill=cell_type)) +
  geom_bar(stat="identity",color='black')+theme_classic()+scale_fill_brewer(palette = "Set3")

