source('Codes/Functions.R')
Initialize()

merged_samples <- readRDS('Results/preproc_rats/merged/merged_rat_samples_2.rds')

### visualize the selected population
cluster_name = 'cluster_3'
merged_samples$final_cluster = as.character(Idents(merged_samples))
umap_df <- data.frame(Embeddings(merged_samples, 'umap_h'))
umap_df$clusters = Idents(merged_samples)
colnames(umap_df)[1:2] = c('UMAP_1', 'UMAP_2')
umap_df$clusters <- ifelse(umap_df$clusters %in% cluster_name,umap_df$clusters , 'other')
ggplot(umap_df,aes(x= UMAP_1, y= UMAP_2, color=clusters))+geom_point()+
  theme_classic()+ggtitle(cluster_name)


cluster3_sub_df$umi <- rownames(cluster3_sub_df)
head(cluster3_sub_df)
umap_df$subclusters <- ifelse(rownames(umap_df) %in% cluster3_sub_df$umi ,cluster3_sub_df$subcluster , 'other')
ggplot(umap_df,aes(x= UMAP_1, y= UMAP_2, color=subclusters))+geom_point()+
  theme_classic()+ggtitle(cluster_name)





### coinstructing a subset the input seurat object
UMIs <- colnames(merged_samples[,Idents(merged_samples) %in% cluster_name])
seur = merged_samples[,colnames(merged_samples) %in% UMIs]

seur@assays$SCT <- NULL
seur <- SetAssayData(
  object = seur,
  assay.type = 'RNA',
  new.data = GetAssayData(seur)[,colnames(seur) %in% UMIs],
  slot = 'data'
)

colnames(seur@meta.data)
seur@meta.data <- seur@meta.data[,colnames(seur@meta.data) %in% c("orig.ident","nCount_RNA","nFeature_RNA","mito_perc",
                                                                  "nCount_SCT","nFeature_SCT","cell_type","sample_type",
                                                                  "strain_type")]
seur@reductions$pca <- NULL
seur@reductions$tsne <- NULL
seur@reductions$umap <- NULL
seur@reductions$umap_h <- NULL

seur <- SCTransform(seur,conserve.memory=F,verbose=T,
                        return.only.var.genes=F,variable.features.n = nrow(seur))

seur <- RunPCA(seur,verbose=T, features=rownames(seur))

## PCA

plot(100 * seur@reductions$pca@stdev^2 / seur@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")

PC_NUMBER = 20

### Harmony
seur <- RunHarmony(seur, "sample_type",assay.use="RNA")
seur <- RunUMAP(seur,dims=1:PC_NUMBER, reduction="harmony")



max_seurat_resolution <- 2.5 ## change this to higher values
FDRthresh <- 0.01 # FDR threshold for statistical tests
min_num_DE <- 5
seurat_resolution <- 0 # Starting resolution is this plus the jump value below.
seurat_resolution_jump <- 0.01

DefaultAssay(seur) = 'RNA'
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
    DRforClust="harmony",
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

dir.create('~/XSpecies/objects/merged_subclusters')
saveRDS(sCVdata_list, '~/XSpecies/objects/merged_subclusters/cluster_3_subclust.rds')

sCVdata_list <- readRDS('~/XSpecies/objects/merged_subclusters/cluster_3_subclust.rds')
lapply(sCVdata_list$res.0.02@Clusters, head)

cluster3_sub_df <- data.frame(subcluster=sCVdata_list$res.0.02@Clusters)
head(cluster3_sub_df)



sCVdata_list <- readRDS('~/XSpecies/objects/merged_subclusters/cluster_5_subclust.rds')
table(sCVdata_list$res.0.1@Clusters)
lapply(sCVdata_list$res.0.1@DEvsRest, head)

#### checking the top markers of sub-clusters 
seur_genes_df <- mapper
Cluster_markers <- sCVdata_list$res.0.1@DEvsRest
Cluster_markers_merged <- sapply(1:length(Cluster_markers), 
                                 function(i){
                                   markers_df <- Cluster_markers[[i]]
                                   markers_df$ensemble_ids = rownames(markers_df)
                                   ## merge the ensemble IDs in the dataframe with the HUGO terms 
                                   markers_df_merged <- merge(markers_df, seur_genes_df, 
                                                              by.x='ensemble_ids', 
                                                              by.y='V1', all.x=T, all.y=F,sort=F)
                                   #markers_df_merged2 <-  markers_df_merged[match(markers_df_merged$ensemble_ids, markers_df$ensemble_ids),]
                                   markers_df_merged2 <- markers_df_merged[order(markers_df_merged$logGER*(-log10(markers_df_merged$FDR)),decreasing = T),]
                                   return(markers_df_merged2)
                                 }, simplify = FALSE)

names(Cluster_markers_merged) = names(sCVdata_list$res.0.1@DEvsRest)         

lapply(Cluster_markers_merged, function(x)head(x, 20))               
