## look at the inflammatory macrophages (Lyz2) in the DA vs Lew and then,
## look at the non-inflammatory clusters in DA vs Lew

## These samples are normalized and scaled before being merged and not afterwards 

source('Codes/Functions.R')
source('~/HumanLiver/Code/PCanalysisUtils.R')
Initialize()
library(plyr)

## pre-processed files with updated QC 
data_file_Lew <- 'objects/rat_LEW_M09_WK_009_3pr_v3/2.seur_dimRed_rat_LEW_M09_WK_009_3pr_v3_mito_50_lib_1000.rds'
data_file_DA <- 'objects/rat_DA_M09_WK_008_3pr_v3/2.seur_dimRed_rat_DA_M09_WK_008_3pr_v3_mito_50_lib_1000.rds'

sample_Lew <- convert_features(readRDS(data_file_Lew))
sample_DA <- convert_features(readRDS(data_file_DA))

all_merged_samples_seur <- readRDS('Results/preproc_rats/merged/newSamples_merged_samples_seur.rds')
sCVdata_list <- readRDS('~/newSamples_merged_rats_sCVdata_list.rds')

lapply(sCVdata_list, head)
getHead(all_merged_samples_seur)

all_merged_samples_seur$clusters = paste0('cluster_', as.character(sCVdata_list$res.0.1@Clusters))


plot_dir = 'plots/newSamples/newSamplesMerged/'
###### UMAP ###### 

ElbowPlot(all_merged_samples_seur, 'pca')
PC_NUMBER = 23
all_merged_samples_seur <- RunUMAP(all_merged_samples_seur,dims=1:PC_NUMBER, reduction="harmony")
df_umap <- data.frame(UMAP_1=getEmb(all_merged_samples_seur, 'umap')[,1], 
                      UMAP_2=getEmb(all_merged_samples_seur, 'umap')[,2], 
                      library_size= all_merged_samples_seur$nCount_RNA, 
                      n_expressed=all_merged_samples_seur$nFeature_RNA,
                      clusters=all_merged_samples_seur$clusters, 
                      sample_name=all_merged_samples_seur$sample_name)

pdf(paste0(plot_dir, 'umap_plots.pdf'))
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=library_size))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=n_expressed))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=sample_name))+geom_point()+theme_classic()+scale_color_manual(values = colorPalatte)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=clusters))+geom_point()+theme_classic()+scale_color_manual(values = colorPalatte)
dev.off()



###### tSNE ###### 

all_merged_samples_seur <- RunTSNE(all_merged_samples_seur,dims=1:PC_NUMBER, reduction="harmony")
df_tsne <- data.frame(tSNE_1=getEmb(all_merged_samples_seur, 'tsne')[,1], 
                      tSNE_2=getEmb(all_merged_samples_seur, 'tsne')[,2], 
                      library_size= all_merged_samples_seur$nCount_RNA, 
                      n_expressed=all_merged_samples_seur$nFeature_RNA,
                      clusters=all_merged_samples_seur$clusters, 
                      sample_name=all_merged_samples_seur$sample_name)

pdf(paste0(plot_dir, 'tsne_plots.pdf'))
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=library_size))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=n_expressed))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=sample_name))+geom_point()+theme_classic()+scale_color_manual(values = colorPalatte)
#ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2))+geom_point(aes(color=lowResCluster),alpha=0.7,size=2)+theme_classic()+
#  scale_color_manual(values = colorPalatte)
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=clusters))+geom_point()+theme_classic()+scale_color_manual(values = colorPalatte)
dev.off()



Marco_expression <- GetAssayData(all_merged_samples_seur)['Marco',]
Lyz2_expression <- GetAssayData(all_merged_samples_seur)['Lyz2',]

ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=Marco_expression))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=Lyz2_expression))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=Marco_expression))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=Lyz2_expression))+geom_point()+theme_classic()+scale_color_viridis(direction = -1)


##### selecting the non-inf and inf macrophage population and
####. finding the DE genes between the two populations
DE_genes_res_0.1 <- lapply(sCVdata_list$res.0.1@DEvsRest, function(x) x[order(x$FDR, decreasing = F),])
lapply(DE_genes_res_0.1, head, 20)

## Non-inf macrophage: Cluster 4
all_merged_samples_seur_NonInfmac <- all_merged_samples_seur[,sCVdata_list$res.0.1@Clusters == '4']
## Inf macrophage: Cluster 5 
all_merged_samples_seur_Infmac <- all_merged_samples_seur[,sCVdata_list$res.0.1@Clusters == '5']


### seems like DA samples has more(x2) non-Inf macrophage population compared to the Lewis samples 
### considering the experimental setting, we can't draw any conclusion about the number of these cells in 
### the rat livers, do we need to make any corrections based on this?

table(all_merged_samples_seur_NonInfmac$sample_name)
table(all_merged_samples_seur_Infmac$sample_name)

Idents(all_merged_samples_seur_Infmac) <- all_merged_samples_seur_Infmac$sample_name
Idents(all_merged_samples_seur_NonInfmac) <- all_merged_samples_seur_NonInfmac$sample_name

DAvsLEW_Infmac <- FindMarkers(all_merged_samples_seur_Infmac, 
                              ident.1 = "rat_DA_M09_WK_008", ident.2 = "rat_LEW_M09_WK_009")

DAvsLEW_NonInfmac <- FindMarkers(all_merged_samples_seur_NonInfmac, 
                              ident.1 = "rat_DA_M09_WK_008", ident.2 = "rat_LEW_M09_WK_009")



write.csv(DAvsLEW_Infmac, file = 'Results/preproc_rats/merged/newSamples_DAvsLEW_Infmac.csv',quote = F)
write.csv(DAvsLEW_NonInfmac, file = 'Results/preproc_rats/merged/newSamples_DAvsLEW_NonInfmac.csv',quote = F)
















