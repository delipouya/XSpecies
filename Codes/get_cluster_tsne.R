## Run this script as: 
# Rscript ~/XSpecies/Codes/get_cluster_tsne.R 'mouse' 50 2000

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

source('Codes/Functions.R')
Initialize()

INPUT_NAME = args[1] 
MIT_CUT_OFF = args[2]  
LIB_SIZE_CUT_OFF = args[3] 

INPUT_FILE = paste0('clusters_', INPUT_NAME, '_mito_', MIT_CUT_OFF , '_lib_', LIB_SIZE_CUT_OFF,'_v2.RData')
OUTPUT_FILE = paste0('clusters_', INPUT_NAME, '_mito_', MIT_CUT_OFF , '_lib_', LIB_SIZE_CUT_OFF,'.pdf')
load(paste0('Results/', INPUT_NAME, '/clusters/', INPUT_FILE))

### manual import
MIT_CUT_OFF = 20
LIB_SIZE_CUT_OFF = 2000
INPUT_NAME =  'rat_DA_M09_WK_008_3pr_v3' # 'rat_LEW_M09_WK_009_3pr_v3' , 'rat_DA_01_reseq'
INPUT_FILE = '2.seur_dimRed_rat_DA_M09_WK_008_3pr_v3_mito_40_lib_1500.rds'

OUTPUT_NAME = gsub('.rds','',gsub('2.seur_dimRed_','',INPUT_FILE ))
INPUT_FILE = paste0('Results/',INPUT_NAME,'/clusters/clusters_',OUTPUT_NAME,'_v2.rds')

output_filename <- paste0('Results/',INPUT_NAME,'/clusters/clusters_',OUTPUT_NAME,'.rds')
seur_output_filename_rds <- paste0('Results/',INPUT_NAME,'/clusters/seur_clustered_',OUTPUT_NAME,'.rds')


INPUT_FILE = '2.seur_dimRed_rat_DA_01_reseq_mito_30_lib_1500.rds'
seur = readRDS(seur_output_filename_rds)
sCVdata_list = readRDS(output_filename)
colnames(seur@meta.data)
i = 2


resolutions <- colnames(seur@meta.data)[grep('SCT_snn_res.', colnames(seur@meta.data))]
pdf(paste0('Results/', INPUT_NAME, '/clusters/', OUTPUT_FILE))
i = length(resolutions)
for( i in 1:length(resolutions)){
  a_resolution <- resolutions[i]
  print(a_resolution)
  title = paste0('resolution= ', gsub('SCT_snn_res.','' , a_resolution), ' (', INPUT_NAME, ' mit: ',MIT_CUT_OFF ,' lib size: ', LIB_SIZE_CUT_OFF,')')
  
  tsne_df <- data.frame(getEmb(seur, 'tsne'), clusters=as.character(seur[[a_resolution]][,1]))
  umap_df <- data.frame(getEmb(seur, 'umap'), clusters=as.character(seur[[a_resolution]][,1]))
  tsne_plot = ggplot(tsne_df, aes(x=tSNE_1, y=tSNE_2, color=clusters))+
    geom_point()+theme_bw()+ggtitle(paste0('tsne ', title))
  umap_plot = ggplot(umap_df, aes(x=UMAP_1, y=UMAP_2, color=clusters))+
    geom_point()+theme_bw()+ggtitle(paste0('umap ', title))
  plot(tsne_plot)
  plot(umap_plot)
}
dev.off()

seurat_resolution = 0.2 # 0.4
seur <- FindNeighbors(seur,reduction="pca",dims=1:PC_NUMBER,verbose=F)
seur <- FindClusters(seur,resolution=seurat_resolution,verbose=F)

saveRDS(seur,output_filename_rds)


a_resolution <- resolutions[length(resolutions)]
title = paste0('resolution= ', gsub('SCT_snn_res.','' , a_resolution), ' (', INPUT_NAME, ' mit: ',MIT_CUT_OFF ,' lib size: ', LIB_SIZE_CUT_OFF,')')

tsne_df <- data.frame(getEmb(seur, 'tsne'), clusters=as.character(seur[[a_resolution]][,1]))
umap_df <- data.frame(getEmb(seur, 'umap'), clusters=as.character(seur[[a_resolution]][,1]))
ggplot(tsne_df, aes(x=tSNE_1, y=tSNE_2, color=clusters))+geom_point()+theme_bw()+ggtitle(paste0('tsne ', title))
ggplot(umap_df, aes(x=UMAP_1, y=UMAP_2, color=clusters))+geom_point()+theme_bw()+ggtitle(paste0('umap ', title))





