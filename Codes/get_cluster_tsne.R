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

INPUT_FILE = paste0('clusters_', INPUT_NAME, '_mito_', MIT_CUT_OFF , '_lib_', LIB_SIZE_CUT_OFF,'.RData')
OUTPUT_FILE = paste0('clusters_', INPUT_NAME, '_mito_', MIT_CUT_OFF , '_lib_', LIB_SIZE_CUT_OFF,'.pdf')
load(paste0('Results/', INPUT_NAME, '/clusters/', INPUT_FILE))


pdf(paste0('Results/', INPUT_NAME, '/clusters/', OUTPUT_FILE))
for( i in 1:length(sCVdata_list)){
  resolutions <- colnames(seur@meta.data)[grep('SCT_snn_res.', colnames(seur@meta.data))]
  a_resolution <- resolutions[i]
  print(a_resolution)
  title = paste0('resolution= ', gsub('SCT_snn_res.','' , a_resolution), ' (', INPUT_NAME, ' mit: ',MIT_CUT_OFF ,' lib size: ', LIB_SIZE_CUT_OFF,')')
  
  tsne_df <- data.frame(getEmb(seur, 'tsne'), clusters=as.character(seur[[a_resolution]][,1]))
  umap_df <- data.frame(getEmb(seur, 'umap'), clusters=as.character(seur[[a_resolution]][,1]))
  tsne_plot = ggplot(tsne_df, aes(x=tSNE_1, y=tSNE_2, color=clusters))+geom_point()+theme_bw()+ggtitle(paste0('tsne ', title))
  umap_plot = ggplot(umap_df, aes(x=UMAP_1, y=UMAP_2, color=clusters))+geom_point()+theme_bw()+ggtitle(paste0('umap ', title))
  plot(tsne_plot)
  plot(umap_plot)
}
dev.off()





