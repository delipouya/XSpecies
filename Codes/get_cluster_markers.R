options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

source('Codes/Functions.R')
Initialize()
INPUT_NAME = args[1] 
# INPUT_NAME = 'mouse'  #'rat_Rnor' 
model_animal_name = args[2]
# model_animal_name = 'mmusculus'   #"rnorvegicus"

Rdata_PATH = paste0('Results/', INPUT_NAME, '/clusters/')
INPUT_FILES = list.files(path = Rdata_PATH , pattern = '.RData', full.names = T, include.dirs = T)
input_version = 1
INPUT_FILE = INPUT_FILES[input_version]
OUTPUT_NAME = gsub('.RData','',gsub(paste0(Rdata_PATH, '/clusters_'),'',INPUT_FILE ))
OUTPUT_NAME = gsub('_v2', '', OUTPUT_NAME)

## load gene symbol and ensemble ids:
input_from_10x <- paste0("Data/", INPUT_NAME,'/')
seur_genes_df <- read.delim(paste0(input_from_10x,'genes.tsv'), header = F)


load(INPUT_FILE)
cluster_names <- levels(seur)
Cluster_markers <- sapply(1:length(cluster_names), 
                          function(i) FindMarkers(seur, ident.1=cluster_names[i], ident.2 = NULL), 
                          simplify = FALSE)

names(Cluster_markers) <- paste0('cluster_', cluster_names)
Cluster_markers_merged <- sapply(1:length(Cluster_markers), 
                                 function(i){
                                   markers_df <- Cluster_markers[[i]]
                                   markers_df$ensemble_ids = rownames(markers_df)
                                   ## merge the ensemble IDs in the dataframe with the HUGO terms 
                                   markers_df_merged <- merge(markers_df, seur_genes_df, 
                                                              by.x='ensemble_ids', 
                                                              by.y='V1', all.x=T, all.y=F)
                                   return(markers_df_merged)
                                 }, simplify = FALSE)

names(Cluster_markers_merged) <- paste0('cluster_', cluster_names)
lapply(Cluster_markers_merged, head)
saveRDS(Cluster_markers_merged, paste0('Results/', INPUT_NAME, '/markers/markers_', OUTPUT_NAME, '.rds'))



