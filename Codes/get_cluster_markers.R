options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

source('Codes/Functions.R')
Initialize()
INPUT_NAME = args[1] 
# INPUT_NAME = 'rat_Rnor' #'mouse'
model_animal_name = args[2]
# model_animal_name = "rnorvegicus" #'mmusculus' 

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


if(INPUT_NAME=='rat_Rnor') Idents(seur) <- as.character(seur$SCT_snn_res.1)
cluster_names <-  levels(seur)
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


### Run this section in order to get the markers in a csv format
RES= ''
if(INPUT_NAME=='rat_Rnor') RES = '_res.1'
# saveRDS(Cluster_markers_merged, 
#          paste0('Results/', INPUT_NAME, '/markers/markers_', OUTPUT_NAME, RES,'.rds'))
length(Cluster_markers_merged)
Cluster_markers_merged <- readRDS(paste0('Results/', INPUT_NAME,
                                              '/markers/markers_', OUTPUT_NAME, RES,'.rds'))



markers_dir = paste0('Results/', INPUT_NAME, '/markers/markers_',OUTPUT_NAME, RES)

dir.create(markers_dir)
for(i in 1:length(Cluster_markers_merged)){
  print(names(Cluster_markers_merged)[i])
  write.csv(Cluster_markers_merged[[i]], 
            paste0(markers_dir,'/',INPUT_NAME,'_', names(Cluster_markers_merged)[i],'.csv'), 
            row.names = T, quote = F)
}

