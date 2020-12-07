options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

source('Codes/Functions.R')
Initialize()
INPUT_NAME = args[1] 
# INPUT_NAME = 'rat_Lew_01' # 'rat_DA 'mouse'
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
seur_genes_df <- read.delim(paste0(input_from_10x,'features.tsv.gz'), header = F)


### manual import
INPUT_NAME = 'rat_DA_M_10WK_003'
INPUT_FILE = 'Results/rat_DA_M_10WK_003/clusters/clusters_rat_DA_M_10WK_003_mito_20_lib_2000_v2.rds'
OUTPUT_NAME = 'rat_DA_M_10WK_003_mito_20_lib_2000_v2'
INPUT_FILE = paste0('Results/',INPUT_NAME,'/clusters/clusters_',OUTPUT_NAME,'.rds')
seur <- readRDS(INPUT_FILE)

# load(INPUT_FILE)


if(INPUT_NAME=='rat_Rnor') Idents(seur) <- paste0('cluster_',as.character(seur$SCT_snn_res.1))
if(INPUT_NAME=='rat_Lew_01') Idents(seur) <- paste0('cluster_',as.character(seur$SCT_snn_res.0.8))
if(INPUT_NAME %in% c('rat_DA', 'rat_Lew_02')) Idents(seur) <- paste0('cluster_',Idents(seur))
if(INPUT_NAME=='rat_DA_M_10WK_003') Idents(seur) <- paste0('cluster_',as.character(seur$SCT_snn_res.0.2))
if(INPUT_NAME %in% c('rat_DA_01_reseq')) Idents(seur) <- paste0('cluster_',as.character(seur@meta.data$SCT_snn_res.0.3))


cluster_names <-  levels(seur)
Cluster_markers <- sapply(1:length(cluster_names), 
                          function(i) FindMarkers(seur, ident.1=cluster_names[i]), 
                          simplify = FALSE)

names(Cluster_markers) <- cluster_names

Cluster_markers_merged <- sapply(1:length(Cluster_markers), 
                                 function(i){
                                   markers_df <- Cluster_markers[[i]]
                                   markers_df$ensemble_ids = rownames(markers_df)
                                   ## merge the ensemble IDs in the dataframe with the HUGO terms 
                                   markers_df_merged <- merge(markers_df, seur_genes_df, 
                                                              by.x='ensemble_ids', 
                                                              by.y='V1', all.x=T, all.y=F,sort=F)
                                   markers_df_merged2 <-  markers_df_merged[match(markers_df_merged$ensemble_ids, markers_df$ensemble_ids),]
                                   return(markers_df_merged)
                                 }, simplify = FALSE)

names(Cluster_markers_merged) <- cluster_names

### checking if the order of DE list has not changed due to merging
sapply(1:length(Cluster_markers), function(i){
  sum(Cluster_markers_merged[[i]]$ensemble_ids != rownames(Cluster_markers[[i]])) }, simplify = F)



# scp delaram@192.168.233.150:~/XSpecies/Results/rat_Lew_02/markers/markers_rat_Lew_02_mito_40_lib_2000.rds .
### Run this section in order to get the markers in a csv format
RES= ''
if(INPUT_NAME=='rat_Rnor') RES = '_res.1'
if(INPUT_NAME=='rat_Lew_01') RES = '_res.0.8'
if(INPUT_NAME=='rat_DA_01_reseq') RES = '_res.0.3'


saveRDS(Cluster_markers_merged, 
          paste0('Results/', INPUT_NAME, '/markers/markers_', OUTPUT_NAME, RES,'.rds'))
length(Cluster_markers_merged)
Cluster_markers_merged <- readRDS(paste0('Results/', INPUT_NAME,
                                              '/markers/markers_', OUTPUT_NAME, RES,'.rds'))


# scp delaram@192.168.233.150:~/XSpecies/Results/rat_DA_M_10WK_003/markers/markers_rat_DA_M_10WK_003_mito_20_lib_2000_v2/*.csv .
markers_dir = paste0('Results/', INPUT_NAME, '/markers/markers_',OUTPUT_NAME, RES)

dir.create(markers_dir)
for(i in 1:length(Cluster_markers_merged)){
  print(names(Cluster_markers_merged)[i])
  write.csv(Cluster_markers_merged[[i]], 
            paste0(markers_dir,'/',INPUT_NAME,'_', names(Cluster_markers_merged)[i],'.csv'), 
            row.names = T, quote = F)
}




