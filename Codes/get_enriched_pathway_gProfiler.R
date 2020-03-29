source('Codes/Functions.R')
Initialize()


get_gprofiler_enrich <- function(markers, model_animal_name){
  gostres <- gost(query = markers,
                  ordered_query = TRUE, exclude_iea =TRUE, 
                  sources=c('GO:BP' ,'REAC'),
                  organism = model_animal_name)
  return(gostres)
}

#### importing the input data
INPUT_NAME = 'rat_Rnor'  # 'rat_Rnor''rat_Lew_01', 'rat_DA' 'mouse' 
model_animal_name = "rnorvegicus" #'mmusculus' 

Rdata_PATH = paste0('Results/', INPUT_NAME, '/clusters/')
INPUT_FILES = list.files(path = Rdata_PATH , pattern = '.RData', full.names = T, include.dirs = T)
input_version = 1
INPUT_FILE = INPUT_FILES[input_version]
OUTPUT_NAME = gsub('.RData','',gsub(paste0(Rdata_PATH, '/clusters_'),'',INPUT_FILE ))
OUTPUT_NAME = gsub('_v2', '', OUTPUT_NAME)

RES= ''
if(INPUT_NAME=='rat_Rnor') RES = '_res.1'
if(INPUT_NAME=='rat_Lew_01')  RES = '_res.0.8' 

Cluster_markers_merged <- readRDS(paste0('Results/', INPUT_NAME, '/markers/markers_', 
                                         OUTPUT_NAME, RES,'.rds'))
cluster_names <- names(Cluster_markers_merged)
lapply(Cluster_markers_merged, dim)



sapply(1:length(Cluster_markers_merged), get_gprofiler_enrich(x[[i]]$ensemble_ids, model_animal_name)
gostres$result = gostres$result[gostres$result$query_size>5 & gostres$result$query_size<350 & 
                                  gostres$result$intersection_size >3, ]
head(gostres$result,20)


p <- gostplot(gostres, capped = T, interactive = FALSE)
# pt2 <- publish_gosttable(gostres, use_colors = TRUE, filename = NULL)


