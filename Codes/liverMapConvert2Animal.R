source('Codes/Functions.R')
source('Codes/Convert_Human2ModelAnimal.R')
Initialize()


### importing the identified markers of liver

### Human Liver Atlas- 2019
all_cluster_files <- list.files('Data/Human_liver_Atlas', pattern = '*.csv', full.names  = T)
all_cluster_files_names <- list.files('Data/Human_liver_Atlas', pattern = '*.csv')
all_cluster_files_names <- as.character(unlist(lapply(all_cluster_files_names, 
                                                      function(x) as.numeric(substr(x, nchar(x)-5, nchar(x)-4)))))
LiverAtlasClusterMarkers <- lapply(all_cluster_files, function(aFile) read.csv(aFile, stringsAsFactors = F))
names(LiverAtlasClusterMarkers) <- paste0('Cluster_',all_cluster_files_names)
lapply(LiverAtlasClusterMarkers , head)
LiverAtlasClusterMarkers_geneNames <- lapply(LiverAtlasClusterMarkers, function(x) x$GeneSymbol)
candidateGenes <- lapply(LiverAtlasClusterMarkers_geneNames, function(x)head(x,10))
candidateGenes <- unique(unlist(candidateGenes))


### McParland liver map- 2018
path_to_markers <- 'Data/McParland_markers'
Known_markers <-list.files(path_to_markers, pattern = '.csv', full.names = T)
file_names <- list.files(path_to_markers, pattern = '.csv')
file_names <- substr(file_names,1, nchar(file_names)-4)
Markers <-lapply(Known_markers, function(input_file) read.csv(input_file, stringsAsFactors = F))
names(Markers) <- file_names
lapply(Markers, head)
# saveRDS(Markers, 'Data/McParland_markers/McParland_Markers.rds')
Markers <- readRDS('Data/McParland_markers/McParland_Markers.rds')
candidateGenes <- unique(unlist(lapply(Markers, function(x)x$Gene)))
candidateGenes <- Markers[['Cluster_1']]$Gene


######################## 
for (i in 1:length(Markers)){
  print('-----------------')
  cluster_id = names(Markers)[i]
  candidateGenes <- Markers[[i]]$Gene
  mappedGenesToOrthologs <- .getMapped_hs2model_df(wanted_attributes, ensembl, candidateGenes)
  # head(mappedGenesToOrthologs)
  write.csv(mappedGenesToOrthologs, 
            file = paste0('Data/McParland_markers/mapped_markers_hs2ratRnor/', cluster_id, '.csv'), row.names = T)
  
  print(paste0(cluster_id, ' is written!'))
}


#### thelist data structure
MappedGenesToOrthologs <- sapply(1:length(Markers), 
                                 function(i){
                                   cluster_id = names(Markers)[i]
                                   candidateGenes <- Markers[[i]]$Gene
                                   mappedGenesToOrthologs <- .getMapped_hs2model_df(wanted_attributes, ensembl, candidateGenes)
                                   return(mappedGenesToOrthologs)}, simplify = F)

names(MappedGenesToOrthologs) <- names(Markers)
lapply(MappedGenesToOrthologs, head)
saveRDS(MappedGenesToOrthologs , 'Data/McParland_markers/mapped_markers_hs2ratRnor/All_mapped_markers_hs2ratRnor_list.rds')



### Check how each of these columns have became
head(mappedGenesToOrthologs)
Orthologs <-  getUnemptyList(mappedGenesToOrthologs$rnorvegicus_homolog_associated_gene_name)
length(Orthologs)


