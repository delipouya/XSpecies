source('Codes/Functions.R')
Initialize()

############################ 
### loading the input data matrix:

INPUT_NAME = 'rat_Rnor'  #'mouse' 
model_animal_name = "rnorvegicus" #'mmusculus'  

Rdata_PATH = paste0('Results/', INPUT_NAME, '/clusters/')
INPUT_FILES = list.files(path = Rdata_PATH , pattern = '.RData', full.names = T, include.dirs = T)
input_version = 1
INPUT_FILE = INPUT_FILES[input_version]
OUTPUT_NAME = gsub('.RData','',gsub(paste0(Rdata_PATH, '/clusters_'),'',INPUT_FILE ))
OUTPUT_NAME = gsub('_v2', '', OUTPUT_NAME)
load(INPUT_FILE)
exprMatrix <- GetAssayData(seur)
getHead(exprMatrix)


### importing the input markers and removing the 
Total_markers_converted_df <- readRDS('Data/Total_markers_converted_df.rds')
Total_markers_names <- names(Total_markers_converted_df)

### markers which aren't present in the expression matrix
Total_markers_converted_df <- sapply(1:length(Total_markers_converted_df),
                                     function(i){
                                       a_mapped_markers_df <- Total_markers_converted_df[[i]]
                                       a_mapped_markers_df <- a_mapped_markers_df[a_mapped_markers_df$rnorvegicus_homolog_ensembl_gene %in% rownames(seur),]
                                       return(a_mapped_markers_df)
                                     }, simplify = F)

names(Total_markers_converted_df) <- Total_markers_names
plasma_cells_cl[plasma_cells_cl %in% Total_markers_converted_df[['B_cells']]$symbol]



#### making the pdf file for all the kown markers,
#### a single pdf file will be made per cell-type
sample_name = 'DA-02'
sapply(1:length(Total_markers_converted_df), 
       function(i){
         cell_type <- gsub(pattern = '_', replacement = ' ', Total_markers_names[i])
         a_mapped_markers_df <- Total_markers_converted_df[[i]]
         
         pdf(paste0('Results/',INPUT_NAME,'/markers/cell_type_markers/',Total_markers_names[i], '_total_markers.pdf'))
         for(j in 1:nrow(a_mapped_markers_df)){
           
           aMarker_symbol = a_mapped_markers_df[j,]$rnorvegicus_homolog_associated_gene_name
           aMarker_ensemble_id_index = a_mapped_markers_df[j,]$rnorvegicus_homolog_ensembl_gene
           aMarker_human_symbol = a_mapped_markers_df[j,]$symbol
           
           aMarker_expression <- exprMatrix[aMarker_ensemble_id_index,]
           tsne_df <- as.data.frame(cbind(getEmb(seur, 'tsne'),gene.expr=aMarker_expression))
           p = Plot.tsne.gene.expr(tsne_df, 
                                   title = paste0('human:', aMarker_human_symbol,'  ', 
                                                  model_animal_name,':',aMarker_symbol),
                                   subtitle = paste0('Cell type: ',cell_type, '   Sample: ', sample_name))
           print(p)
         }
         dev.off()
         
       }, simplify = F)



##################################################
#### Choose a cluster for evaluation
#### get the list of markers for each of the clusters
## table(seur$SCT_snn_res.1)

a_cluster = 'cluster_12'

if(INPUT_NAME=='rat_Rnor') Idents(seur) <- seur$SCT_snn_res.1
cluster_lists <- Idents(seur)

### visualize the particular cluster
tsne_df <- data.frame(getEmb(seur, 'tSNE'), clusters=paste0('cluster_',as.character(cluster_lists)))
#ggplot(tsne_df, aes(x=tSNE_1,y=tSNE_2,color=clusters))+geom_point(alpha=0.8)+
#  theme_classic()+scale_color_manual(values = colorPalatte)

tsne_df$clusters <- ifelse(tsne_df$clusters == a_cluster, a_cluster, 'other')
ggplot(tsne_df, aes(x=tSNE_1,y=tSNE_2,color=clusters))+geom_point()+
  theme_classic()+ggtitle(a_cluster)


### Import the list of DE genes for each of the clusters
RES= ''
if(INPUT_NAME=='rat_Rnor') RES = '_res.1'
Cluster_markers_merged <- readRDS(paste0('Results/', INPUT_NAME, '/markers/markers_', 
                                         OUTPUT_NAME, RES,'.rds'))
lapply(Cluster_markers_merged, head)
head(Cluster_markers_merged[[a_cluster]])
Cluster_markers_merged[[a_cluster]]$V2[1:50]



### Check what percentage of each of DEs are in cell type markers
for(i in 1:length(Total_markers_converted_df)){
  
  a_cluster_DE_genes <- Cluster_markers_merged[[a_cluster]]$ensemble_ids
  cell_type_markers <- Total_markers_converted_df[[i]]$rnorvegicus_homolog_ensembl_gene
  
  print(paste0('percent of DEs in ', names(Total_markers_converted_df)[i], '-markers: ', 
               round(sum(a_cluster_DE_genes %in% cell_type_markers)/length(a_cluster_DE_genes),2)))
  print(paste0('percent of ', names(Total_markers_converted_df)[i], '-markers in DEs: ', 
               round(sum(cell_type_markers %in% a_cluster_DE_genes)/length(cell_type_markers),2)))
  print('-----------------------')
}


### find the DE genes which are present in the most dominant cell-type markers
a_cluster_DE_genes <- Cluster_markers_merged[[a_cluster]]$ensemble_ids
cell_type_markers <- Total_markers_converted_df$KCs$rnorvegicus_homolog_ensembl_gene
Cluster_markers_merged[[a_cluster]]$V2[a_cluster_DE_genes %in% cell_type_markers]




