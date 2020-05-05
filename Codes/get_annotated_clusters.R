source('Codes/Functions.R')
Initialize()

############################ 
### loading the input data matrix:

INPUT_NAME = 'rat_Lew_02'  # rat_Rnor 'rat_Lew_01', 'rat_DA' 'mouse' 
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
Hepatocytes_cl[Hepatocytes_cl %in% Total_markers_converted_df[['Hepatocytes']]$symbol]





# downloading the output pdf files:
# scp delaram@192.168.233.150:~/XSpecies/Results/rat_Lew_01/markers/cell_type_markers/*.pdf 

#### making the pdf file for all the kown markers,
#### a single pdf file will be made per cell-type
dir.create(paste0('Results/',INPUT_NAME,'/markers/cell_type_markers/'))
sample_name = 'Lew-02' #Lew-01
sapply(1:length(Total_markers_converted_df), 
       function(i){
         cell_type <- gsub(pattern = '_', replacement = ' ', Total_markers_names[i])
         a_mapped_markers_df <- Total_markers_converted_df[[i]]
         
         print(cell_type)
         print('################')
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

a_cluster = 'cluster_0'

if(INPUT_NAME=='rat_Rnor') Idents(seur) <- seur$SCT_snn_res.1
if(INPUT_NAME=='rat_Lew_01') Idents(seur) <- seur$SCT_snn_res.0.8
cluster_lists <- Idents(seur)

### visualize the particular cluster
tsne_df <- data.frame(getEmb(seur, 'tSNE'), clusters=paste0('cluster_',as.character(cluster_lists)))
ggplot(tsne_df, aes(x=tSNE_1,y=tSNE_2,color=clusters))+geom_point(alpha=0.8)+
  theme_classic()+scale_color_manual(values = colorPalatte[1:length(table(Idents(seur)))])

umap_df <- data.frame(getEmb(seur, 'umap'), clusters=paste0('cluster_',as.character(cluster_lists)))
ggplot(umap_df, aes(x=UMAP_1,y=UMAP_2,color=clusters))+geom_point(alpha=0.8)+
  theme_classic()+scale_color_manual(values = colorPalatte[1:length(table(Idents(seur)))])



tsne_df$clusters <- ifelse(tsne_df$clusters == a_cluster, a_cluster, 'other')
ggplot(tsne_df, aes(x=tSNE_1,y=tSNE_2,color=clusters))+geom_point()+
  theme_classic()+ggtitle(a_cluster)

umap_df$clusters <-ifelse(umap_df$clusters == a_cluster, a_cluster, 'other')
ggplot(umap_df, aes(x=UMAP_1,y=UMAP_2,color=clusters))+geom_point()+
  theme_classic()+ggtitle(a_cluster)


### Import the list of DE genes for each of the clusters
RES= ''
if(INPUT_NAME=='rat_Rnor') RES = '_res.1'
if(INPUT_NAME=='rat_Lew_01') RES = '_res.0.8'
Cluster_markers_merged <- readRDS(paste0('Results/', INPUT_NAME, '/markers/markers_', 
                                         OUTPUT_NAME, RES,'.rds'))
lapply(Cluster_markers_merged, head)
head(Cluster_markers_merged[[a_cluster]])
Cluster_markers_merged[[a_cluster]]$V2[1:50]


### generate a barplot for this
### Check what percentage of each of DEs are in cell type markers

for(i in 1:length(Total_markers_converted_df)){
  
  a_cluster_DE_genes <- Cluster_markers_merged[[a_cluster]]$ensemble_ids
  cell_type_markers <- Total_markers_converted_df[[i]]$rnorvegicus_homolog_ensembl_gene
  
  DE_in_known_markers <- round(sum(a_cluster_DE_genes %in% cell_type_markers)/length(a_cluster_DE_genes),2)
  print(paste0('percent of DEs in ', names(Total_markers_converted_df)[i], '-markers: ',DE_in_known_markers))
  
  known_markers_in_DE <- round(sum(cell_type_markers %in% a_cluster_DE_genes)/length(cell_type_markers),2)
  print(paste0('percent of ', names(Total_markers_converted_df)[i], '-markers in DEs: ', known_markers_in_DE))
  print('-----------------------')
}


### find the DE genes which are present in the most dominant cell-type markers
a_cluster_DE_genes <- Cluster_markers_merged[[a_cluster]]$ensemble_ids
cell_type_markers <- Total_markers_converted_df$KCs$rnorvegicus_homolog_ensembl_gene
Cluster_markers_merged[[a_cluster]]$V2[a_cluster_DE_genes %in% cell_type_markers]






###########################################################################
#################################### define manual markers to check
# Hba-a2, Hbb, Hba-a3

seur_genes_df <- read.delim(paste0(paste0("Data/", INPUT_NAME,'/'),'genes.tsv'), header = F)
cell_type = 'Hepatocyte'
HEMATOPOIETIC = c('AMPD3'	,'HBB'	,'RHAG'	,'RHD')

### manual tSNE for human marker
manual_marker = 'G6PC'
a_mapped_markers_df <- .getMapped_hs2model_df(ensembl, manual_marker, model_animal_name)
j = 1
aMarker_symbol = a_mapped_markers_df[j,]$rnorvegicus_homolog_associated_gene_name
aMarker_ensemble_id_index = a_mapped_markers_df[j,]$rnorvegicus_homolog_ensembl_gene
aMarker_human_symbol = a_mapped_markers_df[j,]$symbol
aMarker_expression <- exprMatrix[aMarker_ensemble_id_index,]
tsne_df <- as.data.frame(cbind(getEmb(seur, 'tsne'),gene.expr=aMarker_expression))
Plot.tsne.gene.expr(tsne_df, 
                    title = paste0('human:', aMarker_human_symbol,'  ', 
                                   model_animal_name,':',aMarker_symbol),
                    subtitle = paste0('Cell type: ',cell_type, '   Sample: ', sample_name))

### manual tSNE for rat marker
aMarker_symbol = 'Hamp'
aMarker_expression <- exprMatrix[seur_genes_df$V1[seur_genes_df$V2 == aMarker_symbol],]
tsne_df <- as.data.frame(cbind(getEmb(seur, 'tsne'),gene.expr=aMarker_expression))
Plot.tsne.gene.expr(tsne_df, 
                    title = paste0( model_animal_name,':',aMarker_symbol),
                    subtitle = paste0('Cell type: ',cell_type, '   Sample: ', sample_name))





############################################################## checking the DE markers freq in known markers
Cluster_markers_merged <- readRDS(paste0('Results/', INPUT_NAME,'/markers/markers_', OUTPUT_NAME, RES,'.rds'))

##############################   figuring out the general populations
cluster_to_check = 'cluster_14'
Cluster_markers_merged[[cluster_to_check]]$V2[1:20]

for(i in 1:length(Total_markers_converted_df)){
  print(names(Total_markers_converted_df)[i])
  rat_cell_markers <- Total_markers_converted_df[[i]]$rnorvegicus_homolog_associated_gene_name
  DE_genes <- Cluster_markers_merged[[cluster_to_check]]$V2[1:20]
  print(DE_genes[DE_genes%in% rat_cell_markers])
}


##############################   figuring out the sub-populations
################ immune cells - T, B, NK cells 
T_cell_marker_list_rat <- lapply(T_cell_marker_list, function(x) 
  Total_markers_converted_df$T_cells$rnorvegicus_homolog_associated_gene_name[Total_markers_converted_df$T_cells$symbol %in% x])

B_cell_marker_list_rat <- lapply(B_cell_marker_list, function(x) 
  Total_markers_converted_df$B_cells$rnorvegicus_homolog_associated_gene_name[Total_markers_converted_df$B_cells$symbol %in% x])

for(i in 1:length(T_cell_marker_list_rat)){
  print(names(T_cell_marker_list_rat)[i])
  rat_cell_markers <- T_cell_marker_list_rat[[i]]
  DE_genes <- Cluster_markers_merged[[cluster_to_check]]$V2[1:20]
  print(DE_genes[DE_genes%in% rat_cell_markers])
}

for(i in 1:length(B_cell_marker_list_rat)){
  print(names(B_cell_marker_list_rat)[i])
  rat_cell_markers <- B_cell_marker_list_rat[[i]]
  DE_genes <- Cluster_markers_merged[[cluster_to_check]]$V2[1:20]
  print(DE_genes[DE_genes%in% rat_cell_markers])
}








