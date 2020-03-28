## Run this script as: 
# Rscript Codes/get_expression_tsne.R 'rat_Rnor' '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
INPUT_NAME = args[1] 
INPUT_FILE = args[2]

source('Codes/Functions.R')
Initialize()
# rat_Rnor, 2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds
# rat_DA, 2.seur_dimRed_rat_DA_mito_30_lib_1500.rds
# INPUT_NAME = 'rat_Rnor'
# INPUT_FILE = '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'

PATH_TO_FILES = 'Data/McParland_markers/SUPPLEMENTARY_DATA/liver/'
OUTPUT_NAME = gsub('.rds','',gsub('2.seur_dimRed_','',INPUT_FILE ))
model_animal_name = "rnorvegicus"    # 'mmusculus'

seur <- readRDS(paste0('objects/',INPUT_NAME,'/',INPUT_FILE))
exprMatrix <- GetAssayData(seur)



### importing expression matrix and list of markers
#### Sonya markers: 
marker_file <- paste0('candidateGenes_mapped_table_',model_animal_name,'.rds')
candidateGenes_mapped_df <- readRDS(paste0(PATH_TO_FILES, marker_file))
candidateGenes_mapped <- lapply(candidateGenes_mapped_df, function(x) getUnemptyList(x$rnorvegicus_homolog_ensembl_gene))


pdf(paste0("Results/",INPUT_NAME,'/marker/marker_expression_',OUTPUT_NAME,'.pdf'))

for(cell_index in 1:length(candidateGenes_mapped)){
  for(marker_index in 1:length(candidateGenes_mapped[[cell_index]])){
    
    aMarker <- candidateGenes_mapped[[cell_index]][marker_index]
    aMarker_ensemble_id_index <- candidateGenes_mapped_df[[cell_index]]$rnorvegicus_homolog_ensembl_gene == aMarker
    aMarker_symbol = candidateGenes_mapped_df[[cell_index]]$rnorvegicus_homolog_associated_gene_name[aMarker_ensemble_id_index]
    aMarker_expression <- exprMatrix[aMarker,]
    
    tsne_df <- as.data.frame(cbind(getEmb(seur, 'tsne'),gene.expr=aMarker_expression))
    
    p= Plot.tsne.gene.expr(tsne_df, paste0(names(candidateGenes_mapped)[cell_index],' - ', aMarker_symbol))
    print(p)
  }
}

dev.off()



## regenerate the marker distribution

#### Tallulah markers

### importing input markers
general_markers_df <- readRDS('Data/McParland_markers/liver_markers_tallulah/liver_general_markers.rds')
### removing the unmapped markers
general_markers_df$rat_general <- general_markers_df$rat_general[getUnemptyList_bool(general_markers_df$rat_general$ensembl_gene_id) &
                                                                   general_markers_df$rat_general$ensembl_gene_id %in% rownames(seur),]

general_markers_df$human_general <- general_markers_df$human_general[getUnemptyList_bool(general_markers_df$human_general$rnorvegicus_homolog_ensembl_gene) &
                                                                       general_markers_df$human_general$rnorvegicus_homolog_ensembl_gene %in% rownames(seur),]



pdf(paste0("Results/",INPUT_NAME,'/markers/tallulah_general_markers_',INPUT_NAME,'.pdf'))
# >>>>>>>> rat markers
marker_df <- general_markers_df$rat_general
for( i in 1:nrow(marker_df)){
  
  aMarker = marker_df$ensembl_gene_id[i]
  aMarker_expression <- exprMatrix[aMarker,]
  print(paste0('rat: ', i))
  tsne_df <- as.data.frame(cbind(getEmb(seur, 'tsne'),gene.expr=aMarker_expression))
  
  p = Plot.tsne.gene.expr(tsne_df ,
                          marker_df$Gene[i] ,
                          paste0(marker_df$Species[i], ' - ',marker_df$Celltype[i], ' - ', INPUT_NAME))
  print(p)}
# >>>>>>>> human markers
marker_df <- general_markers_df$human_general
for( i in 1:nrow(marker_df)){
  
  aMarker = marker_df$rnorvegicus_homolog_ensembl_gene[i]
  aMarker_expression <- exprMatrix[aMarker,]
  print(paste0('human: ', i))
  tsne_df <- as.data.frame(cbind(getEmb(seur, 'tsne'),gene.expr=aMarker_expression))
  
  p = Plot.tsne.gene.expr(tsne_df,
                          paste0('Human: ', marker_df$Gene[i], '  rat: ', marker_df$rnorvegicus_homolog_associated_gene_name[i]),
                          paste0(marker_df$Species[i], ' - ',marker_df$Celltype[i], ' - ', INPUT_NAME))
  print(p)}

dev.off()

  





