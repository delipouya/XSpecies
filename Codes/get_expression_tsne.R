## Run this script as: 
# Rscript Codes/get_expression_tsne.R 'rat_Rnor' '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'


options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

source('Codes/Functions.R')
Initialize()

# INPUT_NAME = 'rat_Rnor'
# INPUT_FILE = '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'
INPUT_NAME = args[1] 
INPUT_FILE = args[2]
PATH_TO_FILES = 'Data/McParland_markers/SUPPLEMENTARY_DATA/liver/'
OUTPUT_NAME = gsub('.rds','',gsub('2.seur_dimRed_','',INPUT_FILE ))
AUCell_dir = paste0("Results/",INPUT_NAME,"/AUCell/")
model_animal_name = "rnorvegicus"    # 'mmusculus'


### importing expression matrix and list of markers
marker_file <- paste0('candidateGenes_mapped_table_',model_animal_name,'.rds')
candidateGenes_mapped_df <- readRDS(paste0(PATH_TO_FILES, marker_file))
candidateGenes_mapped <- lapply(candidateGenes_mapped_df, function(x) getUnemptyList(x$rnorvegicus_homolog_ensembl_gene))


seur <- readRDS(paste0('objects/',INPUT_NAME,'/',INPUT_FILE))
exprMatrix <- GetAssayData(seur)


pdf(paste0("Results/",INPUT_NAME,'/Marker/marker_expression_',OUTPUT_NAME,'.pdf'))

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










########################################
pdf(paste0(marker_dir, '.pdf'))
for(i in 1:length(human_markers)){
  human_marker = 'VCAN'
  human_marker <- human_markers[i]
  print('------------------')
  print(human_marker)
  df <- .getMapped_hs2model_df(ensembl, human_marker , model_animal_name)
  
  if(is.na(df)){
    print(paste0(human_marker, ' or ', aMarker_ensemble_id_index, ' cant be mapped!'))
    next
  }
  aMarker_ensemble_id_index <- df$rnorvegicus_homolog_ensembl_gene
  
  
  if(!aMarker_ensemble_id_index %in% rownames(seur)){
    print(paste0(human_marker, ' or ', aMarker_ensemble_id_index, ' not available!'))
    next} 
  aMarker_symbol = df$rnorvegicus_homolog_associated_gene_name
  
  
  aMarker_expression <- exprMatrix[aMarker_ensemble_id_index,]
  tsne_df <- as.data.frame(cbind(getEmb(seur, 'tsne'),gene.expr=aMarker_expression))
  p=Plot.tsne.gene.expr(tsne_df, paste0(cell_type,' - ', aMarker_symbol))
  print(p)
}
dev.off()




