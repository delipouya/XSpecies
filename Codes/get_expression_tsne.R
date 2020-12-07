## Run this script as: 
# Rscript Codes/get_expression_tsne.R 'rat_Rnor' '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
INPUT_NAME = args[1] 
INPUT_FILE = args[2]

source('Codes/Functions.R')
source('Codes/convert_human_to_ortholog_functions.R')
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

aMarker_symbol = 'Hp'
seur_data = GetAssayData(seur)
aMarker_expression <- seur_data[aMarker_symbol,]
tsne_df <- as.data.frame(cbind(getEmb(seur, 'tsne'),gene.expr=aMarker_expression))
p= Plot.tsne.gene.expr(tsne_df, aMarker_symbol)
p
