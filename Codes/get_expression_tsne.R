
source('Codes/Functions.R')
Initialize()
PATH_TO_FILES = 'Data/McParland_markers/SUPPLEMENTARY_DATA/liver/'
INPUT_NAME = 'rat_Rnor'

INPUT_FILE = '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'
OUTPUT_NAME = gsub('.rds','',gsub('2.seur_dimRed_','',INPUT_FILE ))



candidateGenes_mapped_df <- readRDS(paste0(PATH_TO_FILES,'candidateGenes_mapped_table.rds'))
candidateGenes_mapped <- lapply(candidateGenes_mapped_df, function(x) getUnemptyList(x$rnorvegicus_homolog_ensembl_gene))

### importing expression matrix and list of markers
candidateGenes_mapped_df <- readRDS(paste0(PATH_TO_FILES,'candidateGenes_mapped_table.rds'))
candidateGenes_mapped <- lapply(candidateGenes_mapped_df, 
                                function(x) getUnemptyList(x$rnorvegicus_homolog_ensembl_gene))

seur <- readRDS(paste0('objects/',INPUT_NAME,'/',INPUT_FILE))
exprMatrix <- as.matrix(seur[['RNA']]@data)

