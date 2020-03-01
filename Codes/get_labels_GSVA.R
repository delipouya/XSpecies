## Run this script as: 
# Rscript Codes/get_labels_AUCell.R 'rat_Rnor' '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)


source('Codes/Functions.R')
Initialize()

INPUT_NAME = args[1] 
INPUT_FILE = args[2]
# INPUT_NAME = 'rat_Rnor'
# INPUT_FILE = '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'
PATH_TO_FILES = 'Data/McParland_markers/SUPPLEMENTARY_DATA/liver/'
OUTPUT_NAME = gsub('.rds','',gsub('2.seur_dimRed_','',INPUT_FILE ))

candidateGenes_mapped_df <- readRDS(paste0(PATH_TO_FILES,'candidateGenes_mapped_table.rds'))
candidateGenes_mapped <- lapply(candidateGenes_mapped_df, 
                                function(x) getUnemptyList(x$rnorvegicus_homolog_ensembl_gene))

#seur <- readRDS(paste0('objects/',INPUT_NAME,'/',INPUT_FILE))
load('Results/rat_Rnor/clusters/clusters_rat_Rnor_mito_50_lib_1500_v2.RData')

#exprMatrix <- as.matrix(seur[['RNA']]@data)
exprMatrix <- as.matrix(GetAssayData(seur, assay.type = "SCT", slot = "data"))
colnames(exprMatrix) <- paste0('cluster_', as.character(seur$SCT_snn_res.1))
rownames(exprMatrix)

#### CHECK THE METHOD AND PARAMETERS OF EACH !!!!

## Gene set variation analysis
gsva_result <- gsva(exprMatrix, candidateGenes_mapped)
getHead(gsva_result)
saveRDS(gsva_result, paste0('Results/',INPUT_NAME,'/GSVA/GSVA_',OUTPUT_NAME,'.rds'))











