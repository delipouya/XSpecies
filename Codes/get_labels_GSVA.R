## Run this script as: 
# Rscript Codes/get_labels_AUCell.R 'rat_Rnor' '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)


source('Codes/Functions.R')
Initialize()

INPUT_NAME = args[1] 
# INPUT_NAME = 'rat_Rnor'

Rdata_PATH = paste0('Results/', INPUT_NAME, '/clusters/')
INPUT_FILES = list.files(path = Rdata_PATH , pattern = '.RData', full.names = T, include.dirs = T)
input_version = 1
INPUT_FILE = INPUT_FILES[input_version]

PATH_TO_FILES = 'Data/McParland_markers/SUPPLEMENTARY_DATA/liver/'
OUTPUT_NAME = gsub('.RData','',gsub(paste0(Rdata_PATH, '/clusters_'),'',INPUT_FILE ))
OUTPUT_NAME = gsub('_v2', '', OUTPUT_NAME)

load(INPUT_FILE)

candidateGenes_mapped_df <- readRDS(paste0(PATH_TO_FILES,'candidateGenes_mapped_table.rds'))
candidateGenes_mapped <- lapply(candidateGenes_mapped_df, 
                                function(x) getUnemptyList(x$rnorvegicus_homolog_ensembl_gene))


#exprMatrix <- as.matrix(seur[['RNA']]@data)
exprMatrix <- as.matrix(GetAssayData(seur, assay.type = "SCT", slot = "data"))
colnames(exprMatrix) <- paste0('cluster_', as.character(seur$SCT_snn_res.1))
rownames(exprMatrix)

#### CHECK THE METHOD AND PARAMETERS OF EACH !!!!

## Gene set variation analysis
gsva_result <- gsva(exprMatrix, candidateGenes_mapped)
getHead(gsva_result)
saveRDS(gsva_result, paste0('Results/',INPUT_NAME,'/GSVA/GSVA_',OUTPUT_NAME,'.rds'))











