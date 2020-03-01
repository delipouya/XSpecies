## Run this script as: 
# Rscript Codes/get_labels_SCINA.R 'rat_Rnor' '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)


source('Codes/Functions.R')
Initialize()

INPUT_NAME = args[1] 
# INPUT_NAME = 'rat_Rnor'
model_animal_name = args[2]
# model_animal_name = 'mmusculus'


Rdata_PATH = paste0('Results/', INPUT_NAME, '/clusters/')
INPUT_FILES = list.files(path = Rdata_PATH , pattern = '.RData', full.names = T, include.dirs = T)
input_version = 1
INPUT_FILE = INPUT_FILES[input_version]

PATH_TO_FILES = 'Data/McParland_markers/SUPPLEMENTARY_DATA/liver/'
OUTPUT_NAME = gsub('.RData','',gsub(paste0(Rdata_PATH, '/clusters_'),'',INPUT_FILE ))
OUTPUT_NAME = gsub('_v2', '', OUTPUT_NAME)
load(INPUT_FILE)

## --------------------------------
##### Check if all the known markers are in the raw dataset
candidateGenes_mapped_df <- readRDS(paste0(PATH_TO_FILES,'candidateGenes_mapped_table_', model_animal_name, '.rds'))
candidateGenes_mapped <- lapply(candidateGenes_mapped_df, function(x) getUnemptyList(x[[paste0(model_animal_name, '_homolog_ensembl_gene')]]))

ensemble_genes <- read.delim(paste0('~/XSpecies/Data/',INPUT_NAME,'/genes.tsv'), header = F)
lapply(candidateGenes_mapped, function(x) sum(x %in% ensemble_genes$V1/length(x)))
lapply(candidateGenes_mapped, function(x) x[!(x %in% ensemble_genes$V1)])



Rdata_PATH = paste0('Results/', INPUT_NAME, '/clusters/')
INPUT_FILES = list.files(path = Rdata_PATH , pattern = '.RData', full.names = T, include.dirs = T)
input_version = 1
INPUT_FILE = INPUT_FILES[input_version]

PATH_TO_FILES = 'Data/McParland_markers/SUPPLEMENTARY_DATA/liver/'
OUTPUT_NAME = gsub('.RData','',gsub(paste0(Rdata_PATH, '/clusters_'),'',INPUT_FILE ))
OUTPUT_NAME = gsub('_v2', '', OUTPUT_NAME)
load(INPUT_FILE)

exprMatrix <- as.matrix(GetAssayData(seur))


SCINA_res = SCINA(exprMatrix, candidateGenes_mapped, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.99, sensitivity_cutoff = 1, 
                rm_overlap=FALSE, allow_unknown=TRUE, log_file='SCINA.log')

rownames(SCINA_res$probabilities)
table(SCINA_res$cell_labels)
saveRDS(SCINA_res, paste0('Results/',INPUT_NAME,'/SCINA/','SCINA_',OUTPUT_NAME,'.rds'))




