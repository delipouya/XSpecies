install.packages('SCINA')
library('SCINA')


source('Codes/Functions.R')
Initialize()
PATH_TO_FILES = 'Data/McParland_markers/SUPPLEMENTARY_DATA/liver/'
SPECIES_NAME = 'rat_Rnor'

candidateGenes_mapped_df <- readRDS(paste0(PATH_TO_FILES,'candidateGenes_mapped_table.rds'))
candidateGenes_mapped <- lapply(candidateGenes_mapped_df, 
                                function(x) getUnemptyList(x$rnorvegicus_homolog_ensembl_gene))

seur <- readRDS(paste0('objects/',SPECIES_NAME,'/2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'))
exprMatrix <- as.matrix(seur[['RNA']]@data)


SCINA_res = SCINA(exprMatrix, candidateGenes_mapped, max_iter = 100, convergence_n = 10, 
                convergence_rate = 0.7, sensitivity_cutoff = 0.7, 
                rm_overlap=FALSE, allow_unknown=FALSE, log_file='SCINA.log')

rownames(SCINA_res$probabilities)
rownames(gsva_result)
View(SCINA_res$cell_labels)
View(SCINA_res$probabilities)


saveRDS(SCINA_res, paste0('Results/',SPECIES_NAME,'/SCINA/SCINA_result.rds'))

