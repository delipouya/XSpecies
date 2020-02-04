
source('Codes/Functions.R')
Initialize()
PATH_TO_FILES = 'Data/McParland_markers/SUPPLEMENTARY_DATA/liver/'
SPECIES_NAME = 'rat_Rnor'


candidateGenes_mapped_df <- readRDS(paste0(PATH_TO_FILES,'candidateGenes_mapped_table.rds'))
candidateGenes_mapped <- lapply(candidateGenes_mapped_df, 
                                function(x) getUnemptyList(x$rnorvegicus_homolog_ensembl_gene))

seur <- readRDS(paste0('objects/',SPECIES_NAME,'/2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'))
exprMatrix <- as.matrix(seur[['RNA']]@data)


#### CHECK THE METHOD AND PARAMETERS OF EACH !!!!

## Gene set variation analysis
gsva_result <- gsva(exprMatrix, candidateGenes_mapped)
getHead(gsva_result)
gsva_result <- readRDS(paste0('Results/',SPECIES_NAME,'/GSE/gsva_result.rds'))


### Gene set enrichment analysis
fgseaRes <- fgseaLabel(candidateGenes_mapped, exprMatrix, 
                       as.numeric(as.factor(colnames(exprMatrix))), nperm = 1000)
saveRDS(fgseaRes, paste0('Results/',SPECIES_NAME,'/GSE/fgsea_result.rds'))



sum(is.na(as.numeric(as.factor(colnames(exprMatrix)))))
