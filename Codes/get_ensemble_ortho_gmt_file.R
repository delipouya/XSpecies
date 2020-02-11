## In this script we'll importing the liver signatures gmt file 
## then, we'll convert them to ensemble id and write them to a new gmt file
##

source('Codes/Functions.R')
Initialize()
PATH_TO_FILES = 'Data/McParland_markers/SUPPLEMENTARY_DATA/liver/'
INPUT_NAME = 'rat_Rnor'


## > DO NOT RUN THIS AGAIN
## --------------------------------
#### generating the gmt file for ensemble id of liver cell markers
### this desction needs to be run once 

source('Codes/convert_human_to_ortholog_functions.R')

files <- list.files(PATH_TO_FILES, include.dirs = T, full.names = T, pattern = '*.tsv')
liver_files <- lapply(files, function(x) read.delim(x))
lapply(liver_files, head)

liver_marker_gene_set_df = read.table(paste0(PATH_TO_FILES,'liver_cell_type_signature_gene_sets.gmt'), fill=T)
rownames(liver_marker_gene_set_df) <- liver_marker_gene_set_df[,1]

liver_marker_gene_set_df <- liver_marker_gene_set_df[,-c(1,2)]
liver_marker_gene_set <- lapply(1:nrow(liver_marker_gene_set_df), 
                                function(i) {
                                  cell_type_markers <- as.character(liver_marker_gene_set_df[i,])
                                  cell_type_markers[cell_type_markers!='']} )

names(liver_marker_gene_set) <- rownames(liver_marker_gene_set_df)

candidateGenes <- getUnemptyList(as.character(unlist(liver_marker_gene_set)))
candidateGenes_mapped_df <- lapply(liver_marker_gene_set, 
                                   function(candidateGenes) 
                                     .getMapped_hs2model_df(wanted_attributes, ensembl, candidateGenes))

saveRDS(candidateGenes_mapped_df, paste0(PATH_TO_FILES,'candidateGenes_mapped_table.rds'))
candidateGenes_mapped <- lapply(candidateGenes_mapped_df, function(x) getUnemptyList(x$rnorvegicus_homolog_ensembl_gene))



sink(paste0(PATH_TO_FILES,"liver_cell_type_signature_gene_sets_ensemble.gmt"))
for(i in 1:length(candidateGenes_mapped)){
  marker <- candidateGenes_mapped[[i]]
  cell_type_name <- names(candidateGenes_mapped)[i]
  
  cat(cell_type_name)
  cat('\t')
  cat(cell_type_name)
  cat('\t')
  for (j in 1:length(marker)){
    cat (marker[j])
    cat('\t')
  }
  cat('\n')
}
sink()



