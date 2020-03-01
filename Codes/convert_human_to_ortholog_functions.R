source('Codes/Functions.R')
Initialize()
# 'rnorvegicus' 
# model_animal_name = 'mmusculus'


.getMapped_hs2model_df <- function(wanted_attributes, ensembl, candidateGenes, model_animal_name){
  
  wanted_attributes <- c(paste0(model_animal_name, '_homolog_ensembl_gene'), 
                         paste0(model_animal_name, '_homolog_associated_gene_name'), 
                         paste0(model_animal_name, '_homolog_orthology_type'),
                         paste0(model_animal_name, '_homolog_perc_id'), 
                         paste0(model_animal_name, '_homolog_perc_id_r1'), 
                         paste0(model_animal_name, '_homolog_dn'), 
                         paste0(model_animal_name, '_homolog_ds') , 
                         paste0(model_animal_name, '_homolog_orthology_confidence'))
  
  
  mappedGenesToOrthologs <- getBM(filters="hgnc_symbol", 
                                  attributes= c("ensembl_gene_id", wanted_attributes),
                                  values=candidateGenes, mart= ensembl)
  ## cleaning the resulting data frame
  hgnc_symbol_to_ensembl <- data.frame(getBM(filters="hgnc_symbol", 
                                             attributes= c("ensembl_gene_id", 'hgnc_symbol'),
                                             values=candidateGenes, mart= ensembl))
  print(dim(mappedGenesToOrthologs))
  print(length(candidateGenes))
  print(dim(hgnc_symbol_to_ensembl))
  
  candidateGenes.df <- data.frame(symbol=candidateGenes)
  candidateGenes.df <- merge(candidateGenes.df, hgnc_symbol_to_ensembl, by.x='symbol', by.y='hgnc_symbol',all.x=T)
  mappedGenesToOrthologs <- merge(candidateGenes.df, mappedGenesToOrthologs, by.x='ensembl_gene_id', by.y='ensembl_gene_id', all.x=T )
  
  return(mappedGenesToOrthologs)
}
# usage:
# .getMapped_hs2model_df(wanted_attributes, ensembl, candidateGenes, )



#### #### #### #### converting the list of genes to their accourding orthologs
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset('hsapiens_gene_ensembl',mart=ensembl)

## checking the appropiate filter and attribute to use
listFilters(ensembl)[grep(listFilters(ensembl)[,1], pattern = 'symbol'),]
listAttributes(ensembl)[grep(listAttributes(ensembl)[,1], pattern = 'rnorvegicus'),] 
listAttributes(ensembl)[grep(listAttributes(ensembl)[,1], pattern = 'mmusculus'),] 
listAttributes(ensembl)[grep(listAttributes(ensembl)[,1], pattern = 'symbol'),] 









