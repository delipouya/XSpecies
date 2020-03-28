## Attension:
## you need to set the animal name here before loading it in another script

source('Codes/Functions.R')
Initialize()

model_animal_name = 'rnorvegicus' # 'mmusculus'


.getMapped_hs2model_df <- function(ensembl, candidateGenes, model_animal_name){
  
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
# .getMapped_hs2model_df(ensembl, candidateGenes, model_animal_name)



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





####  Convert human gene list to ensembl id
get_human_ensembl_ids <- function(human_gene_symbol_list){
  listMarts()
  ensembl <- useMart("ensembl")
  datasets <- listDatasets(ensembl)
  ensembl = useDataset('hsapiens_gene_ensembl',
                       mart=ensembl)
  general_markers_human_ensembl_df <- getBM(filters="hgnc_symbol", 
                                            attributes= c('hgnc_symbol',"ensembl_gene_id"),
                                            values=as.character(human_gene_symbol_list), mart= ensembl)
  return(general_markers_human_ensembl_df)
}

get_rat_ensembl_ids <- function(rat_gene_symbol_list){
  listMarts()
  ensembl <- useMart("ensembl")
  datasets <- listDatasets(ensembl)
  ensembl = useDataset('rnorvegicus_gene_ensembl',
                       mart=ensembl)
  general_markers_rat_ensembl_df <- getBM(filters="rgd_symbol", 
                                            attributes= c('rgd_symbol',"ensembl_gene_id"),
                                            values=as.character(rat_gene_symbol_list), mart= ensembl)
  return(general_markers_rat_ensembl_df)
}

