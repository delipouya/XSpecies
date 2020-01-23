source('Codes/Functions.R')
Initialize()

### importing the identified markers in the human liver atlas project
all_cluster_files <- list.files('Data/Human_liver_Atlas', pattern = '*.csv', full.names  = T)
all_cluster_files_names <- list.files('Data/Human_liver_Atlas', pattern = '*.csv')
all_cluster_files_names <- as.character(unlist(lapply(all_cluster_files_names, 
                                                      function(x) as.numeric(substr(x, nchar(x)-5, nchar(x)-4)))))
LiverAtlasClusterMarkers <- lapply(all_cluster_files, function(aFile) read.csv(aFile, stringsAsFactors = F))
names(LiverAtlasClusterMarkers) <- paste0('Cluster_',all_cluster_files_names)
lapply(LiverAtlasClusterMarkers , head)
LiverAtlasClusterMarkers_geneNames <- lapply(LiverAtlasClusterMarkers, function(x) x$GeneSymbol)
candidateGenes <- lapply(LiverAtlasClusterMarkers_geneNames, function(x)head(x,10))
candidateGenes <- unique(unlist(candidateGenes))

#### converting the list of genes to their accourding orthologs
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset('hsapiens_gene_ensembl',mart=ensembl)

## checking the appropiate filter and attribute to use
listFilters(ensembl)[grep(listFilters(ensembl)[,1], pattern = 'symbol'),]
listAttributes(ensembl)[grep(listAttributes(ensembl)[,1], pattern = 'rnorvegicus'),] 
listAttributes(ensembl)[grep(listAttributes(ensembl)[,1], pattern = 'symbol'),] 

wanted_attributes <- c('rnorvegicus_homolog_ensembl_gene', 'rnorvegicus_homolog_associated_gene_name', 'rnorvegicus_homolog_orthology_type',
                       'rnorvegicus_homolog_perc_id', 'rnorvegicus_homolog_perc_id_r1', 'rnorvegicus_homolog_dn', 'rnorvegicus_homolog_ds' , 
                       'rnorvegicus_homolog_orthology_confidence')

mappedGenesToOrthologs <- getBM(filters="hgnc_symbol", 
                         attributes= c("ensembl_gene_id", wanted_attributes),
                         values=candidateGenes, mart= ensembl)

## cleaning the resulting data frame
hgnc_symbol_to_ensembl <- data.frame(getBM(filters="hgnc_symbol", 
                                attributes= c("ensembl_gene_id", 'hgnc_symbol'),
                                values=candidateGenes, mart= ensembl))

dim(mappedGenesToOrthologs)
length(candidateGenes)
dim(hgnc_symbol_to_ensembl)

candidateGenes.df <- data.frame(symbol=candidateGenes)
candidateGenes.df <- merge(candidateGenes.df, hgnc_symbol_to_ensembl, by.x='symbol', by.y='hgnc_symbol',all.x=T)
mappedGenesToOrthologs <- merge(candidateGenes.df, mappedGenesToOrthologs, by.x='ensembl_gene_id', by.y='ensembl_gene_id', all.x=T )
head(mappedGenesToOrthologs)

## 
Orthologs <-  getUnemptyList(mappedGenesToOrthologs$rnorvegicus_homolog_associated_gene_name)
length(Orthologs)






