source('Codes/Functions.R')
Initialize()

INPUT_NAME = 'rat_DA' #'mouse'
model_animal_name = "rnorvegicus" #'mmusculus' 

Rdata_PATH = paste0('Results/', INPUT_NAME, '/clusters/')
INPUT_FILES = list.files(path = Rdata_PATH , pattern = '.RData', full.names = T, include.dirs = T)
input_version = 1
INPUT_FILE = INPUT_FILES[input_version]
OUTPUT_NAME = gsub('.RData','',gsub(paste0(Rdata_PATH, '/clusters_'),'',INPUT_FILE ))
OUTPUT_NAME = gsub('_v2', '', OUTPUT_NAME)

RES= ''
if(INPUT_NAME=='rat_Rnor') RES = '_res.1'
Cluster_markers_merged <- readRDS(paste0('Results/', INPUT_NAME, '/markers/markers_', 
                                              OUTPUT_NAME, RES,'.rds'))

Cluster_markers_filt <- lapply(Cluster_markers_merged, function(dataframe) 
  dataframe[(dataframe$avg_logFC>0.5 | dataframe$avg_logFC<(-0.5)) & dataframe$p_val_adj<0.01, ])

cluster_names <- names(Cluster_markers_merged)
lapply(Cluster_markers_merged, dim)
lapply(Cluster_markers_filt, dim)

## Pathway enrichment methods
geneSymbolsToMap <- lapply(Cluster_markers_filt, function(aTable)data.frame(ensemble_ids=aTable$ensemble_ids))
lapply(geneSymbolsToMap, head)

#### Query the genes using biomart
## Set up Biomart and choose the dataset
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset('rnorvegicus_gene_ensembl',mart=ensembl)

## checking the appropiate filter and attribute to use
listFilters(ensembl)[grep(listFilters(ensembl)[,1], pattern = 'entrez'),]
listAttributes(ensembl)[grep(listAttributes(ensembl)[,1], pattern = 'ensembl_gene_id'),]

geneSymbolsMapped <- lapply(geneSymbolsToMap, function(aGeneSymbolTable){
  ## mapping each table gene names to mart ensemble 
  MappedEnsemblId <- getBM(filters="ensembl_gene_id", 
                           attributes= c("ensembl_gene_id", 'rgd_symbol','entrezgene_id'),
                           values=aGeneSymbolTable$ensemble_ids, mart= ensembl)
  print('done!')
  ## checking which genes have not been mapped
  geneSymbolsToMap <- merge(aGeneSymbolTable, MappedEnsemblId, by.x='ensemble_ids', by.y='ensembl_gene_id', all.x=T)
  return(geneSymbolsToMap)})

lapply(geneSymbolsMapped, head)
# saveRDS(geneSymbolsMapped, 'objects/5.DEgeneSymbolsMappedEntrez.rds')



## removing the rows(genes symbols) which were not mapped
GenesToEnrich <- lapply(geneSymbolsMapped, function(aGeneTable)
  as.character(aGeneTable$entrezgene_id[!is.na(aGeneTable$entrezgene_id)]))

# https://bioconductor.org/packages/release/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html
# hypergeometric model is implemented to assess whether the number of selected genes associated with reactome pathway is 
# larger than expected. The p values were calculated based the hypergeometric model(Boyle et al. 2004)

EnrichmentResult <-lapply(GenesToEnrich, 
                          function(aGeneVector){print('done!');
                            data.frame(enrichPathway(gene=aGeneVector,organism = "rat", pvalueCutoff=0.05, readable=T))})
 

EnrichmentResult_filt <- lapply(EnrichmentResult, function(anEnrichmentTable){
  anEnrichmentTable <- anEnrichmentTable[anEnrichmentTable$p.adjust<0.01,]
  anEnrichmentTable$GeneRatioVal <- round(as.numeric(str_split_fixed(anEnrichmentTable$GeneRatio, '/', 2)[,1])/ 
                                            as.numeric(str_split_fixed(anEnrichmentTable$GeneRatio, '/', 2)[,2]), 3)
  return(anEnrichmentTable)
})

lapply(EnrichmentResult, dim)
lapply(EnrichmentResult_filt, dim)
head(EnrichmentResult_filt[[1]])


anEnrichmentRes = EnrichmentResult_filt[[1]]

lapply(EnrichmentResult_filt, function(anEnrichmentRes){anEnrichmentRes$Description})



pdf(paste0('Results/', INPUT_NAME, '/markers/pathway_enrichment.pdf'), height = 20, width = 20)
lapply(EnrichmentResult_filt, function(anEnrichmentRes){
  p=ggplot(anEnrichmentRes,aes(y=GeneRatioVal,x=reorder(Description,GeneRatioVal), fill=p.adjust) )+
    geom_bar(stat="identity", color="black")+theme_bw()+coord_flip()
  print(p)
})
dev.off()


## Pathway enrichment using pathfindR:
## for cluster 1 the fisrt few genes are mitochondrial, is this normal? 
test_input <- data.frame(Gene.symbol=rownames(test), logFC=test$avg_logFC, adj.P.Val=test$p_val_adj)
res <- run_pathfindR(test_input)


