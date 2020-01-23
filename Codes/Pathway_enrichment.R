source('Codes/Functions.R')
Initialize()
# can access the clusters this way too: sCVdata_list$res.0.2@Clusters

listOfMarkers <- readRDS('objects/4.listOfMarkers.rds')
lapply(listOfMarkers, head) 

dim(listOfMarkers[['cluster_0']])


listOfMarkers_filt <- lapply(listOfMarkers, function(dataframe) 
  dataframe[(dataframe$avg_logFC>0.5 | dataframe$avg_logFC<(-0.5)) & dataframe$p_val_adj<0.01, ])
lapply(listOfMarkers_filt, dim)

## Pathway enrichment methods
geneSymbolsToMap <- lapply(listOfMarkers_filt, function(aTable)data.frame(gene_symbol=rownames(aTable)))
lapply(geneSymbolsToMap, head)

#### Query the genes using biomart
## Set up Biomart and choose the dataset
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset('rnorvegicus_gene_ensembl',mart=ensembl)

## checking the appropiate filter and attribute to use
listFilters(ensembl)[grep(listFilters(ensembl)[,1], pattern = 'entrez'),]
listAttributes(ensembl)[grep(listAttributes(ensembl)[,1], pattern = 'symbol'),]

geneSymbolsMapped <- lapply(geneSymbolsToMap, function(aGeneSymbolTable){
  ## mapping each table gene names to mart ensemble 
  MappedEnsemblId <- getBM(filters="rgd_symbol", 
                           attributes= c("ensembl_gene_id", 'rgd_symbol','entrezgene_id'),
                           values=aGeneSymbolTable$gene_symbol, mart= ensembl)
  print('done!')
  ## checking which genes have not been mapped
  geneSymbolsToMap <- merge(aGeneSymbolTable, MappedEnsemblId, by.x='gene_symbol', by.y='rgd_symbol', all.x=T)
  return(geneSymbolsToMap)})

lapply(geneSymbolsMapped, head)
saveRDS(geneSymbolsMapped, 'objects/5.DEgeneSymbolsMappedEntrez.rds')



## removing the rows(genes symbols) which were not mapped
GenesToEnrich <- lapply(geneSymbolsMapped, function(aGeneTable)as.character(aGeneTable$entrezgene_id[!is.na(aGeneTable$entrezgene_id)]))

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
pdf('~/Desktop/pathwayEnrichment.pdf', height = 20, width = 20)
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


