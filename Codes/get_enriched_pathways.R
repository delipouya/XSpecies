# https://bioconductor.org/packages/release/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html
# hypergeometric model is implemented to assess whether the number of selected genes associated with reactome pathway is 
# larger than expected. The p values were calculated based the hypergeometric model(Boyle et al. 2004)

source('Codes/Functions.R')
Initialize()

## generate a gene list for gene-set-enrichment analysis
.getGeneList <- function(a_DE_result){
  ## feature 1: numeric vector
  geneList = a_DE_result$avg_logFC
  ## feature 2: named vector
  names(geneList) = as.character(a_DE_result$entrezgene_id)
  ## feature 3: decreasing orde
  return(sort(geneList, decreasing = TRUE))
}

#### importing the input data
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
cluster_names <- names(Cluster_markers_merged)
lapply(Cluster_markers_merged, dim)




###### Pathway enrichment analysis
### DO NOT RUN THIS SECTION
### >>>>>>>>>>>>>>>>>>>>>>>

## getting the list of genes
geneSymbolsToMap <- lapply(Cluster_markers_merged, function(aTable)data.frame(ensemble_ids=aTable$ensemble_ids))

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
  geneSymbolsToMap <- merge(aGeneSymbolTable, MappedEnsemblId, 
                            by.x='ensemble_ids', by.y='ensembl_gene_id', all.x=T)
  return(geneSymbolsToMap)})

### >>>>>>>>>>>>>>>>>>>>>>>

geneSymbolsMapped <- readRDS(paste0('Results/',INPUT_NAME , 
                                    '/pathway_analysis/DE_mapped_entrez_',OUTPUT_NAME, RES,'.rds'))

## removing the rows(genes symbols) which were not mapped
GenesToEnrich <- lapply(geneSymbolsMapped, function(aGeneTable)
  as.character(aGeneTable$entrezgene_id[!is.na(aGeneTable$entrezgene_id)]))

EnrichmentResult <-lapply(GenesToEnrich, 
                          function(aGeneVector){print('done!');
                            enrichPathway(gene=aGeneVector,organism = "rat", pvalueCutoff=0.05, readable=T)})

EnrichmentResult <- readRDS( paste0('Results/',INPUT_NAME , 
                                    '/pathway_analysis/enrich_results_',OUTPUT_NAME, RES,'.rds'))

#### converting to a readable data structure 
EnrichmentResult_filt <- lapply(EnrichmentResult, function(anEnrichmentTable){
  anEnrichmentTable <- data.frame(anEnrichmentTable)
  anEnrichmentTable <- anEnrichmentTable[anEnrichmentTable$p.adjust<0.01,]
  anEnrichmentTable$GeneRatioVal <- round(as.numeric(str_split_fixed(anEnrichmentTable$GeneRatio, '/', 2)[,1])/ 
                                            as.numeric(str_split_fixed(anEnrichmentTable$GeneRatio, '/', 2)[,2]), 3)
  return(anEnrichmentTable)
})


pdf(paste0('Results/',INPUT_NAME , '/pathway_analysis/pathway_enrich_',OUTPUT_NAME, RES,'.pdf'), width = 12, height=5)
for(i in 1:length(EnrichmentResult)){
  print(i)
  anEnrichmentRes = EnrichmentResult[[i]]
  category_number = nrow(anEnrichmentRes)
  if(category_number>15)category_number = 15
  print(dotplot(anEnrichmentRes, showCategory=category_number, title= names(EnrichmentResult)[i]))
}
dev.off()






###### Gene set enrichment analysis
# there were some ties present in the preranked state: 10.02% of the list

Cluster_markers_merged_entrez <- sapply(1:length(Cluster_markers_merged), function(i){
  merged_df = merge(Cluster_markers_merged[[i]], geneSymbolsMapped[[i]], 
                    by.x='ensemble_ids', by.y='ensemble_ids', all.x=T, all.y=F)
  merged_df[!is.na(merged_df$entrezgene_id) & !duplicated(merged_df$entrezgene_id),]
}, simplify = F)


GSA_results <- lapply(Cluster_markers_merged_entrez, function(a_DE_result){
  a_gene_list <- .getGeneList(a_DE_result)
  a_gsa_result <- gsePathway(a_gene_list, nPerm=10000, organism = "rat",
                     pvalueCutoff=0.05, pAdjustMethod="BH", verbose=TRUE)
  return(a_gsa_result)
})

names(GSA_results) <- names(Cluster_markers_merged)
saveRDS(GSA_results, paste0('Results/',INPUT_NAME , '/pathway_analysis/gsea_results_',OUTPUT_NAME, RES,'.rds'))


# rat_DA: no enrichment in cluster_1, 2, 5, 6, 10
pdf(paste0('Results/',INPUT_NAME , '/pathway_analysis/gsea_',OUTPUT_NAME, RES,'.pdf'), width = 10, height=10)
for(i in 1:length(GSA_results)){
  a_GSA_res = as.data.frame(GSA_results[[i]])
  if( nrow(a_GSA_res) == 0 ) { print(paste0('no enrichment in ', names(GSA_results)[i])); next} 
  print(emapplot(GSA_results[[i]], color="pvalue", title=names(GSA_results)[i]))
}
dev.off()






## Pathway enrichment using pathfindR:
## for cluster 1 the fisrt few genes are mitochondrial, is this normal? 
test_input <- data.frame(Gene.symbol=rownames(test), logFC=test$avg_logFC, adj.P.Val=test$p_val_adj)
res <- run_pathfindR(test_input)
