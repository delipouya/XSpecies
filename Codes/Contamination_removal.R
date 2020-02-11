source('Codes/Functions.R')
source('Codes/convert_human_to_ortholog_functions.R')
Initialize()


INPUT_NAME = 'rat_DA'
INPUT_FILE = '2.seur_dimRed_rat_DA_mito_30_lib_1500.rds'
OUTPUT_NAME = gsub('.rds','',gsub('2.seur_dimRed_','',INPUT_FILE ))
PATH_TO_FILES = 'Data/McParland_markers/SUPPLEMENTARY_DATA/liver/'
PC_NUMBER = 18


## Imporing the seur object
seur <- readRDS(paste0('objects/',INPUT_NAME,'/',INPUT_FILE))
seur <- FindNeighbors(seur, dims = 1:PC_NUMBER)
seur <- FindClusters(seur, resolution = 0.5)
seur$seurat_clusters <- paste0('cluster_', seur$seurat_clusters)
table(seur$seurat_clusters)

## Importing dimension reducted data
tsne_df <- as.data.frame(getEmb(seur, 'tsne'))
tsne_df$clusters <- seur$seurat_clusters
ggplot(tsne_df, aes(x=tSNE_1, y=tSNE_2, color=clusters))+geom_point()+theme_bw()



### Known liver markes:
candidateGenes_mapped_df <- readRDS(paste0(PATH_TO_FILES,'candidateGenes_mapped_table.rds'))
candidateGenes_mapped <- lapply(candidateGenes_mapped_df, function(x) getUnemptyList(x$rnorvegicus_homolog_ensembl_gene))
candidateGenes_mapped_symbols <- lapply(candidateGenes_mapped_df, function(x) getUnemptyList(x$rnorvegicus_homolog_associated_gene_name))
names(candidateGenes_mapped_symbols)

### Using SoupX test as a labeling method
pdf('~/Desktop/soupX_plotMarkerMap.pdf')
for(i in  1:length(candidateGenes_mapped_symbols)){
  geneSet2Vis <- candidateGenes_mapped_symbols[[i]]
  geneSet2Vis <- geneSet2Vis[geneSet2Vis %in% rownames(sc$toc) ]
  p= plotMarkerMap(sc,
                   geneSet = geneSet2Vis,
                   tsne_df)+ggtitle(names(candidateGenes_mapped_symbols)[i])+theme_bw()
  print(p)
}
dev.off()


# possible genes to estimate contamination fraction:
# HB genes and erythrocytes, IG genes and B-cells,TPSB2/TPSAB1 and Mast cells, etc. 
### How to choose the gene list???

#####  IG genes in rat 
igGenes = c("IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHD", "IGHE", 
            "IGHM", "IGLC1", "IGLC2", "IGLC3", "IGLC4", "IGLC5", "IGLC6", "IGLC7", "IGKC")
igGenes_mapped_df <- lapply(igGenes, function(candidateGenes) 
                                     .getMapped_hs2model_df(wanted_attributes, ensembl, candidateGenes))
ig_genes <- unlist(lapply(igGenes_mapped_df, function(x) getUnemptyList(x$rnorvegicus_homolog_associated_gene_name)))
## alternative way:
rownames(sc$toc)[grep(c('Igh'), rownames(sc$toc))]


##### hemoglobin genes in rat 
hemoglobin_genes = c(rownames(sc$toc)[grep(c('Hbb'), rownames(sc$toc))], 
                      rownames(sc$toc)[grep('Hba', rownames(sc$toc))])
#####  B cell genes in rat 
b_cell_genes <- candidateGenes_mapped_symbols$B_CELLS

#####  TPSB2/TPSAB1 genes in rat 
tryptase_genes <- rownames(sc$toc)[grep(c('Tpsb'), rownames(sc$toc))]



### Imporing the unfilted input data for SoupX analysis
dataDirs = c("~/Desktop/SoupX_rat_DA/")
sc = load10X(dataDirs, keepDroplets = TRUE, cellIDs = rownames(tsne_df))
sc = estimateSoup(sc)

## saveing the dimension reduction in the meta data
sc = setDR(sc, tsne_df)
sc = setClusters(sc, tsne_df$clusters)
head(sc$metaData)


# visualise the ratio of observed counts for a gene (or set of genes) to this expectation value
plotMarkerMap(sc,geneSet = genesToEstimateContamination,tsne_df)+theme_bw()
## Picking soup specific genes
head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
## I DONT GET THIS PLOT 
marker_dis_plot <- plotMarkerDistribution(sc)




## Finding the cells which actually express these genes and should not be used for contamination estimation
## adding tryptase_genes leads to error in the glm 
NonExpressedGeneList=list(IG=ig_genes,
                          Hemoglobin=hemoglobin_genes,
                          B_Cells=b_cell_genes) 


NonExpressedGeneList <- lapply(NonExpressedGeneList, function(geneset) geneset[geneset%in%rownames(sc$toc)])
# Expressing non-expressing cells
useToEst = estimateNonExpressingCells(sc,
                                      nonExpressedGeneList = NonExpressedGeneList,
                                      clusters=setNames(tsne_df$clusters, rownames(tsne_df)))


Non_expressing_plots <- lapply(names(NonExpressedGeneList), function(set_name){
  plotMarkerMap(sc, DR =tsne_df,
                geneSet = NonExpressedGeneList[[set_name]],
                useToEst = useToEst[,set_name])+theme_bw()+ggtitle(set_name)
})
gridExtra::grid.arrange(grobs=Non_expressing_plots,nrow=2,ncol=2)


lapply(NonExpressedGeneList, class)
## Calculating the comtamination fraction
sc = calculateContaminationFraction(sc, 
                                    nonExpressedGeneList = NonExpressedGeneList, 
                                    useToEst = useToEst)


## Correcting expression profile
out = adjustCounts(sc,clusters = setNames(tsne_df$clusters, rownames(tsne_df)), nCores=detectCores()-2)




## Investigating changes in expression
cntSoggy = rowSums(sc$toc > 0)
cntStrained = rowSums(out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
mostZeroed

tail(sort(rowSums(sc$toc > out)/rowSums(sc$toc > 0)), n = 20)
plotChangeMap(sc, out, "Hba")








