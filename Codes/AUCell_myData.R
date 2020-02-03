## In this script we'll importing the liver signatures gmt file 
## then, we'll convert them to ensemble id and write them to a new gmt file
##

source('Codes/Functions.R')
Initialize()
PATH_TO_FILES = 'Data/McParland_markers/SUPPLEMENTARY_DATA/liver/'
SPECIES_NAME = 'rat_Rnor'


## > DO NOT RUN THIS AGAIN
## --------------------------------
#### generating the gmt file for ensemble id of liver cell markers
### this desction needs to be run once 

source('Codes/Convert_Human2ModelAnimal.R')

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




## --------------------------------
##### Check if all the known markers are in the raw dataset
candidateGenes_mapped_df <- readRDS(paste0(PATH_TO_FILES,'candidateGenes_mapped_table.rds'))
candidateGenes_mapped <- lapply(candidateGenes_mapped_df, function(x) getUnemptyList(x$rnorvegicus_homolog_ensembl_gene))

ensemble_genes <- read.delim(paste0('~/XSpecies/Data/',SPECIES_NAME,'/genes.tsv'), header = F)
lapply(candidateGenes_mapped, function(x) sum(x %in% ensemble_genes$V1/length(x)))
lapply(candidateGenes_mapped, function(x) x[!(x %in% ensemble_genes$V1)])


#### Staring the analysis for signiture finding and labeling using AUCell
AUCell_dir = paste0("Results/",SPECIES_NAME,"/AUCell_myData/")


seur <- readRDS(paste0('objects/',SPECIES_NAME,'/2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'))
exprMatrix <- as.matrix(seur[['RNA']]@data)

gmtFile <- paste0(PATH_TO_FILES,"liver_cell_type_signature_gene_sets_ensemble.gmt")
geneSets <- getGmt(gmtFile)

all_markers <- as.character(unlist(candidateGenes_mapped))
all_markers[!all_markers %in% rownames(seur)]

geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix)) 
cbind(nGenes(geneSets))
geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets), " (", nGenes(geneSets) ,"g)", sep=""))



## build gene expression ranking for each cell
cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=4, plotStats=TRUE)
saveRDS(cells_rankings, paste0(AUCell_dir,"cells_rankings.rds" ))

# Calculate enrichment for the gene signatures (AUC)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) # the top 4 percent are being considered: can change with: aucMaxRank
saveRDS(cells_AUC, file=paste0(AUCell_dir, "cells_AUC.rds"))



cells_rankings <- readRDS(paste0(AUCell_dir,"cells_rankings.rds" ))
cells_AUC <- readRDS(paste0(AUCell_dir, "cells_AUC.rds"))


## converting it to a data frame structure
cells_AUC_df = data.frame(getAUC(cells_AUC))
## Normalizing to consider probabilistic interpretation 
cells_AUC_df_norm <- mapply(`/`, cells_AUC_df,  colSums(cells_AUC_df))  ## normalizing each column to have a probabilistic point of view
rownames(cells_AUC_df_norm) <- rownames(cells_AUC_df)
rough_cell_type_assignment <- rownames(cells_AUC_df_norm)[apply(cells_AUC_df_norm,2,which.max)]
cbind(table(rough_cell_type_assignment)) ## 93% are the hepatocytes!!


cell_tsne <- data.frame(tSNE_1 = getEmb(seur, 'tsne')[,1], tSNE_2 = getEmb(seur, 'tsne')[,2], cell_type=rough_cell_type_assignment)
ggplot(cell_tsne, aes(x=tSNE_1, y=tSNE_2, color=cell_type))+geom_point(alpha=0.5)+theme_bw()

Cell_type_assigned <- sapply(1:length(cells_assignment), function(i) cells_assignment[[i]][['assignment']], simplify = F)
names(Cell_type_assigned) <- names(cells_assignment)


par(mfrow=c(4,3)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
saveRDS(cells_assignment, file=paste0(AUCell_dir, "cells_assignment.rds"))



gsva_result <- gsva(exprMatrix, candidateGenes_mapped)
saveRDS(gsva_result, '~/Desktop/gsva_result.rds')
gsva_result <- readRDS('~/Desktop/gsva_result.rds')
getHead(gsva_result)

pdf(paste0(AUCell_dir,'/marked_AUCell.pdf'))
for (index in 1:length(Cell_type_assigned)){
  
  print(index)
  marked_cells <- Cell_type_assigned[[index]]
  marked_cells_name <- names(Cell_type_assigned)[index]
  
  ####### GSVA results
  cell_tsne_gsva <- data.frame(tSNE_1 = getEmb(seur, 'tsne')[,1],
                              tSNE_2 = getEmb(seur, 'tsne')[,2], 
                              gsva=as.numeric(gsva_result[index,]) ) 
  
  p = ggplot(cell_tsne_gsva, aes(x=tSNE_1, y=tSNE_2, color=gsva))+geom_point(alpha=0.8)+
    theme_bw()+ggtitle(paste0(marked_cells_name, ' gsva results'))+scale_color_viridis(direction=-1)
  print(p)
  
  ####### AUCell results, total expression
  cell_tsne_exp <- data.frame(tSNE_1 = getEmb(seur, 'tsne')[,1],
                              tSNE_2 = getEmb(seur, 'tsne')[,2], 
                              expression=as.numeric(cells_AUC_df[index,]) ) 
  
  p1 = ggplot(cell_tsne_exp, aes(x=tSNE_1, y=tSNE_2, color=expression))+geom_point(alpha=0.8)+
    theme_bw()+ggtitle(marked_cells_name)+scale_color_viridis(direction=-1)
  print(p1)
  
  ####### AUCell results, threshold added
  isMarked <- colnames(seur) %in% marked_cells
  cell_tsne <- data.frame(tSNE_1 = getEmb(seur, 'tsne')[,1], 
                          tSNE_2 = getEmb(seur, 'tsne')[,2], 
                          label=ifelse(isMarked,'Marked','-') ) 
  
  p2 = ggplot(cell_tsne, aes(x=tSNE_1, y=tSNE_2, color=label))+geom_point(alpha=0.8)+theme_bw()+
    ggtitle(marked_cells_name)+scale_color_manual(values = c('cadetblue1','dodgerblue4'))
  print(p2)
}
dev.off()


### Checking how consistent gsva and AUCell results are
pdf('~/Desktop/scoring_correlations.pdf')
for(i in 1:length(Cell_type_assigned)){
  print(i)
  cor_df <- data.frame(gsva=as.numeric(gsva_result[i,]), AUCell=as.numeric(cells_AUC_df[i,]))
  p=ggplot(cor_df, aes(x=gsva, y=AUCell))+geom_point()+theme_bw()+ggtitle(names(Cell_type_assigned)[i])
  print(p)
}
dev.off()



### AUC threshold need to be changes -> need higher AUC threshold for Hepatocytes
getThresholdSelected(cells_assignment)

names(cells_assignment)
cells_assignment[['HEPATOCYTE (6g)']][['aucThr']][['thresholds']]



cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
colnames(assignmentTable)[2] <- "geneSet"
head(assignmentTable)
assignmentMat <- table(assignmentTable[,"geneSet"], assignmentTable[,"cell"])
assignmentMat[,1:2]

#### read through the pipeline and actually think about it!!!
#### check the functions one bu one
#### 


### Check the AUC plots >>> do they make sense 
set.seed(123)
par(mfrow=c(4,3)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 


### plot on tSNE based on the thresholds for each of the cell types 


warningMsg <- sapply(cells_assignment, function(x) x$aucThr$comment)
warningMsg[which(warningMsg!="")]

## deliving into the output object data structures
names(cells_assignment)
names(cells_assignment[['CHOLANGIOCYTES (25g)']])
names(cells_assignment[['CHOLANGIOCYTES (25g)']][['aucThr']])
cells_assignment[['CHOLANGIOCYTES (25g)']][['aucThr']][['selected']]









