## Run this script as: 
# Rscript Codes/get_labels_AUCell.R 'rat_Rnor' '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)



source('Codes/Functions.R')
Initialize()

INPUT_NAME = args[1] 
INPUT_FILE = args[2]
# INPUT_NAME = 'rat_Rnor'
# INPUT_FILE = '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'
PATH_TO_FILES = 'Data/McParland_markers/SUPPLEMENTARY_DATA/liver/'
OUTPUT_NAME = gsub('.rds','',gsub('2.seur_dimRed_','',INPUT_FILE ))


## --------------------------------
##### Check if all the known markers are in the raw dataset
candidateGenes_mapped_df <- readRDS(paste0(PATH_TO_FILES,'candidateGenes_mapped_table.rds'))
candidateGenes_mapped <- lapply(candidateGenes_mapped_df, function(x) getUnemptyList(x$rnorvegicus_homolog_ensembl_gene))

ensemble_genes <- read.delim(paste0('~/XSpecies/Data/',INPUT_NAME,'/genes.tsv'), header = F)
lapply(candidateGenes_mapped, function(x) sum(x %in% ensemble_genes$V1/length(x)))
lapply(candidateGenes_mapped, function(x) x[!(x %in% ensemble_genes$V1)])


#### Staring the analysis for signiture finding and labeling using AUCell
AUCell_dir = paste0("Results/",INPUT_NAME,"/AUCell/")


seur <- readRDS(paste0('objects/',INPUT_NAME,'/',INPUT_FILE))
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
saveRDS(cells_rankings, paste0(AUCell_dir,"cells_rankings_",OUTPUT_NAME,".rds" ))

# Calculate enrichment for the gene signatures (AUC)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings) # the top 4 percent are being considered: can change with: aucMaxRank
saveRDS(cells_AUC, file=paste0(AUCell_dir, "cells_AUC_",OUTPUT_NAME,".rds"))



cells_rankings <- readRDS(paste0(AUCell_dir,"cells_rankings_",OUTPUT_NAME,".rds" ))
cells_AUC <- readRDS(paste0(AUCell_dir,"cells_rankings_",OUTPUT_NAME,".rds" ))



## converting it to a data frame structure
cells_AUC_df = data.frame(getAUC(cells_AUC))

## Normalizing to consider probabilistic interpretation 
cells_AUC_df_norm <- mapply(`/`, cells_AUC_df,  colSums(cells_AUC_df))  ## normalizing each column to have a probabilistic point of view
rownames(cells_AUC_df_norm) <- rownames(cells_AUC_df)
rough_cell_type_assignment <- rownames(cells_AUC_df_norm)[apply(cells_AUC_df_norm,2,which.max)]
cbind(table(rough_cell_type_assignment)) ## 93% are the hepatocytes!!

## visualzing rough estimation
cell_tsne <- data.frame(tSNE_1 = getEmb(seur, 'tsne')[,1], tSNE_2 = getEmb(seur, 'tsne')[,2], cell_type=rough_cell_type_assignment)
ggplot(cell_tsne, aes(x=tSNE_1, y=tSNE_2, color=cell_type))+geom_point(alpha=0.5)+theme_bw()


par(mfrow=c(4,3)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
saveRDS(cells_assignment, file=paste0(AUCell_dir, "cells_assignment_",OUTPUT_NAME,".rds"))

Cell_type_assigned <- sapply(1:length(cells_assignment), function(i) cells_assignment[[i]][['assignment']], simplify = F)
names(Cell_type_assigned) <- names(cells_assignment)




# ### AUC threshold need to be changes -> need higher AUC threshold for Hepatocytes
# getThresholdSelected(cells_assignment)
# 
# names(cells_assignment)
# cells_assignment[['HEPATOCYTE (6g)']][['aucThr']][['thresholds']]
# 
# 
# 
# cellsAssigned <- lapply(cells_assignment, function(x) x$assignment)
# assignmentTable <- reshape2::melt(cellsAssigned, value.name="cell")
# colnames(assignmentTable)[2] <- "geneSet"
# head(assignmentTable)
# assignmentMat <- table(assignmentTable[,"geneSet"], assignmentTable[,"cell"])
# assignmentMat[,1:2]
# 
# #### read through the pipeline and actually think about it!!!
# #### check the functions one bu one
# 
# ### Check the AUC plots >>> do they make sense 
# set.seed(123)
# par(mfrow=c(4,3)) 
# cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 
# 
# 
# ### plot on tSNE based on the thresholds for each of the cell types 
# warningMsg <- sapply(cells_assignment, function(x) x$aucThr$comment)
# warningMsg[which(warningMsg!="")]
# 
# ## deliving into the output object data structures
# names(cells_assignment)
# names(cells_assignment[['CHOLANGIOCYTES (25g)']])
# names(cells_assignment[['CHOLANGIOCYTES (25g)']][['aucThr']])
# cells_assignment[['CHOLANGIOCYTES (25g)']][['aucThr']][['selected']]








