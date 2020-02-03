
## https://bioconductor.org/packages/release/bioc/vignettes/AUCell/inst/doc/AUCell.html
# Build the rankings
# Calculate the Area Under the Curve (AUC)
# Set the assignment thresholds
# To support paralell execution:
BiocManager::install(c("doMC", "doRNG"))
# For the main example:
BiocManager::install(c("mixtools", "GEOquery", "SummarizedExperiment"))
# For the examples in the follow-up section of the tutorial:
BiocManager::install(c("DT", "plotly", "NMF", "d3heatmap", "shiny", "rbokeh"))

#dir.create("AUCell_tutorial")
setwd("~/XSpecies/AUCell_tutorial")

library(GEOquery)
attemptsLeft <- 20
while(attemptsLeft>0)
{
  geoFile <- tryCatch(getGEOSuppFiles("GSE60361", makeDirectory=FALSE), error=identity) 
  if(methods::is(geoFile,"error")) 
  {
    attemptsLeft <- attemptsLeft-1
    Sys.sleep(5)
  }
  else
    attemptsLeft <- 0
}

gzFile <- grep(".txt.gz", basename(rownames(geoFile)), fixed=TRUE, value=TRUE)
message(gzFile)
txtFile <- gsub(".gz", "", gzFile, fixed=TRUE)
message(txtFile)
gunzip(filename=gzFile, destname=txtFile, remove=TRUE)

library(data.table)
geoData <- fread(txtFile, sep="\t")
geneNames <- unname(unlist(geoData[,1, with=FALSE]))
exprMatrix <- as.matrix(geoData[,-1, with=FALSE])
rm(geoData)
dim(exprMatrix)
rownames(exprMatrix) <- geneNames
exprMatrix[1:5,1:4]

# Remove file
file.remove(txtFile)

# Save for future use
mouseBrainExprMatrix <- exprMatrix
save(mouseBrainExprMatrix, file="exprMatrix_AUCellVignette_MouseBrain.RData")
load("exprMatrix_AUCellVignette_MouseBrain.RData")


set.seed(333)
exprMatrix <- mouseBrainExprMatrix[sample(rownames(mouseBrainExprMatrix), 5000),]
getHead(exprMatrix)
library(AUCell)
library(GSEABase)
gmtFile <- paste(file.path(system.file('examples', package='AUCell')), "geneSignatures.gmt", sep="/")
geneSets <- getGmt(gmtFile)

geneSets <- subsetGeneSets(geneSets, rownames(exprMatrix)) 
cbind(nGenes(geneSets))
geneSets <- setGeneSetNames(geneSets, newNames=paste(names(geneSets), " (", nGenes(geneSets) ,"g)", sep=""))

# Random
set.seed(321)
extraGeneSets <- c(
  GeneSet(sample(rownames(exprMatrix), 50), setName="Random (50g)"),
  GeneSet(sample(rownames(exprMatrix), 500), setName="Random (500g)"))

countsPerGene <- apply(exprMatrix, 1, function(x) sum(x>0))
# Housekeeping-like
extraGeneSets <- c(extraGeneSets,
                   GeneSet(sample(names(countsPerGene)[which(countsPerGene>quantile(countsPerGene, probs=.95))], 100), setName="HK-like (100g)"))

geneSets <- GeneSetCollection(c(geneSets,extraGeneSets))
names(geneSets)

### adding a few random genes and also many house keeping genes
# Random
set.seed(321)
extraGeneSets <- c(
  GeneSet(sample(rownames(exprMatrix), 50), setName="Random (50g)"),
  GeneSet(sample(rownames(exprMatrix), 500), setName="Random (500g)"))

countsPerGene <- apply(exprMatrix, 1, function(x) sum(x>0))
# Housekeeping-like
extraGeneSets <- c(extraGeneSets,
                   GeneSet(sample(names(countsPerGene)[which(countsPerGene>quantile(countsPerGene, probs=.95))], 100), setName="HK-like (100g)"))

geneSets <- GeneSetCollection(c(geneSets,extraGeneSets))
names(geneSets)

cells_rankings <- AUCell_buildRankings(exprMatrix, nCores=1, plotStats=TRUE)
save(cells_rankings, file="cells_rankings.RData")

cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
save(cells_AUC, file="cells_AUC.RData")




# In order to calculate the AUC, by default only the top 5% of the genes in the 
# ranking are used (i.e. checks whether the genes in the gene-set or signature 
# are within the top 5%). This allows faster execution on bigger datasets,
# and reduce the effect of the noise at the bottom of the ranking (e.g. where many 
# genes might be tied at 0 counts). The percentage to take into account can be modified with 
# the argument aucMaxRank. For datasets where most cells express many genes (e.g. a filtered 
# dataset), or these have high expression values, it might be good to increase this threshold. 
# Check the histogram provided by AUCell_buildRankings to get an estimation on 
# where this threshold lies within the dataset.

set.seed(123)
par(mfrow=c(3,3)) 
cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE) 

# The thresholds calcuated for each gene set are stored in the $aucThr slot. 
## For example, the thresholds suggested for the oligodendrocyte gene-set:
sapply(1:length(cells_assignment), function(i)cells_assignment[[i]]$aucThr$thresholds, simplify = F)
cells_assignment[[1]]$aucThr$thresholds
cells_assignment$Oligodendrocyte_Cahoy$aucThr$thresholds
cells_assignment$Oligodendrocyte_Cahoy$aucThr$selected

warningMsg <- sapply(cells_assignment, function(x) x$aucThr$comment)
warningMsg[which(warningMsg!="")]


par(mfrow=c(1,1))
geneSetName <- rownames(cells_AUC)[grep("Oligodendrocyte_Cahoy", rownames(cells_AUC))]
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=0.25)
abline(v=0.25)

# Assigning cells to this new threshold:
newSelectedCells <- names(which(getAUC(cells_AUC)[geneSetName,]>0.08))
length(newSelectedCells)
head(newSelectedCells)


