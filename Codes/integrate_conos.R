## integration using Conos integration method

source('Codes/Functions.R')
Initialize()
library(conos)
library(dplyr)
library(pagoda2)

### importing data
dir = 'Results/preproc_rats'
Samples <- readRDS(paste0(dir,'/merged/rat_samples_list.rds'))
Samples_2 = sapply(1:length(Samples), function(i){ RenameCells(Samples[[i]], 
                                                             add.cell.id = names(Samples)[i], 
                                                             for.merge = FALSE)}, simplify = F)
names(Samples_2) = names(Samples)
Samples_data <- lapply(Samples_2, GetAssayData)
any(duplicated(unlist(lapply(Samples_data,colnames))))

Samples_preproc <- lapply(Samples_data, basicP2proc, n.cores=detectCores()-2, 
                             min.cells.per.gene=0, n.odgenes=2e3, 
                             get.largevis=FALSE, make.geneknn=FALSE)

str(Samples_preproc,1)
con <- Conos$new(Samples_preproc, n.cores=detectCores()-2)
str(con$samples,1)
con$plotPanel(clustering="multilevel", use.local.clusters=T, title.size=6)

## if you want to change the parameters: con$pairs$PCA <- NULL
con$buildGraph(k=30, k.self=5, space='PCA', ncomps=30, 
               n.odgenes=2000, matching.method='mNN', metric='angular', 
               score.component.variance=TRUE, verbose=TRUE)
plotComponentVariance(con, space='PCA')

con$findCommunities(method=leiden.community, resolution=1)
con$plotPanel(font.size=4)
plotClusterBarplots(con, legend.height = 0.1)



