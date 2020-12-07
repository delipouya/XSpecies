source('Codes/Functions.R')
Initialize()

AUCell_dir <- '~/XSpecies/Results/preproc_rats/merged/AUCell/'
gene_set_name = 'PPARG'

#######################################
#### Calculating the enrichment of genesets using AUCell  ####


PPARG_geneSets_converted2rat <- readRDS('~/PPARG_geneSets_converted2rat.rds')
### checking PPARG instead of all gene-sets
geneSets_converted2rat <- PPARG_geneSets_converted2rat
a_geneSet_names <- names(geneSets_converted2rat)

rat_geneSets <- sapply(1:length(geneSets_converted2rat), 
                       function(i) GeneSet(geneSets_converted2rat[[i]], 
                                           setName=names(geneSets_converted2rat)[i]), simplify =F)

names(rat_geneSets) = names(geneSets_converted2rat)


#### loading the imputated count matrix
merged_samples_imputed <- readRDS('Results/preproc_rats/merged/merged_samples_imputed.rds')
exprMatrix <- as.matrix(merged_samples_imputed$estimate)


cells_rankings <- AUCell_buildRankings(exprMatrix, 
                                       nCores=detectCores()-2, 
                                       plotStats=TRUE)
cells_rankings <- readRDS('~/XSpecies/cells_rankings_imputed.rds')
dim(cells_rankings)



###
cells_AUC <- sapply(1:length(rat_geneSets), function(i){ # 645 gene sets
  print('>>>>>>>>>>>>>')
  print(i)
  a_gene_set <- rat_geneSets[[i]]
  print(sum(geneIds(a_gene_set) %in% rownames(exprMatrix))/length(geneIds(a_gene_set)))
  a_gene_set_name <- gsub(' ', '_' ,names(rat_geneSets)[i])
  
  # the top 4 percent are being considered: can change with: aucMaxRank
  cells_AUC <- AUCell_calcAUC(a_gene_set, cells_rankings) 
  
  return(cells_AUC)
}, simplify = F )

cells_AUC <- readRDS(file=paste0(AUCell_dir, gene_set_name ,"_cells_AUC_imputed.rds"))
names(cells_AUC) <-  a_geneSet_names

cells_AUC_df <- lapply(cells_AUC, function(x) data.frame(t(data.frame(getAUC(x)))))
cells_AUC_df <- do.call(cbind, cells_AUC_df)
cells_AUC_df$UMI = rownames(cells_AUC_df)

cells_AUC_df <- readRDS(paste0(AUCell_dir, gene_set_name ,"_cells_AUC_df_imputed.rds"))




