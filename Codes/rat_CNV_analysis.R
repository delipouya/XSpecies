
######################
source('Codes/Functions.R')
Initialize()
library(plyr)


#####
mapper <- read.delim('Data/rat_DA/features.tsv.gz', header = F)
rat_genome_pos <- read.delim('Data/mart_export.txt')
table(rat_genome_pos$Chromosome.scaffold.name[rat_genome_pos$Karyotype.band==''])

rat_genome_pos$chrom_brand = paste0(rat_genome_pos$Chromosome.scaffold.name, '_',rat_genome_pos$Karyotype.band)
tail(rat_genome_pos)
df_sub <- rat_genome_pos[,colnames(rat_genome_pos) %in% c('Gene.stable.ID','chrom_brand')]
karotype_lists <- split.data.frame(df_sub,df_sub$chrom_brand )
names(karotype_lists)[1] <- 'unconvensional'
karotype_lists$unconvensional$Karyotype.band <- 'None'
lapply(karotype_lists, head)

karotype_gmt <- lapply(karotype_lists, function(x) x$Gene.stable.ID)
lapply(karotype_gmt, head)


f3 = '1'
if (f3=="2") {rseq="Poisson"} else {rseq="Gaussian"}

### importing data
dir = 'Results/preproc_rats'
Samples <- readRDS('Results/preproc_rats/merged/merged_rat_samples_all_features.rds')

######################
strain_names_types <- unique(Samples$strain_type)
### dividing the expression matrix based on the clusters
strain_specific_expression <- sapply(1:length(strain_names_types), function(i){
  a_strain_name = strain_names_types[i]
  Samples[,Samples$strain_type == a_strain_name]
}, simplify = F)

names(strain_specific_expression) = strain_names_types
lapply(strain_specific_expression, dim)

strain_average_exp <- lapply(strain_specific_expression, function(x){
  ## calculate the average expression of each gene in each cluster
  df = data.frame(average=rowSums(x)/ncol(x))
  return(df)
})
strain_average_exp_df = do.call(cbind,strain_average_exp)
colnames(strain_average_exp_df) = names(strain_average_exp)

strain_average_exp_df$genes <- rownames(strain_average_exp_df)
head(strain_average_exp_df)


## maybe try this using harmony embeddings as well
loadings_total <- Loadings(Samples, 'pca')

######### Varimax rotation
initial_data <- t(GetAssayData(Samples)[rownames(Samples)%in%rownames(loadings_total),])
## apply varimax rotation on the loadings
varimax_res <- varimax(loadings_total)
rotatedLoadings <- as.matrix(varimax_res$loadings)
PC_loadings <- data.frame(genes=rownames(rotatedLoadings),PC_loading=rotatedLoadings[,1:10])
colnames(PC_loadings)[2:ncol(PC_loadings)] = paste0('PC_', 1:10)
head(PC_loadings)


input_df <- merge(strain_average_exp_df, PC_loadings, by.x='genes',by.y='genes',all.x=F, all.y=F, sort=F)
rownames(input_df) <- input_df$genes
input_df <- input_df[,-1]
dim(input_df)

PC_loadings_gsva <- gsva(as.matrix(input_df), karotype_gmt, min.sz=5, max.sz=2000, mx.diff=TRUE, 
                           verbose=T, kcdf=rseq, parallel.sz=0)# ,,method='ssgsea'

colnames(PC_loadings_gsva)[1:2] <- paste0('CNV_',colnames(PC_loadings_gsva)[1:2])
p1=pheatmap::pheatmap(PC_5_loadings_gsva, main='CNV inference based on ave-exp',fontsize_row = 4, fontsize_col =10)
p2=pheatmap::pheatmap(cor(PC_loadings_gsva), main='CNV inference based on ave-exp',fontsize=4)


pdf('plots/strain_cnv_averageExp_pc5.pdf', height = 20, width=10)
print(p1)
dev.off()

# CNV_rat_DA CNV_rat_Lew PC5_loading
# CNV_rat_DA   1.00000000  -0.6737886 -0.02722777
# CNV_rat_Lew -0.67378859   1.0000000 -0.45483014
# PC5_loading -0.02722777  -0.4548301  1.00000000

#GSVA and output
ratStrain_cnv_gsva <- gsva(as.matrix(strain_average_exp_df), karotype_gmt, min.sz=0, max.sz=2000, mx.diff=TRUE, 
               verbose=T, kcdf=rseq, parallel.sz=0)# ,,method='ssgsea'
ratStrain_cnv_gsva <- as.matrix(ratStrain_cnv_gsva)
ratStrain_cnv_gsva








############################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("infercnv")



###
#####
mapper$V2 <- make.unique(mapper$V2)
datam_2 <- merge(datam, mapper, by.x='X', by.y='V2', all.x=T, sort=F)
rownames(datam_2) <- datam_2$V1
datam_2 <- as.matrix(datam_2)
datam_3 <- datam_2[,c(2:(ncol(datam_2)-2))]

datam_2 <- datam[,-1]
rownames(datam_2) <- datam[,1]
datam_3 = as.matrix(datam_2)

