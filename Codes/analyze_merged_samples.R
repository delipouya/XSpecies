source('Codes/Functions.R')
Initialize()

PC_NUMBER = 18

dir = 'Results/preproc_rats'
seur_genes_df <- read.delim('Data/rat_DA_M_10WK_003/features.tsv.gz', header = F)

  
merged_samples <- readRDS(paste0(dir,'/merged/merged_rat_samples.rds'))
merged_samples$strain = ifelse(merged_samples$sample_type =='rat_DA_M_10WK_003' | merged_samples$sample_type=='rat_DA_01', 
                               'rat_DA', 'rat_Lew')
merged_samples$merged_clusters <- as.character(Idents(merged_samples))
#saveRDS(merged_samples, paste0(dir,'/merged/merged_rat_samples.rds'))


#### in each of the clusters, finding the DE genes between DA and lewis samples
strain_markers = sapply(0:(length(levels(merged_samples))-1), function(i){
  
  print(i)
  merged_samples_sub = merged_samples[,merged_samples$merged_clusters == as.character(i)]
  
  Idents(merged_samples_sub) = as.factor(merged_samples_sub$strain)
  markers_DA_vs_Lew <- FindMarkers(merged_samples_sub, ident.1 = 'rat_DA', ident.2 = 'rat_Lew')
  markers_Lew_vs_DA <- FindMarkers(merged_samples_sub, ident.1 = 'rat_Lew', ident.2 = 'rat_DA')
  a_cluster_strain_markers <- list(DA_vs_Lew=markers_DA_vs_Lew, Lew_vs_DA=markers_Lew_vs_DA)
  
  a_cluster_strain_markers = sapply(a_cluster_strain_markers, 
                                    function(markers_df) {
                                      markers_df$ensemble_ids = rownames(markers_df)
                                      merge(markers_df, seur_genes_df, 
                                            by.x='ensemble_ids', 
                                            by.y='V1', all.x=T, all.y=F, sort=F)}
                                    ,simplify = F )
  
}, simplify = F)

names(strain_markers) = as.character(0:(length(levels(merged_samples))-1))
saveRDS(strain_markers, 'Results/preproc_rats/strain_markers.rds')



### checking the sample (DA vs lewis) composition of each clusters
table(merged_samples$cell_type, merged_samples$strain)
total = colSums(table(merged_samples$cell_type, merged_samples$strain))
cell_type_strain = table(merged_samples$cell_type, merged_samples$strain)/total
cell_type_strain_melt = data.frame(reshape2::melt(cell_type_strain))
colnames(cell_type_strain_melt)[1:2] = c('cell_type','strain')

ggplot(data=cell_type_strain_melt, aes(x=cell_type, y=value, fill=strain)) +
  geom_bar(stat="identity",color='black')+theme_classic()+scale_fill_brewer(palette = "Paired")

total = rowSums(table(merged_samples$cell_type, merged_samples$strain))
cell_type_strain = table(merged_samples$cell_type, merged_samples$strain)
for( i in 1:length(total)){
  cell_type_strain[i,] = cell_type_strain[i,]/total[i]
}

cell_type_strain_melt = data.frame(reshape2::melt(cell_type_strain))
colnames(cell_type_strain_melt)[1:2] = c('cell_type','strain')

ggplot(data=cell_type_strain_melt, aes(x=cell_type, y=value, fill=strain)) +
  geom_bar(stat="identity",color='black')+theme_classic()+scale_fill_brewer(palette = "Paired")




### finding markers for each of the clusters:
clusters_markers <- sapply(0:(length(levels(merged_samples))-1), function(i){
  markers_df = FindMarkers(merged_samples, ident.1 = as.character(i), ident.2 = NULL)
  markers_df$ensemble_ids = rownames(markers_df)
  merge(markers_df, seur_genes_df, 
        by.x='ensemble_ids', 
        by.y='V1', all.x=T, all.y=F, sort=F)
}, simplify = F)

names(clusters_markers) = paste0('cluster_', 0:(length(levels(merged_samples))-1) )
lapply(clusters_markers, head)
saveRDS(clusters_markers, 'Results/preproc_rats/merged_samples_clusters_markers.rds')
clusters_markers <- readRDS('Results/preproc_rats/merged_samples_clusters_markers.rds')





#####################
pc <- data.frame(Embeddings(merged_samples, 'pca'))
pc_df = data.frame(PC_1=pc$PC_1, PC_2=pc$PC_2, strain=merged_samples$strain, 
                   sample_type=merged_samples$sample_type, cell_type=merged_samples$cell_type)
ggplot(pc_df, aes(x=PC_1, y=PC_2, color=cell_type))+geom_point()+theme_classic()


harmony_c <- data.frame(Embeddings(merged_samples, 'harmony'))
head(harmony_c)
harmony_df = data.frame(harmony_1=harmony_c$harmony_1, harmony_2=harmony_c$harmony_2, 
                        strain=merged_samples$strain, sample_type=merged_samples$sample_type, 
                        cell_type=merged_samples$cell_type)
ggplot(harmony_df, aes(x=harmony_1, y=harmony_2, color=cell_type))+geom_point()+theme_classic()

diff = data.frame(diff_1=pc_df$PC_1 - harmony_df$harmony_1, 
                  diff_2=pc_df$PC_2 - harmony_df$harmony_2,
                  strain=merged_samples$strain, sample_type=merged_samples$sample_type, 
                  cell_type=merged_samples$cell_type)

ggplot(diff, aes(x=diff_1, y=diff_2, color=cell_type))+geom_point()+theme_classic()




##### checking the loadings of the total merged samples
loading_df = data.frame(readRDS('~/XSpecies/Results/pc_res_rotation_merged_samples.rds'))

threshold = sqrt(1/ncol(data)) 
ggplot(loading_df, aes(PC1, PC2))+geom_point()+
  geom_vline(xintercept=threshold, linetype="dashed", color = "red")+
  geom_vline(xintercept=-threshold, linetype="dashed", color = "red")+theme_classic()

loading_df$PC1 = abs(loading_df$PC1)
loadings_PC1 <- data.frame(gene=rownames(loading_df)[loading_df$PC1 > threshold], 
                           loading = loading_df$PC1[loading_df$PC1 > threshold])

loadings_PC1 <- loadings_PC1[order(loadings_PC1$loading, decreasing = T),]

model_animal_name = "rnorvegicus"

get_gprofiler_enrich <- function(markers, model_animal_name){
  gostres <- gost(query = markers,
                  ordered_query = TRUE, exclude_iea =TRUE, 
                  sources=c('GO:BP' ,'REAC'),evcodes = T,
                  organism = model_animal_name, correction_method='fdr')
  return(gostres)
}
test = get_gprofiler_enrich(loadings_PC1$gene, model_animal_name)
head(test$result)




### split the gene expression matrix based on the cell-type
cell_type_names = names(table(merged_samples$cell_type))
merged_samples_split <- sapply(1:length(cell_type_names), 
       function(i){merged_samples[,merged_samples$cell_type == cell_type_names[i]]}, simplify = F)
names(merged_samples_split) = cell_type_names


### run PCA on each sample
split_samples_pcs = lapply(merged_samples_split, RunPCA)
split_samples_pcs = lapply(split_samples_pcs, function(x) Embeddings(x, 'pca'))
names(split_samples_pcs) = cell_type_names
lapply(split_samples_pcs, head)

pdf('~/XSpecies/Results/merged_samples_splitted_pca.pdf', width = 6, height = 5)
for( i in 1:length(cell_type_names)){
  pc = data.frame(split_samples_pcs[[i]])
  pc_df = data.frame(PC_1=pc$PC_1, PC_2=pc$PC_2, 
                     strain=merged_samples_split[[i]]$strain, 
                     sample_type=merged_samples_split[[i]]$sample_type, 
                     cell_type=merged_samples_split[[i]]$cell_type)
  
  p1=ggplot(pc_df, aes(x=PC_1, y=PC_2, color=sample_type))+geom_point()+theme_classic()+ggtitle(cell_type_names[i])
  p2=ggplot(pc_df, aes(x=PC_1, y=PC_2, color=strain))+geom_point()+theme_classic()+ggtitle(cell_type_names[i])
  print(p1)
  print(p2)
}
dev.off()


### run harmony on each cell-type
split_samples_harmony = lapply(merged_samples_split, function(x) RunHarmony(x, "sample_type",assay.use="SCT"))
split_samples_harmony = lapply(split_samples_harmony, function(x) Embeddings(x, 'harmony'))
names(split_samples_harmony) = cell_type_names

pdf('~/XSpecies/Results/merged_samples_splitted_harmony.pdf', width = 6, height = 5)
for( i in 1:length(cell_type_names)){
  harmony = data.frame(split_samples_harmony[[i]])
  harmony_df = data.frame(harmony_1=harmony$harmony_1, harmony_2=harmony$harmony_2, 
                     strain=merged_samples_split[[i]]$strain, 
                     sample_type=merged_samples_split[[i]]$sample_type, 
                     cell_type=merged_samples_split[[i]]$cell_type)
  
  p1=ggplot(harmony_df, aes(x=harmony_1, y=harmony_2, color=sample_type))+geom_point()+theme_classic()+ggtitle(cell_type_names[i])
  p2=ggplot(harmony_df, aes(x=harmony_1, y=harmony_2, color=strain))+geom_point()+theme_classic()+ggtitle(cell_type_names[i])
  print(p1)
  print(p2)
}
dev.off()



#### checking markers expression distribution
## Klrb1 gene as an NK marker
x = GetAssayData(merged_samples)
sum(rownames(x) == 'ENSRNOG00000057410')
## ApoE gene 
sum(rownames(x) == 'ENSRNOG00000018454')
ApoE_expression = x[rownames(x) == 'ENSRNOG00000018454',]

##### Running UMAP 
umap_emb <- data.frame(Embeddings(merged_samples, 'umap'))
umap_emb$sample_type = merged_samples$sample_type
umap_emb$cell_type = merged_samples$cell_type
umap_emb$ApoE_expression = ApoE_expression

ggplot(umap_emb, aes(x=UMAP_1, y=UMAP_2))+
  geom_point(aes(color=ApoE_expression),alpha=0.7,size=2)+
  theme_classic()+scale_color_viridis(direction = -1,option = "plasma")


