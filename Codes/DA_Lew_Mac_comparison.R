## import data 
mapper = read.delim('~/XSpecies/Data/rat_DA/features.tsv.gz',header=F)
#merged_samples <- readRDS('Results/preproc_rats/merged/merged_rat_samples_all_features.rds')
merged_samples <- readRDS('Results/preproc_rats/merged/merged_rat_samples_2.rds')

merged_samples$strain = ifelse(merged_samples$sample_type =='rat_DA_M_10WK_003' | merged_samples$sample_type=='rat_DA_01', 
                               'rat_DA', 'rat_Lew')
merged_samples$merged_clusters <- as.character(Idents(merged_samples))
merged_samples$merged_clusters <- as.character(merged_samples$RNA_snn_res.0.05)

DA_cluster_4_UMIs <- colnames(merged_samples[,merged_samples$merged_clusters == '4' & merged_samples$strain_type=='RAT_DA'])
LEW_cluster_4_UMIs <- colnames(merged_samples[,merged_samples$merged_clusters == '4' & merged_samples$strain_type=='RAT_LEW'])
merged_samples$Mac_ident <- ifelse(colnames(merged_samples) %in% DA_cluster_4_UMIs, 'DA_Mac', 
                                   ifelse(colnames(merged_samples) %in% LEW_cluster_4_UMIs, 'LEW_Mac', 'None'))

markers_df <- FindMarkers(object = merged_samples, assay='SCT', 
                          features = rownames(GetAssayData(merged_samples, assay = 'SCT')),
                          group.by = 'Mac_ident', 
                          ident.1='DA_Mac', ident.2='LEW_Mac', 
                          logfc.threshold=0,
                          min.pct=0)

markers_df$ensemble_ids = rownames(markers_df)
markers_df_mapped <- merge(markers_df, mapper, 
      by.x='ensemble_ids', by.y='V1', 
      all.x=T, all.y=F, sort=F)

head(markers_df_mapped, 30)
dim(markers_df)


markers_df_mapped <- readRDS('Results/DA_Lew_Mac_compareDF.rds')
ranked_list <- data.frame(gene_symbol= markers_df_mapped$V2, 
                          gene_ensemble = markers_df_mapped$ensemble_ids,  
                          score=-sign(markers_df_mapped$avg_logFC) * log10(markers_df_mapped$p_val_adj))

ranked_list <- ranked_list[order(ranked_list$score, decreasing = T), ]
head(ranked_list,20)
tail(ranked_list)

write.table(ranked_list, 'Results/DA_Lew_Mac_compareDF.rnk', quote = F, sep = '\t', row.names = F, col.names = F)



i = 1
gene_name = markers_df_mapped$V2[i]
gene_expression <- GetAssayData(merged_samples)[mapper$V1[mapper$V2==gene_name],]


#### Visualizing the UMAP plot ####
umap_df <- data.frame(Embeddings(merged_samples, reduction = 'umap_h'),
                      cell_type=merged_samples$cell_type,
                      cluster=merged_samples$merged_clusters,
                      sample_type=merged_samples$sample_type,
                      strain=merged_samples$strain,
                      mito_perc=merged_samples$mito_perc,
                      libSize=merged_samples$nCount_RNA,
                      numDetectedGenes= merged_samples$nFeature_RNA,
                      a_cluster=merged_samples$a_cluster)

#Marco_EXP = Marco_expression)
ggplot(umap_df, aes(x=umap_h_1, y=umap_h_2, color=a_cluster))+geom_point()+theme_classic()

ggplot(umap_df, aes(x=umap_h_1, y=umap_h_2, color=gene_expression))+geom_point(alpha=0.5)+
  theme_classic()+ylab(paste0('PC_',i))+scale_color_viridis(direction = -1)+ggtitle(gene_name)

ggplot(umap_df, aes(x=umap_h_1, y=umap_h_2, color=strain))+geom_point()+
  theme_classic()+scale_color_manual(values = colorPalatte)

ggplot(umap_df, aes(x=umap_h_1, y=umap_h_2, color=mito_perc))+geom_point(alpha=0.5)+
  theme_classic()+ylab(paste0('PC_',i))+scale_color_viridis(direction = -1)

cluster_to_check <- 'cluster_7'
umap_df$a_cluster <- umap_df$cluster == cluster_to_check
ggplot(umap_df, aes(x=umap_h_1, y=umap_h_2, color=a_cluster))+geom_point()+
  theme_classic()+ggtitle(cluster_to_check)+xlab('UMAP_1')+ylab('UMAP_2')

umap_df$subclusters = all_UMIs_2$cluster
umap_df$subclusters = ifelse(merged_samples$merged_clusters %in% c('4', '8'), 
                             merged_samples$merged_clusters, 'other')
ggplot(umap_df, aes(x=umap_h_1, y=umap_h_2, color=subclusters))+geom_point()+
  theme_classic()+xlab('UMAP_1')+ylab('UMAP_2')




#### running GSEA ####

gsea_jar <- '~/GSEA_Linux_4.1.0/gsea-cli.sh'
dest_gmt_file <- '~/Rat_GOBP_AllPathways_no_GO_iea_May_01_2020_symbol.gmt'

#### making the directory #### 
working_dir = '~/XSpecies/Results/DA_Lew_Mac_compare/'
dir.create(working_dir)

rnk_file = '~/XSpecies/Results/DA_Lew_Mac_compare/DA_Lew_Mac_compareDF.rnk'
print(rnk_file)
analysis_name = 'DA_Lew_Mac_compare'
working_dir = working_dir

# max: 200 min: 15
GSEA_command = paste("",gsea_jar,  "GSEAPreRanked -gmx", dest_gmt_file, "-rnk" ,
                     rnk_file , "-collapse false -nperm 1000 -scoring_scheme weighted -rpt_label ",
                     analysis_name,"  -plot_top_x 20 -rnd_seed 12345  -set_max 200 -set_min 15 -zip_report false -out" , 
                     working_dir, " > gsea_output.txt") ## removing the filter on gene set size
system(GSEA_command)



####