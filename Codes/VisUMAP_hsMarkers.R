source('Codes/Functions.R')
Initialize()
## import data 
mapper = read.delim('~/XSpecies/Data/rat_DA/features.tsv.gz',header=F)
#merged_samples <- readRDS('Results/preproc_rats/merged/merged_rat_samples_all_features.rds')
merged_samples <- readRDS('Results/preproc_rats/merged/merged_rat_samples_2.rds')

merged_samples$strain = ifelse(merged_samples$sample_type =='rat_DA_M_10WK_003' | merged_samples$sample_type=='rat_DA_01', 
                               'rat_DA', 'rat_Lew')
merged_samples$merged_clusters <- as.character(Idents(merged_samples))
merged_samples$merged_clusters <- as.character(merged_samples$RNA_snn_res.0.05)

hsMarkers = list(Heps = c('ALB', 'HAMP', 'ARG1', 'PCK1', 'AFP', 'BCHE'), 
     LSECs = c('CALCRL', 'CD32B', 'VWF'), 
     cholangiocytes = c('KRT19', 'EPCAM', 'FXYD2', 'CLDN4', 'CLDN10', 'SOX9', 'MMP7', 'CXCL1', 'CFTR', 
                        'TFF2', 'KRT7', 'CD24'),
     hepatic_stellate_cell = c('ACTA2', 'COL1A1', 'TAGLN', 'COL1A2', 'COL3A1', 'SPARC', 'RBP1', 'DCN', 'MYL9'),
     nonInfMac =c('CD68', 'MARCO'),
     abTcell = c('PTPRC', 'CD45', 'CD3D', 'TRAC'),
     gdTcell = c('IL7R', 'CD127', 'KLRB1', 'CD161', 'NKG7', 'FCGR3A', 'CD16'), 
     NKcell = c('GZMA', 'GZMB', 'GZMK', 'PRF1'),
     Plasmacell = c('CD79A', 'CD79B', 'CD27', 'IGHG1'),
     matureBcell= c('MS4A1', 'CD20', 'LTB', 'CD52', 'IGHD')
)

normed_gene_names_ens <- rownames(merged_samples@assays$SCT)
normed_gene_names_symbol <- mapper$V2[mapper$V1 %in% normed_gene_names_ens]

pdf('plots/hsLiverPaperUmaps.pdf')
for(i in 1:length(hsMarkers)){
  a_cell_type_markers = hsMarkers[[i]]
  for(j in 1:length(a_cell_type_markers)){
    hs_gene_name = a_cell_type_markers[j]
    gene_name <- .getMapped_hs2model_df(ensembl=ensembl, 
                                        candidateGenes=hs_gene_name, 
                                        model_animal_name= 'rnorvegicus' )$rnorvegicus_homolog_associated_gene_name[1]
    
    print(gene_name)
    if(is.na(gene_name)) next
    if(!gene_name %in% normed_gene_names_symbol) {print('not in matrix');next}
    feature <- mapper$V1[mapper$V2==gene_name]
    feature <- feature[feature %in% rownames(GetAssayData(merged_samples, assay='SCT'))]
    gene_expression <- GetAssayData(merged_samples, assay='SCT')[feature,]
    #gene_expression <- GetAssayData(merged_samples, assay='RNA')[mapper$V1[mapper$V2==gene_name],]
    
    #### Visualizing the UMAP plot ####
    umap_df <- data.frame(Embeddings(merged_samples, reduction = 'umap_h'),
                          gene_expression=gene_expression)
    
    p=ggplot(umap_df, aes(x=umap_h_1, y=umap_h_2, color=gene_expression))+geom_point(alpha=0.5)+
      theme_classic()+ylab(paste0('PC_',i))+scale_color_viridis(direction = -1)+
      ggtitle(paste0(gene_name, ' - ', names(hsMarkers)[i]))
    print(p)
    print('done')
  }
}
dev.off()
