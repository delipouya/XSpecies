source('Codes/Functions.R')
Initialize()

path_to_files <- list.files('objects/rat_Rnor/mit_thresh/', full.names = T, include.dirs = T)
Seur_mt <- lapply(path_to_files, readRDS)
names(Seur_mt) <- c('MAD12','thresh_40', 'thresh50_lowLibSize' )


index = 3
seur = Seur_mt[[index]]
seur_index_name = names(Seur_mt)[index]

.get_Included_mapped_Ortholog_markers <- function(vector_of_orthologs, seur){
  vector_of_orthologs <- getUnemptyList(vector_of_orthologs)
  return(vector_of_orthologs[vector_of_orthologs %in% rownames(seur) ])
}

MappedGenesToOrthologs <- readRDS('Data/McParland_markers/mapped_markers_hs2ratRnor/All_mapped_markers_hs2ratRnor_list.rds')
Marker_orthologs_cleaned <-   lapply(MappedGenesToOrthologs, 
                                     function(x) 
                                       .get_Included_mapped_Ortholog_markers(x$rnorvegicus_homolog_associated_gene_name, seur))

pdf('~/Desktop/mit_tsne_mitDis.pdf')
for( i in 1:length(Seur_mt)){
  seur = Seur_mt[[i]]
  tsne_df <- data.frame(tSNE_1 = getEmb(seur, 'tsne')[,1],
                        tSNE_2 = getEmb(seur, 'tsne')[,2],
                        gene.expr = seur$percent.mt)
  print(Plot.tsne.gene.expr(tsne_df,  names(Seur_mt)[i]))
}
dev.off()

pdf('~/Desktop/mit_tsne_LibSizeDis.pdf')
for( i in 1:length(Seur_mt)){
  seur = Seur_mt[[i]]
  tsne_df <- data.frame(tSNE_1 = getEmb(seur, 'tsne')[,1],
                        tSNE_2 = getEmb(seur, 'tsne')[,2],
                        gene.expr = seur$nCount_RNA)
  print(Plot.tsne.gene.expr(tsne_df,  names(Seur_mt)[i]))
}
dev.off()


pdf('~/Desktop/mit_tsne_expGeneDis.pdf')
for( i in 1:length(Seur_mt)){
  seur = Seur_mt[[i]]
  tsne_df <- data.frame(tSNE_1 = getEmb(seur, 'tsne')[,1],
                        tSNE_2 = getEmb(seur, 'tsne')[,2],
                        gene.expr = seur$nFeature_RNA)
  print(Plot.tsne.gene.expr(tsne_df,  names(Seur_mt)[i]))
}
dev.off()


## importing the seur object 
# load('objects/rat_Rnor/3.seur_clustered.RData')

for(i in 1:length(Marker_orthologs_cleaned)){
  cluster_id = names(Marker_orthologs_cleaned)[i]
  a_cluster_marker_orthologs = Marker_orthologs_cleaned[[i]]
  print(cluster_id)
  pdf(paste0('~/Desktop/mit/',seur_index_name,'/', cluster_id, '.pdf'))
  for(j in 1:length(a_cluster_marker_orthologs)){
    GENE_NAME = a_cluster_marker_orthologs[j]
    gene.expr <- GetAssayData(object= seur, assay='SCT')[rownames(seur)==GENE_NAME,]
    tsne_df <- data.frame(tSNE_1 = getEmb(seur, 'tsne')[,1],
                          tSNE_2 = getEmb(seur, 'tsne')[,2],
                          gene.expr = gene.expr)
    
    print(Plot.tsne.gene.expr(tsne_df, GENE_NAME))
  }
  dev.off()
}


