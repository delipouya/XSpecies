
source('Codes/Functions.R')
source('Codes/convert_human_to_ortholog_functions.R')
Initialize()
# rat_Rnor, 2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds
# rat_DA, 2.seur_dimRed_rat_DA_mito_30_lib_1500.rds
# INPUT_NAME = 'rat_DA' #rat_Rnor
# INPUT_FILE = '2.seur_dimRed_rat_DA_mito_30_lib_1500.rds'
PATH_TO_FILES = 'Data/McParland_markers/SUPPLEMENTARY_DATA/liver/'
OUTPUT_NAME = gsub('.rds','',gsub('2.seur_dimRed_','',INPUT_FILE ))
model_animal_name = "rnorvegicus"    # 'mmusculus'

seur <- readRDS(paste0('objects/',INPUT_NAME,'/',INPUT_FILE))
exprMatrix <- GetAssayData(seur)

RES= ''
if(INPUT_NAME=='rat_Rnor') RES = '_res.1'
if(INPUT_NAME=='rat_Lew_01') RES = '_res.0.8'
Cluster_markers_merged <- readRDS(paste0('Results/', INPUT_NAME,'/markers/markers_', OUTPUT_NAME, RES,'.rds'))

##############################################################
#### Tallulah markers

### importing input markers
general_markers_df <- readRDS('Data/McParland_markers/liver_markers_tallulah/liver_general_markers.rds')
### removing the unmapped markers
general_markers_df$rat_general <- general_markers_df$rat_general[getUnemptyList_bool(general_markers_df$rat_general$ensembl_gene_id) &
                                                                   general_markers_df$rat_general$ensembl_gene_id %in% rownames(seur),]

general_markers_df$human_general <- general_markers_df$human_general[getUnemptyList_bool(general_markers_df$human_general$rnorvegicus_homolog_ensembl_gene) &
                                                                       general_markers_df$human_general$rnorvegicus_homolog_ensembl_gene %in% rownames(seur),]


# check Tallulah hepatocyte populations cell markers
hep_populations = c('centralHep', 'Hepatocytes', 'midHep', 'periportHep','pericentralHep','periportalHep') #last 2 for rat
df1 = general_markers_df$rat_general[general_markers_df$rat_general$Celltype %in% hep_populations, c('Gene','Celltype')]
df2 = general_markers_df$human_general[general_markers_df$human_general$Celltype %in% hep_populations, c('rnorvegicus_homolog_associated_gene_name', 'Celltype')]
colnames(df2)[1] = 'Gene'
hep_markers <- rbind(df1, df2)
hepatocytes_split <- split( hep_markers , f = hep_markers$Celltype )
hepatocytes_split_genes <- lapply(hepatocytes_split, function(x) x$Gene)

cluster_to_check = 'cluster_1'

for(i in 1:length(hepatocytes_split_genes)){
  print('----------------')
  print(names(hepatocytes_split_genes)[i])
  print(length(hepatocytes_split_genes[[i]]))
  rat_cell_markers <- hepatocytes_split_genes[[i]]
  DE_genes <- Cluster_markers_merged[[cluster_to_check]]$V2[1:50]
  print(DE_genes[DE_genes%in% rat_cell_markers])
}


# check Tallulah immune cell markers

###### DONT RUN THIS SECTION
immune_markers <- read.csv('Data/McParland_markers/liver_markers_tallulah/liver_immune_markers.csv')
immune_markers_human <- immune_markers[immune_markers$Species == 'Human',]
data.frame(table(immune_markers_human$Celltype))

mapped_markers_df = .getMapped_hs2model_df( ensembl, immune_markers_human$Gene , model_animal_name)
mapped_markers_df_2 = mapped_markers_df[mapped_markers_df$rnorvegicus_homolog_ensembl_gene != '' & 
                                          !is.na(mapped_markers_df$rnorvegicus_homolog_ensembl_gene),]

mapped_markers_df_2 <- merge(mapped_markers_df_2,immune_markers_human ,by.x='symbol', by.y='Gene',sort=F,all.x=F,all.y=F)
######
mapped_markers_df_2 <- readRDS('Data/McParland_markers/liver_markers_tallulah/liver_immune_markers_mapped.rds')
immune_split <- split( mapped_markers_df_2 , f = mapped_markers_df_2$Celltype )
immune_split <- lapply(immune_split, function(x) x[x$rnorvegicus_homolog_ensembl_gene %in% rownames(seur),] )
immune_split_genes <- lapply(immune_split, function(x) x$rnorvegicus_homolog_associated_gene_name)


cluster_to_check = 'cluster_4'

for(i in 1:length(immune_split_genes)){
  print('----------------')
  print(names(immune_split_genes)[i])
  print(length(immune_split_genes[[i]]))
  rat_cell_markers <- immune_split_genes[[i]]
  DE_genes <- Cluster_markers_merged[[cluster_to_check]]$V2[1:20]
  print(DE_genes[DE_genes%in% rat_cell_markers])
}
##############################################################









############################################################## Visualization
pdf(paste0("Results/",INPUT_NAME,'/markers/tallulah_general_markers_',INPUT_NAME,'.pdf'))
# >>>>>>>> rat markers
marker_df <- general_markers_df$rat_general
for( i in 1:nrow(marker_df)){
  
  aMarker = marker_df$ensembl_gene_id[i]
  aMarker_expression <- exprMatrix[aMarker,]
  print(paste0('rat: ', i))
  tsne_df <- as.data.frame(cbind(getEmb(seur, 'tsne'),gene.expr=aMarker_expression))
  
  p = Plot.tsne.gene.expr(tsne_df ,
                          marker_df$Gene[i] ,
                          paste0(marker_df$Species[i], ' - ',marker_df$Celltype[i], ' - ', INPUT_NAME))
  print(p)}
# >>>>>>>> human markers
marker_df <- general_markers_df$human_general
for( i in 1:nrow(marker_df)){
  
  aMarker = marker_df$rnorvegicus_homolog_ensembl_gene[i]
  aMarker_expression <- exprMatrix[aMarker,]
  print(paste0('human: ', i))
  tsne_df <- as.data.frame(cbind(getEmb(seur, 'tsne'),gene.expr=aMarker_expression))
  
  p = Plot.tsne.gene.expr(tsne_df,
                          paste0('Human: ', marker_df$Gene[i], '  rat: ', marker_df$rnorvegicus_homolog_associated_gene_name[i]),
                          paste0(marker_df$Species[i], ' - ',marker_df$Celltype[i], ' - ', INPUT_NAME))
  print(p)}

dev.off()


######### visualizing known immune markers 
pdf(paste0('Results/',INPUT_NAME,'/markers/tallulah_immune_markers_',INPUT_NAME,'.pdf'))
for(i in 1:nrow(mapped_markers_df_2)){
  aMarker = mapped_markers_df_2$rnorvegicus_homolog_ensembl_gene[i]
  
  aMarker_expression <- exprMatrix[aMarker,]
  tsne_df <- as.data.frame(cbind(getEmb(seur, 'tsne'),gene.expr=aMarker_expression))
  p=Plot.tsne.gene.expr(tsne_df,
                        paste0('Human: ', mapped_markers_df_2$symbol[i], 
                               '  rat: ', mapped_markers_df_2$rnorvegicus_homolog_associated_gene_name[i]),
                        paste0(mapped_markers_df_2$Species[i], ' - ',
                               mapped_markers_df_2$Celltype[i], ' - ', INPUT_NAME))
  print(p)
  print('--------------------')
  print(paste0(aMarker, ' - ', 
               mapped_markers_df_2$rnorvegicus_homolog_associated_gene_name[i], ' - ',
               mapped_markers_df_2$Celltype[i]))
}
dev.off()
##############################################################





