## Run this script as: 
# Rscript Codes/get_labels_SCINA.R 'rat_Rnor' '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

source('Codes/Functions.R')
Initialize()
# INPUT_NAME = 'rat_Rnor'
# INPUT_FILE = '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'
INPUT_NAME = args[1] 
INPUT_FILE = args[2]
PATH_TO_FILES = 'Data/McParland_markers/SUPPLEMENTARY_DATA/liver/'
OUTPUT_NAME = gsub('.rds','',gsub('2.seur_dimRed_','',INPUT_FILE ))
AUCell_dir = paste0("Results/",INPUT_NAME,"/AUCell/")


gsva_result <- readRDS(paste0('Results/',INPUT_NAME,'/GSVA/GSVA_',OUTPUT_NAME,'.rds'))
SCINA_res <- readRDS(paste0('Results/',INPUT_NAME,'/SCINA/','SCINA_',OUTPUT_NAME,'.rds'))

cells_assignment <- readRDS(file=paste0(AUCell_dir, "cells_assignment_",OUTPUT_NAME,".rds"))
Cell_type_assigned <- sapply(1:length(cells_assignment), function(i) cells_assignment[[i]][['assignment']], simplify = F)
names(Cell_type_assigned) <- names(cells_assignment)

cells_AUC <- readRDS(paste0(AUCell_dir,"cells_rankings_",OUTPUT_NAME,".rds" ))
cells_AUC_df = data.frame(getAUC(cells_AUC))




### importing expression matrix and list of markers
candidateGenes_mapped_df <- readRDS(paste0(PATH_TO_FILES,'candidateGenes_mapped_table.rds'))
candidateGenes_mapped <- lapply(candidateGenes_mapped_df, 
                                function(x) getUnemptyList(x$rnorvegicus_homolog_ensembl_gene))

seur <- readRDS(paste0('objects/',INPUT_NAME,'/',INPUT_FILE))
exprMatrix <- as.matrix(seur[['RNA']]@data)


### Checking how consistent gsva and AUCell results are

pdf(paste0("Results/",INPUT_NAME,'/labeling_',OUTPUT_NAME,'.pdf'), height = 14, width = 18)

for (index in 1:length(Cell_type_assigned)){ 
  
  print(index)
  marked_cells <- Cell_type_assigned[[index]]
  marked_cells_name <- names(Cell_type_assigned)[index]
  
  labeling_df<- data.frame(gsva=as.numeric(gsva_result[index,]), 
                           AUCell=as.numeric(cells_AUC_df[index,]),
                           AUCell_label=ifelse(colnames(seur) %in% marked_cells,'Marked','-'),
                           SCINA= as.numeric(SCINA_res$probabilities[index,]),
                           tSNE_1 = getEmb(seur, 'tsne')[,1],
                           tSNE_2 = getEmb(seur, 'tsne')[,2])
  
  p_hist_1=ggplot(labeling_df, aes(x=gsva))+geom_histogram(bins=40)+theme_bw()+ggtitle(paste0('GSVA - ', names(Cell_type_assigned)[index]))
  
  p_hist_2=ggplot(labeling_df, aes(x=AUCell))+geom_histogram(bins=40)+theme_bw()+
    ggtitle(paste0('AUCell - ', names(Cell_type_assigned)[index]))+
    geom_vline(xintercept=getThresholdSelected(cells_assignment)[index], linetype="dashed", color = "red")
  
  p_hist_3=ggplot(labeling_df, aes(x=SCINA))+geom_histogram(bins=40)+theme_bw()+ggtitle(paste0('SCINA - ', names(Cell_type_assigned)[index]))
  
  
  p_cor_1=ggplot(labeling_df, aes(x=gsva, y=AUCell))+geom_point()+theme_bw()+
    ggtitle(paste0('GSVA - AUCell - ', names(Cell_type_assigned)[index]))+
    geom_hline(yintercept=getThresholdSelected(cells_assignment)[index], linetype="dashed", color = "red")
  
  p_cor_2=ggplot(labeling_df, aes(x=gsva, y=SCINA))+geom_point()+theme_bw()+
    ggtitle(paste0('GSVA - SCINA - ', names(Cell_type_assigned)[index]))
  
  p_cor_3=ggplot(labeling_df, aes(x=SCINA, y=AUCell))+geom_point()+theme_bw()+
    ggtitle(paste0('SCINA - AUCell - ', names(Cell_type_assigned)[index]))+
    geom_hline(yintercept=getThresholdSelected(cells_assignment)[index], linetype="dashed", color = "red")
  
  ####### GSVA results
  p_tsne_1 = ggplot(labeling_df, aes(x=tSNE_1, y=tSNE_2, color=gsva))+geom_point(alpha=0.8)+
    theme_bw()+ggtitle(paste0(marked_cells_name, ' gsva results'))+scale_color_viridis(direction=-1)
  
  ####### AUCell results, total expression
  p_tsne_2 = ggplot(labeling_df, aes(x=tSNE_1, y=tSNE_2, color=AUCell))+geom_point(alpha=0.8)+
    theme_bw()+ggtitle(paste0(marked_cells_name, ' AUCell results'))+scale_color_viridis(direction=-1)
  
  ####### AUCell results, threshold added
  p_tsne_3 = ggplot(labeling_df, aes(x=tSNE_1, y=tSNE_2, color=AUCell_label))+geom_point(alpha=0.8)+theme_bw()+
    ggtitle(paste0(marked_cells_name, ' AUCell results'))+scale_color_manual(values = c('cadetblue1','dodgerblue4'))
  
  #####   SCINA results 
  p_tsne_4 = ggplot(labeling_df, aes(x=tSNE_1, y=tSNE_2, color=SCINA))+geom_point(alpha=0.8)+
      theme_bw()+ggtitle(paste0(marked_cells_name, ' SCINA results'))+scale_color_viridis(direction=-1)
    
    
  gridExtra::grid.arrange(p_hist_1, p_hist_2, p_hist_3,
                          p_cor_1, p_cor_2, p_cor_3,ncol=3, nrow=2)
  
  gridExtra::grid.arrange(p_tsne_1, p_tsne_2, p_tsne_3, p_tsne_4, ncol=2, nrow=2)
}
dev.off()






