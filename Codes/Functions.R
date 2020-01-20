## This script includes functions needed in the analysis 
## and need to be loaded in the beginning of scripts

Initialize <- function(){
  options(stringsAsFactors = FALSE)
  library("Seurat")
  library("BiocManager")
  library("devtools")
  library("Matrix")
  library(scran)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(scClustViz)
  library(presto)
  library(ggplot2)
  library('pathfindR')
  library("biomaRt")
  library("ReactomePA")
  library(stringr)
  library('org.Rn.eg.db')
}


colorPalatte <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "brown"
)

getHead <- function(dataframe){
  print(dataframe[1:5, 1:5])
}

Plot.tsne.gene.expr <- function(tsne.gene.df, GENE_NAME){
  MAX_VAL = 7.35
  MIN_VAL = 0  
  ggplot(tsne.gene.df, aes(x=tSNE_1, y=tSNE_2, color=gene.expr))+
    geom_point(alpha=0.8)+theme_bw()+ggtitle(GENE_NAME) + 
    theme(plot.title = element_text(hjust = 0.5))+ xlab('tSNE1')+ylab('tSNE2')+
    scale_colour_gradientn(name='Expression',
                           colours = c("gray","orange","red","red2","red4"),
                           values = c(0, 1, 2, 4, 7.35), 
                           breaks=c(0, 1, 2, 4, 7.35), 
                           limits=c(0, 7.35),
                           labels=c('Min', '', '','' ,'Max')) 
  
}

Plot.tsne.gene.expr_2 <- function(tsne.gene.df, GENE_NAME){
  ggplot(tsne.gene.df, aes(x=tSNE_1, y=tSNE_2, color=gene.expr))+
    geom_point(alpha=0.6)+theme_bw()+ggtitle(GENE_NAME) + 
    theme(plot.title = element_text(hjust = 0.5))+ xlab('tSNE1')+ylab('tSNE2')+
    scale_colour_gradientn(name='Expression',colours = c("gray",'orange',"tomato",'orangered',"red", "red2","red4"),
                           values = c(min(gene.expr), max(gene.expr)/6,max(gene.expr)/5,
                                      max(gene.expr)/4,max(gene.expr)/3,max(gene.expr)/2,
                                      max(gene.expr))) 
  
  
}

getTsneDF <- function(seur){
  listOfClusters <- as.character(Idents(seur))
  DataFrame <- as.data.frame(cbind(getEmb(seur, 'tsne'),clusters=listOfClusters))
  DataFrame$tSNE_1 <- as.numeric(DataFrame$tSNE_1)
  DataFrame$tSNE_2 <- as.numeric(DataFrame$tSNE_2)
  return(DataFrame)
}

getUmapDF <- function(seur){
  listOfClusters <- as.character(Idents(seur))
  DataFrame <- as.data.frame(cbind(getEmb(seur, 'umap'),clusters=listOfClusters))
  DataFrame$UMAP_1 <- as.numeric(DataFrame$UMAP_1)
  DataFrame$UMAP_2 <- as.numeric(DataFrame$UMAP_2)
  return(DataFrame)
}

getPcaDF <- function(seur){
  listOfClusters <- as.character(Idents(seur))
  DataFrame <- as.data.frame(cbind(getEmb(seur, 'pca')[,1:2],clusters=listOfClusters))
  DataFrame$PC_1 <- as.numeric(DataFrame$PC_1)
  DataFrame$PC_2 <- as.numeric(DataFrame$PC_2)
  return(DataFrame) 
}


#breaks=c(0,0.5,1),labels=c("Minimum",0.5,"Maximum"),
# limits=c(0,1)


#plot_tsne(cell_coord=getEmb(seur,"tsne"),
#          md=gene.expr,
#          md_title=GENE_NAME,
#          md_log=F)
