## This script includes functions needed in the analysis 
## and need to be loaded in the beginning of scripts

# function to check to see if packages are installed. 
#  Install them if they are not, then load them into the R session.
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    BiocManager::install(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


## loading required packages    
Initialize <- function(){
  options(stringsAsFactors = FALSE)
  
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  options(stringsAsFactors = FALSE)
  listOfPackages <- c('BiocManager', 'Seurat', 'viridis', 'org.Rn.eg.db', 'stringr', 'ReactomePA',  'biomaRt', 
                      'pathfindR', 'ggplot2', 'presto', 'scClustViz', 'org.Hs.eg.db', 'AnnotationDbi', 
                      'scran', 'Matrix', 'devtools', 'AUCell', 'GSEABase','GSVA', 'fgsea','limma', 'SCINA',
                      'DelayedArray', 'DelayedMatrixStats','garnett','monocle')
  ipak(listOfPackages)
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

getHead <- function(dataframe){print(dataframe[1:5, 1:5])}

getUnemptyList <- function(chrList){ chrList[!is.na(chrList) & chrList != '' ]}



Plot.tsne.gene.expr <- function(tsne.gene.df, GENE_NAME){
  ggplot(tsne.gene.df, aes(x=tSNE_1, y=tSNE_2, color=gene.expr))+
    geom_point(alpha=0.6)+theme_bw()+ggtitle(GENE_NAME) + 
    theme(plot.title = element_text(hjust = 0.5))+ xlab('tSNE1')+ylab('tSNE2')+scale_color_viridis(direction = -1)
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
