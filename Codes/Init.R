source('Codes/Functions.R')
Initialize()

## ToDO ## 
# have the subsequent code ready in order to do the rest of the analysis
# import the resolution list and start the analysis

## -------------------------------------------- import input Data 
input_from_10x <- "Data/rat_Rnor/"
seur <- CreateSeuratObject(counts=Read10X(input_from_10x),
                           min.cells=1,min.features=1, 
                           project = "snRNAseq")
summary(seur$nCount_RNA)
row.names(seur)[grep( pattern = 'Mt-', row.names(seur))]


seur[["percent.mt"]] <- PercentageFeatureSet(seur, pattern = '^Mt-')
show(seur)
getHead(GetAssayData(seur@assays$RNA))

libSize <- colSums(GetAssayData(seur@assays$RNA))


pdf('init_results.pdf')
par(mfrow=c(1,1))
hist(libSize, main = 'Histogram of library size')


## -------------------------------------------- QC and filtering the cells
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
## visualizing metrics
p1=ggplot(data.frame(seur$nCount_RNA), aes(seur.nCount_RNA))+
  geom_histogram(bins = 60,color='black',fill='pink',alpha=0.5)+
  theme_bw()+ggtitle('library size for all cells')+xlab('Library sizes')+
  ylab('Number of cells')
p1
summary(seur$nCount_RNA)

## distribuiton of number of genes with 
## expression value higher than zero for each cell
p2=ggplot(data.frame(seur$nFeature_RNA), aes(seur.nFeature_RNA))+
  geom_histogram(bins = 60,color='black',fill='blue',alpha=0.3)+
  theme_bw()+ggtitle('# expressed genes for all cells')+xlab('Number of expressed genes')+
  ylab('Number of cells')
p2
## proportion of reads mapped to the mitochondrian genes for each cell
## -> high proportion means low quality cells
p3=ggplot(data.frame(seur$percent.mt), aes(seur.percent.mt))+
  geom_histogram(bins = 60,color='black',fill='green',alpha=0.3)+
  theme_bw()+ggtitle('proportion of reads mapped to Mt genes')+xlab('Mitochondrial proportion (%)')+
  ylab('Number of cells')
p3
gridExtra::grid.arrange(p1,p2,p3,nrow=1,ncol=3)


### inspecting mitochondrial expression
mito_gene_identifier <- "^Mt-" 
mads_thresh <- 12
hard_thresh <- 80

mito_thresh <- median(seur$percent.mt) + mad(seur$percent.mt) * mads_thresh
drop_mito <- seur$percent.mt > mito_thresh | seur$percent.mt > hard_thresh
summary(seur$percent.mt[drop_mito])
sum(drop_mito)
mean(seur$percent.mt)
median(seur$percent.mt)

par(mar=c(3,3,2,1),mgp=2:0)
hist(seur$percent.mt,breaks=50,xlab="% mitochondrial mRNA")
abline(v=mito_thresh,col="red",lwd=2)
mtext(paste(paste0(round(mean(seur$percent.mt,2)),"% mean mitochondrial mRNA"),
            paste0(min(mito_thresh, hard_thresh), '% threshold'),
            paste0(sum(drop_mito)," cells removed"),
            sep="\n"),
      side=3,line=-3,at=mito_thresh,adj=-0.05)


temp_col <- colorspace::sequential_hcl(100,palette="Viridis",alpha=0.5,rev=T)
par(mfrow=c(1,2),mar=c(3,3,2,1),mgp=2:0)

plot(seur$nCount_RNA,seur$nFeature_RNA,log="xy",pch=20,
     xlab="nCount_RNA",ylab="nFeature_RNA",
     col=temp_col[cut(c(0,1,seur$percent.mt),100,labels=F)[c(-1,-2)]])

legend("topleft",bty="n",title="Mito %",
       legend=c(0,50,100),pch=20,col=temp_col[c(1,50,100)])

plot(seur$nCount_RNA,seur$nFeature_RNA,log="xy",pch=20,
     xlab="nCount_RNA",ylab="total_features",
     col=temp_col[cut(c(0,1,seur$percent.mt),100,labels=F)[c(-1,-2)]])
points(seur$nCount_RNA[drop_mito],seur$nFeature_RNA[drop_mito],
       pch=4,col="red")
legend("topleft",bty="n",pch=4,col="red",
       title=paste("Mito % >",round(mito_thresh,2)),
       legend=paste(sum(drop_mito),"cells"))


## filtering based on mitochondrial expression
seur <- seur[,!drop_mito]
show(seur)


# It is important to manually inspect the relationship between library size and 
# gene detection rates per cell to identify obvious outliers. In this case, we’ve identified a 
# population of cells with a different relationship between library size and complexity, as well as one
# cell with a clearly outlying library size.

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
feature.df <- data.frame(nCount_RNA=seur$nCount_RNA,percent.mt=seur$percent.mt)
p1=ggplot(feature.df, aes(x=nCount_RNA, y=percent.mt))+geom_point()+theme_bw()
feature.df <- data.frame(nCount_RNA=seur$nCount_RNA,nFeature_RNA=seur$nFeature_RNA)
p2=ggplot(feature.df, aes(x=nCount_RNA, y=nFeature_RNA))+geom_point()+theme_bw()
gridExtra::grid.arrange(p1,p2,nrow=1,ncol=2)


filt_intercept <- 60
filt_slope <- .1
to_inspect <- seur$nFeature_RNA < (seur$nCount_RNA * filt_slope + filt_intercept)
sum(to_inspect)
temp_col <- colorspace::sequential_hcl(100,palette="Viridis",alpha=0.5,rev=T)

par(mfrow=c(1,2),mar=c(3,3,2,1),mgp=2:0)
plot(seur$nCount_RNA,seur$nFeature_RNA,log="",pch=20,
     xlab="nCount_RNA",ylab="total_features",
     main="Select outliers to inspect",
     col=temp_col[cut(c(0,1,seur$percent.mt),100,labels=F)[c(-1,-2)]])
legend("topleft",bty="n",title="Mito %",legend=c(0,50,100),pch=20,col=temp_col[c(1,50,100)])
abline(filt_intercept,filt_slope,lwd=2,col="red")


plot(seur$nCount_RNA,seur$nFeature_RNA,log="xy",pch=20,
     xlab="nCount_RNA",ylab="total_features",
     main="Select outliers to inspect",
     col=temp_col[cut(c(0,1,seur$percent.mt),100,labels=F)[c(-1,-2)]])
points(seur$nCount_RNA[to_inspect],seur$nFeature_RNA[to_inspect],pch=1,col="red")
legend("topleft",bty="n",pch=1,col="red",legend="Outliers")

summary(seur$nCount_RNA)
summary(seur$nFeature_RNA)

seur <- seur[,!to_inspect]
dim(seur)


# Filtering cells based on the proportion of mitochondrial gene transcripts per cell.
## A high proportion of mitochondrial gene transcripts are indicative of poor quality cells, 
# probably due to compromised cell membranes.

## --------------------------------------  Normalization
seur <- SCTransform(seur,conserve.memory=T,verbose=F)
# "iteration limit reached" warning can be safely ignored
show(seur)
saveRDS(seur, 'objects/1.seur_normed.rds')

## -------------------------------------- PCA
seur <- FindVariableFeatures(seur, assay = 'RNA')
seur <- RunPCA(seur,assay="RNA",verbose=F) # Data has not been scaled. Please run ScaleData and retry

seur <- RunPCA(seur,assay="SCT",verbose=F)
par(mfrow=c(1,1))
plot(100 * seur@reductions$pca@stdev^2 / seur@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")

sum(seur@reductions$pca@stdev^2/seur@reductions$pca@misc$total.variance) ## check this!!?? 


## -------------------------------------- tSNE
# Select the number of principle components to use in downstream analysis, and set n_pc accordingly.
n_pc <- 20
seur <- RunTSNE(seur,dims=1:n_pc,reduction="pca",perplexity=30)

plot_tsne(cell_coord=getEmb(seur,"tsne"),
          md=getMD(seur)$nFeature_SCT,
          md_title=("nFeature_SCT"),
          md_log=F)
plot_tsne(cell_coord=getEmb(seur,"tsne"),
          md=getMD(seur)$percent.mt,
          md_title=("percent.mt: "),
          md_log=F)

## -------------------------------------- UMAP
# Playing with the perplexity parameter can improve the visualization. Perplexity can be interpretted as the number of nearby cells to consider when trying to minimize distance between neighbouring cells.

# only run if you've installed UMAP - see ?RunUMAP
seur <- RunUMAP(seur,dims=1:n_pc,reduction="pca")

plot_tsne(cell_coord=getEmb(seur,"umap"),
          md=getMD(seur)$nFeature_SCT,
          md_title="nFeature_SCT",
          md_log=F)
plot_tsne(cell_coord=getEmb(seur,"umap"),
          md=getMD(seur)$percent.mt,
          md_title="percent.mt",
          md_log=F)

saveRDS(seur, 'objects/2.seur_dimRed.rds')
seur <- readRDS('objects/2.seur_dimRed.rds')
dev.off()
## -------------------------------------- Iterative clustering with scClustViz

# Seurat implements an interpretation of SNN-Cliq (https://doi.org/10.1093/bioinformatics/btv088) 
# for clustering of single-cell expression data. They use PCs to define the distance metric, then embed the 
# cells in a graph where edges between cells (nodes) are weighted based on their similarity (euclidean distance in PCA space). 
# These edge weights are refined based on Jaccard distance (overlap in local neighbourhoods), and then communities (“quasi-cliques”) are identified in the graph
# using a smart local moving algorithm (SLM, http://dx.doi.org/10.1088/1742-5468/2008/10/P10008) to optimize the modularity measure of the defined communities in the graph.
# This code block iterates through “resolutions” of the Seurat clustering method, testing each for overfitting. Overfitting is determined by testing differential 
# expression between all pairs of clusters using a wilcoxon rank-sum test. If there are no significantly differentially expressed genes between nearest 
# neighbouring clusters, iterative clustering is stopped. The output is saved as an sCVdata object for use in scClustViz.

source('Codes/Functions.R')
Initialize()
seur <- readRDS('objects/2.seur_dimRed.rds')
n_pc = 20
max_seurat_resolution <- 1.8
## ^ change this to something large (5?) to ensure iterations stop eventually.
output_filename <- "objects/3.seur_clustered.RData"
FDRthresh <- 0.01 # FDR threshold for statistical tests
min_num_DE <- 10
seurat_resolution <- 0 # Starting resolution is this plus the jump value below.
seurat_resolution_jump <- 0.05

seur <- FindNeighbors(seur,reduction="pca",dims=1:n_pc,verbose=F)

sCVdata_list <- list()
DE_bw_clust <- TRUE
while(DE_bw_clust) {
  if (seurat_resolution >= max_seurat_resolution) { break }
  seurat_resolution <- seurat_resolution + seurat_resolution_jump 
  # ^ iteratively incrementing resolution parameter 
  
  seur <- FindClusters(seur,resolution=seurat_resolution,verbose=F)
  
  message(" ")
  message("------------------------------------------------------")
  message(paste0("--------  res.",seurat_resolution," with ",
                 length(levels(Idents(seur)))," clusters --------"))
  message("------------------------------------------------------")
  
  if (length(levels(Idents(seur))) <= 1) { 
    message("Only one cluster found, skipping analysis.")
    next 
  } 
  # ^ Only one cluster was found, need to bump up the resolution!
  
  if (length(sCVdata_list) >= 1) {
    temp_cl <- length(levels(Clusters(sCVdata_list[[length(sCVdata_list)]])))
    if (temp_cl == length(levels(Idents(seur)))) { 
      temp_cli <- length(levels(interaction(
        Clusters(sCVdata_list[[length(sCVdata_list)]]),
        Idents(seur),
        drop=T
      )))
      if (temp_cli == length(levels(Idents(seur)))) { 
        message("Clusters unchanged from previous, skipping analysis.")
        next 
      }
    }
  }
  
  curr_sCVdata <- CalcSCV(
    inD=seur,
    assayType="SCT",
    assaySlot="counts",
    cl=Idents(seur), 
    # ^ your most recent clustering results get stored in the Seurat "ident" slot
    exponent=NA, 
    # ^ going to use the corrected counts from SCTransform
    pseudocount=NA,
    DRthresh=0.1,
    DRforClust="pca",
    calcSil=T,
    calcDEvsRest=T,
    calcDEcombn=T
  )
  
  DE_bw_NN <- sapply(DEneighb(curr_sCVdata,FDRthresh),nrow)
  # ^ counts # of DE genes between neighbouring clusters at your selected FDR threshold
  message(paste("Number of DE genes between nearest neighbours:",min(DE_bw_NN)))
  
  if (min(DE_bw_NN) < min_num_DE) { DE_bw_clust <- FALSE }
  # ^ If no DE genes between nearest neighbours, don't loop again.
  
  sCVdata_list[[paste0("res.",seurat_resolution)]] <- curr_sCVdata
}


# cleaning redundant metadata
# seur@meta.data <- seur@meta.data[,colnames(seur@meta.data) != "seurat_clusters"]
# seur@meta.data <- seur@meta.data[,!grepl("^SCT_snn_res",colnames(seur@meta.data))]

# shrinks the size of the Seurat object by removing the scaled matrix
seur <- DietSeurat(seur,dimreducs=Reductions(seur))
save(sCVdata_list,seur,file=output_filename)


