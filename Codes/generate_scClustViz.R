source('Codes/Functions.R')
Initialize()

### importing data and changing the features to gene-symbols
mapper = read.delim('~/XSpecies/Data/rat_DA/features.tsv.gz',header=F)
merged_samples <- readRDS('Results/preproc_rats/merged/merged_rat_samples_2.rds')
assay_data <- GetAssayData(merged_samples, 'data')
rownames(assay_data) <- mapper$V2


## creating the object
your_scRNAseq_data_object <- CreateSeuratObject(
  counts=assay_data,
  project = "MergedRats",
  assay = "RNA",
  names.field = 1,
  names.delim = "_",
  meta.data = NULL)


## adding the scaled data matrix 
assay_scaled.data <- GetAssayData(merged_samples, 'scale.data')
assay_scaled.data[2000:2010,2000:2010] ## check if the data is normalized
rownames(assay_scaled.data) <- mapper$V2

your_scRNAseq_data_object <- SetAssayData(
  object = your_scRNAseq_data_object,
  assay.type = 'SCT',
  new.data = assay_scaled.data,
  slot = 'scale.data'
)


### adding dimension reduction embeddings
your_scRNAseq_data_object[["pca"]] <- CreateDimReducObject(embeddings = Embeddings(merged_samples, 'pca'), 
                                                               key = "PCA_", 
                                                               assay = DefaultAssay(your_scRNAseq_data_object))

your_scRNAseq_data_object[["harmony"]] <- CreateDimReducObject(embeddings = Embeddings(merged_samples, 'harmony'), 
                                                               key = "Harmony_", 
                                                               assay = DefaultAssay(your_scRNAseq_data_object))

your_scRNAseq_data_object[["umap"]] <- CreateDimReducObject(embeddings = Embeddings(merged_samples, 'umap_h'), 
                                                               key = "UMAP_", 
                                                               assay = DefaultAssay(your_scRNAseq_data_object))

your_scRNAseq_data_object[["tsne"]] <- CreateDimReducObject(embeddings = Embeddings(merged_samples, 'tsne_h'), 
                                                               key = "tSNE_", 
                                                               assay = DefaultAssay(your_scRNAseq_data_object))

## creating the meta.data dataframe
your_cluster_results <- data.frame(RNA_snn_res.0.05=merged_samples@meta.data$RNA_snn_res.0.05,
                                   RNA_snn_res.0.1=merged_samples@meta.data$RNA_snn_res.0.1)
rownames(your_cluster_results)  = rownames(merged_samples@meta.data)

### calculating the differentially expressed marker genes
sCVdata_list <- CalcAllSCV(
  inD=your_scRNAseq_data_object,
  clusterDF=your_cluster_results,
  assayType='RNA', #specify assay slot of data
  DRforClust="harmony",#reduced dimensions for silhouette calc
  exponent=exp(1), #log base of normalized data
  pseudocount=1,
  DRthresh=0.1, #gene filter - minimum detection rate
  testAll=F, #stop testing clusterings when no DE between clusters
  FDRthresh=0.05,
  calcSil=T, #use cluster::silhouette to calc silhouette widths
  calcDEvsRest=T,
  calcDEcombn=T
)

save(your_scRNAseq_data_object,sCVdata_list,
     file="for_scClustViz.RData")
# This file can now be shared so anyone 
# can view your results with the Shiny app!
load("for_scClustViz.RData")


runShiny(
  filePath="for_scClustViz.RData",
  
  outPath="./",
  # Save any further analysis performed in the app to the
  # working directory rather than library directory.
  
  annotationDB="org.Rn.eg.db",
  # This is an optional argument, but will add annotations.
  
  imageFileType="png"
  #Set the file format of any saved figures from the app.
)




