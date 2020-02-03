source('Codes/Functions.R')
Initialize()

input_from_10x <- "Data/rat_Rnor/"
seur <- CreateSeuratObject(counts=Read10X(input_from_10x),
                           min.cells=1,min.features=1, 
                           project = "snRNAseq")
seur[["percent.mt"]] <- PercentageFeatureSet(seur, pattern = '^Mt-')


## making a list of boolian filters
Drop_mito = list(0,0,0)

mads_thresh <- 12
mito_thresh <- median(seur$percent.mt) + mad(seur$percent.mt) * mads_thresh
Drop_mito[[1]] <- seur$percent.mt > mito_thresh

mito_thresh = 40
Drop_mito[[2]] <- seur$percent.mt > mito_thresh

mito_thresh = 50
Drop_mito[[3]] <- seur$percent.mt > mito_thresh & (seur$nCount_RNA < median(seur$nCount_RNA))


names(Drop_mito) <- c('MAD12','thresh_40', 'thresh50_lowLibSize' )
lapply(Drop_mito, sum)
length(Drop_mito[[1]])



## going through each filtered seurat object and 
## precessing and visualization, in order to find the best 
## filteration cutofff for mitochondria percentage

for(i in 1:length(Drop_mito)){
  iteration_name = names(Drop_mito)[i]
  print(iteration_name)
  drop_mito <- Drop_mito[[i]]
  seur.2 <- seur[,!drop_mito]
  
  
  pdf(paste0('Results/rat_Rnor/', iteration_name, '.pdf'))
  print('start normalizarion')
  seur.2 <- SCTransform(seur.2,conserve.memory=T,verbose=F)
  print('normalization is done!')
  seur.2 <- FindVariableFeatures(seur.2, assay = 'RNA')
  print('starting pca...')
  seur.2 <- RunPCA(seur.2,assay="SCT",verbose=F)
  n_pc <- 20
  print('pca is done.')
  seur.2 <- RunTSNE(seur.2,dims=1:n_pc,reduction="pca",perplexity=30)
  plot_tsne(cell_coord=getEmb(seur.2,"tsne"),
            md=getMD(seur.2)$nFeature_SCT,
            md_title=("nFeature_SCT"),
            md_log=F)
  plot_tsne(cell_coord=getEmb(seur.2,"tsne"),
            md=getMD(seur.2)$percent.mt,
            md_title=("percent.mt: "),
            md_log=F)
  print('tsne is done')
  seur.2 <- RunUMAP(seur.2,dims=1:n_pc,reduction="pca")
  saveRDS(seur.2, paste0('objects/rat_Rnor/mit_thresh/', iteration_name,'.rds'))
  plot_tsne(cell_coord=getEmb(seur.2,"umap"),
            md=getMD(seur.2)$nFeature_SCT,
            md_title="nFeature_SCT",
            md_log=F)
  plot_tsne(cell_coord=getEmb(seur.2,"umap"),
            md=getMD(seur.2)$percent.mt,
            md_title="percent.mt",
            md_log=F)
  print('umap is done.')
  dev.off()
  print('finished.')
}




