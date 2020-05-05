# tutorial
# https://satijalab.org/seurat/v3.1/immune_alignment.html

source('Codes/Functions.R')
Initialize()

### importing data
dir = 'Results/preproc_rats'
Samples <- readRDS(paste0(dir,'/merged/rat_samples_list.rds'))
lapply(Samples, head)

anchors <- FindIntegrationAnchors(object.list = Samples, dims = 1:20)
combined <- IntegrateData(anchorset = anchors, dims = 1:20)
head(combined)

### add sample name to the integrated object
combined$ds_index = str_split_fixed(colnames(combined),pattern = '_',2)[,2]
map = data.frame(ds_index=c('1', '2', '3', '4', '5'),sample=names(Samples))
ds_index_map = data.frame(ds_index=combined$ds_index, umi=colnames(combined),index=1:ncol(combined))
ds_index_map = merge(ds_index_map, map,by.x='ds_index',by.y='ds_index',all.x=T,sort=F)
combined$sample = ds_index_map$sample


DefaultAssay(combined) <- "integrated"
# Run the standard workflow for visualization and clustering
combined <- ScaleData(combined, verbose = T)
combined <- RunPCA(combined, npcs = 30, verbose = T)

# umap
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20, verbose = T)

umap_emb <- data.frame(Embeddings(combined, 'umap'))
umap_emb$sample = combined$sample
umap_emb$cell_type = combined$cell_type

ggplot(umap_emb, aes(x=UMAP_1, y=UMAP_2))+geom_point(aes(color=sample),alpha=0.7,size=2)+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5))#+ggtitle('Seurat-CCA integration')

ggplot(umap_emb, aes(x=UMAP_1, y=UMAP_2))+
  geom_point(aes(shape=sample, color=cell_type),alpha=0.7,size=2)+
  theme_classic()+scale_colour_brewer(palette = "Paired")+
  theme(plot.title = element_text(hjust = 0.5))#+ggtitle('Seurat-CCA integration')


# tsne
combined <- RunTSNE(combined, reduction = "pca", dims = 1:20, verbose = T)

tsne_emb <- data.frame(Embeddings(combined, 'tsne'))
tsne_emb$sample = combined$sample
tsne_emb$cell_type = combined$cell_type

ggplot(tsne_emb, aes(x=tSNE_1, y=tSNE_2))+geom_point(aes(color=sample),alpha=0.7,size=2)+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5))+ggtitle('Seurat-CCA integration')

ggplot(tsne_emb, aes(x=tSNE_1, y=tSNE_2))+
  geom_point(aes(shape=sample, color=cell_type),alpha=0.7,size=2)+
  theme_classic()+scale_colour_brewer(palette = "Paired")+
  theme(plot.title = element_text(hjust = 0.5))+ggtitle('Seurat-CCA integration')


saveRDS(combined, paste0(dir,'/merged/CCA_combined.rds'))

### clustering
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)




### FindIntegrationAnchors steps:
# Computing 2000 integration features
# Scaling features for provided objects
# Finding all pairwise anchors
# Merging objects
# Finding neighborhoods
# Finding anchors
# Found 14166 anchors
# Filtering anchors
# Retained 4375 anchors
# Extracting within-dataset neighbors

