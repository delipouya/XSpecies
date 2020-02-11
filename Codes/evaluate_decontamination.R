
source('Codes/Functions.R')
Initialize()


INPUT_NAME = 'rat_DA'
INPUT_FILE = '2.seur_dimRed_rat_DA_mito_30_lib_1500.rds'
OUTPUT_NAME = gsub('.rds','',gsub('2.seur_dimRed_','',INPUT_FILE ))
PATH_TO_FILES = 'Data/McParland_markers/SUPPLEMENTARY_DATA/liver/'
PC_NUMBER = 18
MIT_CUT_OFF = 30
LIB_SIZE_CUT_OFF = 1500
  
TITLE = paste0('mito threshold: ', MIT_CUT_OFF,' , library size threshold: ', LIB_SIZE_CUT_OFF)

input_from_SoupX <- paste0('Results/rat_DA/Xsoup_refined_',OUTPUT_NAME,'/')
input_from_10x <- paste0("Data/", INPUT_NAME,'/')

genes_df_10x <- read.delim(paste0(input_from_10x,'features.tsv.gz'), header = F)
genes_df_soupx <- read.delim(paste0(input_from_SoupX,'genes.tsv'), header = F)
genes_df_soupx_2 = merge(genes_df_soupx, genes_df_10x, by.x='V2', by.y='V2', all.x=T, all.y=F)
genes_df_soupx_3 <- genes_df_soupx_2[!duplicated(genes_df_soupx_2$V2),]
genes_df_soupx_3 <- genes_df_soupx_3[,c('V2','V1.y')]
write.table(genes_df_soupx_3, file=paste0(input_from_SoupX,'genes.tsv'), quote=FALSE, sep='\t', col.names = FALSE)


seur <- CreateSeuratObject(counts=Read10X(input_from_SoupX, gene.column = 1),
                           min.cells=0,min.features=1, 
                           project = "snRNAseq")
### converting the fist colum row names to ensemble IDs
seur[['RNA']] <- AddMetaData(seur[['RNA']], seur_genes_df$V2, col.name = 'symbol')

mito_genes_index <- grep(pattern = '^Mt-', seur[['RNA']]@meta.features$symbol )
seur[["mito_perc"]] <- PercentageFeatureSet(seur, features = mito_genes_index)

## adding ensemble id as a meta data to the object

getHead(GetAssayData(seur@assays$RNA))
libSize <- colSums(GetAssayData(seur@assays$RNA))

## Normalization
### SCTransform
seur <- NormalizeData(seur, normalization.method = "LogNormalize", scale.factor = 10000)
dim(seur[['RNA']]@data)

## Finding variable genes
seur <- FindVariableFeatures(seur, selection.method = "vst", nfeatures = 2000)
head(seur[['RNA']]@var.features)

## scaling data
seur <- ScaleData(seur, features = rownames(seur))

##  PCA
seur <- RunPCA(seur,verbose=F)
plot(100 * seur@reductions$pca@stdev^2 / seur@reductions$pca@misc$total.variance,
     pch=20,xlab="Principal Component",ylab="% variance explained",log="y")

## tSNE
seur <- RunTSNE(seur,dims=1:PC_NUMBER,reduction="pca",perplexity=30)

TITLE_tsne = paste0('tSNE (',TITLE,')')
df_tsne <- data.frame(tSNE_1=getEmb(seur, 'tsne')[,1], tSNE_2=getEmb(seur, 'tsne')[,2], 
                      library_size= seur$nCount_RNA, mito_perc=seur$mito_perc , n_expressed=seur$nFeature_RNA) 

ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=library_size))+geom_point()+theme_bw()+scale_color_viridis(direction= -1)+ggtitle(TITLE_tsne)+labs(caption = INPUT_NAME)
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=mito_perc))+geom_point()+theme_bw()+scale_color_viridis(direction= -1)+ggtitle(TITLE_tsne)+labs(caption = INPUT_NAME)
ggplot(df_tsne, aes(x=tSNE_1, y=tSNE_2, color=n_expressed))+geom_point()+theme_bw()+scale_color_viridis(direction= -1)+ggtitle(TITLE_tsne)+labs(caption = INPUT_NAME)


##UMAP
seur <- RunUMAP(seur,dims=1:PC_NUMBER, reduction="pca")

TITLE_umap = paste0('UMAP (',TITLE,')')
df_umap <- data.frame(UMAP_1=getEmb(seur, 'umap')[,1], UMAP_2=getEmb(seur, 'umap')[,2], 
                      library_size= seur$nCount_RNA, mito_perc=seur$mito_perc , n_expressed=seur$nFeature_RNA )

ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=library_size))+geom_point()+theme_bw()+scale_color_viridis(direction = -1)+ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=mito_perc))+geom_point()+theme_bw()+scale_color_viridis(direction = -1)+ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)
ggplot(df_umap, aes(x=UMAP_1, y=UMAP_2, color=n_expressed))+geom_point()+theme_bw()+scale_color_viridis(direction = -1)+ggtitle(TITLE_umap)+labs(caption = INPUT_NAME)



#### checking marker genes:
candidateGenes_mapped_df <- readRDS(paste0(PATH_TO_FILES,'candidateGenes_mapped_table.rds'))
candidateGenes_mapped <- lapply(candidateGenes_mapped_df, 
                                function(x) getUnemptyList(x$rnorvegicus_homolog_associated_gene_name))

exprMatrix <- as.matrix(seur[['RNA']]@data)
SCINA_res = SCINA(exprMatrix, candidateGenes_mapped, max_iter = 100, convergence_n = 10, 
                  convergence_rate = 0.99, sensitivity_cutoff = 1, 
                  rm_overlap=FALSE, allow_unknown=TRUE, log_file='SCINA.log')


pdf('~/Desktop/Check_soupX_result.pdf')
Cell_type_assigned <- names(candidateGenes_mapped)
for (index in 1:length(Cell_type_assigned)){ 
  
  print(index)
  marked_cells_name <- Cell_type_assigned[index]
  
  labeling_df<- data.frame(#gsva=as.numeric(gsva_result[index,]), 
    SCINA= as.numeric(SCINA_res$probabilities[index,]),
    tSNE_1 = getEmb(seur, 'tsne')[,1],
    tSNE_2 = getEmb(seur, 'tsne')[,2])
  
  #####   SCINA results 
  p = ggplot(labeling_df, aes(x=tSNE_1, y=tSNE_2, color=SCINA))+geom_point(alpha=0.8)+
    theme_bw()+ggtitle(paste0(marked_cells_name, ' SCINA results'))+scale_color_viridis(direction=-1)
  print(p)
}
dev.off()










