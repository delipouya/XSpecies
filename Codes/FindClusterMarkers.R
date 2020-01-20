source('Codes/Functions.R')
Initialize()
# can access the clusters this way too: sCVdata_list$res.0.2@Clusters

output_filename <- "objects/3.seur_clustered.RData"
load(output_filename)
RES = 0.1
Idents(seur) <- paste0('cluster_',as.character(seur@meta.data$SCT_snn_res.0.1))

TsneDF <- getTsneDF(seur)
PcaDF <- getPcaDF(seur)
UmapDF <- getUmapDF(seur)

ggplot(TsneDF, aes(tSNE_1, tSNE_2, color=clusters))+geom_point()+theme_bw()+ggtitle(paste0('Resolution: ', RES))
ggplot(PcaDF, aes(PC_1, PC_2, color=clusters))+geom_point()+theme_bw()+ggtitle(paste0('Resolution: ', RES))
ggplot(UmapDF, aes(UMAP_1, UMAP_2, color=clusters))+geom_point()+theme_bw()+ggtitle(paste0('Resolution: ', RES))


## we're going with 1.2 resolution
table(seur@meta.data$SCT_snn_res.0.1)

Idents(seur) <- paste0('cluster_',as.character(seur@meta.data$SCT_snn_res.0.1))
table(Idents(seur))
TsneDF <- getTsneDF(seur)
head(TsneDF)
ggplot(TsneDF, aes(x=tSNE_1, y=tSNE_2, color= clusters))+
  geom_point(alpha=0.5)+scale_color_manual(values=colorPalatte)+theme_bw()


listOfIdentityCLasses <- as.character(unique(Idents(seur)))
listOfMarkers <- sapply(1:length(listOfIdentityCLasses), 
                        function(i) {print(i);FindMarkers(seur, ident.1 = listOfIdentityCLasses[i], ident.2 = NULL)},
                        simplify = F)

names(listOfMarkers) <- listOfIdentityCLasses
saveRDS(listOfMarkers, 'objects/4.listOfMarkers.rds')
lapply(listOfMarkers, head) 


