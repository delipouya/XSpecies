source('Codes/Functions.R')
Initialize()

dir = 'Results/preproc_rats'
seur <- readRDS(paste0(dir,'/merged/merged_rat_samples.rds'))

## PCA
seur <- RunPCA(seur,verbose=T)
PC_NUMBER = 18

### Harmony
seur <- RunHarmony(seur, "sample_type",assay.use="SCT")

max_seurat_resolution <- 2.4 ## change this to higher values
FDRthresh <- 0.01 # FDR threshold for statistical tests
min_num_DE <- 5
seurat_resolution <- 0 # Starting resolution is this plus the jump value below.
seurat_resolution_jump <- 0.1

seur <- FindNeighbors(seur,reduction="harmony",dims=1:PC_NUMBER,verbose=F)

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
    DRforClust="harmony", #pca
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

saveRDS(sCVdata_list, '~/sCVdata_list_merged_rat.rds') 

# shrinks the size of the Seurat object by removing the scaled matrix
#seur <- DietSeurat(seur,dimreducs=Reductions(seur))

