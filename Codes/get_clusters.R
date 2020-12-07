
source('Codes/Functions.R')
Initialize()


INPUT_NAME =  'rat_DA_M09_WK_008_3pr_v3' # 'rat_LEW_M09_WK_009_3pr_v3' , 'rat_DA_01_reseq'
INPUT_FILE = '2.seur_dimRed_rat_DA_M09_WK_008_3pr_v3_mito_40_lib_1500.rds'
# 2.seur_dimRed_rat_LEW_M09_WK_009_3pr_v3_mito_40_lib_1500.rds
# '2.seur_dimRed_rat_DA_M09_WK_008_3pr_v3_mito_40_lib_1500.rds'
# 'seur_dimRed_rat_DA_01_reseq_mito_30_lib_1500.rds' 
# '2.seur_dimRed_rat_DA_01_reseq_mito_30_lib_1500.rds'


OUTPUT_NAME = gsub('.rds','',gsub('2.seur_dimRed_','',INPUT_FILE ))

seur <- readRDS(paste0('objects/',INPUT_NAME,'/',INPUT_FILE))
output_filename <- paste0('Results/',INPUT_NAME,'/clusters/clusters_',OUTPUT_NAME,'.rds')
seur_output_filename_rds <- paste0('Results/',INPUT_NAME,'/clusters/seur_clustered_',OUTPUT_NAME,'.rds')
dir.create('Results/rat_LEW_M09_WK_009_3pr_v3/clusters')

PC_NUMBER = 22

max_seurat_resolution <- 2.0 ## change this to higher values
FDRthresh <- 0.01 # FDR threshold for statistical tests
min_num_DE <- 10
seurat_resolution <- 0 # Starting resolution is this plus the jump value below.
seurat_resolution_jump <- 0.1

seur <- FindNeighbors(seur,reduction="pca",dims=1:PC_NUMBER,verbose=F)
#seur <- FindNeighbors(seur,reduction="harmony",dims=1:PC_NUMBER,verbose=F)

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
    DRforClust="pca", #pca
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


# shrinks the size of the Seurat object by removing the scaled matrix
#seur <- DietSeurat(seur,dimreducs=Reductions(seur))
#save(sCVdata_list,seur,file=output_filename)

saveRDS(seur,seur_output_filename_rds)
saveRDS(sCVdata_list, output_filename)

seur = readRDS(seur_output_filename_rds)
sCVdata_list = readRDS(output_filename)

