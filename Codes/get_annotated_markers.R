
### list of markers:

##########################################################################  
################################## Macrophages
## M1 > inflammatory macrophage
## M2 > non-inflammatory macrophage

## 2 distinct population of CD68 cells

inflammatory_KC <-  c('CD68','LYZ','CSTA','VCAN','CD74', 'IL18','S100A8', 'S100A9')
inflammatory_KC_DE <- c('S100A8', 'LYZ', 'S100A9', 'HLA-DPB1', 'S100A12', 'RP11-1143G9.4', 'EVI2A', 
                        'HLA-DPA1', 'VCAN', 'S100A6', 'CXCL8', 'HLA-DRA', 'MNDA', 'TYR-OBP', 
                        'HLA-DRB1', 'FCN1', 'HLA-DQA1', 'IL18', 'C1QC', 'CD74', 'HLA-DRB5')

non_inflammatory_KC <- c("CD68", 'CD5L', 'MARCO', 'VSIG4' ,'CPVL','CD163', 'VCAM1')
non_inflammatory_KC_DE <- c('CD5L', 'MARCO', 'VSIG4', 'CPVL', 'CD163', 'CCDC88A', 'C5AR1',
                            'LIPA', 'LILRB5', 'MAF', 'CTSB', 'MS4A7', 'VMO1', 'RAB31', 'SLC31A2',
                            'TTYH3', 'VCAM1', 'KLF4', 'HMOX1', 'AIF1l', 'TMIGD3')

M1_KC_cl <- c('CD14','CD68','CD86','IL1B','IL2RA','TNF')
M2_KC_cl <- c('ACTB', 'ARG1', 'AXL', 'BASP1', 'CCL18', 'CCL22', 'CD163', 'CD209', 'CD5L', 'CETP',
              'CLEC7A', 'CXCL12', 'EPB41L2',	'FCGR3A','FOLR2','GPNMB','HMOX1','IL10',
              'ITLN1','LILRB5','LIPA','MARCO','MRC1','MS4A7','PDK4','RBP7','SDC3','SLC40A1','VCAM1')

KC_total <- unique(c(inflammatory_KC, inflammatory_KC_DE, non_inflammatory_KC,
                     non_inflammatory_KC_DE, M1_KC_cl, M2_KC_cl))

##########################################################################  
################################## Stellate cells

## upon activation, express actin, collagen, RBP1: retinol binding protein:a Vit-A associated protein
stellate_cell <- c('ACTA2', 'COL1A1', 'COL1A2', 'TAGLN', 'COL3A1','SPARC', 'RBP1') 
stellate_cells_cl <- c('CRP2', 'DES', 'FAP', 'GFAP', 'MYB', 'NCX', 'PRNP', 'RELN', 'VIM')
stellate_cell_DE <- c('DCN', 'MYL9', 'TPM2', 'MEG3', 'BGN', 'IGFBP7', 'IGFBP3', 
                      'CYR61', 'OLFML3', 'IGFBP6', 'CCL2', 'COLEC11', 'CTGF', 'HGF')

Stellates_total <- unique(c(stellate_cell, stellate_cells_cl, stellate_cell_DE))


##########################################################################  
################################## Cholangiocyte

cholangiocyte <- c('KRT7', 'KRT19','TFF2','SCGB3A1', 'FXYD2','DEFB1','CD24','EPCAM')
cholangiocytes_cl <- c('AFP', 'ALB', 'ANXA4', 'AQP1', 'CD24', 'CFTR', 'CHST4', 'CLDN10', 
                         'CLDN4', 'CXCL1', 'CXCL6',	'DEFB1', 'EPCAM', 'ERICH5','FAM150B', 'FGFR2',
                         'FXYD2', 'GGT1', 'HNF1B', 'HNF4A', 'KRT7', 'MMP7', 'PLPP2',	'S100A14',
                         'SCTR', 'SLC17A4', 'SLC4A2', 'SOX9','SPP1','SSTR2')

cholangiocytes_DE <- c('TFF2', 'SCGB3A1', 'FXYD2', 'KRT7', 'DEFB1', 'CD24', 
                       'TFF3', 'LCN2', 'KRT19', 'CXCL1', 'PIGR', 'TFF1', 'CXCL6', 
                       'LGALS2', 'TACSTD2', 'ELF3', 'SPP1', 'MUC5B', 'LGALS4')

Cholangiocyte_total <- c(cholangiocyte, cholangiocytes_cl, cholangiocytes_DE)


##########################################################################  
################################## Hepatocytes
### 6 populations, less proliferative, enriched expression of albumin(ALB)
### ALB, HAMP, ARG1, PKC1, 

Hepatocytes_cl	<- c('AFP', 'ALB', 'CRP', 'EPHB6','PROX1','REX3', 'TAT', 'TDO2')
interzonal_hepatocytes <- c('HPR','GSTA2HPR','GSTA2', 'AKR1C1')
interzonal_hepatocytes_DE <- c('HPR', 'GSTA2', 'AKR1C1', 'MASP2', 'NNT', 'SAA4', 'MRPS18C', 'OCIAD1', 'APOA5', 'TTR')

midzonal_hepatocytes_DE <- c('RPP25L', 'HSD11B1', 'HAMP', 'GHR', 'APOM', 'APOC4-APOC2',
                             'TKFC', 'G6PC', 'G0S2', 'PON3', 'C1orf53', 'TTC36', 'GOLT1A',
                             'RCAN1', 'RP4-710M16.2', 'FST', 'MCC', 'AQP9', 'PLIN1')

## periportal functions: complement activation, immune activation, fibrinolysis, and triglyceride biosynthesis
periportal_hepatocytes_DE <- c('CYP2A7', 'ENTPD5', 'CYP3A7', 'CYP2A6', 'C4B', 'EID2',
                               'TP53INP2', 'SULT1A1', 'ATIC', 'SERPINH1', 'SAMD5', 
                               'GRB14', 'SEC16B', 'SLBP', 'RND3', 'ABCD3', 'RHOB', 
                               'EPB41L4B', 'GPAT4', 'SPTBN1', 'SDC2', 'PHLDA1', 'WTAP', 'ACADM')

# zone 1: gluconeogenesis and beta-oxidation, genes involved in lipid and cholesterol synthesis
zone1_hepatocytes <- c('SCD','HMGCS1','ACSS2', 'TM7SF2', 'SEC16B', 
                       'SLBP', 'RND3', 'ARG1', 'PCK1', 'CPS1') # last 2 in mice
zone1_hepatocytes_DE <- c('SCD', 'HMGCS1', 'ACSS2', 'TM7SF2', 'TMEM97', 'CP', 'CRP', 
              'SLPI', 'C2orf82', 'ACAT2', 'TM4SF5', 'MSMO1', 'LEPR')

# Zone 3 or central venous hepatocytes play a role in drug metabolism and detoxification, high cytochrome P450 enzymes
zone3_hepatocytes <- c('BCHE', 'G6PC', 'GHR', 'ALDH6A1', 'RCAN1', 'AR',
                       'RPP25L', 'HSD11B1', 'HAMP', 'GHR', 'APOM')

zone3_hepatocytes_DE <- c('BCHE', 'G6PC', 'GHR', 'ALDH6A1', 'RCAN1', 'AR', 
                          'RP4-710M16.2', 'LINC00261', 'PLIN1', 'RP11-390F4.3')


Hepatocytes_total <- unique(c(Hepatocytes_cl, interzonal_hepatocytes, interzonal_hepatocytes_DE, 
                       midzonal_hepatocytes_DE,periportal_hepatocytes_DE, zone1_hepatocytes,
                       zone1_hepatocytes_DE, zone3_hepatocytes, zone3_hepatocytes_DE))


##########################################################################  
################################## Endothelial cells

## 3 population , less proliferative than immune cells,  express: 'CALCRL', 'RAMP2'
### liver sinusoid cells
LSECs_cl <-	c('ADIRF', 'C7', 'CALCRL', 'CD34', 'CD36', 'CLEC14A', 'CTGF', 'F8',
              'FCGR2B', 'FLT1', 'GATA4', 'GJA4','HYAL2', 'ID3', 'LYVE1', 'MGP',
              'NOSTRIN', 'PCAT19', 'PECAM1', 'PRSS23', 'PTGDS', 'RAMP2', 'RAMP3',
              'SPARCL1', 'SRPX', 'STAB2', 'TSPAN7', 'VWF')

# zone1=periportal: high expression of F8, PECAM1, with little expression of CD32B, LYVE-1, STAB2, and CD14 
zone1_LSECs <- c('MGP', 'SPARCL1', 'TMSF1', 'CLEC14A', 'IDI', 'IGFBP7','CTGF','VWF') 
zone1_LSECs_DE <- c('MGP', 'SPARCL1', 'TM4SF1', 'CLEC14A', 'ID1','IGFBP7', 
                         'ADIRF', 'CTGF', 'VWF','CD9', 'C7', 'SRPX','ID3', 
                         'CAV1', 'GNG11', 'AQP1', 'HSPG2', 'EMP1', 'SOX18', 'CLDN5')

# second abundant LSECs: zone 2-3, central venous- 
# enriched expression of CD32B, LYVE1, STAB2, with little expression of VWF
zone2.3_LSECs <- c('CD32B', 'LYVE1', 'STAB2','CCL14', 'CLEC1B','FCN2','S100A13')
zone2.3_LSECs_DE <- c('CCL14', 'CLEC1B', 'FCN2', 'S100A13', 'FCN3', 'CRHBP', 'STAB1', 'GNG11', 
                      'IFI27', 'CLEC4G', 'CLDN5', 'CCL23', 'OIT3', 'RAMP3', 'SGK1', 'DNASE1L3', 
                      'LIFR', 'SPARC', 'ADGRL4', 'EGFL7', 'PCAT19', 'CDKN1C')

# least abundant endothelial cells> low or no expression of LSEC markers (LYVE1, STAB2, CD32B)
# least abundant: These cells are likely non-LSEC endothelial cells including central vein and portal
# arterial and venous endothelial cells based on the expression of ENG (protein alias CD105) and PECAM1 (protein alias CD31) 
portal_endothelial_cells <- c('CD32B', 'CD105', 'VWF', 'PECAM', 'PRSS23', 'RAMP3', 'INMT')
non_LSECs_DE <- c('RAMP3', 'INMT', 'DNASE1L3', 'LIFR', 'PTGDS', 'C7', 'CTGF', 'TIMP3', 
               'RNASE1', 'ID3', 'ENG', 'MGP', 'PCAT19', 'HSPG2', 'GPM6A', 'PTPRB', 
               'VWF', 'FAM167B', 'SRPX', 'LTC4S', 'IFI27')

endothelial_total <- unique(c(LSECs_cl, zone1_LSECs, zone1_LSECs_DE,
                              zone2.3_LSECs, zone2.3_LSECs_DE, 
                              portal_endothelial_cells, non_LSECs_DE))


##########################################################################  
################################## T cells
### 3 clusters of CD3+ T-cells
## conventional T cells: CD4+ and CD8+ > most of CD3+ T cells in the liver
## unconventional T cells: can exoress each either alpha-beta T cells and gamma-delta T cells

T_cells_cl <- c('CD160','CD247', 'CD3E', 'CD4','CD7', 'CD8A')
alpha_beta_T_cells <- c('CD2', 'CD3D', 'TRAC', 'GZMK')
alpha_beta_T_cells_DE <- c('CD2', 'CD3D', 'TRAC', 'GZMK', 'CCL5', 'CCL4L2', 'PYHIN1', 'TRBC1', 'TRBC2',
                           'GZMA', 'CD3E', 'JUNB', 'CD69', 'IL7R', 'DUSP2', 'IFNG', 'LTB', 'IL32', 'CD52')

gamma_delta_T_cells <- c('CD3D', 'CD247', 'NKG7', 'TRDC', 'STMN1', 'MKI67', 'NKG2A', 'TRDC')
gamma_delta_T_cells_DE <- c('GNLY', 'PTGDS', 'GZMB', 'S100B', 'FGFBP2', 'NKG7', 'PRF1', 
                            'KLRF1', 'HOPX', 'CST7', 'KLRD1' ,'CTSW', 'SPON2', 'IFITM1', 
                            'GZMA', 'CD247', 'CLIC3', 'CD7', 'ADGRG1', 'CCL5', 'TRDC')

gamma_delta_T_cells_HE <- c('TBX21', 'KLRB1','CD161', 'FCGR3A', 'CD16', 
                            'NKG7','GNLY', 'NKG5', 'TRDC', 'TRGC1') 


### third population >> enriched expression of GNLY, NKG2A, TYMS, and TOP2A >> more complete in the paper
cytotoxic_cells_cl <-	c('GNLY', 'GZMB', 'PRF1')
Cytotoxic_Cells 



T_Cells <- unique(c(T_cells_cl, alpha_beta_T_cells, alpha_beta_T_cells_DE, 
             gamma_delta_T_cells, gamma_delta_T_cells_DE, gamma_delta_T_cells_HE))


##########################################################################  
################################## NK cells
# low expression of FCGR3A (CD16) or ITGA1 (CD49a) 
NK_cells <- c('CD7', 'KLRB1', 'KLRD1', 'NKG7', 'CD69', 'NCAM1')
NK_cells_HE <- c('CD7', 'KLRD1' ,'CD94', 'GZMK' , 'NCR1' ,'NKp46', 'NCAM1', 'CD56', 'EOMES')
NK_cells_DE <- c('CD7', 'CMC1', 'XCL2', 'KLRB1', 'XCL1', 'KLRC1', 'KLRF1', 'IL2RB', 
                 'CD160', 'CCL3', 'KLRD1', 'NKG7', 'TXK', 'ALOX5AP', 'TRDC', 'CD69', 
                 'TMIGD2', 'CLIC3', 'GZMK', 'DUSP2', 'MATK', 'IFITM1', 'CCL4', 'CD247')


##########################################################################  
##################################  B cells
## 2 populations of liver resident B cells
B_cells_cl <- c('CD19','CD79A','CD79B')

## mature B cells >> without expression of CD27 and CD138
mature_B_cell <- c('CD20', 'LTB', 'CD37', 'CD79B')
mature_B_cell_HE <- c('IGHD', 'CD19', 'MS4A1', 'CD20', 'CD22' ,'CD52')
mature_B_cell_DE <- c('MS4A1', 'LTB', 'CD37', 'CD79B', 'CD52', 'HLA-DQB1', 'TNFRSF13C', 
                      'TCL1A', 'LINC00926', 'STAG3', 'IGHD', 'BANK1', 'IRF8', 'BIRC3', 
                      'P2RX5', 'RP11-693J15.5', 'RP5-887A10.1', 'VPREB3', 'CD22', 'CD74', 'SELL')
## Plasma cells
plasma_cells_cl <- c('IGLC2', 'IGHG1', 'IGKC', 'CD79A', 'CD38','CD27') # and low expression of MS4A1(CD20)
plasma_cells_DE <- c('IGLC2', 'IGHG1', 'IGKC', 'IGHG2', 'IGHG3', 'IGHGP', 'IGLC3', 'JCHAIN', 'IGHA1',
                     'IGHG4', 'IGHA2', 'IGHM', 'IGLV3-1', 'IGLC7', 'MZB1', 'CD79A', 'SSR4', 'IL16')

B_Cells <- c(B_cells_cl, mature_B_cell, mature_B_cell_HE, mature_B_cell_DE, plasma_cells_cl, plasma_cells_DE)








### loading the input data:
INPUT_NAME = 'rat_Rnor'  #'mouse' 
model_animal_name = "rnorvegicus" #'mmusculus'  

Rdata_PATH = paste0('Results/', INPUT_NAME, '/clusters/')
INPUT_FILES = list.files(path = Rdata_PATH , pattern = '.RData', full.names = T, include.dirs = T)
input_version = 1
INPUT_FILE = INPUT_FILES[input_version]
OUTPUT_NAME = gsub('.RData','',gsub(paste0(Rdata_PATH, '/clusters_'),'',INPUT_FILE ))
OUTPUT_NAME = gsub('_v2', '', OUTPUT_NAME)
load(INPUT_FILE)


#### Choose a cluster for evaluation
#### get the list of markers for each of the clusters
a_cluster = 'cluster_0'
tsne_df <- data.frame(getEmb(seur, 'tSNE'), clusters=paste0('cluster_',as.character(Idents(seur))))
tsne_df$clusters <- ifelse(tsne_df$clusters == a_cluster, a_cluster, 'other')
ggplot(tsne_df, aes(x=tSNE_1,y=tSNE_2,color=clusters))+geom_point()+theme_classic()+ggtitle(a_cluster)
head(tsne_df)

Cluster_markers_merged <- readRDS(paste0('Results/', INPUT_NAME, '/markers/markers_', OUTPUT_NAME, '.rds'))
lapply(Cluster_markers_merged, head)
head(Cluster_markers_merged[[a_cluster]])
Cluster_markers_merged[[a_cluster]]$V2[1:20]


#### 
source('Codes/convert_human_to_ortholog_functions.R')
seur_genes_df <- read.delim(paste0(paste0("Data/", INPUT_NAME,'/'),'genes.tsv'), header = F)


sample_name = 'DA-02'
cell_type <- 'inflammatory KC'
human_markers <- inflammatory_KC
mapped_markers_df <- .getMapped_hs2model_df(ensembl, human_markers , model_animal_name)
mapped_markers_df_2 <- mapped_markers_df[mapped_markers_df$rnorvegicus_homolog_ensembl_gene != '' & 
                                           !is.na(df$rnorvegicus_homolog_ensembl_gene),]
mapped_markers_df_2 <- mapped_markers_df_2[mapped_markers_df_2$rnorvegicus_homolog_ensembl_gene %in% rownames(seur),]

for(i in 1:nrow(mapped_markers_df_2)){
  aMarker_symbol = mapped_markers_df_2[i,]$rnorvegicus_homolog_associated_gene_name
  aMarker_ensemble_id_index = mapped_markers_df_2[i,]$rnorvegicus_homolog_ensembl_gene
  aMarker_expression <- exprMatrix[aMarker_ensemble_id_index,]
  tsne_df <- as.data.frame(cbind(getEmb(seur, 'tsne'),gene.expr=aMarker_expression))
  p = Plot.tsne.gene.expr(tsne_df, 
                          title = paste0('human:', human_marker,'  ', 
                                       model_animal_name,':',aMarker_symbol),
                          subtitle = paste0('Cell type: ',cell_type, '   Sample: ', sample_name))
  print(p)
}






