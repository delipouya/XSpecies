### getting the cleaned list of Tallulah markers

source('Codes/convert_human_to_ortholog_functions.R')
INPUT_NAME = 'rat_DA' #'mouse'
model_animal_name = "rnorvegicus" #'mmusculus' 

### import input data
general_markers <- read.csv('Data/McParland_markers/liver_markers_tallulah/liver_general_markers.csv')


### adding to human attributes
general_markers_human <- general_markers[general_markers$Species=='Human',]
#### convert the human ensemble ids to rat orthologs
general_markers_human_mapped = .getMapped_hs2model_df(ensembl, 
                                                      candidateGenes=general_markers_human$Gene, 
                                                      model_animal_name)
general_markers_human_mapped <- merge(general_markers_human, general_markers_human_mapped, by.x='Gene', by.y='symbol') 


### adding to rat attributes
general_markers_rat <- general_markers[general_markers$Species=='Rat',]
general_markers_rat_df <- get_rat_ensembl_ids(general_markers_rat$Gene)
general_markers_rat_df <- merge(general_markers_rat, general_markers_rat_df, by.x='Gene', by.y='rgd_symbol', all.x=T)



general_markers_df <- list(general_markers_human_mapped, general_markers_rat_df)
names(general_markers_df) <- c('human_general', 'rat_general')
saveRDS(general_markers_df, 'Data/McParland_markers/liver_markers_tallulah/liver_general_markers.rds')




