

# mouse
Rscript Codes/get_labels_AUCell.R 'mouse' ${input_file##*/}
Rscript Codes/get_labels_GSVA.R 'mouse' ${input_file##*/}
Rscript Codes/get_labels_SCINA.R 'mouse' ${input_file##*/}
Rscript Codes/get_integrated_labels.R 'mouse'  ${input_file##*/}

# rat_Rnor   
Rscript Codes/get_labels_AUCell.R 'rat_Rnor' '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'
Rscript Codes/get_labels_GSVA.R 'rat_Rnor' '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'
Rscript Codes/get_labels_SCINA.R 'rat_Rnor' '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'
Rscript Codes/get_integrated_labels.R 'rat_Rnor' '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'

# rat_DA
Rscript Codes/get_labels_AUCell.R 'rat_DA' '2.seur_dimRed_rat_DA_mito_30_lib_1500.rds'
Rscript Codes/get_labels_GSVA.R 'rat_DA' '2.seur_dimRed_rat_DA_mito_30_lib_1500.rds'
Rscript Codes/get_labels_SCINA.R 'rat_DA' '2.seur_dimRed_rat_DA_mito_30_lib_1500.rds'
Rscript Codes/get_integrated_labels.R 'rat_DA' '2.seur_dimRed_rat_DA_mito_30_lib_1500.rds'

# rat_Lew_01
Rscript Codes/get_labels_AUCell.R 'rat_Lew_01' ${input_file##*/}
Rscript Codes/get_labels_GSVA.R 'rat_Lew_01' ${input_file##*/}
Rscript Codes/get_labels_SCINA.R 'rat_Lew_01' ${input_file##*/}
Rscript Codes/get_integrated_labels.R 'rat_Lew_01' ${input_file##*/}

# rat_Lew_02
Rscript Codes/get_labels_AUCell.R 'rat_Lew_02' ${input_file##*/}
Rscript Codes/get_labels_GSVA.R 'rat_Lew_02' ${input_file##*/}
Rscript Codes/get_labels_SCINA.R 'rat_Lew_02' ${input_file##*/}
Rscript Codes/get_integrated_labels.R 'rat_Lew_02' ${input_file##*/}