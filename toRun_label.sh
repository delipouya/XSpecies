#!/bin/bash

# INPUT_NAME = 'rat_Rnor'
# INPUT_FILE = '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'
for input_name in rat_Rnor rat_DA rat_Lew_01 rat_Lew_02; do
  for input_file in objects/${input_name}/2.seur_dimRed_* ; do 
    echo ${input_name}
    echo ${input_file##*/}
    Rscript Codes/get_labels_AUCell.R ${input_name} ${input_file##*/}
    #Rscript Codes/get_labels_GSVA.R ${input_name} ${input_file##*/}
    Rscript Codes/get_labels_SCINA.R ${input_name} ${input_file##*/}
    Rscript Codes/get_integrated_labels.R ${input_name} ${input_file##*/}
  done
done


