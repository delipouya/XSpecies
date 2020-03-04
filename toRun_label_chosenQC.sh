
for input_name in rat_Rnor rat_DA rat_Lew_01 rat_Lew_02; do
  echo ${input_name}
  Rscript Codes/get_labels_AUCell.R ${input_name} 'rnorvegicus'
  #Rscript Codes/get_labels_GSVA.R ${input_name} 'rnorvegicus'
  Rscript Codes/get_labels_SCINA.R ${input_name} 'rnorvegicus'
  Rscript Codes/get_integrated_labels.R ${input_name} 'rnorvegicus'
done


Rscript Codes/get_labels_AUCell.R 'mouse' 'mmusculus'
Rscript Codes/get_labels_SCINA.R 'mouse' 'mmusculus'
Rscript Codes/get_integrated_labels.R 'mouse' 'mmusculus'