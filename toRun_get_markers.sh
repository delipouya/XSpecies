for input_name in rat_Rnor rat_DA rat_Lew_01 rat_Lew_02; do
  echo ${input_name}
  Rscript Codes/get_cluster_markers.R ${input_name} 'rnorvegicus'
done

Rscript Codes/get_cluster_markers.R 'mouse' 'mmusculus'
