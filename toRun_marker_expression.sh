#!/bin/bash

Rscript Codes/get_expression_tsne.R 'rat_Rnor' '2.seur_dimRed_rat_Rnor_mito_50_lib_1500.rds'
Rscript Codes/get_expression_tsne.R 'rat_DA' '2.seur_dimRed_rat_DA_mito_30_lib_1500.rds'
Rscript Codes/get_expression_tsne.R 'rat_Lew_01' '2.seur_dimRed_rat_Lew_01_mito_40_lib_2000.rds'
Rscript Codes/get_expression_tsne.R 'rat_Lew_02' '2.seur_dimRed_rat_Lew_02_mito_40_lib_2000.rds'