#!/bin/bash
read input_name
echo ${input_name}
# mkdir Data/${input_name}
mkdir objects/${input_name}
mkdir Results/${input_name}
mkdir Results/${input_name}/QC
mkdir Results/${input_name}/AUCell
mkdir Results/${input_name}/SCINA

