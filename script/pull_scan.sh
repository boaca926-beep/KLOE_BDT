#!/bin/bash

header=../header_bdt/pull_scan.h

sample_type=norm
data_type=bdt
tuning_type=tuning #raw: kinematic fitted; tuning: kinematic fitted + pi0 decay photon (pull bias correction + scale correction)
main_folder="/home/bo/Desktop/${data_type}_${tuning_type}_TDATA_${sample_type}"

echo -e "\nPull scan ... using ${main_folder}"

TREE_NM=("TISR3PI_SIG_PEAK" "TDATA")
SAMPLE_TYPE=("Signal" "Data")

output_folder="../pull_scan"

#check output folder and update output files
if [[ -d $output_folder ]]; then
    
    echo updating $output_folder
    rm $output_folder/*.pdf
    rm $output_folder/*.png
    rm $output_folder/*.root
    
else
    
    echo root file $output_folder does not exsit;
    mkdir $output_folder
    
fi
