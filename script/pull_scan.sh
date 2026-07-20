#!/bin/bash

sample_type=norm
data_type=bdt
tuning_type=tuning   # raw: kinematic fitted; tuning: kinematic fitted + pi0 decay photon (pull bias correction + scale correction)
input_file_nm="/home/bo/Desktop/${data_type}_${tuning_type}_TDATA_${sample_type}/cut/tree_pre.root"

echo -e "\nPull scan ... using ${main_folder}"

# Define trees and sample types
TREES=("TISR3PI_SIG_PEAK" "TDATA")
SAMPLES=("Signal" "Data")

#TREES=("TDATA")
#SAMPLES=("Data")

output_folder="../pull_scan"

# Create output folder if it doesn't exist, and clean it
if [[ -d $output_folder ]]; then
    echo "Updating $output_folder"
    rm -f $output_folder/*.pdf $output_folder/*.png $output_folder/*.root
else
    echo "Output folder $output_folder does not exist; creating it."
    mkdir -p $output_folder
fi

# Loop over trees and run the macro
for ((i=0; i<${#TREES[@]}; ++i)); do
    tree=${TREES[i]}
    sample=${SAMPLES[i]}
    echo "Processing $tree ($sample)"
    
    # Run ROOT in batch mode, passing arguments to pull_scan.C
    root -l -q -b "../run_bdt/pull_scan.C(\"$tree\", \"$sample\", \"pull\", true, \"$input_file_nm\", \"crystalBall\")"
done

root -l -q -b "../run_bdt/pull_compr.C()"

echo "All scans completed."
