#!/bin/bash

sample_type=chain
data_type=bdt
tuning_type=tuning   # raw: kinematic fitted; tuning: kinematic fitted + pi0 decay photon (pull bias correction + scale correction)
input_file_nm="/home/bo/Desktop/${data_type}_${tuning_type}_TDATA_${sample_type}/cut/tree_pre.root"

echo -e "\nPull scan ... using ${main_folder}"

# Define trees and sample types
TREES_DATA=("TDATA")
SAMPLES_DATA=("Data")

TREES_MC=("TISR3PI_SIG_PEAK")
SAMPLES_MC=("Signal")

output_folder="../pull_scan"

# Create output folder if it doesn't exist, and clean it
if [[ -d $output_folder ]]; then
    echo "Updating $output_folder"
    rm -f $output_folder/*.pdf $output_folder/*.png $output_folder/*.root
else
    echo "Output folder $output_folder does not exist; creating it."
    mkdir -p $output_folder
fi

# Loop over trees and run the macro for MC
for ((i=0; i<${#TREES_MC[@]}; ++i)); do
    tree=${TREES_MC[i]}
    sample=${SAMPLES_MC[i]}
    echo "Processing $tree ($sample)"
    
    # Run ROOT in batch mode, passing arguments to pull_scan.C
    root -l -q -b "../run_bdt/pull_scan.C(\"$tree\", \"$sample\", \"pull\", true, \"$input_file_nm\", \"crystalBall\", \"old\")" # old pull

    root -l -q -b "../run_bdt/pull_scan.C(\"$tree\", \"$sample\", \"pull\", true, \"$input_file_nm\", \"crystalBall\", \"new\")" # new pull
    
done

# Loop over trees and run the macro for DATA
for ((i=0; i<${#TREES_DATA[@]}; ++i)); do
    tree=${TREES_DATA[i]}
    sample=${SAMPLES_DATA[i]}
    echo "Processing $tree ($sample)"
    
    # Run ROOT in batch mode, passing arguments to pull_scan.C
    root -l -q -b "../run_bdt/pull_scan.C(\"$tree\", \"$sample\", \"pull\", true, \"$input_file_nm\", \"crystalBall\", \"old\")" # old pull

done

# Before pull tuning
root -l -q -b "../run_bdt/bias_compr.C(false)" 
root -l -q -b "../run_bdt/resol_compr.C(false)"

# After pull tuning
root -l -q -b "../run_bdt/bias_compr.C(true)" 
root -l -q -b "../run_bdt/resol_compr.C(true)"

echo "All scans completed."
