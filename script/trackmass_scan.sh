#!/bin/bash

sample_type=norm
data_type=bdt
tuning_type=raw   # raw: kinematic fitted; tuning: kinematic fitted + pi0 decay photon (pull bias correction + scale correction)
pull_status=false

main_folder="/home/bo/Desktop/${data_type}_${tuning_type}_TDATA_${sample_type}_${pull_status}"
input_file_nm="${main_folder}/cut/tree_pre.root";
#input_file_nm="../../${data_type}_${tuning_type}_TDATA_${sample_type}_${pull_status}/cut/tree_pre.root"
outputSfw2D="${main_folder}/sfw2d/";

echo "========================================"
echo "PULL SCAN"
echo "========================================"

echo -e "\nPull scan ... using ${main_folder}"

# Define trees and sample types
output_folder="../trackmass_scan"

TREES_MC=("TISR3PI_SIG_PEAK")
SAMPLES_MC=("Signal")

# Create output folder if it doesn't exist, and clean it
if [[ -d $output_folder ]]; then
    echo "Updating $output_folder"
    rm -f $output_folder/*.pdf $output_folder/*.h $output_folder/*.png $output_folder/*.root $output_folder/*.zip $output_folder/*#
else
    echo "Output folder $output_folder does not exist; creating it."
    mkdir -p $output_folder
fi

TREES=("TISR3PI_SIG_PEAK")
SAMPLES=("Signal")

# Loop over trees and run the macro for DATA
for ((i=0; i<${#TREES[@]}; ++i)); do
    tree=${TREES[i]}
    sample=${SAMPLES[i]}
    echo "Processing $tree ($sample)"
    
    # Run ROOT in batch mode, passing arguments to pull_scan.C
    root -l -q -b "../run_bdt/trackmass_scan.C(\"$tree\", \"$sample\", \"pull\", true, \"$input_file_nm\", \"gausPoly\", \"old\")" # old pull

done

#================================================
# Track bias correction
#================================================
input_file="../trackmass_scan/pull_scan_TISR3PI_SIG_PEAK.root"
root -l -q -b "../run_bdt/fit_trk_bias_mc.C(\"${input_file}\")"
