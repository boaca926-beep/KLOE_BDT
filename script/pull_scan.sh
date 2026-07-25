#!/bin/bash

sample_type=chain
data_type=bdt
tuning_type=tuning   # raw: kinematic fitted; tuning: kinematic fitted + pi0 decay photon (pull bias correction + scale correction)
pull_status=true
main_folder="/home/bo/Desktop/${data_type}_${tuning_type}_TDATA_${sample_type}_${pull_status}"
input_file_nm="${main_folder}/cut/tree_pre.root";
#input_file_nm="../../${data_type}_${tuning_type}_TDATA_${sample_type}_${pull_status}/cut/tree_pre.root"
outputSfw2D="${main_folder}/sfw2d/";

echo "========================================"
echo "PULL SCAN"
echo "========================================"

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
    #echo "Processing $tree ($sample)"
    
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
root -l -q -b "../run_bdt/scale_compr.C(false)"

# After pull tuning
root -l -q -b "../run_bdt/bias_compr.C(true)" 
root -l -q -b "../run_bdt/scale_compr.C(true)"

echo "All scans completed!"

# Normalize MC mpi0 distribution to data
echo "Get data mpi0"

output_folder="../tuning_${pull_status}_m_gg_bdt"

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

root -l -q -b "../run_bdt/compr_bdt.C(\"$input_file_nm\", \"$output_folder\", \"${outputSfw2D}\", \"m_gg_bdt\", \"M_{#gamma#gamma}\", \"MeV/c^{2}\", 120, 120, 150)"

root -l -q -b "../run_bdt/BiasMgg.C(\"tuning_${pull_status}\", \"m_gg_bdt\", \"M_{#gamma#gamma} [MeV/c^{2}]\")"

MASSBIAS_BDT="../pull_scan/massbias_bdt.txt"
TARGET_FILE="../header_bdt/energy_shift_tuning_sum.h"

if [[ -f "$MASSBIAS_BDT" ]]; then
    mpi0_data=$(grep -oP '(?:const\s+double\s+)?mpi0_data\s*=\s*\K[0-9.eE+-]+' "$MASSBIAS_BDT" | head -1)
    mpi0_data_err=$(grep -oP '(?:const\s+double\s+)?mpi0_data_err\s*=\s*\K[0-9.eE+-]+' "$MASSBIAS_BDT" | head -1)
    mpi0_mc=$(grep -oP '(?:const\s+double\s+)?mpi0_mc\s*=\s*\K[0-9.eE+-]+' "$MASSBIAS_BDT" | head -1)
    mpi0_mc_err=$(grep -oP '(?:const\s+double\s+)?mpi0_mc_err\s*=\s*\K[0-9.eE+-]+' "$MASSBIAS_BDT" | head -1)
    echo "MPI0_DATA = $mpi0_data +/- $mpi0_data_err"
    echo "MPI0_MC = $mpi0_mc +/- $mpi0_mc_err"
fi

echo "Updating $TARGET_FILE with mpi0_data=$mpi0_data, mpi0_data_err=$mpi0_data_err"
sed -i "s/\(const double mpi0_data\s*=\s*\)[0-9.eE+-]*;/\1$mpi0_data;/" "$TARGET_FILE"
sed -i "s/\(const double mpi0_data_err\s*=\s*\)[0-9.eE+-]*;/\1$mpi0_data_err;/" "$TARGET_FILE"
        
echo "Updating $TARGET_FILE with mpi0_mc=$mpi0_mc, mpi0_mc_err=$mpi0_mc_err"
sed -i "s/\(const double mpi0_mc\s*=\s*\)[0-9.eE+-]*;/\1$mpi0_mc;/" "$TARGET_FILE"
sed -i "s/\(const double mpi0_mc_err\s*=\s*\)[0-9.eE+-]*;/\1$mpi0_mc_err;/" "$TARGET_FILE"

        


