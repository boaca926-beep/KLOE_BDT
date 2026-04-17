#!/bin/bash

#FILE_TYPE=("vertex" "norm" "vmd_norm")
#FILE_TYPE=("vmd_mini" "mini")
#FILE_TYPE=("norm" "vmd_norm" "isrlumi_norm")
FILE_TYPE=("chain")
#FILE_TYPE=("norm")
#FILE_TYPE=("mini")
#FILE_TYPE=("vertex" "efficy_norm" "norm" "vmd_norm" "vmd_efficy_corr_norm")
#FILE_TYPE=("efficy_norm")
#FILE_TYPE=("efficy_norm")
#FILE_TYPE=("norm")
#FILE_TYPE=("isrlumi_norm" "norm")
#FILE_TYPE=("vmd_norm" "norm" "efficy_norm" "isrlumi_norm")

#sample_size="mini"
file_type1=${FILE_TYPE[0]}
file_type2=${FILE_TYPE[1]}

# Plot main folder
plot_folder=../plots_omega/
if [[ -d "$plot_folder" ]]; then
    ls $plot_folder
    echo "Remove plot folder"
    rm -rf $plot_folder
else
    echo "Create ${plot_folder}"
fi
mkdir $plot_folder


header="../header_bdt/plot_omega.h"

for ((i=0; i<${#FILE_TYPE[@]}; ++i)); do

    file_type=${FILE_TYPE[i]}
    
    file_indx=$(echo "$i  " | bc)
    echo $file_indx $file_type
    
    #input_file="/home/bo/Desktop/analysis_root_v6/analysis_${file_type}_TDATA/omega_fit/omega_fit.root"
    input_file="../../input_bdt_TDATA_${file_type}/omega_fit/omega_fit.root"
    #input_file="../../analysis_${file_type}_TDATA/omega_fit/omega_fit.root"
    #input_file="../../input_${file_type}_TDATA_bkgrej/omega_fit/omega_fit.root"
    
    sed -i 's|\(const TString input_file =\)\(.*\)|\1 "'"${input_file}"'";|' "$header"
    sed -i 's|\(const TString file_type =\)\(.*\)|\1 "'"${file_type}"'";|' "$header"
    sed -i 's|\(const TString outputPlot =\)\(.*\)|\1 "'"${plot_folder}"'";|' "$header"
    
    plot_script=plot_script.C
    echo '#include <iostream>' > $plot_script
    echo "void plot_script() {" >> $plot_script
    echo '  gROOT->ProcessLine(".L ../run_bdt/plot_omega.C");' >> $plot_script
    echo '  gROOT->ProcessLine("plot_omega()");' >> $plot_script
    echo '}' >> $plot_script
    root -l -n -q -b $plot_script

done

#header_compr="../header_bdt/plot_omega_compr.h"

#echo $file_type1 $file_type2

#sed -i 's|\(const TString file_type1 =\)\(.*\)|\1 "'"${file_type1}"'";|' "$header_compr"
#sed -i 's|\(const TString file_type2 =\)\(.*\)|\1 "'"${file_type2}"'";|' "$header_compr"
#sed -i 's|\(const TString outputFolder =\)\(.*\)|\1 "'"${plot_folder}"'";|' "$header_compr"

#plot_compr_script=plot_compr_script.C
#echo '#include <iostream>' > $plot_compr_script
#echo "void plot_compr_script() {" >> $plot_compr_script
#echo '  gROOT->ProcessLine(".L ../run_plot/plot_omega_compr.C");' >> $plot_compr_script
#echo '  gROOT->ProcessLine("plot_omega_compr()");' >> $plot_compr_script
#echo '}' >> $plot_compr_script
#root -l -n -q -b $plot_compr_script


rm $plot_script
#rm $plot_compr_script
#ls $plot_folder
