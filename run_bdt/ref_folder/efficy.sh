#!/bin/bash

# input folder
sample_size=small # norm_TUFO small_TUFO
syst_type=filfo # evtcls, filfo
exp_type=TUFO # UFO
input_path=../../input_${sample_size}_${syst_type}_${exp_type}
#echo $input_path

if [[ -d "$input_path" ]]; then
    echo "${input_path} exist!"
    ls ${input_path}
else
    echo "Folder ${input_path} doesn't exist!"
    exit 1
fi

# Main folder
efficy_path=../../efficy_${syst_type}
if [[ -d "$efficy_path" ]]; then
    ls $efficy_path
    echo "Remove syst. folder"
    rm -rf $efficy_path
else
    echo "Create ${efficy_path}"
fi
mkdir ${efficy_path}
echo "Syst. folder is created at ${efficy_path}"
ls ${efficy_path}

# Log files 
log_efficy=${efficy_path}/log_efficy.txt
echo "" > ${log_efficy}
#touch $log_efficy

#TREE_TYPE=("TISR3PI_SIG" "ksl" "exp" "eeg" "TUFO")
TREE_TYPE=("TISR3PI_SIG" "TUFO")

# loop over sample types
path_header=../header/efficy.h

sed -i 's|\(const TString inputFile =\)\(.*\)|\1 "'"${input_path}/cut/tree_pre"'";|' "$path_header"
    
for ((i=0;i<${#TREE_TYPE[@]};++i)); do

    treeType=${TREE_TYPE[i]}
    
    sed -i 's|\(const TString outputCut =\)\(.*\)|\1 "'"${efficy_path}/"'";|' "$path_header"
    sed -i 's|\(const TString treeType =\)\(.*\)|\1 "'"${treeType}"'";|' "$path_header"

    efficy_script=efficy_script.C
    echo "void efficy_script() {" > $efficy_script
    echo '  gROOT->ProcessLine(".L ../run_ufo/efficy_'${syst_type}'.C");' >> $efficy_script
    echo '  gROOT->ProcessLine("efficy_'${syst_type}'()");' >> $efficy_script
    echo '}' >> $efficy_script
    root -l -n -q -b $efficy_script >> ${log_efficy}
    
done
echo "Efficiency is calculated!"
#rm $efficy_script

input_header="../header/plot_efficy.h"

sed -i 's|\(const TString input_folder =\)\(.*\)|\1 "'"${efficy_path}"'";|' "$input_header"

plot_efficy_script=plot_efficy_script.C
echo "void plot_efficy_script() {" > $plot_efficy_script
echo '  gROOT->ProcessLine(".L ../run_plot/plot_efficy.C");' >> $plot_efficy_script
echo '  gROOT->ProcessLine("plot_efficy()");' >> $plot_efficy_script
echo '}' >> $plot_efficy_script
root -l -n -q -b $plot_efficy_script

rm $plot_efficy_script
