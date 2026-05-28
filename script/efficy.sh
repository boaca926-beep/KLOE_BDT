#!/bin/bash
main_path=/home/kloe/Desktop/
work_space=KLOE_BDT/
output_folder=${main_path}${work_space}plots_efficy/
if [[ -d "$output_folder" ]]; then
    ls $output_folder
    echo "Remove syst. folder"
    rm -rf $output_folder
else
    echo "Create ${output_folder}"
fi
mkdir ${output_folder}
echo "Syst. folder is created at ${output_folder}"
ls ${output_folder}

SAMPLE_TYPE=("bdt" "kloe")

path_header=${main_path}${work_space}header_bdt/efficy_plot.h

for ((i=0;i<${#SAMPLE_TYPE[@]};++i)); do
    sample_type=${SAMPLE_TYPE[i]}

    echo ${sample_type}
    sed -i 's|\(const TString infile =\)\(.*\)|\1 "'"${main_path}input_${sample_type}_TDATA_chain/hist/hist.root"'";|' "$path_header"
    sed -i 's|\(const TString sample_type =\)\(.*\)|\1 "'"${sample_type}"'";|' "$path_header"
    sed -i 's|\(const TString output_folder =\)\(.*\)|\1 "'"${output_folder}"'";|' "$path_header"

    efficy_script=efficy_script.C
    echo "void efficy_script() {" > $efficy_script
    echo '  gROOT->ProcessLine(".L ../run_bdt/get_efficy.C");' >> $efficy_script
    echo '  gROOT->ProcessLine("get_efficy()");' >> $efficy_script
    echo '}' >> $efficy_script
    root -l -n -q -b $efficy_script 
done

echo "Efficiency plots are created!"
rm $efficy_script
