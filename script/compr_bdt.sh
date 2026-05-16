#!/bin/bash

echo -e "\nPlotting histo comparison ..."

VAR_NM="lagvalue_min_7C"
VAR_SYMB="#chi^{2}_{7C}"
UNIT=""

XMIN=0
XMAX=100
BINS=200

output_folder="../output_"${VAR_NM[0]}

#check output folder and update output files
if [[ -d $output_folder ]]; then
    
    echo updating $output_folder
    #rm $output_folder/*.pdf
    #rm $output_folder/*.root
    
else
    
    echo root file $output_folder does not exsit;
    mkdir $output_folder
    
fi

compr=../header_bdt/compr.h

for ((i=0;i<${#VAR_NM[@]};++i)); do

    echo ${VAR_NM[i]}
    
    sed -i 's/\(const int binsize =\)\(.*\)/\1 '${BINS[i]}';/' $compr
    sed -i 's/\(const double var_min =\)\(.*\)/\1 '${XMIN[i]}';/' $compr
    sed -i 's/\(const double var_max =\)\(.*\)/\1 '${XMAX[i]}';/' $compr
    
    sed -i 's/\(const TString var_nm =\)\(.*\)/\1 "'${VAR_NM[i]}'";/' $compr
    sed -i 's/\(const TString unit =\)\(.*\)/\1 "'${UNIT[i]}'";/' $compr
    sed -i 's/\(const TString var_symb =\)\(.*\)/\1 "'${VAR_SYMB[i]}'";/' $compr

    compr_script=compr_script.C
    echo '#include <iostream>' > $compr_script
    echo "void compr_script() {" >> $compr_script
    echo '  gROOT->ProcessLine(".L ../run_bdt/compr_bdt.C");' >> $compr_script
    echo '  gROOT->ProcessLine("compr_bdt()");' >> $compr_script
    echo '}' >> $compr_script
    root -l -n -q -b $compr_script 

    plot_script=plot_script.C
    echo '#include <iostream>' > $plot_script
    echo "void plot_script() {" >> $plot_script
    echo '  gROOT->ProcessLine(".L ../run_bdt/plot_compr.C");' >> $plot_script
    echo '  gROOT->ProcessLine("plot_compr()");' >> $plot_script
    echo '}' >> $plot_script
    root -l -n -q -b $plot_script >> output.txt
    rm $plot_script
done

rm $compr_script
rm $plot_script

