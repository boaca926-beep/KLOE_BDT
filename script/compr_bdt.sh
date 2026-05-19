#!/bin/bash

echo -e "\nPlotting histo comparison ..."

#VAR_NM=("IM3pi_7C" "angle_isr_7C" "angle_pi0gam12" "ppIM" "lagvalue_min_7C" "pvalue" "deltaE" "trkmass" "Eisr" "Emax_clust" "betapi0")
#VAR_SYMB=("M_{3#pi}" "cos#theta_{#gamma_{3}}" "#angle(#gamma_{1},#gamma_{2})" "M_{2#pi}" "#chi^{2}_{7C}" "" "E_{miss}" "M_{trk}" "Eisr" "Emax^{max}_clust" "")
#UNIT=("[MeV\/c^{2}]" "" "[#circ]" "[MeV\/c^{2}]" "" "" "[MeV]" "[MeV\/c^{2}]" "[MeV]" "" "")

#XMIN=(300 -1 0 250 0 0 -800 100 0.3)
#XMAX=(1020 1 180 750 50 1 100 600 1)
#BINS=(200 100 180 200 100 100 200 200 200)

##################################################################
#VAR_NM="betapi0"
#VAR_SYMB="#beta_{#pi}"
#UNIT=""

#XMIN=0.3
#XMAX=1
#BINS=150

##################################################################
#VAR_NM="lagvalue_min_7C"
#VAR_SYMB="#chi^{2}_{7C}"
#UNIT=""

#XMIN=0
#XMAX=45
#BINS=150

##################################################################
#VAR_NM=("pvalue")
#VAR_SYMB=("p-value")
#UNIT=("")

#XMIN=(0)
#XMAX=(1)
#BINS=(200)

##################################################################
VAR_NM=("m3pi_bdt")
VAR_SYMB=("M_{3#pi}")
UNIT=("[MeV\/c^{2}]")

XMIN=(770) #300 600
XMAX=(800) #1020 1050
BINS=(100)

##################################################################
#VAR_NM=("Eisr")
#VAR_SYMB=("E_{#gamma_{3}}")
#UNIT=("[MeV]")

#XMIN=(0)
#XMAX=(500)
#BINS=(200)

##################################################################
#VAR_NM=("Eprompt_max")
#VAR_SYMB=("E_{#gamma}^{max}")
#UNIT=("[MeV]")

#XMIN=(50)
#XMAX=(500)
#BINS=(450)

##################################################################
#VAR_NM=("angle_pi0gam12")
#VAR_SYMB=("#angle_{#gamma#gamma}")
#UNIT=("[#circ]")

#XMIN=(0) #20
#XMAX=(180) #140
#BINS=(180) #120

##################################################################
#VAR_NM=("deltaE")
#VAR_SYMB=("E_{diff}")
#UNIT=("[MeV]")

#XMIN=(-700) #-700, -500
#XMAX=(-150) #-200, 50
#BINS=(200) #150, 550 

##################################################################
#VAR_NM=("IM_pi0_7C")
#VAR_SYMB=("M_{#gamma#gamma}")
#UNIT=("[MeV\/c^{2}]")

#XMIN=(100)
#XMAX=(180)
#BINS=(180)

##################################################################
#VAR_NM="trkmass"
#VAR_SYMB="M_{trk}"
#UNIT="[MeV\/c^{2}]"

#XMIN=100
#XMAX=450
#BINS=550


output_folder="../output_"${VAR_NM[0]}

#check output folder and update output files
if [[ -d $output_folder ]]; then
    
    echo updating $output_folder
    rm $output_folder/*.pdf
    rm $output_folder/*.root
    
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

