#!/bin/bash

compr=../header_bdt/compr.h
sample_type=norm
data_type=bdt
tuning_type=tuning #raw: kinematic fitted; tuning: kinematic fitted + pi0 decay photon (pull bias correction + scale correction)
tuning_status=false
main_folder="/home/bo/Desktop/${data_type}_${tuning_type}_TDATA_${sample_type}_${tuning_status}"
tree_file_nm="${main_folder}/cut/tree_pre.root";
outputSfw2D="${main_folder}/sfw2d/";

echo -e "\nPlotting histo comparison ... from ${main_folder}"

#VAR_NM=("IM3pi_7C" "angle_isr_7C" "angle_pi0gam12" "ppIM" "lagvalue_min_7C" "pvalue" "deltaE" "trkmass" "Eisr" "Emax_clust" "betapi0")
#VAR_SYMB=("M_{3#pi}" "cos#theta_{#gamma_{3}}" "#angle(#gamma_{1},#gamma_{2})" "M_{2#pi}" "#chi^{2}_{7C}" "" "E_{miss}" "M_{trk}" "Eisr" "Emax^{max}_clust" "")
#UNIT=("[MeV\/c^{2}]" "" "[#circ]" "[MeV\/c^{2}]" "" "" "[MeV]" "[MeV\/c^{2}]" "[MeV]" "" "")

#XMIN=(300 -1 0 250 0 0 -800 100 0.3)
#XMAX=(1020 1 180 750 50 1 100 600 1)
#BINS=(200 100 180 200 100 100 200 200 200)

##################################################################
#VAR_NM="bdt_score_max"
#VAR_SYMB="BDT Score"
#UNIT=""

#XMIN=0
#XMAX=1
#BINS=200

##################################################################
#name_tmp="E3"
#VAR_NM="pull_"${name_tmp}
#VAR_SYMB="Pull "${name_tmp}

#UNIT="[MeV]"
#XMIN=-10
#XMAX=10
#BINS=200

##################################################################
#VAR_NM=("m3pi_bdt") #"m3pi_bdt", "IM3pi_7C"
#VAR_SYMB=("M_{3#pi} [MeV\/c^{2}]")
#UNIT=("[MeV\/c^{2}]")

#XMIN=(760) #300 600 760 (analysis)
#XMAX=(800) #1020 1050 800 (analysis)
#BINS=(100)

##################################################################
VAR_NM=("m_gg_bdt") # "IM_pi0_7C", "m_gg_bdt"
VAR_SYMB=("M_{#gamma#gamma}")
UNIT=("[MeV\/c^{2}]")

XMIN=(100)
XMAX=(170)
BINS=(200)

##################################################################
#VAR_NM="Eprompt_max"
#VAR_SYMB="E^{max}_{#gamma}"
#UNIT=""

#XMIN=150
#XMAX=500
#BINS=200

##################################################################
#VAR_NM="betapi0_bdt"
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
#XMAX=20
#BINS=100

##################################################################
#VAR_NM=("pvalue")
#VAR_SYMB=("p-value")
#UNIT=("")

#XMIN=(0)
#XMAX=(1)
#BINS=(200)

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

#XMIN=(100)
#XMAX=(350)
#BINS=(100)

##################################################################
#VAR_NM=("angle_pi0gam12_bdt")
#VAR_SYMB=("#angle_{#gamma#gamma}")
#UNIT=("[#circ]")

#XMIN=(20) #20
#XMAX=(140) #140
#BINS=(180) #120

##################################################################
#VAR_NM=("angle_ppl_pmi")
#VAR_SYMB=("#angle_{+-}")
#UNIT=("[#circ]")

#XMIN=(0) #20
#XMAX=(180) #140
#BINS=(100) #120

##################################################################
#VAR_NM=("angle_trk_neutral")
#VAR_SYMB=("#angle_{trk_neutral}")
#UNIT=("[#circ]")

#XMIN=(150) #20
#XMAX=(180) #140
#BINS=(100) #120

##################################################################
#VAR_NM=("deltaE")
#VAR_SYMB=("E_{diff}")
#UNIT=("[MeV]")

#XMIN=(-460) #-700, -500
#XMAX=(-220) #-200, 50
#BINS=(150) #150, 550 

##################################################################
#VAR_NM="ppIM"
#VAR_SYMB="M_{trk}"
#UNIT="[MeV\/c^{2}]"

#XMIN=300
#XMAX=650
#BINS=100

##################################################################
#VAR_NM="bdt_score"
#VAR_SYMB="BDT value"
#UNIT=""

#XMIN=0
#XMAX=1
#BINS=150

##################################################################
#VAR_NM="e_asym"
#VAR_SYMB=""

#UNIT="[]"
#XMIN=0
#XMAX=1
#BINS=200

##################################################################
#VAR_NM="beta_3pi"
#VAR_SYMB="#beta_{3#pi}"

#UNIT=""
#XMIN=0
#XMAX=1
#BINS=200

output_folder="../${tuning_type}_${tuning_status}_"${VAR_NM[0]}

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

sed -i "s|\(const TString out_dir =\).*|\1 \"${output_folder}\";|" $compr

for ((i=0;i<${#VAR_NM[@]};++i)); do

    echo ${VAR_NM[i]}

    # Use | as delimiter instead of / to avoid conflicts with paths
    sed -i "s|\(const TString tree_file_nm =\).*|\1 \"${tree_file_nm}\";|" $compr
    sed -i "s|\(const TString outputSfw2D =\).*|\1 \"${outputSfw2D}\";|" $compr
    sed -i "s|\(const int binsize =\).*|\1 ${BINS[i]};|" $compr
    sed -i "s|\(const double var_min =\).*|\1 ${XMIN[i]};|" $compr
    sed -i "s|\(const double var_max =\).*|\1 ${XMAX[i]};|" $compr

    sed -i "s|\(const TString var_nm =\).*|\1 \"${VAR_NM[i]}\";|" $compr
    sed -i "s|\(const TString unit =\).*|\1 \"${UNIT[i]}\";|" $compr
    sed -i "s|\(const TString var_symb =\).*|\1 \"${VAR_SYMB[i]}\";|" $compr

    compr_script=compr_script.C
    echo '#include <iostream>' > $compr_script
    echo "void compr_script() {" >> $compr_script
    echo '  gROOT->ProcessLine(".L ../run_bdt/compr_bdt.C");' >> $compr_script
    echo '  gROOT->ProcessLine("compr_bdt()");' >> $compr_script

    #echo '  gROOT->ProcessLine(".L ../run_bdt/compr_bdt_merged.C");' >> $compr_script
    #echo '  gROOT->ProcessLine("compr_bdt_merged()");' >> $compr_script
    echo '}' >> $compr_script
    root -l -n -q -b $compr_script 

    plot_script=plot_script.C
    echo '#include <iostream>' > $plot_script
    echo "void plot_script() {" >> $plot_script
    #echo '  gROOT->ProcessLine(".L ../run_bdt/plot_compr.C");' >> $plot_script
    #echo '  gROOT->ProcessLine("plot_compr()");' >> $plot_script
    
    echo '  gROOT->ProcessLine(".L ../run_bdt/plot_compr_roofit.C");' >> $plot_script
    echo '  gROOT->ProcessLine("plot_compr_roofit()");' >> $plot_script
    echo '}' >> $plot_script
    #root -l -n -q -b $plot_script >> ${output_folder}/output.txt
done

rm $compr_script
rm $plot_script

