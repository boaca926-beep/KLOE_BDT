#!/bin/bash

##################################################################
#VAR_TYPE=("Br_betapi0_bdt"
#	  "Br_e1_bdt"
#	  "Br_e2_bdt"
#	  "Br_e3_bdt"
#	  "Br_m3pi_bdt"
#	  "Br_angle_pi0gam12_bdt"
#	  "Br_m_gg_bdt")

#VAR_TYPE_TRUE=("Br_betapi0_bdt_true"
#	       "Br_e1_bdt_true"
#	       "Br_e2_bdt_true"
#	       "Br_e3_bdt_true"
#	       "Br_m3pi_true_bdt"
#	       "Br_angle_pi0gam12_bdt_true"
#	       "Br_m_gg_true_bdt")

#UNIT=(""
#      "[MeV]"
#      "[MeV]"
#      "[MeV]"
#      "[MeV/c^{2}]"
#      "[#circ]"
#      "[MeV/c^{2}]")

#XTITLE=("#beta^{rec}_{#pi}-#beta^{true}_{#pi}"
#	"E^{rec}_{1}-E^{true}_{1}"
#	"E^{rec}_{2}-E^{true}_{2}"
#	"E^{rec}_{3}-E^{true}_{3}"
#	"M^{rec}_{3#pi}-M^{true}_{3#pi}"
#	"#angle^{rec}_{#gamma#gamma}-#angle^{true}_{#gamma#gamma}"
#	"M^{rec}_{#gamma#gamma}-M^{true}_{#gamma#gamma}")

#BIN_SIZE=(400
#	  1000
#	  1000
#	  1000
#	  1000
#	  4000
#	  1000)

#FIT_FACTOR=(0.5
#	    1.0
#	    1.0
#	    .5
#	    1.0
#	    1.0
#	    1.0)

# UPDATED: More reasonable range factors
#RANGE_FACTOR=(2
#	      5.0
#	      5.0
#	      5.0
#	      1.5
#	      5.0
#	      1.5)

#XMIN=(-0.05
#      -100
#      -100
#      -100
#      -100
#      -100
#      -100)

#XMAX=(0.05
#      100
#      100
#      100
#      100
#      100
#      100)

#==============================================================
# Resolution of kinematic variables (KLOE Raw)
#==============================================================
#VAR_TYPE=("Br_E1" 
#	  "Br_E2" 
#	  "Br_E3")

#VAR_TYPE_TRUE=("Br_E1_true" 
#	       "Br_E2_true" 
#	       "Br_E3_true")

VAR_TYPE=("Br_e1_bdt"
	  "Br_e2_bdt"
	  "Br_e3_bdt")

VAR_TYPE_TRUE=("Br_e1_bdt_true"
	       "Br_e2_bdt_true"
	       "Br_e3_bdt_true")

UNIT=("[MeV]" 
      "[MeV]" 
      "[MeV]")

XTITLE=("E^{rec}_{1}-E^{true}_{1}" 
	"E^{rec}_{2}-E^{true}_{2}" 
	"E^{rec}_{3}-E^{true}_{3}")

BIN_SIZE=(500 
	  500 
	  500)

FIT_FACTOR=(1.0 
	    1.0 
	    0.5)

RANGE_FACTOR=(2 
	      2 
	      2)

XMIN=(-50 
      -50 
      -20)

XMAX=(50 
      50 
      20)


output_folder="../plots_resol"

#check output folder and update output files
if [[ -d $output_folder ]]; then
    
    echo updating $output_folder
    rm $output_folder/*.png
else
    
    echo root file $output_folder does not exsit;
    mkdir $output_folder
    
fi

header=../header_bdt/plot_resol.h
sample_type=norm
data_type=bdt
tuning_type=raw
main_folder="/home/bo/Desktop/${data_type}_${tuning_type}_TDATA_${sample_type}"
#main_folder=/home/bo/Desktop/bdt_raw_TDATA_norm
treeFile="${main_folder}/cut/tree_pre.root";

for ((i=0;i<${#VAR_TYPE[@]};++i)); do

    echo ${VAR_TYPE[i]} ${VAR_TYPE_TRUE[i]}
    xtitle=${XTITLE[i]}" "${UNIT[i]}
    echo $xtitle
    
    # Use | as delimiter instead of / to avoid conflicts with paths
    sed -i "s|\(const TString treeFile =\).*|\1 \"${treeFile}\";|" $header
    sed -i "s|\(const TString var_type =\).*|\1 \"${VAR_TYPE[i]}\";|" $header
    sed -i "s|\(const TString var_type_true =\).*|\1 \"${VAR_TYPE_TRUE[i]}\";|" $header
    sed -i "s|\(const TString unit =\).*|\1 \"${UNIT[i]}\";|" $header
    sed -i "s|\(const TString x_title =\).*|\1 \"${xtitle}\";|" $header
    sed -i "s|\(const int bin_size =\).*|\1 ${BIN_SIZE[i]};|" $header
    sed -i "s|\(const double fit_factor =\).*|\1 ${FIT_FACTOR[i]};|" $header
    sed -i "s|\(const double range_factor =\).*|\1 ${RANGE_FACTOR[i]};|" $header
    sed -i "s|\(const double XMIN =\).*|\1 ${XMIN[i]};|" $header
    sed -i "s|\(const double XMAX =\).*|\1 ${XMAX[i]};|" $header

    plot_script=plot_script.C
    echo '#include <iostream>' > $plot_script
    echo "void plot_script() {" >> $plot_script
    echo '  gROOT->ProcessLine(".L ../run_bdt/plot_resol.C");' >> $plot_script
    echo '  gROOT->ProcessLine("plot_resol()");' >> $plot_script
    echo '}' >> $plot_script
    root -l -n -q -b $plot_script 
done

rm $plot_script
