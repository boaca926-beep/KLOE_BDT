#!/bin/bash

header=../header_bdt/plot2d.h
infile_nm="/home/kloe/Desktop/input_kloe_TDATA_chain/hist/hist.root"
output_path=../run_bdt/

echo -e "\nPlotting sfw2d ..."

#rm plots/hist.root 
#rm plots/*.pdf 

echo "Get hist.root"

gethist_script=gethist_script.C
echo '#include <iostream>' > $gethist_script
echo "void gethist_script() {" >> $gethist_script
echo '  gROOT->ProcessLine(".L ../run_bdt/get2Dhist.C");' >> $gethist_script
echo '  gROOT->ProcessLine("get2Dhist()");' >> $gethist_script
echo '}' >> $gethist_script
root -l -n -q -b $gethist_script
rm $gethist_script

# plot sfw2d
#hist_type=("h2d_sfw_TDATA" "h2d_sfw_TISR3PI_SIG" "h2d_sfw_TETAGAM" "hbkgsum_noeta")
#cv_text=("Data" "Signal" "#eta#gamma" "Others" )
#pt1_x0=(0.5 0.5 0.5 0.5)
#pt1_x1=(0.8 0.8 0.8 0.8)

hist_type=("h2d_sfw_TDATA")
cv_text=("Data")
pt1_x0=(0.5)
pt1_x1=(0.8)

#You need to use IFS to stop space as element delimiter.
IFS=""
for ((i=0;i<${#hist_type[@]};++i)); do

    #echo $i
    
    sed -i 's/\(const double pt1_x0 =\)\(.*\)/\1 '${pt1_x0[i]}';/' $header
    sed -i 's/\(const double pt1_x1 =\)\(.*\)/\1 '${pt1_x1[i]}';/' $header
    
    sed -i 's/\(const TString hist_type =\)\(.*\)/\1 "'${hist_type[i]}'";/' $header
    sed -i 's|\(const TString infile_nm =\)\(.*\)|\1 "'"${infile_nm}"'";|' $header
    sed -i 's|\(const TString output_path =\)\(.*\)|\1 "'"${output_path}"'";|' $header
    sed -i 's/\(const TString cv_nm =\)\(.*\)/\1 "'${cv_nm[i]}'";/' $header
    sed -i 's/\(const TString cv_text =\)\(.*\)/\1 "'${cv_text[i]}'";/' $header
    #sed -i 's/\(sprintf(display,\)\(.*\)/\1 "'${cv_text[i]}'");/' sfw2d.C

    sfw2d_script=sfw2d_script.C
    echo '#include <iostream>' > $sfw2d_script
    echo "void sfw2d_script() {" >> $sfw2d_script
    echo '  gROOT->ProcessLine(".L sfw2d.C");' >> $sfw2d_script
    echo '  gROOT->ProcessLine("sfw2d()");' >> $sfw2d_script
    echo '}' >> $sfw2d_script
    #root -l -n -q -b $sfw2d_script
    rm $sfw2d_script
done

# plot mcsum
#hist_type=("hmcsum" "hbkgsum" "hmcsum_noeta")
#cv_nm=("mcsum" "bkgsum" "mcsum_noeta")
#cv_text=("MC_Sum" "Bkg_Sum" "MC_without_#eta#gamma")
#pt1_x0=(0.6 0.6 0.45)
#pt1_x1=(0.7 0.7 0.7)

#hist_type=("hmcsum_noeta")
#cv_nm=("mcsum_noeta")
##cv_text=(" ")
#pt1_x0=(0.45)
#pt1_x1=(0.7)

#for ((i=0;i<${#hist_type[@]};++i)); do

#    #echo $i
    
#    sed -i 's/\(const double pt1_x0 =\)\(.*\)/\1 '${pt1_x0[i]}';/' header
#    sed -i 's/\(const double pt1_x1 =\)\(.*\)/\1 '${pt1_x1[i]}';/' header

#    sed -i 's/\(const TString hist_type =\)\(.*\)/\1 "'${hist_type[i]}'";/' header
#    sed -i 's/\(const TString cv_nm =\)\(.*\)/\1 "'${cv_nm[i]}'";/' header
#    sed -i 's/\(sprintf(display,\)\(.*\)/\1 "'${cv_text[i]}'");/' mcsum.C

#    mcsum_script=mcsum_script.C
#    echo '#include <iostream>' > $mcsum_script
#    echo "void mcsum_script() {" >> $mcsum_script
#    echo '  gROOT->ProcessLine(".L mcsum.C");' >> $mcsum_script
#    echo '  gROOT->ProcessLine("mcsum()");' >> $mcsum_script
#    echo '}' >> $mcsum_script
#    root -l -n -q -b $mcsum_script
#    rm $mcsum_script
    
#done
