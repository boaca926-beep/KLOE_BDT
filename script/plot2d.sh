#!/bin/bash

# Quick start
# ./plot2d.sh sfw2d             # for the original sideband plots
# ./plot2d.sh ppIM_vs_betapi0   # for the new dipion mass vs π⁰ β correlation

hist_root="/home/kloe/Desktop/input_bdt_TDATA_chain/hist/hist.root"

# Usage: ./plot2d.sh [sfw2d|ppIM_vs_betapi0]
plot_type=${1:-sfw2d}
header=../header_bdt/plot2d.h
output_path=../run_bdt/
plots_dir=../plots2d_${plot_type}/

# ----------------------------------------------------------------------
# Create required directories
# ----------------------------------------------------------------------
mkdir -p $output_path
mkdir -p $plots_dir

echo -e "\nPlotting $plot_type ..."

sed -i 's|\(const TString hist_root =\)\(.*\)|\1 "'"${hist_root}"'";|' $header

# ----------------------------------------------------------------------
# Get hist.root and produce combined 2D outputs (sfw2d or ppIM_vs_betapi0)
# ----------------------------------------------------------------------
gethist_script=gethist_script.C
cat > $gethist_script <<EOF
#include <iostream>
void gethist_script() {
    gROOT->ProcessLine(".L ../run_bdt/get2Dhist.C");
    gROOT->ProcessLine("get2Dhist(\"$plot_type\")");
}
EOF
root -l -n -q -b $gethist_script
rm $gethist_script

# ----------------------------------------------------------------------
# Configure arrays and macro based on plot type
# ----------------------------------------------------------------------
if [ "$plot_type" == "sfw2d" ]; then
    outfile_name="sfw2d_output.root"
    hist_type=("h2d_sfw_TDATA" "h2d_sfw_TISR3PI_SIG" "h2d_sfw_TETAGAM" "hbkgsum_noeta")
    cv_nm=("data" "signal" "etagam" "bkgsum_noeta")
    cv_text=("Data" "Signal" "#eta#gamma" "Others")
    pt1_x0=(0.5 0.5 0.5 0.5)
    pt1_x1=(0.8 0.8 0.8 0.8)
    macro="plot2d_sfw.C"
    func="plot2d_sfw"
elif [ "$plot_type" == "ppIM_vs_betapi0" ]; then
    outfile_name="ppIM_vs_betapi0_output.root"
    hist_type=("h2d_ppIM_vs_betapi0_TDATA" "h2d_ppIM_vs_betapi0_TISR3PI_SIG" "h2d_ppIM_vs_betapi0_TETAGAM")
    cv_nm=("data_ppbeta" "signal_ppbeta" "etagam_ppbeta")
    cv_text=("Data" "Signal" "#eta#gamma")
    pt1_x0=(0.5 0.5 0.5)
    pt1_x1=(0.8 0.8 0.8)
    macro="plot2d_ppbeta.C"
    func="plot2d_ppbeta"
else
    echo "Unknown plot type: $plot_type. Use 'sfw2d' or 'ppIM_vs_betapi0'."
    exit 1
fi

# ----------------------------------------------------------------------
# Loop over histograms
# ----------------------------------------------------------------------
IFS=""
for ((i=0; i<${#hist_type[@]}; ++i)); do
    # Update header with current histogram settings
    sed -i 's/\(const double pt1_x0 =\)\(.*\)/\1 '${pt1_x0[i]}';/' $header
    sed -i 's/\(const double pt1_x1 =\)\(.*\)/\1 '${pt1_x1[i]}';/' $header
    sed -i 's/\(const TString hist_type =\)\(.*\)/\1 "'${hist_type[i]}'";/' $header
    sed -i 's|\(const TString infile_nm =\)\(.*\)|\1 "'"${output_path}${outfile_name}"'";|' $header
    sed -i 's|\(const TString output_path =\)\(.*\)|\1 "'"${output_path}"'";|' $header
    sed -i 's/\(const TString cv_nm =\)\(.*\)/\1 "'${cv_nm[i]}'";/' $header
    sed -i 's/\(const TString cv_text =\)\(.*\)/\1 "'${cv_text[i]}'";/' $header

    # Create and run the ROOT macro
    cat > plot_script.C <<EOF
#include <iostream>
void plot_script() {
    gROOT->ProcessLine(".L ../run_bdt/${macro}");
    gROOT->ProcessLine("${func}(\"$plots_dir\")");
}
EOF
    root -l -n -q -b plot_script.C
    rm plot_script.C
done

echo "All $plot_type plots finished."
