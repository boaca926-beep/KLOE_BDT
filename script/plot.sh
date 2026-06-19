#!/bin/bash

# Quick start
# ./plot.sh sfw2d             # for the original sideband plots
# ./plot.sh ppIM_vs_betapi0   # for the new dipion mass vs π⁰ β correlation

sample_type=norm
main_folder="/home/bo/Desktop/input_bdt_TDATA_${sample_type}"
hist_root="${main_folder}/hist/hist.root"

# Usage: ./plot.sh [sfw2d|ppIM_vs_betapi0]
plot_type=${1:-sfw2d}
header=../header_bdt/plot.h
output_path=../run_bdt/
plots_dir=../plots_${plot_type}/

# ----------------------------------------------------------------------
# Create required directories
# ----------------------------------------------------------------------
mkdir -p "$output_path"
mkdir -p "$plots_dir"

echo -e "\nPlotting $plot_type ..."

# Check if histogram file exists
if [ ! -f "$hist_root" ]; then
    echo "ERROR: Histogram file not found: $hist_root"
    exit 1
fi

# Update header with hist_root path
sed -i 's|\(const TString hist_root =\)\(.*\)|\1 "'"${hist_root}"'";|' "$header"

# ----------------------------------------------------------------------
# Get hist.root and produce histos for plotting e.g. combined 2D outputs (sfw2d or ppIM_vs_betapi0)
# ----------------------------------------------------------------------
gethist_script=gethist_script.C
cat > "$gethist_script" <<EOF
#include <iostream>
void gethist_script() {
    gROOT->ProcessLine(".L ../run_bdt/getplothist.C");
    gROOT->ProcessLine("getplothist(\"$plot_type\")");
}
EOF
root -l -n -q -b "$gethist_script"
rm -f "$gethist_script"

# ----------------------------------------------------------------------
# Configure arrays and macro based on plot type
# ----------------------------------------------------------------------
if [ "$plot_type" == "sfw2d" ]; then
    outfile_name="sfw2d_output.root"
    hist_type=("h2d_sfw_TDATA" "h2d_sfw_TISR3PI_SIG_peak" "h2d_sfw_TISR3PI_SIG_non_reson" "hbkgsum_nosig")
    cv_nm=("data" "signal_peak" "signal_non_reson" "bkgsum_nosig")
    cv_text=("Data" "Signal (#omega peak)" "Signal (Non-resonant)" "Others")
    #pt1_x0=(0.3 0.3 0.3 0.3) #0.5
    #pt1_x1=(0.4 0.4 0.4 0.4) #0.8
    macro="plot_sfw.C"
    func="plot_sfw"
elif [ "$plot_type" == "ppIM_vs_betapi0" ]; then
    outfile_name="ppIM_vs_betapi0_output.root"
    # FIXED: removed stray commas, corrected histogram name "hbkgsum_nosig" (was "hbkgrum_nosig")
    hist_type=("h2d_ppIM_vs_betapi0_TDATA" "h2d_ppIM_vs_betapi0_TISR3PI_SIG_peak" "h2d_ppIM_vs_betapi0_TISR3PI_SIG_non_reson" "hbkgsum_nosig")
    cv_nm=("data_ppbeta" "signal_ppbeta_peak" "signal_ppbeta_non_reson" "bkgsum_nosig")
    cv_text=("Data" "Signal (#omega peak)" "Signal (Non-resonant)" "Others")
    #pt1_x0=(0.6 0.6 0.6 0.6)
    #pt1_x1=(0.85 0.85 0.85 0.85)
    macro="plot_ppbeta.C"
    func="plot_ppbeta"
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
    #sed -i 's/\(const double pt1_x0 =\)\(.*\)/\1 '"${pt1_x0[i]}"';/' "$header"
    #sed -i 's/\(const double pt1_x1 =\)\(.*\)/\1 '"${pt1_x1[i]}"';/' "$header"
    sed -i 's/\(const TString hist_type =\)\(.*\)/\1 "'"${hist_type[i]}"'";/' "$header"
    sed -i 's|\(const TString infile_nm =\)\(.*\)|\1 "'"${output_path}${outfile_name}"'";|' "$header"
    sed -i 's|\(const TString output_path =\)\(.*\)|\1 "'"${output_path}"'";|' "$header"
    sed -i 's/\(const TString cv_nm =\)\(.*\)/\1 "'"${cv_nm[i]}"'";/' "$header"
    sed -i 's/\(const TString cv_text =\)\(.*\)/\1 "'"${cv_text[i]}"'";/' "$header"

    # Create and run the ROOT macro
    cat > plot_script.C <<EOF
#include <iostream>
void plot_script() {
    gROOT->ProcessLine(".L ../run_bdt/${macro}");
    gROOT->ProcessLine("${func}(\"$plots_dir\")");
}
EOF
    root -l -n -q -b plot_script.C
    rm -f plot_script.C
done

echo "All $plot_type plots finished."
