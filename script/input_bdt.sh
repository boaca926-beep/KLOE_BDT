#!/bin/bash
set -e   # exit immediately if any command fails

# ============================================================
# Configuration
# ============================================================
sample_size=norm          # norm; small; mini; chain
sample_path=../path_${sample_size}/
exp_type=TDATA             # DATA
tuning_type=tuning         # raw: no tuning; tuning: tuned + scale
gsf=1                      # DATA
pull_status=false          # main switch: false = no corrections, true = apply corrections

# ------------------- Correction flags -------------------
# Set these to true/false to independently enable/disable each correction
APPLY_PULL=true            # apply bias shift (mean) + scale ratio (width)
APPLY_MASS_SCALE=true      # apply mass scale (MASS_SCALE_PI0)
# ---------------------------------------------------------

## Initialize tuning corrections conditions
# step1: set pull_status to false to get pull correction parameters, input for pull_scan (run only once!)
# step2: set pull_status to true to apply pull corrections (always on for analysis)

SOURCE_BIAS="../pull_scan/bias_shift.txt"
SOURCE_SCALE="../pull_scan/scale_ratio.txt"   # adjust if your file is scale_ratio_corr.txt
TARGET_FILE="../header_bdt/energy_shift_tuning_sum.h"

# ============================================================
# Update header based on pull_status and individual flags
# ============================================================
if [ "$pull_status" = true ]; then
    echo "Pull corrections are applied!"

    # ---------- PULL (BIAS + SCALE) ----------
    if [ "$APPLY_PULL" = true ]; then
	echo "Energy Pull Tuning is Applied!"

        # ---- BIAS ----
        if [[ -f "$SOURCE_BIAS" ]]; then
            bias_shift=$(grep -oP '(?:const\s+double\s+)?bias_shift\s*=\s*\K[0-9.eE+-]+' "$SOURCE_BIAS" | head -1)
            bias_shift_err=$(grep -oP '(?:const\s+double\s+)?bias_shift_err\s*=\s*\K[0-9.eE+-]+' "$SOURCE_BIAS" | head -1)
        else
            echo "Warning: $SOURCE_BIAS not found. Using defaults."
            bias_shift=0.0
            bias_shift_err=0.0
        fi
        if [[ -n "$bias_shift" && -n "$bias_shift_err" ]]; then
            echo "Updating $TARGET_FILE with bias_shift=$bias_shift, bias_shift_err=$bias_shift_err"
            sed -i "s/\(const double bias_shift\s*=\s*\)[0-9.eE+-]*;/\1$bias_shift;/" "$TARGET_FILE"
            sed -i "s/\(const double bias_shift_err\s*=\s*\)[0-9.eE+-]*;/\1$bias_shift_err;/" "$TARGET_FILE"
        else
            echo "Error: Could not extract bias shift values from $SOURCE_BIAS"
        fi

        # ---- SCALE ----
        if [[ -f "$SOURCE_SCALE" ]]; then
            scale_ratio=$(grep -oP '(?:const\s+double\s+)?scale_ratio\s*=\s*\K[0-9.eE+-]+' "$SOURCE_SCALE" | head -1)
            scale_ratio_err=$(grep -oP '(?:const\s+double\s+)?scale_ratio_err\s*=\s*\K[0-9.eE+-]+' "$SOURCE_SCALE" | head -1)
        else
            echo "Warning: $SOURCE_SCALE not found. Using defaults."
            scale_ratio=1.0
            scale_ratio_err=0.0
        fi
        if [[ -n "$scale_ratio" && -n "$scale_ratio_err" ]]; then
            echo "Updating $TARGET_FILE with scale_ratio=$scale_ratio, scale_ratio_err=$scale_ratio_err"
            sed -i "s/\(const double scale_ratio\s*=\s*\)[0-9.eE+-]*;/\1$scale_ratio;/" "$TARGET_FILE"
            sed -i "s/\(const double scale_ratio_err\s*=\s*\)[0-9.eE+-]*;/\1$scale_ratio_err;/" "$TARGET_FILE"
        else
            echo "Error: Could not extract scale ratio values from $SOURCE_SCALE"
        fi
    else
        echo "Pull correction disabled. Setting bias_shift=0, scale_ratio=1"
        sed -i "s/\(const double bias_shift\s*=\s*\)[0-9.eE+-]*;/\10.0;/" "$TARGET_FILE"
        sed -i "s/\(const double bias_shift_err\s*=\s*\)[0-9.eE+-]*;/\10.0;/" "$TARGET_FILE"
        sed -i "s/\(const double scale_ratio\s*=\s*\)[0-9.eE+-]*;/\11.0;/" "$TARGET_FILE"
        sed -i "s/\(const double scale_ratio_err\s*=\s*\)[0-9.eE+-]*;/\10.0;/" "$TARGET_FILE"
    fi

    # ---------- MASS SCALE ----------
    if [ "$APPLY_MASS_SCALE" = true ]; then
	echo "Energy Scaling is Applied!"

        if [[ -f "../pull_scan/massbias_bdt.txt" ]]; then
            mpi0_data=$(grep -oP 'mpi0_data\s*=\s*\K[0-9.eE+-]+' "../pull_scan/massbias_bdt.txt" | head -1)
            mpi0_data_err=$(grep -oP 'mpi0_data_err\s*=\s*\K[0-9.eE+-]+' "../pull_scan/massbias_bdt.txt" | head -1)
            mpi0_mc=$(grep -oP 'mpi0_mc\s*=\s*\K[0-9.eE+-]+' "../pull_scan/massbias_bdt.txt" | head -1)
            mpi0_mc_err=$(grep -oP 'mpi0_mc_err\s*=\s*\K[0-9.eE+-]+' "../pull_scan/massbias_bdt.txt" | head -1)
            sed -i "s/\(const double mpi0_data\s*=\s*\)[0-9.eE+-]*;/\1$mpi0_data;/" "$TARGET_FILE"
            sed -i "s/\(const double mpi0_data_err\s*=\s*\)[0-9.eE+-]*;/\1$mpi0_data_err;/" "$TARGET_FILE"
            sed -i "s/\(const double mpi0_mc\s*=\s*\)[0-9.eE+-]*;/\1$mpi0_mc;/" "$TARGET_FILE"
            sed -i "s/\(const double mpi0_mc_err\s*=\s*\)[0-9.eE+-]*;/\1$mpi0_mc_err;/" "$TARGET_FILE"
            if [[ -n "$mpi0_data" && -n "$mpi0_mc" ]]; then
                mass_scale=$(echo "scale=8; $mpi0_data / $mpi0_mc" | bc)
                echo "Updating MASS_SCALE_PI0 = $mass_scale"
                sed -i "s/\(const double MASS_SCALE_PI0\s*=\s*\)[0-9.eE+-]*;/\1$mass_scale;/" "$TARGET_FILE"
            fi
        else
            echo "Warning: ../pull_scan/massbias_bdt.txt not found. MASS_SCALE_PI0 unchanged."
        fi
    else
        echo "Mass scale correction disabled. Setting MASS_SCALE_PI0=1"
        sed -i "s/\(const double MASS_SCALE_PI0\s*=\s*\)[0-9.eE+-]*;/\11.0;/" "$TARGET_FILE"
    fi

    echo "Pull correction parameters updated."

else
    echo "No pull corrections are applied!"

    # Reset all corrections to neutral values
    echo "Setting bias_shift=0, scale_ratio=1, and MASS_SCALE_PI0=1"
    sed -i "s/\(const double bias_shift\s*=\s*\)[0-9.eE+-]*;/\10.0;/" "$TARGET_FILE"
    sed -i "s/\(const double bias_shift_err\s*=\s*\)[0-9.eE+-]*;/\10.0;/" "$TARGET_FILE"
    sed -i "s/\(const double scale_ratio\s*=\s*\)[0-9.eE+-]*;/\11.0;/" "$TARGET_FILE"
    sed -i "s/\(const double scale_ratio_err\s*=\s*\)[0-9.eE+-]*;/\10.0;/" "$TARGET_FILE"
    sed -i "s/\(const double mpi0_mc\s*=\s*\)[0-9.eE+-]*;/\10.0;/" "$TARGET_FILE"
    sed -i "s/\(const double mpi0_mc_err\s*=\s*\)[0-9.eE+-]*;/\10.0;/" "$TARGET_FILE"
    sed -i "s/\(const double mpi0_data\s*=\s*\)[0-9.eE+-]*;/\10.0;/" "$TARGET_FILE"
    sed -i "s/\(const double mpi0_data_err\s*=\s*\)[0-9.eE+-]*;/\10.0;/" "$TARGET_FILE"
    sed -i "s/\(const double MASS_SCALE_PI0\s*=\s*\)[0-9.eE+-]*;/\11.0;/" "$TARGET_FILE"
fi

# ============================================================
# Result path and directories (unchanged)
# ============================================================
result_path=../../bdt_${tuning_type}_${exp_type}_${sample_size}_${pull_status}
#result_path=/media/bo/Analysis_Disk/


## Initialize the normal conditions
# Pre-selection
egammamin=15 # 15, 20
Rhovmax=4
Zvmax=10
nb_sigma_T_clust=3 # 3, 4

class_header=../header/MyClass.h
sed -i 's/\(egammamin =\)\(.*\)/\1 '$egammamin';/' $class_header
sed -i 's/\(Rhovmax =\)\(.*\)/\1 '$Rhovmax';/' $class_header
sed -i 's/\(Zvmax =\)\(.*\)/\1 '$Zvmax';/' $class_header
sed -i 's/\(nb_sigma_T_clust =\)\(.*\)/\1 '$nb_sigma_T_clust';/' $class_header

# Selection cuts
Eprompt_max_cut=300
chi2_cut=20 #43 20
angle_cut=138 #138 66
beta_cut=1.98
deltaE_min=-440
deltaE_max=-240 #-150

beta_3pi_min=0.23
beta_3pi_max=0.28

bdt_cut=0.4
c0=0.11
c1=0.8
cut_nm=""
cut_value=0

cut_header=../header_bdt/cut_para.h
cat > $cut_header <<EOF
const double chi2_cut = $chi2_cut;
const double angle_cut = $angle_cut;
const double beta_cut = $beta_cut;
const double Eprompt_max_cut = $Eprompt_max_cut;
const double deltaE_min = $deltaE_min;
const double deltaE_max = $deltaE_max;
const double beta_3pi_min = $beta_3pi_min;
const double beta_3pi_max = $beta_3pi_max;
const double bdt_cut = $bdt_cut;
const double c0 = $c0;
const double c1 = $c1;
double cut_value = -1;
const TString cut_nm = "";
EOF

# histo
mass_sigma_nb=1
sfw2d_sigma_nb=1

hist_header=../header_bdt/hist.h
sed -i 's/\(const double mass_sigma_nb =\)\(.*\)/\1 '$mass_sigma_nb';/' $hist_header
sed -i 's/\(const double sfw2d_sigma_nb =\)\(.*\)/\1 '$sfw2d_sigma_nb';/' $hist_header

# omega fit
fit_min=760
fit_max=800

omega_header=../header_bdt/omega_fit.h
sed -i 's/\(const double fit_min =\)\(.*\)/\1 '$fit_min';/' $omega_header
sed -i 's/\(const double fit_max =\)\(.*\)/\1 '$fit_max';/' $omega_header

# sm_para
Lumi_tot=1724470

sm_header=../header_bdt/sm_para.h
sed -i 's/\(const double Lumi_tot =\)\(.*\)/\1 '$Lumi_tot';/' $sm_header

## Samples
DATA_TYPE=("sig" "ksl" "exp" "eeg" "ufo")
#DATA_TYPE=("sig")

## Folders
input_path=${result_path}/input/
cut_path=${result_path}/cut/
gen_path=${result_path}/gen/
hist_path=${result_path}/hist/
sfw2d_path=${result_path}/sfw2d/
sfw1d_path=${result_path}/sfw1d/
omega_path=${result_path}/omega_fit/
log_path=${result_path}/log/

if [[ -d "$result_path" ]]; then
    echo "${result_path} ..."
    rm -rf $result_path
fi

mkdir ${result_path} # result folder
mkdir ${input_path} # input root files
mkdir ${cut_path} # trees: after all cuts
mkdir ${gen_path} # tree: signal MC generated
mkdir ${hist_path} # histos
mkdir ${sfw2d_path} # mc normalization
mkdir ${sfw1d_path} # mc signal tuning
mkdir ${omega_path} # omega parameters
mkdir ${log_path} # log files
echo "Results folder is created at ${result_path}"

## Initializing path_header and log files
path_header=../header_bdt/path.h
echo "Initializing $path_header and log files!"

echo -e 'const TString rootFile = "";' > $path_header
echo -e 'const TString sampleFile = "";' >> $path_header
echo -e 'const TString outputCut = "";' >> $path_header
echo -e 'const TString sig_path = "";' >> $path_header
echo -e 'const TString outputGen = "";' >> $path_header
echo -e 'const TString outputHist = "";' >> $path_header
echo -e 'const TString outputSfw2D = "";' >> $path_header
echo -e 'const TString outputSfw1D = "";' >> $path_header
echo -e 'const TString outputOmega = "";' >> $path_header
echo -e 'const TString data_type = "";' >> $path_header
echo -e 'const TString exp_type = "'$exp_type'";' >> $path_header
echo -e 'const TString tuning_type = "'$tuning_type'";' >> $path_header
echo -e "double gsf = $gsf;" >> $path_header

sed -i 's|\(const TString sig_path =\)\(.*\)|\1 "'"${input_path}"'";|' "$path_header"
sed -i 's|\(const TString outputGen =\)\(.*\)|\1 "'"${gen_path}"'";|' "$path_header"
sed -i 's|\(const TString outputHist =\)\(.*\)|\1 "'"${hist_path}"'";|' "$path_header"
hist_script=hist_script.C
sed -i 's|\(const TString outputSfw2D =\)\(.*\)|\1 "'"${sfw2d_path}"'";|' "$path_header"
sfw2d_script=sfw2d_script.C
sed -i 's|\(const TString outputSfw1D =\)\(.*\)|\1 "'"${sfw1d_path}"'";|' "$path_header"
sfw1d_script=sfw1d_script.C
sed -i 's|\(const TString outputOmega =\)\(.*\)|\1 "'"${omega_path}"'";|' "$path_header"

log_input=${log_path}log_input.txt
touch $log_input

log_cut=${log_path}log_cut.txt
touch $log_cut

log_hist=${log_path}log_hist.txt
touch $log_hist

log_sfw2d=${log_path}log_sfw2d.txt
touch $log_sfw2d

log_sfw1d=${log_path}log_sfw1d.txt
touch $log_sfw1d

log_omega_fit=${log_path}log_omega_fit.txt
touch $log_omega_fit

echo "Looping over data samples ..."
# Loop over data samples
for ((i=0;i<${#DATA_TYPE[@]};++i)); do

    data_type=${DATA_TYPE[i]}

    # Ensure the sampleFile path matches data_type
    sample_file_name="${input_path}${data_type}"

    INPUT_FILE=${sample_path}${data_type}_path
    ROOT_FILE=${input_path}${data_type}
    echo $INPUT_FILE

    cat > "$path_header" <<EOF
const TString rootFile = "${INPUT_FILE}";
const TString sampleFile = "${ROOT_FILE}";
const TString outputCut = "${cut_path}";
const TString sig_path = "${input_path}";
const TString outputGen = "${gen_path}";
const TString outputHist = "${hist_path}";
const TString outputSfw2D = "${sfw2d_path}";
const TString outputSfw1D = "${sfw1d_path}";
const TString outputOmega = "${omega_path}";
const TString data_type = "${data_type}";
const TString exp_type = "${exp_type}";
const TString tuning_type = "${tuning_type}";
double gsf = ${gsf};
EOF

    # Verify consistency
    if ! grep -q "data_type = \"${data_type}\"" "$path_header"; then
        echo "ERROR: path.h data_type mismatch!" >> ${log_cut}
        exit 1
    fi

    ## Input trees
    run_script=run_script.C

    echo "void run_script() {" > $run_script
    echo '  gROOT->ProcessLine(".L ../run/MyClass.C");' >> $run_script
    echo '  gROOT->ProcessLine(".L ../run_bdt/Analys_class.C");' >> $run_script
    echo '  gROOT->ProcessLine("Analys_class(rootFile, sampleFile)");' >> $run_script
    echo '}' >> $run_script
    root -l -n -q -b $run_script >> ${log_input} 2>&1 || { echo "ROOT failed at run_script for $data_type"; exit 1; }

    ## Selection cuts
    tree_cut_script=tree_cut_script.C
    echo '#include <iostream>' > $tree_cut_script
    echo "void tree_cut_script() {" >> $tree_cut_script
    echo "gROOT->ProcessLine(\".L ../run_bdt/tree_cut_bdt_${tuning_type}.C\");" >> $tree_cut_script
    echo "gROOT->ProcessLine(\"tree_cut_bdt_${tuning_type}()\");" >> $tree_cut_script
    echo '}' >> $tree_cut_script
    root -l -n -q -b $tree_cut_script >> ${log_cut} 2>&1 || { echo "ROOT failed at tree_cut_script for $data_type"; exit 1; }
done
echo "Selection cuts applied!"

# ----------------------------------------------------------------------
# Reset path.h for aggregated steps (no data_type)
# ----------------------------------------------------------------------
cat > "$path_header" <<EOF
const TString rootFile = "${INPUT_FILE}";
const TString sampleFile = "${ROOT_FILE}";
const TString outputCut = "${cut_path}";
const TString sig_path = "${input_path}";
const TString outputGen = "${gen_path}";
const TString outputHist = "${hist_path}";
const TString outputSfw2D = "${sfw2d_path}";
const TString outputSfw1D = "${sfw1d_path}";
const TString outputOmega = "${omega_path}";
const TString data_type = "";
const TString exp_type = "${exp_type}";
double gsf = ${gsf};
EOF

## Signal MC generated
tree_gen_script=tree_gen_script.C
echo '#include <iostream>' > $tree_gen_script
echo "void tree_gen_script() {" >> $tree_gen_script
echo 'gROOT->ProcessLine(".L ../run_bdt/tree_gen.C");' >> $tree_gen_script
echo 'gROOT->ProcessLine("tree_gen()");' >> $tree_gen_script
echo '}' >> $tree_gen_script
root -l -n -q -b $tree_gen_script 2>&1 || { echo "ROOT failed at tree_gen_script"; exit 1; }
echo "Signal MC is generated!"

## Histos
echo '#include <iostream>' > $hist_script
echo "void hist_script() {" >> $hist_script
echo 'gROOT->ProcessLine(".L ../run_bdt/gethist.C");' >> $hist_script
echo 'gROOT->ProcessLine("gethist()");' >> $hist_script
echo '}' >> $hist_script
root -l -n -q -b $hist_script >> ${log_hist} 2>&1 || { echo "ROOT failed at hist_script"; exit 1; }
echo "Histos are created!"

## Normalization
echo '#include <iostream>' > $sfw2d_script
echo "void sfw2d_script() {" >> $sfw2d_script
#echo 'gROOT->ProcessLine(".L ../run_bdt/sfw2d.C");' >> $sfw2d_script
#echo 'gROOT->ProcessLine("sfw2d()");' >> $sfw2d_script

echo 'gROOT->ProcessLine(".L ../run_bdt/sfw2d_noeta.C");' >> $sfw2d_script
echo 'gROOT->ProcessLine("sfw2d_noeta()");' >> $sfw2d_script

#echo 'gROOT->ProcessLine(".L ../run_bdt/sfw2d_merged.C");' >> $sfw2d_script
#echo 'gROOT->ProcessLine("sfw2d_merged()");' >> $sfw2d_script
echo '}' >> $sfw2d_script
root -l -n -q -b $sfw2d_script >> ${log_sfw2d} 2>&1 || { echo "ROOT failed at sfw2d_script"; exit 1; }
echo "MC normalization!"

## MC signal tuning (ω peak correction)
echo '#include <iostream>' > $sfw1d_script
echo "void sfw1d_script() {" >> $sfw1d_script
echo 'gROOT->ProcessLine(".L ../run_bdt/sfw1d.C");' >> $sfw1d_script
echo 'gROOT->ProcessLine("sfw1d()");' >> $sfw1d_script
echo '}' >> $sfw1d_script
#root -l -n -q -b $sfw1d_script >> ${log_sfw1d} 2>&1 || { echo "ROOT failed at sfw1d_script"; exit 1; }
echo "omega peak correction!"

## Omega parameters
omega_fit_script=omega_fit_script.C
echo '#include <iostream>' > $omega_fit_script
echo "void omega_fit_script() {" >> $omega_fit_script
echo 'gROOT->ProcessLine(".L ../run_bdt/omega_fit.C");' >> $omega_fit_script
echo 'gROOT->ProcessLine("omega_fit()");' >> $omega_fit_script
echo '}' >> $omega_fit_script
#root -l -n -q -b $omega_fit_script >> ${log_omega_fit}
echo "Omega parameters are extracted!"

# Optional: keep temporary scripts for debugging
if [ -z "$DEBUG" ]; then
    rm -f $run_script $tree_cut_script $tree_gen_script $hist_script $sfw2d_script $sfw1d_script $omega_fit_script
else
    echo "Debug mode: temporary scripts retained in current directory."
fi
