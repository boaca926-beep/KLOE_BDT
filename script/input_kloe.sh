#!/bin/bash

sample_size=chain # norm; small; mini; chain
sample_path=../path_${sample_size}/ 
exp_type=TDATA # DATA
tuning_type=raw #raw: no tuning; scaled: pi0 photon energy scale
gsf=1 # DATA

result_path=../../input_kloe_${tuning_type}_${exp_type}_${sample_size}
#result_path=/media/bo/Analysis_Disk/

## Initialize the normial conditions
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

c0=0.11
c1=0.8
cut_nm=""
cut_value=0

cut_header=../header/cut_para.h
cat > $cut_header <<EOF
const double chi2_cut = $chi2_cut;
const double angle_cut = $angle_cut;
const double beta_cut = $beta_cut;
const double Eprompt_max_cut = $Eprompt_max_cut;
const double deltaE_min = $deltaE_min;
const double deltaE_max = $deltaE_max;
const double beta_3pi_min = $beta_3pi_min;
const double beta_3pi_max = $beta_3pi_max;
const double c0 = $c0;
const double c1 = $c1;
double cut_value = -1;
const TString cut_nm = "";
EOF


# histo
mass_sigma_nb=1
sfw2d_sigma_nb=1

hist_header=../header/hist.h
sed -i 's/\(const double mass_sigma_nb =\)\(.*\)/\1 '$mass_sigma_nb';/' $hist_header
sed -i 's/\(const double sfw2d_sigma_nb =\)\(.*\)/\1 '$sfw2d_sigma_nb';/' $hist_header

# omega fit
fit_min=760
fit_max=800

omega_header=../header/omega_fit.h
sed -i 's/\(const double fit_min =\)\(.*\)/\1 '$fit_min';/' $omega_header
sed -i 's/\(const double fit_max =\)\(.*\)/\1 '$fit_max';/' $omega_header

# sm_para
Lumi_tot=1724470

sm_header=../header/sm_para.h
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
massbias_path=${result_path}/massbias/
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

echo "Initializing $path_header and log files!"
## Initializing $path_header and log files!"
path_header=../header/path.h
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
omega_fit_script=omega_fit_script.C


log_input=${log_path}log_input.txt
#echo "" > ${log_input}
touch $log_input

log_cut=${log_path}log_cut.txt
#echo "" > ${log_cut}
touch $log_cut

log_hist=${log_path}log_hist.txt
#echo "" > ${log_hist}
touch $log_hist

log_sfw2d=${log_path}log_sfw2d.txt
#echo "" > ${log_sfw2d}
touch $log_sfw2d

log_sfw1d=${log_path}log_sfw1d.txt
#echo "" > ${log_sfw1d}
touch $log_sfw1d

log_omega_fit=${log_path}log_omega_fit.txt
#echo "" > ${log_omega_fit}
touch $log_omega_fit

echo "Looping over data samples ..."
# Loop over data smaples
for ((i=0;i<${#DATA_TYPE[@]};++i)); do

    data_type=${DATA_TYPE[i]}
    #echo ${data_type}

    # Ensure the sampleFile path matches data_type
    sample_file_name="${input_path}${data_type}"
    
    INPUT_FILE=${sample_path}${data_type}_path
    ROOT_FILE=${input_path}${data_type}
    echo $INPUT_FILE
    #echo $ROOT_FILE

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
    echo '  gROOT->ProcessLine(".L ../run/Analys_class.C");' >> $run_script
    echo '  gROOT->ProcessLine("Analys_class(rootFile, sampleFile)");' >> $run_script
    echo '}' >> $run_script
    root -l -n -q -b $run_script >> ${log_input}
    ## Selection cuts
    tree_cut_script=tree_cut_script.C
    echo '#include <iostream>' > $tree_cut_script
    echo "void tree_cut_script() {" >> $tree_cut_script
    
    echo "gROOT->ProcessLine(\".L ../run/tree_cut_${tuning_type}.C\");" >> $tree_cut_script
    echo "gROOT->ProcessLine(\"tree_cut_${tuning_type}()\");" >> $tree_cut_script
    echo '}' >> $tree_cut_script
    root -l -n -q -b $tree_cut_script >> ${log_cut}
done
echo "Selection cuts applied!"

## Signal MC generated
tree_gen_script=tree_gen_script.C
echo '#include <iostream>' > $tree_gen_script
echo "void tree_gen_script() {" >> $tree_gen_script
echo 'gROOT->ProcessLine(".L ../run/tree_gen.C");' >> $tree_gen_script
echo 'gROOT->ProcessLine("tree_gen()");' >> $tree_gen_script
echo '}' >> $tree_gen_script
root -l -n -q -b $tree_gen_script
echo "Signal MC is generated!"

## Histos
echo '#include <iostream>' > $hist_script
echo "void hist_script() {" >> $hist_script
echo 'gROOT->ProcessLine(".L ../run/gethist.C");' >> $hist_script
echo 'gROOT->ProcessLine("gethist()");' >> $hist_script
echo '}' >> $hist_script
root -l -n -q -b $hist_script >> ${log_hist}
echo "Histos are created!"

## Normalization
echo '#include <iostream>' > $sfw2d_script
echo "void sfw2d_script() {" >> $sfw2d_script
echo 'gROOT->ProcessLine(".L ../run/sfw2d.C");' >> $sfw2d_script
echo 'gROOT->ProcessLine("sfw2d()");' >> $sfw2d_script
echo '}' >> $sfw2d_script
root -l -n -q -b $sfw2d_script >> ${log_sfw2d}
#cp ../header/sfw2d.txt ${outputSfw2D}
#ls ${outputSfw2D}
echo "MC normalization!"

## MC signal tuning
echo '#include <iostream>' > $sfw1d_script
echo "void sfw1d_script() {" >> $sfw1d_script
echo 'gROOT->ProcessLine(".L ../run/sfw1d.C");' >> $sfw1d_script
echo 'gROOT->ProcessLine("sfw1d()");' >> $sfw1d_script
echo '}' >> $sfw1d_script
root -l -n -q -b $sfw1d_script >> ${log_sfw1d}
#cp ../header/sfw1d.txt ${outputSfw1D}
#ls ${outputSfw1D}
echo "MC signal tuning!"


## Omega parameters
echo '#include <iostream>' > $omega_fit_script
echo "void omega_fit_script() {" >> $omega_fit_script
echo 'gROOT->ProcessLine(".L ../run/omega_fit.C");' >> $omega_fit_script
echo 'gROOT->ProcessLine("omega_fit()");' >> $omega_fit_script
echo '}' >> $omega_fit_script
root -l -n -q -b $omega_fit_script >> ${log_omega_fit}
echo "Omega parameters are extracted!"

rm $run_script
rm $tree_cut_script
rm $tree_gen_script
rm $hist_script
rm $sfw2d_script
rm $sfw1d_script
rm $omega_fit_script
