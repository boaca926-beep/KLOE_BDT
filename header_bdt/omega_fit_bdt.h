// modify to omega_fit folder, input_bdt.sh
//const TString mainFolder = "/media/bo/Backup/bdt_output/bdt_raw_TDATA_norm_false"; // kinematic fit
//const TString mainFolder = "/media/bo/Backup/bdt_output/bdt_tuning_TDATA_chain_false_trk"; // track correction
//const TString mainFolder = "/media/bo/Backup/bdt_output/bdt_tuning_TDATA_norm_true_pull"; // pull correction
//const TString mainFolder = "/media/bo/Backup/bdt_output/bdt_tuning_TDATA_norm_true_corr"; // final correction

const TString mainFolder = "/home/bo/Desktop/bdt_tuning_post_fit_TDATA_chain_true"; 
const TString treeFile = mainFolder + "/cut/tree_pre.root";
const TString sfw2dFile = mainFolder + "/sfw2d/sfw2d.root";
const TString output_path = "../plots_m3pi_corr/";

// binning.h
#ifndef BINNING_H
#define BINNING_H
const double MASS_MIN = 760.; //600.0;
const double MASS_MAX = 800.; //900.0;
const double IM3pi_sigma = 1.66;
const double mass_sigma_nb = 0.5;
const int NBINS = TMath::Nint((MASS_MAX - MASS_MIN) / mass_sigma_nb / IM3pi_sigma);

#endif
