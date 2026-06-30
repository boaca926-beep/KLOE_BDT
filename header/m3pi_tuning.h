const TString mainFolder = "/home/bo/Desktop/input_kloe_TDATA_chain";
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
