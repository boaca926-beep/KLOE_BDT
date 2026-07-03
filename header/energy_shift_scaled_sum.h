const double energy_shift0 = -0.28989; // iteration 0. run tree_cut_raw.C
const double energy_shift0_err = 0.0223559;

const double energy_shift1 = -0.160623; // iteration 1. run tree_cut_scaled.C
const double energy_shift1_err = 0.022352;

const double energy_shift2 = -0.084715; // iteration 2. run tree_cut_scaled.C
const double energy_shift2_err = 0.0223685;

const double energy_shift_total = energy_shift0 + energy_shift1 + energy_shift2;
const double energy_shift_total_err = TMath::Sqrt(TMath::Power(energy_shift0_err, 2) + TMath::Power(energy_shift1_err, 2) + TMath::Power(energy_shift2_err, 2));

// from tuning.h
const double bias_E12 = -0.113633;
const double sigma_scale_E12 = 1.01361;
const double bias_E3 = -0.162934;
const double sigma_scale_E3 = 1.03341;
