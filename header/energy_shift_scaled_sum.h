const double energy_shift0 = -0.28989; // iteration 1; 
const double energy_shift0_err = 0.0223559;

const double energy_shift1 = -0.160623; // iteration 2; 
const double energy_shift1_err = 0.022352;

const double energy_shift2 = -0.084715; // iteration 3; 
const double energy_shift2_err = 0.0223685;

const double energy_shift3 = -0.0449349; // iteration 4; 
const double energy_shift3_err = 0.022372;

const double energy_shift_total = energy_shift0 + energy_shift1 + energy_shift2 + energy_shift3;
const double energy_shift_total_err = TMath::Sqrt(TMath::Power(energy_shift0_err, 2) + TMath::Power(energy_shift1_err, 2) + TMath::Power(energy_shift2_err, 2) + + TMath::Power(energy_shift3_err, 2));

double alpha_delta = energy_shift_total / 353.36;   // Sum of all iterations
//const double MASS_SCALE_PI0 = 1 + alpha_delta;
const double MASS_SCALE_PI0 = 1.; // Set to determine mass bias from pull tuning

// from tuning_raw.h, data pulls
const double bias_E12 = -0.0794267;
const double sigma_scale_E12 = 1.02669;
const double bias_E3 = -0.0632655;
const double sigma_scale_E3 = 1.03991;
