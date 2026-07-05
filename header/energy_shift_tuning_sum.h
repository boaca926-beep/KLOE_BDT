const double energy_shift0 = 3.0679; // iteration 1; massbias_tuning.h
const double energy_shift0_err = 0.0224552;

const double energy_shift1 = 1.94332;
const double energy_shift1_err = 0.0224152;

const double energy_shift_total = energy_shift0 + energy_shift1;
const double energy_shift_total_err = TMath::Sqrt(TMath::Power(energy_shift0_err, 2) + TMath::Power(energy_shift1_err, 2));

double alpha_delta = energy_shift_total / 353.36;   // Sum of all iterations
const double MASS_SCALE_PI0 = 1 + alpha_delta;
//const double MASS_SCALE_PI0 = 1.; // Set to determine mass bias from pull tuning

// from tuning_raw.h, data pulls
const double bias_E12 = -0.0794267;
const double sigma_scale_E12 = 1.02669;
const double bias_E3 = -0.0632655;
const double sigma_scale_E3 = 1.03991;
