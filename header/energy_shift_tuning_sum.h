const double energy_shift1 = 3.0679; // iteration 1; massbias_tuning.h
const double energy_shift1_err = 0.0224552;

const double energy_shift2 = 1.94332;
const double energy_shift2_err = 0.0224152;

const double energy_shift3 = 1.23598;
const double energy_shift3_err = 0.0223806;

const double energy_shift4 = 0.787253;
const double energy_shift4_err = 0.0223592;

const double energy_shift5 = 0.505824;
const double energy_shift5_err = 0.0223641;

const double energy_shift6 = 0.323055;
const double energy_shift6_err = 0.022359;

const double energy_shift7 = 0.204365;
const double energy_shift7_err = 0.0223529;

const double energy_shift8 = 0.133591;
const double energy_shift8_err = 0.0223555;

const double energy_shift_total = energy_shift1 + energy_shift2 + energy_shift3 + energy_shift4 + energy_shift5 + energy_shift6 + energy_shift7 + energy_shift8;
const double energy_shift_total_err = TMath::Sqrt(TMath::Power(energy_shift1_err, 2) + TMath::Power(energy_shift2_err, 2) + TMath::Power(energy_shift3_err, 2) + TMath::Power(energy_shift4_err, 2) + TMath::Power(energy_shift5_err, 2) + TMath::Power(energy_shift6_err, 2) + TMath::Power(energy_shift7_err, 2) + TMath::Power(energy_shift8_err, 2));


double alpha_delta = energy_shift_total / 353.36;   // Sum of all iterations
const double MASS_SCALE_PI0 = 1 + alpha_delta;
//const double MASS_SCALE_PI0 = 1.; // Set to determine mass bias from pull tuning

// from tuning_raw.h, data pulls
const double bias_E12 = -0.0794267;
const double sigma_scale_E12 = 1.02669;
const double bias_E3 = -0.0632655;
const double sigma_scale_E3 = 1.03991;
