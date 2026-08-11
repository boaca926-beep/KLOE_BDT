// MC resolution parameters from modified logistic fit
// sigma(M) = low + (high(M)-low)/(1+exp(-k*(M-M0)))
// high(M) = h0 + h1*(M - Mref)  with Mref = 500 MeV
const double SIGMA_LOW  = 0.8;
const double SIGMA_LOW_ERR  = 0.504575;
const double H0         = 1.10666;
const double H0_ERR     = 0.0924445;
const double H1         = -0.000738332;
const double H1_ERR     = 0.00011882;
const double K          = 0.17161;
const double K_ERR      = 0.215728;
const double M0         = 314.359;
const double M0_ERR     = 3.55056;
const double MREF       = 500;
const double CHI2_NDF   = 4.16075;

// Standard logistic plateau (constant) for comparison
const double SIGMA_HIGH_CONST = 1.12491;
const double SIGMA_HIGH_CONST_ERR = 0.00636273;
const double CHI2_NDF_STD = 5.59376;

// Omega width ratio (to be filled after running omega_resolution.C)
// const double R_OMEGA = ...; // sigma_data_omega / sigma_mc_omega
