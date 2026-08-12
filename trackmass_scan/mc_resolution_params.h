// MC resolution parameters from modified logistic fit
// sigma(M) = low + (high(M)-low)/(1+exp(-k*(M-M0)))
// high(M) = h0 + h1*(M - Mref)  with Mref = 500 MeV
const double SIGMA_LOW  = 0.274517;
const double SIGMA_LOW_ERR  = 0.0912288;
const double H0         = 1.12367;
const double H0_ERR     = 0.0249805;
const double H1         = -0.000881756;
const double H1_ERR     = 3.3147e-05;
const double K          = 0.0794985;
const double K_ERR      = 0.00629665;
const double M0         = 300.014;
const double M0_ERR     = 2.75775;
const double MREF       = 500;
const double CHI2_NDF   = 37.8923;

// Standard logistic plateau (constant) for comparison
const double SIGMA_HIGH_CONST = 1.13652;
const double SIGMA_HIGH_CONST_ERR = 0.00174418;
const double CHI2_NDF_STD = 64.7093;

// Omega width ratio (to be filled after running omega_resolution.C)
// const double R_OMEGA = ...; // sigma_data_omega / sigma_mc_omega
