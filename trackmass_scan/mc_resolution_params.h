// MC resolution parameters from modified logistic fit
// sigma(M) = low + (high(M)-low)/(1+exp(-k*(M-M0)))
// high(M) = h0 + h1*(M - Mref)  with Mref = 500 MeV
const double SIGMA_LOW  = 0.178101;
const double SIGMA_LOW_ERR  = 0.0378855;
const double H0         = 1.07536;
const double H0_ERR     = 0.00303821;
const double H1         = -0.000921395;
const double H1_ERR     = 5.91642e-05;
const double K          = 0.0860827;
const double K_ERR      = 0.00788239;
const double M0         = 300;
const double M0_ERR     = 1.19054;
const double MREF       = 500;
const double CHI2_NDF   = 12.9497;

// Standard logistic plateau (constant) for comparison
const double SIGMA_HIGH_CONST = 1.09344;
const double SIGMA_HIGH_CONST_ERR = 0.0031034;
const double CHI2_NDF_STD = 21.9372;

// Omega width ratio (to be filled after running omega_resolution.C)
// const double R_OMEGA = ...; // sigma_data_omega / sigma_mc_omega
