// MC resolution parameters from modified logistic fit
// sigma(M) = low + (high(M)-low)/(1+exp(-k*(M-M0)))
// high(M) = h0 + h1*(M - Mref)  with Mref = 500 MeV
const double SIGMA_LOW  = 0.8;
const double SIGMA_LOW_ERR  = 0.128942;
const double H0         = 1.10757;
const double H0_ERR     = 0.0922852;
const double H1         = -0.000709056;
const double H1_ERR     = 0.00011823;
const double K          = 0.182909;
const double K_ERR      = 0.227312;
const double M0         = 314.133;
const double M0_ERR     = 3.42022;
const double MREF       = 500;
const double CHI2_NDF   = 4.2698;

// Standard logistic plateau (constant) for comparison
const double SIGMA_HIGH_CONST = 1.12508;
const double SIGMA_HIGH_CONST_ERR = 0.00634788;
const double CHI2_NDF_STD = 5.56962;

// Omega width ratio (to be filled after running omega_resolution.C)
// const double R_OMEGA = ...; // sigma_data_omega / sigma_mc_omega
