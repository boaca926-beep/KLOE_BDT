const double me = 0.511 * 1e-3; // mass of electron in GeV
const double mass_omega = 782.65; // omega mass in MeV/c^2
const double Gam_omega = 8.49; // omega mass width in MeV
const double BB_omega = 6.38e-5; // branching ratio product

const double Lumi_tot = 1724470;
double Lumi_int = Lumi_tot;

const double crx_scale = 0.3894e6; // 1 [GeV]^{-2}=0.3894 mb=0.3894e6nb

const double pi = TMath::Pi();
const double alpha = 1. / 137; // fine structure constant 
const double alphapi = alpha / pi;
const double sqrtS = 1.019; // in GeV 
const double theta0_degr = 23.;
const double theta0 = theta0_degr * pi / 180.; // in radius

const double frac = 0.2;
const double wideBias = 0.5;
const double wideSigma = 13;
const double smallBias = 0.34;
const double smallSigma = 3.2;

const double a3pi = 0.; //0.1; // [GeV]^{-2}
const double s1 = 0.; //1.
const double gpipi2 = 35.3098; // 37
const double phase_angle = 155. / 180. * pi; // 155 or 163 phi-omega interference angle [degr] -> [rad]

// Masses in MeV/c²
const double mpi0 = 134.9768; // neutral pi mass
const double mpi = 139.57039; // charge pi mass
const double mrho0 = 775.26; // rho0, rho+,  rho- 
const double mK = 493.677; // K+/K- mass, MeV/c²
const double mK0 = 497.611; // K0 mass, MeV/c^2
const double meta = 547.862; // eta mass

const double M_phi = 1019.461; // phi mass
const double M_omega = 782.65; // omega mass

// Width in MeV
const double Gamma_rho0 = 147.4; // rho0
const double Gamma_omega = 8.68; // omega
const double Gamma_phi = 4.249; // phi

// Branching ratio
const double Br_2pi = 1.53e-2; // branching fraction (omega -> pi+ pi-) = 0.0153 + 0.0011 - 0.0013
const double Br_pigam = 8.35e-2; // branching fraction (omega -> pi0 + gamma) = 0.0834 +/- 0.0026
const double Br_3pi_omega = 89.2e-2; // branching fraction (omega -> pi+ pi- pi0) = 0.892 +/- 0.007
const double Br_KpKm = 49.2e-2; // phi -> K+ K- (49.2 +/- 0.5)%
const double Br_KK = 34.0e-2; // phi -> KL KS (34.0 +/- 0.4)%
const double Br_etagam = 1.303e-2; // phi -> eta gamma (1.303 +/- 0.025)%
const double Br_3pi_phi = 15.24e-2; // phi -> rho pi -> 3pi (15.24 +/- 0.33)%
const double BB_phi = 4.51e-5; // branching ratio product B(e+ e- -> phi)xB(phi->pi+ pi- pi0)
//const double BB_phi = 4.35e-5; // branching ratio product B(e+ e- -> phi)xB(phi->pi+ pi- pi0)

// Cross section in nb
const double Sigma_phi = 637.034; // Calculate phi peak cross section
const double Sigma_omega = 1457.; // omega peak cross section nb
const double Sigma_bg = 12.; // back ground cross section nb, refer to https://arxiv.org/pdf/hep-ex/0002017.pdf; 12 +/- 5 nb

