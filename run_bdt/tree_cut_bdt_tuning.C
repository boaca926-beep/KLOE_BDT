#include "../header_bdt/cut_para.h"
#include "../trackmass_scan/mc_bias_params.h"
#include "../trackmass_scan/mc_resolution_params.h"
#include "../trackmass_scan/omega_params.h"
#include "../header_bdt/energy_shift_tuning_sum.h" 
//#include "../header_bdt/tuning.h"
#include "../header_bdt/sm_para.h"
#include "../header_bdt/path.h"
#include "../header_bdt/method.h"
#include "../header_bdt/resol.h"
#include <TStopwatch.h>
#include <TMVA/RBDT.hxx>
#include <TMVA/RTensor.hxx>
#include <TMath.h>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace TMVA::Experimental;

// ----------------------------------------------------------------------
// Configuration
// ----------------------------------------------------------------------
#define BDT_MODEL_PATH "/home/bo/Desktop/KLOE_BDT/models/bdt_pi0_TCOMB.root"

constexpr double ENERGY_THRESHOLD = 5.0;   // MeV
const double M_PION = 139.57061;           // MeV/c^2

// ----------------------------------------------------------------------
// Structures and helper prototypes
// ----------------------------------------------------------------------
struct EventData {
  double photons[3][4];
  double tracks[2][4];
  double lagvalue_min_7C;
  double deltaE;
  double betapi0;
  double angle_pi0gam12;
  double ppIM;
  int bkg_indx;
  int recon_indx;
};

struct BDTResult {
  double max_score;
  double mean_scores;
  double scores[3];
  double score;
  int best_pair_index;
  int pi0_indices[2];
  int prompt_index;
  bool is_valid;
};

// Helper function prototypes
double compute_invariant_mass(int i, int j, const double photons[3][4]);
double compute_3pi_mass(int pi0_idx1, int pi0_idx2, const double photons[3][4], const double tracks[2][4]);
double compute_dipion_mass(const double tracks[2][4]);
double compute_cos_theta(int i, int j, const double photons[3][4]);
double compute_invariant_mass_from_values(double e1, double e2,
                                          double px1, double py1, double pz1,
                                          double px2, double py2, double pz2);
std::vector<float> extract_features(int i_idx, int j_idx, int unpaired_idx,
                                    const double photons[3][4], double energy_threshold);
BDTResult find_best_pion_pair(const EventData& event, TMVA::Experimental::RBDT& bdt);

// ----------------------------------------------------------------------
// Helper to recompute deltaE (same formula as in MyClass.C)
// ----------------------------------------------------------------------
double compute_deltaE(const TLorentzVector &ppl, const TLorentzVector &pmi,
                      const TLorentzVector &beam, double m_pion) {
  TVector3 boost = -beam.BoostVector();
  TLorentzVector ppl_b = ppl; ppl_b.Boost(boost);
  TLorentzVector pmi_b = pmi; pmi_b.Boost(boost);
  double ppl_mag = ppl_b.Vect().Mag();
  double pmi_mag = pmi_b.Vect().Mag();
  double ppl_E = sqrt(m_pion*m_pion + ppl_mag*ppl_mag);
  double pmi_E = sqrt(m_pion*m_pion + pmi_mag*pmi_mag);
  double delta = (ppl_b.Vect() + pmi_b.Vect()).Mag() - (beam.M() - ppl_E - pmi_E);
  return delta;
}

// ----------------------------------------------------------------------
int tree_cut_bdt_tuning() {
  TStopwatch timer;
  timer.Start();
  
  cout << "Input path: " << sampleFile << endl;
  
  TFile *f_input = new TFile(sampleFile + ".root");
  TTree *ALLCHAIN_CUT = (TTree*)f_input->Get("ALLCHAIN_CUT");
  if (!ALLCHAIN_CUT) {
    cerr << "ERROR: Cannot find tree ALLCHAIN_CUT" << endl;
    return 1;
  }
  
  // Load BDT model
  if (gSystem->AccessPathName(BDT_MODEL_PATH)) {
    cerr << "ERROR: BDT model file not found: " << BDT_MODEL_PATH << endl;
    return 1;
  }
  TMVA::Experimental::RBDT bdt("BDT_pi0", BDT_MODEL_PATH);
  cout << "✓ BDT model loaded from " << BDT_MODEL_PATH << endl;
  
  gRandom->SetSeed(12345);
  
  // ---------- Variables ----------
  double lagvalue_min_7C = 0., deltaE = 0., betapi0 = 0., angle_pi0gam12 = 0.;
  double m02 = 0., mplus2 = 0.;
  double m3pi = 0.;
  double ppIM = 0., ppIM_true = 0.;
  double IM3pi_7C = 0., IM3pi_true = 0.;
  double IM_pi0_7C = 0., IM_pi0_nofit = 0.;
  double Eisr = 0., Epi0_pho1 = 0., Epi0_pho2 = 0.;
  double pull_E1 = 0., pull_x1 = 0., pull_y1 = 0., pull_z1 = 0., pull_t1 = 0.;
  double pull_E2 = 0., pull_x2 = 0., pull_y2 = 0., pull_z2 = 0., pull_t2 = 0.;
  double pull_E3 = 0., pull_x3 = 0., pull_y3 = 0., pull_z3 = 0., pull_t3 = 0.;
  double ppl_E = 0., ppl_px = 0., ppl_py = 0., ppl_pz = 0.;
  double pmi_E = 0., pmi_px = 0., pmi_py = 0., pmi_pz = 0.;
  double pho_E1 = 0., pho_px1 = 0., pho_py1 = 0., pho_pz1 = 0.;
  double pho_E2 = 0., pho_px2 = 0., pho_py2 = 0., pho_pz2 = 0.;
  double pho_E3 = 0., pho_px3 = 0., pho_py3 = 0., pho_pz3 = 0.;
  
  // Unfitted 4‑vectors in fit order (read from input)
  double pho_E1_nofit = 0., pho_px1_nofit = 0., pho_py1_nofit = 0., pho_pz1_nofit = 0.;
  double pho_E2_nofit = 0., pho_px2_nofit = 0., pho_py2_nofit = 0., pho_pz2_nofit = 0.;
  double pho_E3_nofit = 0., pho_px3_nofit = 0., pho_py3_nofit = 0., pho_pz3_nofit = 0.;
  
  double ppl_E_true = 0., ppl_px_true = 0., ppl_py_true = 0., ppl_pz_true = 0.;
  double pmi_E_true = 0., pmi_px_true = 0., pmi_py_true = 0., pmi_pz_true = 0.;
  double pho_E1_true = 0., pho_px1_true = 0., pho_py1_true = 0., pho_pz1_true = 0.;
  double pho_E2_true = 0., pho_px2_true = 0., pho_py2_true = 0., pho_pz2_true = 0.;
  double pho_E3_true = 0., pho_px3_true = 0., pho_py3_true = 0., pho_pz3_true = 0.;
  int phid = 0, sig_type = 0;
  int bkg_indx = 0, recon_indx = 0;
  double evnt_tot = 0;
  double Eprompt_max = 0.;
  
  int pho_indx[3], EPI0NTMC[4];
  
  // Arrays for GetLeaf (DO NOT use SetBranchAddress for these!)
  double ENERGYLIST[4];
  double MASSLIST[12];
  double ANGLELIST[1];
  
  // BDT‑specific variables (reco)
  double bdt_score = 0.;
  double e1_bdt = 0., e2_bdt = 0., e3_bdt = 0.;
  double px1_bdt = 0., py1_bdt = 0., pz1_bdt = 0.;
  double px2_bdt = 0., py2_bdt = 0., pz2_bdt = 0.;
  double px3_bdt = 0., py3_bdt = 0., pz3_bdt = 0.;
  
  // BDT‑selected unfitted values (saved to output)
  double e1_nofit_bdt = 0., e2_nofit_bdt = 0., e3_nofit_bdt = 0.;
  double px1_nofit_bdt = 0., py1_nofit_bdt = 0., pz1_nofit_bdt = 0.;
  double px2_nofit_bdt = 0., py2_nofit_bdt = 0., pz2_nofit_bdt = 0.;
  double px3_nofit_bdt = 0., py3_nofit_bdt = 0., pz3_nofit_bdt = 0.;
  double sigma_denom_e1 = 0., sigma_denom_e2 = 0., sigma_denom_e3 = 0.;   
  
  // BDT-selected energy sigmas from SIGMA_FIT_LIST
  double sigma_e1_bdt = 0., sigma_e2_bdt = 0., sigma_e3_bdt = 0.;
  
  double e1_fit = 0., e2_fit = 0., e3_fit = 0.;
  double e1_bdt_true = 0., e2_bdt_true = 0., e3_bdt_true = 0.;
  double m_gg_bdt = 0., m3pi_bdt = 0., m_gg_true = 0.;
  double m2pi_true = 0., m3pi_true = 0.;
  double angle_pi0gam12_bdt = 0., betapi0_bdt = 0.;
  double angle_pi0gam12_bdt_true = 0., betapi0_bdt_true = 0.;
  double true_px_piplus, true_py_piplus, true_pz_piplus, true_E_piplus;
  double true_px_piminus, true_py_piminus, true_pz_piminus, true_E_piminus;
  double true_px_pi0, true_py_pi0, true_pz_pi0, true_E_pi0;
  double true_m3pi;
  double angle_ppl_pmi, angle_trk_neutral;
  double e_asym;
  double beta_3pi;
  
  // Correction
  double delta_e1_bdt = 0., delta_e2_bdt = 0.;
  
  // Pull variables
  double e1_pull_bdt = 0., e2_pull_bdt = 0., e3_pull_bdt = 0.;
  double e1_pull_bdt_new = 0., e2_pull_bdt_new = 0.;
  double e1_pull = 0., e2_pull = 0., e3_pull = 0.;
  double px1_pull = 0., py1_pull = 0., pz1_pull = 0.;
  double px2_pull = 0., py2_pull = 0., pz2_pull = 0.;
  double px3_pull = 0., py3_pull = 0., pz3_pull = 0.;
  double m_gg_pull = 0., m3pi_pull = 0.;
  double m2pi_pull = 0.;
  
  int recon_indx_bdt = 0;
  int isr_recon_quality = 0;
  int total_recon_quality = 0;
  
  // ---- Mapping indices from fit order to original cluster order ----
  int isrgam_indx = 0, pi0gam1_indx = 0, pi0gam2_indx = 0;
  
  // ---------- Beam variables ----------
  float bpx = 0., bpy = 0., bpz = 0., bene = 0.;   // <-- ADDED
  
  // ---------- Output trees ----------
  const int list_size = 13;
  const TString TNM[list_size] = {"TDATA", "TOMEGAPI", "TKPM", "TKSL", "T3PIGAM", "TRHOPI", "TETAGAM", "TBKGREST", "TUFO", "TEEG", "TISR3PI_SIG", "TISR3PI_SIG_PEAK", "TISR3PI_SIG_NON_RESON"};
  TTree *TTList[list_size];
  TCollection* tree_list = new TList;
  
  for (int i = 0; i < list_size; i++) {
    TTList[i] = new TTree(TNM[i], "recreate");
    TTList[i]->SetAutoSave(0);
    tree_list->Add(TTList[i]);
  }
  
  // Add branches to all trees
  TObject* treeout = 0;
  TIter treeliter(tree_list);
  while ((treeout = treeliter.Next()) != 0) {
    TTree* tree_tmp = dynamic_cast<TTree*>(treeout);
    tree_tmp->Branch("Br_ppl_E", &ppl_E, "Br_ppl_E/D");
    tree_tmp->Branch("Br_ppl_px", &ppl_px, "Br_ppl_px/D");
    tree_tmp->Branch("Br_ppl_py", &ppl_py, "Br_ppl_py/D");
    tree_tmp->Branch("Br_ppl_pz", &ppl_pz, "Br_ppl_pz/D");
    tree_tmp->Branch("Br_pmi_E", &pmi_E, "Br_pmi_E/D");
    tree_tmp->Branch("Br_pmi_px", &pmi_px, "Br_pmi_px/D");
    tree_tmp->Branch("Br_pmi_py", &pmi_py, "Br_pmi_py/D");
    tree_tmp->Branch("Br_pmi_pz", &pmi_pz, "Br_pmi_pz/D");
    tree_tmp->Branch("Br_ppl_E_true", &ppl_E_true, "Br_ppl_E_true/D");
    tree_tmp->Branch("Br_ppl_px_true", &ppl_px_true, "Br_ppl_px_true/D");
    tree_tmp->Branch("Br_ppl_py_true", &ppl_py_true, "Br_ppl_py_true/D");
    tree_tmp->Branch("Br_ppl_pz_true", &ppl_pz_true, "Br_ppl_pz_true/D");
    tree_tmp->Branch("Br_pmi_E_true", &pmi_E_true, "Br_pmi_E_true/D");
    tree_tmp->Branch("Br_pmi_px_true", &pmi_px_true, "Br_pmi_px_true/D");
    tree_tmp->Branch("Br_pmi_py_true", &pmi_py_true, "Br_pmi_py_true/D");
    tree_tmp->Branch("Br_pmi_pz_true", &pmi_pz_true, "Br_pmi_pz_true/D");
    
    tree_tmp->Branch("Br_E1", &pho_E1, "Br_pho_E1/D");
    tree_tmp->Branch("Br_px1", &pho_px1, "Br_pho_px1/D");
    tree_tmp->Branch("Br_py1", &pho_py1, "Br_pho_py1/D");
    tree_tmp->Branch("Br_pz1", &pho_pz1, "Br_pho_pz1/D");
    tree_tmp->Branch("Br_E2", &pho_E2, "Br_pho_E2/D");
    tree_tmp->Branch("Br_px2", &pho_px2, "Br_pho_px2/D");
    tree_tmp->Branch("Br_py2", &pho_py2, "Br_pho_py2/D");
    tree_tmp->Branch("Br_pz2", &pho_pz2, "Br_pho_pz2/D");
    tree_tmp->Branch("Br_E3", &pho_E3, "Br_pho_E3/D");
    tree_tmp->Branch("Br_E1_true", &pho_E1_true, "Br_pho_E1_true/D");
    tree_tmp->Branch("Br_px1_true", &pho_px1_true, "Br_pho_px1_true/D");
    tree_tmp->Branch("Br_py1_true", &pho_py1_true, "Br_pho_py1_true/D");
    tree_tmp->Branch("Br_pz1_true", &pho_pz1_true, "Br_pho_pz1_true/D");
    tree_tmp->Branch("Br_E2_true", &pho_E2_true, "Br_pho_E2_true/D");
    tree_tmp->Branch("Br_px2_true", &pho_px2_true, "Br_pho_px2_true/D");
    tree_tmp->Branch("Br_py2_true", &pho_py2_true, "Br_pho_py2_true/D");
    tree_tmp->Branch("Br_pz2_true", &pho_pz2_true, "Br_pho_pz2_true/D");
    tree_tmp->Branch("Br_E3_true", &pho_E3_true, "Br_pho_E3_true/D");
    tree_tmp->Branch("Br_px3_true", &pho_px3_true, "Br_pho_px3_true/D");
    tree_tmp->Branch("Br_py3_true", &pho_py3_true, "Br_pho_py3_true/D");
    tree_tmp->Branch("Br_pz3_true", &pho_pz3_true, "Br_pho_pz3_true/D");
    tree_tmp->Branch("Br_px3", &pho_px3, "Br_pho_px3/D");
    tree_tmp->Branch("Br_py3", &pho_py3, "Br_pho_py3/D");
    tree_tmp->Branch("Br_pz3", &pho_pz3, "Br_pho_pz3/D");
    tree_tmp->Branch("Br_pull_E1", &pull_E1, "Br_pull_E1/D");
    tree_tmp->Branch("Br_pull_x1", &pull_x1, "Br_pull_x1/D");
    tree_tmp->Branch("Br_pull_y1", &pull_y1, "Br_pull_y1/D");
    tree_tmp->Branch("Br_pull_z1", &pull_z1, "Br_pull_z1/D");
    tree_tmp->Branch("Br_pull_t1", &pull_t1, "Br_pull_t1/D");
    tree_tmp->Branch("Br_pull_E2", &pull_E2, "Br_pull_E2/D");
    tree_tmp->Branch("Br_pull_x2", &pull_x2, "Br_pull_x2/D");
    tree_tmp->Branch("Br_pull_y2", &pull_y2, "Br_pull_y2/D");
    tree_tmp->Branch("Br_pull_z2", &pull_z2, "Br_pull_z2/D");
    tree_tmp->Branch("Br_pull_t2", &pull_t2, "Br_pull_t2/D");
    tree_tmp->Branch("Br_pull_E3", &pull_E3, "Br_pull_E3/D");
    tree_tmp->Branch("Br_pull_x3", &pull_x3, "Br_pull_x3/D");
    tree_tmp->Branch("Br_pull_y3", &pull_y3, "Br_pull_y3/D");
    tree_tmp->Branch("Br_pull_z3", &pull_z3, "Br_pull_z3/D");
    tree_tmp->Branch("Br_pull_t3", &pull_t3, "Br_pull_t3/D");
    tree_tmp->Branch("Br_e1_pull_bdt", &e1_pull_bdt, "Br_e1_pull_bdt/D");
    tree_tmp->Branch("Br_e2_pull_bdt", &e2_pull_bdt, "Br_e2_pull_bdt/D");
    tree_tmp->Branch("Br_e3_pull_bdt", &e3_pull_bdt, "Br_e3_pull_bdt/D");
    tree_tmp->Branch("Br_e1_pull_bdt_new", &e1_pull_bdt_new, "Br_e1_pull_bdt_new/D");
    tree_tmp->Branch("Br_e2_pull_bdt_new", &e2_pull_bdt_new, "Br_e2_pull_bdt_new/D");
    tree_tmp->Branch("Br_sig_type", &sig_type, "Br_sig_type/I");
    tree_tmp->Branch("Br_bkg_indx", &bkg_indx, "Br_bkg_indx/I");
    tree_tmp->Branch("Br_recon_indx", &recon_indx, "Br_recon_indx/I");
    tree_tmp->Branch("Br_IM3pi_7C", &IM3pi_7C, "Br_IM3pi_7C/D");
    tree_tmp->Branch("Br_IM3pi_true", &IM3pi_true, "Br_IM3pi_true/D");
    tree_tmp->Branch("Br_IM_pi0_7C", &IM_pi0_7C, "Br_IM_pi0_7C/D");
    tree_tmp->Branch("Br_IM_pi0_nofit", &IM_pi0_nofit, "Br_IM_pi0_nofit/D");
    tree_tmp->Branch("Br_mplus2", &mplus2, "Br_mplus2/D");
    tree_tmp->Branch("Br_m02", &m02, "Br_m02/D");
    tree_tmp->Branch("Br_ppIM", &ppIM, "Br_ppIM/D");
    tree_tmp->Branch("Br_ppIM_true", &ppIM_true, "Br_ppIM_true/D");
    tree_tmp->Branch("Br_Eisr", &Eisr, "Br_Eisr/D");
    tree_tmp->Branch("Br_Epi0_pho1", &Epi0_pho1, "Br_Epi0_pho1/D");
    tree_tmp->Branch("Br_Epi0_pho2", &Epi0_pho2, "Br_Epi0_pho2/D");
    tree_tmp->Branch("Br_angle_pi0gam12", &angle_pi0gam12, "Br_angle_pi0gam12/D");
    tree_tmp->Branch("Br_betapi0", &betapi0, "Br_betapi0/D");
    tree_tmp->Branch("Br_Eprompt_max", &Eprompt_max, "Br_Eprompt_max/D");
    tree_tmp->Branch("Br_lagvalue_min_7C", &lagvalue_min_7C, "Br_lagvalue_min_7C/D");
    tree_tmp->Branch("Br_deltaE", &deltaE, "Br_deltaE/D");
    tree_tmp->Branch("Br_m3pi", &m3pi, "Br_m3pi/D");
    
    // BDT branches
    tree_tmp->Branch("Br_bdt_score", &bdt_score, "Br_bdt_score/D");
    tree_tmp->Branch("Br_e1_bdt", &e1_bdt, "Br_e1_bdt/D");
    tree_tmp->Branch("Br_e2_bdt", &e2_bdt, "Br_e2_bdt/D");
    tree_tmp->Branch("Br_e3_bdt", &e3_bdt, "Br_e3_bdt/D");
    tree_tmp->Branch("Br_px1_bdt", &px1_bdt, "Br_px1_bdt/D");
    tree_tmp->Branch("Br_py1_bdt", &py1_bdt, "Br_py1_bdt/D");
    tree_tmp->Branch("Br_pz1_bdt", &pz1_bdt, "Br_pz1_bdt/D");
    tree_tmp->Branch("Br_px2_bdt", &px2_bdt, "Br_px2_bdt/D");
    tree_tmp->Branch("Br_py2_bdt", &py2_bdt, "Br_py2_bdt/D");
    tree_tmp->Branch("Br_pz2_bdt", &pz2_bdt, "Br_pz2_bdt/D");
    tree_tmp->Branch("Br_px3_bdt", &px3_bdt, "Br_px3_bdt/D");
    tree_tmp->Branch("Br_py3_bdt", &py3_bdt, "Br_py3_bdt/D");
    tree_tmp->Branch("Br_pz3_bdt", &pz3_bdt, "Br_pz3_bdt/D");
    
    // nofitted values
    tree_tmp->Branch("Br_e1_nofit_bdt", &e1_nofit_bdt, "Br_e1_nofit_bdt/D");
    tree_tmp->Branch("Br_e2_nofit_bdt", &e2_nofit_bdt, "Br_e2_nofit_bdt/D");
    tree_tmp->Branch("Br_e3_nofit_bdt", &e3_nofit_bdt, "Br_e3_nofit_bdt/D");
    
    tree_tmp->Branch("Br_px1_nofit_bdt", &px1_nofit_bdt, "Br_px1_nofit_bdt/D");
    tree_tmp->Branch("Br_py1_nofit_bdt", &py1_nofit_bdt, "Br_py1_nofit_bdt/D");
    tree_tmp->Branch("Br_pz1_nofit_bdt", &pz1_nofit_bdt, "Br_pz1_nofit_bdt/D");
    
    tree_tmp->Branch("Br_px2_nofit_bdt", &px2_nofit_bdt, "Br_px2_nofit_bdt/D");
    tree_tmp->Branch("Br_py2_nofit_bdt", &py2_nofit_bdt, "Br_py2_nofit_bdt/D");
    tree_tmp->Branch("Br_pz2_nofit_bdt", &pz2_nofit_bdt, "Br_pz2_nofit_bdt/D");
    
    tree_tmp->Branch("Br_px3_nofit_bdt", &px3_nofit_bdt, "Br_px3_nofit_bdt/D");
    tree_tmp->Branch("Br_py3_nofit_bdt", &py3_nofit_bdt, "Br_py3_nofit_bdt/D");
    tree_tmp->Branch("Br_pz3_nofit_bdt", &pz3_nofit_bdt, "Br_pz3_nofit_bdt/D");
    
    // BDT-selected sigma raw branches
    tree_tmp->Branch("Br_sigma_denom_e1", &sigma_denom_e1, "Br_sigma_denom_e1/D");
    tree_tmp->Branch("Br_sigma_denom_e2", &sigma_denom_e2, "Br_sigma_denom_e2/D");
    tree_tmp->Branch("Br_sigma_denom_e3", &sigma_denom_e3, "Br_sigma_denom_e3/D");
    
    // BDT-selected sigma kinematic fitted branches
    tree_tmp->Branch("Br_sigma_e1_bdt", &sigma_e1_bdt, "Br_sigma_e1_bdt/D");
    tree_tmp->Branch("Br_sigma_e2_bdt", &sigma_e2_bdt, "Br_sigma_e2_bdt/D");
    tree_tmp->Branch("Br_sigma_e3_bdt", &sigma_e3_bdt, "Br_sigma_e3_bdt/D");
    // Corrections
    tree_tmp->Branch("Br_delta_e1_bdt", &delta_e1_bdt, "Br_delta_e1_bdt/D");
    tree_tmp->Branch("Br_delta_e2_bdt", &delta_e2_bdt, "Br_delta_e2_bdt/D");
    
    // kinematic fitted + BDT selected (uncorrected)
    tree_tmp->Branch("Br_e1_fit", &e1_fit, "Br_e1_fit/D"); 
    tree_tmp->Branch("Br_e2_fit", &e2_fit, "Br_e2_fit/D");
    tree_tmp->Branch("Br_e3_fit", &e3_fit, "Br_e3_fit/D");
    
    tree_tmp->Branch("Br_e1_bdt_true", &e1_bdt_true, "Br_e1_bdt_true/D");
    tree_tmp->Branch("Br_e2_bdt_true", &e2_bdt_true, "Br_e2_bdt_true/D");
    tree_tmp->Branch("Br_e3_bdt_true", &e3_bdt_true, "Br_e3_bdt_true/D");
    tree_tmp->Branch("Br_m_gg_bdt", &m_gg_bdt, "Br_m_gg_bdt/D");
    tree_tmp->Branch("Br_m_gg_true_bdt", &m_gg_true, "Br_m_gg_true_bdt/D");
    tree_tmp->Branch("Br_m2pi_true", &m2pi_true, "Br_m2pi_true/D");
    tree_tmp->Branch("Br_m3pi_bdt", &m3pi_bdt, "Br_m3pi_bdt/D");
    tree_tmp->Branch("Br_m3pi_true_bdt", &m3pi_true, "Br_m3pi_true_bdt/D");
    tree_tmp->Branch("Br_angle_pi0gam12_bdt", &angle_pi0gam12_bdt, "Br_angle_pi0gam12_bdt/D");
    tree_tmp->Branch("Br_angle_pi0gam12_bdt_true", &angle_pi0gam12_bdt_true, "Br_angle_pi0gam12_bdt_true/D");
    tree_tmp->Branch("Br_betapi0_bdt", &betapi0_bdt, "Br_betapi0_bdt/D");
    tree_tmp->Branch("Br_betapi0_bdt_true", &betapi0_bdt_true, "Br_betapi0_bdt_true/D");
    tree_tmp->Branch("Br_angle_ppl_pmi", &angle_ppl_pmi, "Br_angle_ppl_pmi/D");
    tree_tmp->Branch("Br_angle_trk_neutral", &angle_trk_neutral, "Br_angle_trk_neutral/D");
    tree_tmp->Branch("Br_beta_3pi", &beta_3pi, "Br_beta_3pi/D");
    tree_tmp->Branch("Br_e_asym", &e_asym, "Br_e_asym/D");
    // generator MC true
    tree_tmp->Branch("Br_true_px_piminus", &true_px_piminus, "Br_true_px_piminus/D");
    tree_tmp->Branch("Br_true_py_piminus", &true_py_piminus, "Br_true_py_piminus/D");
    tree_tmp->Branch("Br_true_pz_piminus", &true_pz_piminus, "Br_true_pz_piminus/D");
    tree_tmp->Branch("Br_true_E_piminus", &true_E_piminus, "Br_true_E_piminus/D");
    tree_tmp->Branch("Br_true_px_piplus", &true_px_piplus, "Br_true_px_piplus/D");
    tree_tmp->Branch("Br_true_py_piplus", &true_py_piplus, "Br_true_py_piplus/D");
    tree_tmp->Branch("Br_true_pz_piplus", &true_pz_piplus, "Br_true_pz_piplus/D");
    tree_tmp->Branch("Br_true_E_piplus", &true_E_piplus, "Br_true_E_piplus/D");
    tree_tmp->Branch("Br_true_px_pi0", &true_px_pi0, "Br_true_px_pi0/D");
    tree_tmp->Branch("Br_true_py_pi0", &true_py_pi0, "Br_true_py_pi0/D");
    tree_tmp->Branch("Br_true_pz_pi0", &true_pz_pi0, "Br_true_pz_pi0/D");
    tree_tmp->Branch("Br_true_E_pi0", &true_E_pi0, "Br_true_E_pi0/D");
    tree_tmp->Branch("Br_true_m3pi", &true_m3pi, "Br_true_m3pi/D");

    // Pulls
    tree_tmp->Branch("Br_e1_pull", &e1_pull, "Br_e1_pull/D");
    tree_tmp->Branch("Br_e2_pull", &e2_pull, "Br_e2_pull/D");
    tree_tmp->Branch("Br_e3_pull", &e3_pull, "Br_e3_pull/D");
    tree_tmp->Branch("Br_px1_pull", &px1_pull, "Br_px1_pull/D");
    tree_tmp->Branch("Br_py1_pull", &py1_pull, "Br_py1_pull/D");
    tree_tmp->Branch("Br_pz1_pull", &pz1_pull, "Br_pz1_pull/D");
    tree_tmp->Branch("Br_px2_pull", &px2_pull, "Br_px2_pull/D");
    tree_tmp->Branch("Br_py2_pull", &py2_pull, "Br_py2_pull/D");
    tree_tmp->Branch("Br_pz2_pull", &pz2_pull, "Br_pz2_pull/D");
    tree_tmp->Branch("Br_px3_pull", &px3_pull, "Br_px3_pull/D");
    tree_tmp->Branch("Br_py3_pull", &py3_pull, "Br_py3_pull/D");
    tree_tmp->Branch("Br_pz3_pull", &pz3_pull, "Br_pz3_pull/D");
    tree_tmp->Branch("Br_m_gg_pull", &m_gg_pull, "Br_m_gg_pull/D");
    tree_tmp->Branch("Br_m3pi_pull", &m3pi_pull, "Br_m3pi_pull/D");
    tree_tmp->Branch("Br_m2pi_pull", &m2pi_pull, "Br_m2pi_pull/D");
    tree_tmp->Branch("Br_recon_indx_bdt", &recon_indx_bdt, "Br_recon_indx_bdt/I");
    tree_tmp->Branch("Br_isr_recon_quality", &isr_recon_quality, "Br_isr_recon_quality/I");
    tree_tmp->Branch("Br_total_recon_quality", &total_recon_quality, "Br_total_recon_quality/I");
  }
  
  TLorentzVector pi0gam1, pi0gam2, isrgam, trkplus, trkmin;
  TLorentzVector trksum, neutralsum, system_3pi;
  
  // ============================================================
  // SET BRANCH ADDRESSES
  // ============================================================
  
  // Reco tracks
  ALLCHAIN_CUT->SetBranchAddress("Br_ppl_E", &ppl_E);
  ALLCHAIN_CUT->SetBranchAddress("Br_ppl_px", &ppl_px);
  ALLCHAIN_CUT->SetBranchAddress("Br_ppl_py", &ppl_py);
  ALLCHAIN_CUT->SetBranchAddress("Br_ppl_pz", &ppl_pz);
  ALLCHAIN_CUT->SetBranchAddress("Br_pmi_E", &pmi_E);
  ALLCHAIN_CUT->SetBranchAddress("Br_pmi_px", &pmi_px);
  ALLCHAIN_CUT->SetBranchAddress("Br_pmi_py", &pmi_py);
  ALLCHAIN_CUT->SetBranchAddress("Br_pmi_pz", &pmi_pz);
  
  
  ALLCHAIN_CUT->SetBranchAddress("Br_E1", &pho_E1);
  ALLCHAIN_CUT->SetBranchAddress("Br_px1", &pho_px1);
  ALLCHAIN_CUT->SetBranchAddress("Br_py1", &pho_py1);
  ALLCHAIN_CUT->SetBranchAddress("Br_pz1", &pho_pz1);
  ALLCHAIN_CUT->SetBranchAddress("Br_E2", &pho_E2);
  ALLCHAIN_CUT->SetBranchAddress("Br_px2", &pho_px2);
  ALLCHAIN_CUT->SetBranchAddress("Br_py2", &pho_py2);
  ALLCHAIN_CUT->SetBranchAddress("Br_pz2", &pho_pz2);
  ALLCHAIN_CUT->SetBranchAddress("Br_E3", &pho_E3);
  ALLCHAIN_CUT->SetBranchAddress("Br_px3", &pho_px3);
  ALLCHAIN_CUT->SetBranchAddress("Br_py3", &pho_py3);
  ALLCHAIN_CUT->SetBranchAddress("Br_pz3", &pho_pz3);
  
  // Unfitted (pre‑kinematic‑fit) 4‑vectors in fit order
  ALLCHAIN_CUT->SetBranchAddress("Br_E1_nofit", &pho_E1_nofit);
  ALLCHAIN_CUT->SetBranchAddress("Br_px1_nofit", &pho_px1_nofit);
  ALLCHAIN_CUT->SetBranchAddress("Br_py1_nofit", &pho_py1_nofit);
  ALLCHAIN_CUT->SetBranchAddress("Br_pz1_nofit", &pho_pz1_nofit);
  
  ALLCHAIN_CUT->SetBranchAddress("Br_E2_nofit", &pho_E2_nofit);
  ALLCHAIN_CUT->SetBranchAddress("Br_px2_nofit", &pho_px2_nofit);
  ALLCHAIN_CUT->SetBranchAddress("Br_py2_nofit", &pho_py2_nofit);
  ALLCHAIN_CUT->SetBranchAddress("Br_pz2_nofit", &pho_pz2_nofit);
  
  ALLCHAIN_CUT->SetBranchAddress("Br_E3_nofit", &pho_E3_nofit);
  ALLCHAIN_CUT->SetBranchAddress("Br_px3_nofit", &pho_px3_nofit);
  ALLCHAIN_CUT->SetBranchAddress("Br_py3_nofit", &pho_py3_nofit);
  ALLCHAIN_CUT->SetBranchAddress("Br_pz3_nofit", &pho_pz3_nofit);
  
  // Mapping indices from fit order to original cluster order
  ALLCHAIN_CUT->SetBranchAddress("Br_isrgam_indx", &isrgam_indx);
  ALLCHAIN_CUT->SetBranchAddress("Br_pi0gam1_indx", &pi0gam1_indx);
  ALLCHAIN_CUT->SetBranchAddress("Br_pi0gam2_indx", &pi0gam2_indx);
  
  // Reconstructed true tracks
  ALLCHAIN_CUT->SetBranchAddress("Br_ppl_E_true", &ppl_E_true);
  ALLCHAIN_CUT->SetBranchAddress("Br_ppl_px_true", &ppl_px_true);
  ALLCHAIN_CUT->SetBranchAddress("Br_ppl_py_true", &ppl_py_true);
  ALLCHAIN_CUT->SetBranchAddress("Br_ppl_pz_true", &ppl_pz_true);
  ALLCHAIN_CUT->SetBranchAddress("Br_pmi_E_true", &pmi_E_true);
  ALLCHAIN_CUT->SetBranchAddress("Br_pmi_px_true", &pmi_px_true);
  ALLCHAIN_CUT->SetBranchAddress("Br_pmi_py_true", &pmi_py_true);
  ALLCHAIN_CUT->SetBranchAddress("Br_pmi_pz_true", &pmi_pz_true);
  
  // Reconstructed true photons
  ALLCHAIN_CUT->SetBranchAddress("Br_E1_true", &pho_E1_true);
  ALLCHAIN_CUT->SetBranchAddress("Br_px1_true", &pho_px1_true);
  ALLCHAIN_CUT->SetBranchAddress("Br_py1_true", &pho_py1_true);
  ALLCHAIN_CUT->SetBranchAddress("Br_pz1_true", &pho_pz1_true);
  ALLCHAIN_CUT->SetBranchAddress("Br_E2_true", &pho_E2_true);
  ALLCHAIN_CUT->SetBranchAddress("Br_px2_true", &pho_px2_true);
  ALLCHAIN_CUT->SetBranchAddress("Br_py2_true", &pho_py2_true);
  ALLCHAIN_CUT->SetBranchAddress("Br_pz2_true", &pho_pz2_true);
  ALLCHAIN_CUT->SetBranchAddress("Br_E3_true", &pho_E3_true);
  ALLCHAIN_CUT->SetBranchAddress("Br_px3_true", &pho_px3_true);
  ALLCHAIN_CUT->SetBranchAddress("Br_py3_true", &pho_py3_true);
  ALLCHAIN_CUT->SetBranchAddress("Br_pz3_true", &pho_pz3_true);
  
  // Generator-level true particles
  ALLCHAIN_CUT->SetBranchAddress("true_px_piplus", &true_px_piplus);
  ALLCHAIN_CUT->SetBranchAddress("true_py_piplus", &true_py_piplus);
  ALLCHAIN_CUT->SetBranchAddress("true_pz_piplus", &true_pz_piplus);
  ALLCHAIN_CUT->SetBranchAddress("true_E_piplus", &true_E_piplus);
  ALLCHAIN_CUT->SetBranchAddress("true_px_piminus", &true_px_piminus);
  ALLCHAIN_CUT->SetBranchAddress("true_py_piminus", &true_py_piminus);
  ALLCHAIN_CUT->SetBranchAddress("true_pz_piminus", &true_pz_piminus);
  ALLCHAIN_CUT->SetBranchAddress("true_E_piminus", &true_E_piminus);
  ALLCHAIN_CUT->SetBranchAddress("true_px_pi0", &true_px_pi0);
  ALLCHAIN_CUT->SetBranchAddress("true_py_pi0", &true_py_pi0);
  ALLCHAIN_CUT->SetBranchAddress("true_pz_pi0", &true_pz_pi0);
  ALLCHAIN_CUT->SetBranchAddress("true_E_pi0", &true_E_pi0);
  
  // Other single variables
  ALLCHAIN_CUT->SetBranchAddress("Br_bkg_indx", &bkg_indx);
  ALLCHAIN_CUT->SetBranchAddress("Br_recon_indx", &recon_indx);
  ALLCHAIN_CUT->SetBranchAddress("Br_phid", &phid);
  ALLCHAIN_CUT->SetBranchAddress("Br_sig_type", &sig_type);
  ALLCHAIN_CUT->SetBranchAddress("Br_lagvalue_min_7C", &lagvalue_min_7C);
  ALLCHAIN_CUT->SetBranchAddress("Br_betapi0", &betapi0);
  ALLCHAIN_CUT->SetBranchAddress("Br_IM3pi_7C", &IM3pi_7C);
  ALLCHAIN_CUT->SetBranchAddress("Br_IM_pi0_7C", &IM_pi0_7C);
  ALLCHAIN_CUT->SetBranchAddress("Br_IM_pi0_nofit", &IM_pi0_nofit);
  ALLCHAIN_CUT->SetBranchAddress("Br_IM3pi_true", &IM3pi_true);
  
  // ----- Beam branch addresses -----
  ALLCHAIN_CUT->SetBranchAddress("Br_bpx", &bpx);   // <-- ADDED
  ALLCHAIN_CUT->SetBranchAddress("Br_bpy", &bpy);   // <-- ADDED
  ALLCHAIN_CUT->SetBranchAddress("Br_bpz", &bpz);   // <-- ADDED
  ALLCHAIN_CUT->SetBranchAddress("Br_bene", &bene); // <-- ADDED
  
  // ============================================================
  // ARRAYS - Use GetLeaf (SetBranchAddress causes crash!)
  // ============================================================
  
  cout << "MASS_SCALE_PI0 = " << MASS_SCALE_PI0  << endl;
  
  // ---------- Event loop ----------
  Long64_t nentries = ALLCHAIN_CUT->GetEntries();
  cout << "Processing " << nentries << " events" << endl;
  
  for (Long64_t irow = 0; irow < nentries; irow++) {
    ALLCHAIN_CUT->GetEntry(irow);
    if (irow % 100000 == 0) cout << "Event " << irow << "/ " << nentries << endl;
    
    // ============================================================
    // GetLeaf for arrays (these work!)
    // ============================================================
    
    pho_indx[0] = ALLCHAIN_CUT->GetLeaf("Br_pho_indx")->GetValue(0);
    pho_indx[1] = ALLCHAIN_CUT->GetLeaf("Br_pho_indx")->GetValue(1);
    pho_indx[2] = ALLCHAIN_CUT->GetLeaf("Br_pho_indx")->GetValue(2);
    
    EPI0NTMC[0] = ALLCHAIN_CUT->GetLeaf("Br_EPI0NTMC_save")->GetValue(0);
    EPI0NTMC[1] = ALLCHAIN_CUT->GetLeaf("Br_EPI0NTMC_save")->GetValue(1);
    EPI0NTMC[2] = ALLCHAIN_CUT->GetLeaf("Br_EPI0NTMC_save")->GetValue(2);
    EPI0NTMC[3] = ALLCHAIN_CUT->GetLeaf("Br_EPI0NTMC_save")->GetValue(3);
    
    pull_E1 = ALLCHAIN_CUT->GetLeaf("Br_PULLIST")->GetValue(0);
    pull_x1 = ALLCHAIN_CUT->GetLeaf("Br_PULLIST")->GetValue(1);
    pull_y1 = ALLCHAIN_CUT->GetLeaf("Br_PULLIST")->GetValue(2);
    pull_z1 = ALLCHAIN_CUT->GetLeaf("Br_PULLIST")->GetValue(3);
    pull_t1 = ALLCHAIN_CUT->GetLeaf("Br_PULLIST")->GetValue(4);
    pull_E2 = ALLCHAIN_CUT->GetLeaf("Br_PULLIST")->GetValue(5);
    pull_x2 = ALLCHAIN_CUT->GetLeaf("Br_PULLIST")->GetValue(6);
    pull_y2 = ALLCHAIN_CUT->GetLeaf("Br_PULLIST")->GetValue(7);
    pull_z2 = ALLCHAIN_CUT->GetLeaf("Br_PULLIST")->GetValue(8);
    pull_t2 = ALLCHAIN_CUT->GetLeaf("Br_PULLIST")->GetValue(9);
    pull_E3 = ALLCHAIN_CUT->GetLeaf("Br_PULLIST")->GetValue(10);
    pull_x3 = ALLCHAIN_CUT->GetLeaf("Br_PULLIST")->GetValue(11);
    pull_y3 = ALLCHAIN_CUT->GetLeaf("Br_PULLIST")->GetValue(12);
    pull_z3 = ALLCHAIN_CUT->GetLeaf("Br_PULLIST")->GetValue(13);
    pull_t3 = ALLCHAIN_CUT->GetLeaf("Br_PULLIST")->GetValue(14);
    
    ppIM_true = ALLCHAIN_CUT->GetLeaf("Br_ppIM_true")->GetValue(0);
    ppIM = ALLCHAIN_CUT->GetLeaf("Br_MASSLIST")->GetValue(5);
    m02 = ALLCHAIN_CUT->GetLeaf("Br_MASSLIST")->GetValue(10);
    mplus2 = ALLCHAIN_CUT->GetLeaf("Br_MASSLIST")->GetValue(11);
    
    Eisr = ALLCHAIN_CUT->GetLeaf("Br_ENERGYLIST")->GetValue(0);
    Epi0_pho1 = ALLCHAIN_CUT->GetLeaf("Br_ENERGYLIST")->GetValue(1);
    deltaE = ALLCHAIN_CUT->GetLeaf("Br_ENERGYLIST")->GetValue(2);
    Epi0_pho2 = ALLCHAIN_CUT->GetLeaf("Br_ENERGYLIST")->GetValue(3);
    
    angle_pi0gam12 = ALLCHAIN_CUT->GetLeaf("Br_ANGLELIST")->GetValue(0);
    
    // Read SIGMA_FIT_LIST and SIGMA_DENOM_LIST (15 values)
    double SIGMA_FIT_LIST[15];
    double SIGMA_DENOM_LIST[15];
    
    for (int idx = 0; idx < 15; ++idx) {
      SIGMA_FIT_LIST[idx] = ALLCHAIN_CUT->GetLeaf("Br_sigma_fit")->GetValue(idx);
      SIGMA_DENOM_LIST[idx] = ALLCHAIN_CUT->GetLeaf("Br_sigma_denom")->GetValue(idx);
    }
    
    // ============================================================
    // ppl_E, pho_E1, etc. are filled by SetBranchAddress + GetEntry!
    // Do NOT use GetLeaf for these!
    // ============================================================
    
    // Build Lorentz vectors
    pi0gam1.SetPxPyPzE(pho_px1, pho_py1, pho_pz1, pho_E1);
    pi0gam2.SetPxPyPzE(pho_px2, pho_py2, pho_pz2, pho_E2);
    isrgam.SetPxPyPzE(pho_px3, pho_py3, pho_pz3, pho_E3);
    trkplus.SetPxPyPzE(ppl_px, ppl_py, ppl_pz, ppl_E);
    trkmin.SetPxPyPzE(pmi_px, pmi_py, pmi_pz, pmi_E);
    trksum = trkplus + trkmin;
    neutralsum = pi0gam1 + pi0gam2 + isrgam;
    angle_ppl_pmi = trkplus.Angle(trkmin.Vect())*TMath::RadToDeg();
    angle_trk_neutral = trksum.Angle(neutralsum.Vect())*TMath::RadToDeg();
    m3pi = (pi0gam1 + pi0gam2 + trkplus + trkmin).M();
    
    // Maximum prompt photon energy before bdt and pull tuning
    // Effective of etagam suppression
    Eprompt_max = std::max({Eisr, Epi0_pho1, Epi0_pho2});
    
    // Photons are sorted by energies
    TLorentzVector photon[3] = {pi0gam1, pi0gam2, isrgam};
    std::sort(photon, photon+3, [](const TLorentzVector& a, const TLorentzVector& b) {
      return a.E() > b.E();
    });
    
    evnt_tot++;
    
    // ---------- BDT evaluation ----------
    EventData event;
    event.photons[0][0] = pho_E1; event.photons[0][1] = pho_px1; event.photons[0][2] = pho_py1; event.photons[0][3] = pho_pz1;
    event.photons[1][0] = pho_E2; event.photons[1][1] = pho_px2; event.photons[1][2] = pho_py2; event.photons[1][3] = pho_pz2;
    event.photons[2][0] = pho_E3; event.photons[2][1] = pho_px3; event.photons[2][2] = pho_py3; event.photons[2][3] = pho_pz3;
    event.tracks[0][0] = ppl_E; event.tracks[0][1] = ppl_px; event.tracks[0][2] = ppl_py; event.tracks[0][3] = ppl_pz;
    event.tracks[1][0] = pmi_E; event.tracks[1][1] = pmi_px; event.tracks[1][2] = pmi_py; event.tracks[1][3] = pmi_pz;
    event.lagvalue_min_7C = lagvalue_min_7C;
    event.deltaE = deltaE;
    event.betapi0 = betapi0;
    event.angle_pi0gam12 = angle_pi0gam12;
    event.ppIM = ppIM;
    event.bkg_indx = bkg_indx;
    event.recon_indx = recon_indx;
    
    // Prepare event true data structure
    double true_photons[3][4] = {
      {pho_E1_true, pho_px1_true, pho_py1_true, pho_pz1_true},
      {pho_E2_true, pho_px2_true, pho_py2_true, pho_pz2_true},
      {pho_E3_true, pho_px3_true, pho_py3_true, pho_pz3_true}
    };
    double true_tracks[2][4] = {
      {ppl_E_true, ppl_px_true, ppl_py_true, ppl_pz_true},
      {pmi_E_true, pmi_px_true, pmi_py_true, pmi_pz_true}
    };
    
    // Prepare event unfitted data structure
    double unfit_photons[3][4] = {
      {pho_E1_nofit, pho_px1_nofit, pho_py1_nofit, pho_pz1_nofit},
      {pho_E2_nofit, pho_px2_nofit, pho_py2_nofit, pho_pz2_nofit},
      {pho_E3_nofit, pho_px3_nofit, pho_py3_nofit, pho_pz3_nofit}
    };
    
    
    BDTResult result = find_best_pion_pair(event, bdt);
    if (!result.is_valid) continue;
    
    // Mapping from event index (0,1,2) to block index in sigma arrays:
    // block 0 = ISR, block 1 = pi0 photon1, block 2 = pi0 photon2
    int block_index[3] = {1, 2, 0};   // event[0]->block1, event[1]->block2, event[2]->block0
    
    // Mapping from event index to natural index for pulls
    int natural_index[3] = {pi0gam1_indx, pi0gam2_indx, isrgam_indx};
    
    int idx0 = result.pi0_indices[0];
    int idx1 = result.pi0_indices[1];
    int idx_prompt = result.prompt_index;
    
    // Sigma denominator (pull denominator)
    sigma_denom_e1 = SIGMA_DENOM_LIST[5 * block_index[idx0]];
    sigma_denom_e2 = SIGMA_DENOM_LIST[5 * block_index[idx1]];
    sigma_denom_e3 = SIGMA_DENOM_LIST[5 * block_index[idx_prompt]];
    
    // Fitted sigma
    sigma_e1_bdt = SIGMA_FIT_LIST[5 * block_index[idx0]];
    sigma_e2_bdt = SIGMA_FIT_LIST[5 * block_index[idx1]];
    sigma_e3_bdt = SIGMA_FIT_LIST[5 * block_index[idx_prompt]];
    
    // Pulls
    double pull_array[3] = {pull_E1, pull_E2, pull_E3};  // natural order
    e1_pull_bdt = pull_array[natural_index[idx0]];
    e2_pull_bdt = pull_array[natural_index[idx1]];
    e3_pull_bdt = pull_array[natural_index[idx_prompt]];
    
    // Save kinematic fitted values (uncorrected)
    e1_fit = event.photons[result.pi0_indices[0]][0];
    e2_fit = event.photons[result.pi0_indices[1]][0];
    e3_fit = event.photons[result.prompt_index][0];
    
    double px1_fit = event.photons[result.pi0_indices[0]][1];
    double py1_fit = event.photons[result.pi0_indices[0]][2];
    double pz1_fit = event.photons[result.pi0_indices[0]][3];
    
    double px2_fit = event.photons[result.pi0_indices[1]][1];
    double py2_fit = event.photons[result.pi0_indices[1]][2];
    double pz2_fit = event.photons[result.pi0_indices[1]][3];
    
    double px3_fit = event.photons[result.prompt_index][1];
    double py3_fit = event.photons[result.prompt_index][2];
    double pz3_fit = event.photons[result.prompt_index][3];
    
    // ---- Assign unfitted BDT‑selected values ----
    e1_nofit_bdt = unfit_photons[result.pi0_indices[0]][0];
    e2_nofit_bdt = unfit_photons[result.pi0_indices[1]][0];
    e3_nofit_bdt = unfit_photons[result.prompt_index][0];
    px1_nofit_bdt = unfit_photons[result.pi0_indices[0]][1];
    py1_nofit_bdt = unfit_photons[result.pi0_indices[0]][2];
    pz1_nofit_bdt = unfit_photons[result.pi0_indices[0]][3];
    px2_nofit_bdt = unfit_photons[result.pi0_indices[1]][1];
    py2_nofit_bdt = unfit_photons[result.pi0_indices[1]][2];
    pz2_nofit_bdt = unfit_photons[result.pi0_indices[1]][3];
    px3_nofit_bdt = unfit_photons[result.prompt_index][1];
    py3_nofit_bdt = unfit_photons[result.prompt_index][2];
    pz3_nofit_bdt = unfit_photons[result.prompt_index][3];
    
    
    delta_e1_bdt = bias_shift * sigma_denom_e1;
    delta_e2_bdt = bias_shift * sigma_denom_e2;
    
    //cout << "bias_shift = " << bias_shift << ", sigma_denom_e1 = " << sigma_denom_e1 << ", delta_e1_bdt = " << delta_e1_bdt << "\n";
    
    e1_pull_bdt_new = 0.;
    e2_pull_bdt_new = 0.;
    
    // Full transformation of the pull distribution to match data exactly
    // Pull_corr = alpha * Pull + beta
    // where alpha = scale_ratio; beta = bias_data - alpha * bias_sig
    // Assume bias_sig close to zero, beta = bias_shift
    // Solving for E_corr: E_corr = (1 - alpha)E_raw + alpha*E_fit - beta*sigma_denom

    const double M_PION_CONST = 139.57061;  // MeV/c^2

    // ---- Apply corrections ----
    if (data_type == "sig") {
      double M_true = ppIM_true; // true dipion mass

      /*
      if (M_true > 0) {
	double s = M_true / (M_true + BIAS_P0 + BIAS_P1 * M_true);
	
	// Safety: clamp to reasonable range (prevents pathological low-mass extrapolation)
	if (s < 0.90) s = 0.90;
	if (s > 1.10) s = 1.10;
	
	ppl_px *= s;
	ppl_py *= s;
	ppl_pz *= s;
	ppl_E = sqrt(M_PION_CONST*M_PION_CONST + ppl_px*ppl_px + ppl_py*ppl_py + ppl_pz*ppl_pz);
	
	pmi_px *= s;
	pmi_py *= s;
	pmi_pz *= s;
	pmi_E = sqrt(M_PION_CONST*M_PION_CONST + pmi_px*pmi_px + pmi_py*pmi_py + pmi_pz*pmi_pz);
      }
      else {// Fallback: if true mass is unphysical (should not happen for signal), apply constant scale
	// 1. Bias correction: scale momenta by (1 - BIAS_P1)
	const double scale = 1.0 - BIAS_P1;
	ppl_px *= scale;
	ppl_py *= scale;
	ppl_pz *= scale;
	ppl_E = sqrt(M_PION_CONST*M_PION_CONST + ppl_px*ppl_px + ppl_py*ppl_py + ppl_pz*ppl_pz);
	
	pmi_px *= scale;
	pmi_py *= scale;
	pmi_pz *= scale;
	pmi_E = sqrt(M_PION_CONST*M_PION_CONST + pmi_px*pmi_px + pmi_py*pmi_py + pmi_pz*pmi_pz);
      }
      
      // 2. Mass-dependent smearing

      // MC resolution from logistic (linear plateau variant)
      double high = H0 + H1 * (M_true - MREF);
      double sigma_mc = SIGMA_LOW + (high - SIGMA_LOW) / (1.0 + exp(-K * (M_true - M0)));
      
      // Target data resolution (scaled by omega width ratio)
      double sigma_data = R_OMEGA * sigma_mc;
      
      // Relative momentum smearing width
      double delta = (M_true > 0) ? sigma_data / M_true : 0.0;
      
      // Smear each track independently
      double smear_plus  = 1.0 + gRandom->Gaus(0.0, delta);
      double smear_minus = 1.0 + gRandom->Gaus(0.0, delta);
      
      ppl_px *= smear_plus; ppl_py *= smear_plus; ppl_pz *= smear_plus;
      ppl_E  = sqrt(M_PION_CONST*M_PION_CONST + ppl_px*ppl_px + ppl_py*ppl_py + ppl_pz*ppl_pz);
      
      pmi_px *= smear_minus; pmi_py *= smear_minus; pmi_pz *= smear_minus;
      pmi_E  = sqrt(M_PION_CONST*M_PION_CONST + pmi_px*pmi_px + pmi_py*pmi_py + pmi_pz*pmi_pz);
      
      // ---- Update the event structure ----
      event.tracks[0][0] = ppl_E;
      event.tracks[0][1] = ppl_px;
      event.tracks[0][2] = ppl_py;
      event.tracks[0][3] = ppl_pz;
          
      event.tracks[1][0] = pmi_E;
      event.tracks[1][1] = pmi_px;
      event.tracks[1][2] = pmi_py;
      event.tracks[1][3] = pmi_pz;
      
      // ---- Update Lorentz vectors and ppIM with corrected tracks ----
      trkplus.SetPxPyPzE(ppl_px, ppl_py, ppl_pz, ppl_E);
      trkmin.SetPxPyPzE(pmi_px, pmi_py, pmi_pz, pmi_E);
      ppIM = compute_dipion_mass(event.tracks);   // recompute corrected dipion mass
      */
      
      // ---- Recompute deltaE using beam vector ----
      TLorentzVector beamVec; beamVec.SetPxPyPzE(bpx, bpy, bpz, bene);
      deltaE = compute_deltaE(trkplus, trkmin, beamVec, M_PION_CONST);  // <-- MODIFIED
      //cout << "(bpx, bpy, bpz, bene) = (" << bpx << ", " << bpy << ", " << bpz << ", " << bene << ")\n";
      
      // The rest of the photon corrections (e1_bdt, etc.) remain unchanged.
      // ... (existing code for photon corrections)
      
      //const double delta_e1_bdt = bias_E12 * resol_E12;  
      //const double bias_MeV_E3 = bias_E3 * resol_E3;
      
      double sigma_e1_bdt_corr = scale_ratio * sigma_e1_bdt;
      double sigma_e2_bdt_corr = scale_ratio * sigma_e2_bdt;
      
      // New pulls
      e1_pull_bdt_new = e1_pull_bdt * scale_ratio + bias_shift;
      e2_pull_bdt_new = e2_pull_bdt * scale_ratio + bias_shift;
      
      // Correct MC energies (with bias shift and mass scale)
      //e1_bdt = (e1_fit - delta_e1_bdt) * MASS_SCALE_PI0;
      //e2_bdt = (e2_fit - delta_e2_bdt) * MASS_SCALE_PI0;
      e1_bdt = (1 - scale_ratio) * e1_nofit_bdt + scale_ratio * e1_fit - delta_e1_bdt;
      e1_bdt *= MASS_SCALE_PI0;
      
      e2_bdt = (1 - scale_ratio) * e2_nofit_bdt + scale_ratio * e2_fit - delta_e2_bdt;
      e2_bdt *= MASS_SCALE_PI0;
      
      e3_bdt = e3_fit; // ISR unchanged
      
      // Correct momenta: preserve photon direction from fitted momentum, set magnitude = energy
      double p1_mag = sqrt(px1_fit*px1_fit + py1_fit*py1_fit + pz1_fit*pz1_fit);
      if (p1_mag > 0) {
        px1_bdt = e1_bdt * (px1_fit / p1_mag);
        py1_bdt = e1_bdt * (py1_fit / p1_mag);
        pz1_bdt = e1_bdt * (pz1_fit / p1_mag);
      } else {
        // fallback: just apply mass scale
        px1_bdt = px1_fit * MASS_SCALE_PI0;
        py1_bdt = py1_fit * MASS_SCALE_PI0;
        pz1_bdt = pz1_fit * MASS_SCALE_PI0;
      }
      
      // Same for photon 2
      double p2_mag = sqrt(px2_fit*px2_fit + py2_fit*py2_fit + pz2_fit*pz2_fit);
      if (p2_mag > 0) {
        px2_bdt = e2_bdt * (px2_fit / p2_mag);
        py2_bdt = e2_bdt * (py2_fit / p2_mag);
        pz2_bdt = e2_bdt * (pz2_fit / p2_mag);
      } else {
        px2_bdt = px2_fit * MASS_SCALE_PI0;
        py2_bdt = py2_fit * MASS_SCALE_PI0;
        pz2_bdt = pz2_fit * MASS_SCALE_PI0;
      }
      
      // ISR photon: no energy correction, so momentum unchanged
      px3_bdt = px3_fit;
      py3_bdt = py3_fit;
      pz3_bdt = pz3_fit;
      
      event.photons[result.pi0_indices[0]][0] = e1_bdt;
      event.photons[result.pi0_indices[0]][1] = px1_bdt;
      event.photons[result.pi0_indices[0]][2] = py1_bdt;
      event.photons[result.pi0_indices[0]][3] = pz1_bdt;
          
      event.photons[result.pi0_indices[1]][0] = e2_bdt;
      event.photons[result.pi0_indices[1]][1] = px2_bdt;
      event.photons[result.pi0_indices[1]][2] = py2_bdt;
      event.photons[result.pi0_indices[1]][3] = pz2_bdt;
          
      event.photons[result.prompt_index][0] = e3_bdt;
      event.photons[result.prompt_index][1] = px3_bdt;
      event.photons[result.prompt_index][2] = py3_bdt;
      event.photons[result.prompt_index][3] = pz3_bdt;
    } else {
      e1_bdt = e1_fit;
      e2_bdt = e2_fit;
      e3_bdt = e3_fit;
      
      px1_bdt = px1_fit;
      py1_bdt = py1_fit;
      pz1_bdt = pz1_fit;
      
      px2_bdt = px2_fit;
      py2_bdt = py2_fit;
      pz2_bdt = pz2_fit;
      
      px3_bdt = px3_fit;
      py3_bdt = py3_fit;
      pz3_bdt = pz3_fit;
    }
    
    m_gg_bdt = compute_invariant_mass(result.pi0_indices[0], result.pi0_indices[1], event.photons);
    m3pi_bdt = compute_3pi_mass(result.pi0_indices[0], result.pi0_indices[1], event.photons, event.tracks);
    
    TLorentzVector pi0gam1_bdt, pi0gam2_bdt;
    pi0gam1_bdt.SetPxPyPzE(px1_bdt, py1_bdt, pz1_bdt, e1_bdt);
    pi0gam2_bdt.SetPxPyPzE(px2_bdt, py2_bdt, pz2_bdt, e2_bdt);
    angle_pi0gam12_bdt = pi0gam1_bdt.Angle(pi0gam2_bdt.Vect()) * TMath::RadToDeg();
    
    TLorentzVector pi0_bdt = pi0gam1_bdt + pi0gam2_bdt;
    betapi0_bdt = (pi0_bdt.Vect()).Mag() / pi0_bdt.E();
    bdt_score = result.score;
    
    // Generator true
    TLorentzVector piplus_gen, piminus_gen, pi0_gen;
    piplus_gen.SetPxPyPzE(true_px_piplus, true_py_piplus, true_pz_piplus, true_E_piplus);
    piminus_gen.SetPxPyPzE(true_px_piminus, true_py_piminus, true_pz_piminus, true_E_piminus);
    pi0_gen.SetPxPyPzE(true_px_pi0, true_py_pi0, true_pz_pi0, true_E_pi0);
    true_m3pi = (piplus_gen + piminus_gen + pi0_gen).M();
        
    // True quantities
    double e1_true = true_photons[result.pi0_indices[0]][0];
    double e2_true = true_photons[result.pi0_indices[1]][0];
    double e3_true = true_photons[result.prompt_index][0];
    double px1_true = true_photons[result.pi0_indices[0]][1];
    double py1_true = true_photons[result.pi0_indices[0]][2];
    double pz1_true = true_photons[result.pi0_indices[0]][3];
    double px2_true = true_photons[result.pi0_indices[1]][1];
    double py2_true = true_photons[result.pi0_indices[1]][2];
    double pz2_true = true_photons[result.pi0_indices[1]][3];
    double px3_true = true_photons[result.prompt_index][1];
    double py3_true = true_photons[result.prompt_index][2];
    double pz3_true = true_photons[result.prompt_index][3];
    
    e1_bdt_true = e1_true;
    e2_bdt_true = e2_true;
    e3_bdt_true = e3_true;
    
    m_gg_true = compute_invariant_mass(result.pi0_indices[0], result.pi0_indices[1], true_photons);
    m3pi_true = compute_3pi_mass(result.pi0_indices[0], result.pi0_indices[1], true_photons, true_tracks);
    m2pi_true = compute_dipion_mass(true_tracks);
    
    TLorentzVector pi0gam1_bdt_true, pi0gam2_bdt_true;
    pi0gam1_bdt_true.SetPxPyPzE(px1_true, py1_true, pz1_true, e1_true);
    pi0gam2_bdt_true.SetPxPyPzE(px2_true, py2_true, pz2_true, e2_true);
    angle_pi0gam12_bdt_true = pi0gam1_bdt_true.Angle(pi0gam2_bdt_true.Vect()) * TMath::RadToDeg();
    
    TLorentzVector pi0_bdt_true = pi0gam1_bdt_true + pi0gam2_bdt_true;
    betapi0_bdt_true = (pi0_bdt_true.Vect()).Mag() / pi0_bdt_true.E();
    
    // Compute pulls (corrected - true)
    e1_pull = e1_bdt - e1_true;
    e2_pull = e2_bdt - e2_true;
    e3_pull = e3_bdt - e3_true;
    
    px1_pull = px1_bdt - px1_true;
    py1_pull = py1_bdt - py1_true;
    pz1_pull = pz1_bdt - pz1_true;
    
    px2_pull = px2_bdt - px2_true;
    py2_pull = py2_bdt - py2_true;
    pz2_pull = pz2_bdt - pz2_true;
    
    px3_pull = px3_bdt - px3_true;
    py3_pull = py3_bdt - py3_true;
    pz3_pull = pz3_bdt - pz3_true;

    m_gg_pull = m_gg_bdt - m_gg_true;
    m3pi_pull = m3pi_bdt - m3pi_true;
    m2pi_pull = compute_dipion_mass(event.tracks) - m2pi_true;
    
    e_asym = abs(e1_bdt - e2_bdt) / (e1_bdt + e2_bdt);
    
    system_3pi = pi0gam1_bdt + pi0gam2_bdt + trkplus + trkmin;
    beta_3pi = system_3pi.P() / system_3pi.E();
    
    // Reconstruction quality
    int correct_bdt = 0;
    bool checked[2] = {false, false};
    for (int i = 0; i < 2; ++i) {
      if (pho_indx[result.pi0_indices[0]] == EPI0NTMC[i]) {
	correct_bdt++;
	checked[i] = true;
	break;
      }
    }
    for (int i = 0; i < 2; ++i) {
      if (pho_indx[result.pi0_indices[1]] == EPI0NTMC[i] && !checked[i]) {
	correct_bdt++;
	break;
      }
    }
    recon_indx_bdt = correct_bdt;
    
    bool isr_correct = (pho_indx[result.prompt_index] == EPI0NTMC[2] || pho_indx[result.prompt_index] == EPI0NTMC[3]);
    isr_recon_quality = isr_correct ? 1 : 0;
    total_recon_quality = recon_indx_bdt + isr_recon_quality;
    
    // After assigning e1_bdt, e2_bdt, e3_bdt (or after event.photons is updated)
    //Eprompt_max = std::max({e1_bdt, e2_bdt, e3_bdt});
    
    // Selection cuts (deltaE now updated for signal)
    if (lagvalue_min_7C > chi2_cut) continue;
    if (angle_pi0gam12_bdt > angle_cut) continue;
    if (betapi0_bdt > GetFBeta(beta_cut, c0, c1, ppIM)) continue;
    if (deltaE < deltaE_min || deltaE > deltaE_max) continue;
    if (beta_3pi < beta_3pi_min || beta_3pi > beta_3pi_max) continue;
    if (Eprompt_max > Eprompt_max_cut) continue;
    if (bdt_score <= bdt_cut) continue;
    
    //cout << "Eisr = " << Eisr << ", Epi0_pho1 = " << Epi0_pho1 << ", Epi0_pho2 = " << Epi0_pho2 << ", e1_bdt = " << e1_bdt << ", e2_bdt = " << e2_bdt << ", e3_bdt = " << e3_bdt << endl;
    
    
    // Fill output trees
    if (data_type == "exp") {
      TTList[0]->Fill();
    } else if (data_type == "ufo") {
      TTList[8]->Fill();
    } else if (data_type == "eeg") {
      TTList[9]->Fill();
    } else if (data_type == "sig") {
      TTList[10]->Fill();
      if (total_recon_quality == 3) {
	TTList[11]->Fill();
      } else {
	TTList[12]->Fill();
      }
    } else if (data_type == "ksl") {
      if (phid == 0) {
	TTList[1]->Fill();
      } else if (phid == 1) {
	TTList[2]->Fill();
      } else if (phid == 2) {
	TTList[3]->Fill();
      } else if (phid == 3) {
	if (sig_type == 1) TTList[4]->Fill();
	else TTList[5]->Fill();
      } else if (phid == 5) {
	if (sig_type == 1) TTList[6]->Fill();
	else TTList[7]->Fill();
      } else {
	TTList[7]->Fill();
      }
    }
  }
  
  // Write output
  TFile *f_output = new TFile(outputCut + "tree_pre.root", "update");
  f_output->cd();
  
  if (data_type == "exp") {
    TTList[0]->Write();
    cout << "TDATA saved with " << TTList[0]->GetEntries() << " entries" << endl;
  } else if (data_type == "ufo") {
    TTList[8]->Write();
    cout << "TUFO saved with " << TTList[8]->GetEntries() << " entries" << endl;
  } else if (data_type == "eeg") {
    TTList[9]->Write();
    cout << "TEEG saved with " << TTList[9]->GetEntries() << " entries" << endl;
  } else if (data_type == "sig") {
    TTList[10]->Write();
    TTList[11]->Write();
    TTList[12]->Write();
    cout << "TISR3PI_SIG trees saved: "
	 << TTList[10]->GetEntries() << ", "
	 << TTList[11]->GetEntries() << ", "
	 << TTList[12]->GetEntries() << " entries" << endl;
  } else if (data_type == "ksl") {
    TTList[1]->Write();
    TTList[2]->Write();
    TTList[3]->Write();
    TTList[4]->Write();
    TTList[5]->Write();
    TTList[6]->Write();
    TTList[7]->Write();
    cout << "KSL trees saved: "
	 << TTList[1]->GetEntries() << ", "
	 << TTList[2]->GetEntries() << ", "
	 << TTList[3]->GetEntries() << ", "
	 << TTList[4]->GetEntries() << ", "
	 << TTList[5]->GetEntries() << ", "
	 << TTList[6]->GetEntries() << ", "
	 << TTList[7]->GetEntries() << " entries" << endl;
  }
  
  f_output->Close();
  f_input->Close();
  delete tree_list;

  cout << "\n=========================================\n"
       << "Output file: " << outputCut << "tree_pre_bdt.root" << endl;
  cout << "Total events processed: " << evnt_tot << endl;
  
  timer.Stop();
  timer.Print();
  return 0;
}

// ----------------------------------------------------------------------
// Helper function implementations (unchanged)
// ----------------------------------------------------------------------
double compute_invariant_mass(int i, int j, const double photons[3][4]) {
  double E_sum = photons[i][0] + photons[j][0];
  double px_sum = photons[i][1] + photons[j][1];
  double py_sum = photons[i][2] + photons[j][2];
  double pz_sum = photons[i][3] + photons[j][3];
  double mass2 = E_sum*E_sum - (px_sum*px_sum + py_sum*py_sum + pz_sum*pz_sum);
  return (mass2 > 0) ? sqrt(mass2) : 0.0;
}

double compute_3pi_mass(int pi0_idx1, int pi0_idx2, const double photons[3][4], const double tracks[2][4]) {
  double E_sum = photons[pi0_idx1][0] + photons[pi0_idx2][0] + tracks[0][0] + tracks[1][0];
  double px_sum = photons[pi0_idx1][1] + photons[pi0_idx2][1] + tracks[0][1] + tracks[1][1];
  double py_sum = photons[pi0_idx1][2] + photons[pi0_idx2][2] + tracks[0][2] + tracks[1][2];
  double pz_sum = photons[pi0_idx1][3] + photons[pi0_idx2][3] + tracks[0][3] + tracks[1][3];
  double mass2 = E_sum*E_sum - (px_sum*px_sum + py_sum*py_sum + pz_sum*pz_sum);
  return (mass2 > 0) ? sqrt(mass2) : 0.0;
}

double compute_dipion_mass(const double tracks[2][4]) {
  double E_sum = tracks[0][0] + tracks[1][0];
  double px_sum = tracks[0][1] + tracks[1][1];
  double py_sum = tracks[0][2] + tracks[1][2];
  double pz_sum = tracks[0][3] + tracks[1][3];
  double mass2 = E_sum*E_sum - (px_sum*px_sum + py_sum*py_sum + pz_sum*pz_sum);
  return (mass2 > 0) ? sqrt(mass2) : 0.0;
}

double compute_cos_theta(int i, int j, const double photons[3][4]) {
  double px_sum = photons[i][1] + photons[j][1];
  double py_sum = photons[i][2] + photons[j][2];
  double pz_sum = photons[i][3] + photons[j][3];
  double p_mag = sqrt(px_sum*px_sum + py_sum*py_sum + pz_sum*pz_sum);
  if (p_mag < 1e-10) return 0.0;
  return pz_sum / p_mag;
}

double compute_invariant_mass_from_values(double e1, double e2,
                                          double px1, double py1, double pz1,
                                          double px2, double py2, double pz2) {
  double E_sum = e1 + e2;
  double px_sum = px1 + px2;
  double py_sum = py1 + py2;
  double pz_sum = pz1 + pz2;
  double mass2 = E_sum*E_sum - (px_sum*px_sum + py_sum*py_sum + pz_sum*pz_sum);
  return (mass2 > 0) ? sqrt(mass2) : 0.0;
}

std::vector<float> extract_features(int i_idx, int j_idx, int unpaired_idx,
                                    const double photons[3][4], double energy_threshold) {
  std::vector<float> features(10, 0.0f);
  double e1 = photons[i_idx][0];
  double e2 = photons[j_idx][0];
  double e3 = photons[unpaired_idx][0];
  features[5] = (float)e1;
  features[6] = (float)e2;
  features[7] = (float)e3;
  if (e1 >= energy_threshold && e2 >= energy_threshold) {
    double m_gg = compute_invariant_mass(i_idx, j_idx, photons);
    double cos_theta = compute_cos_theta(i_idx, j_idx, photons);
    double opening_angle = acos(std::max(-1.0, std::min(1.0, cos_theta)));
    double denominator = e1 + e2;
    double E_asym = (denominator > 1e-10) ? fabs(e1 - e2) / denominator : 0.0;
    E_asym = std::max(0.0, std::min(1.0, E_asym));
    double e_min_x_angle = std::min(e1, e2) * opening_angle;
    double E_diff = fabs(e1 - e2);
    double asym_x_angle = E_asym * opening_angle;
    features[0] = (float)m_gg;
    features[1] = (float)opening_angle;
    features[2] = (float)cos_theta;
    features[3] = (float)E_asym;
    features[4] = (float)e_min_x_angle;
    features[8] = (float)asym_x_angle;
    features[9] = (float)E_diff;
  }
  return features;
}

BDTResult find_best_pion_pair(const EventData& event, TMVA::Experimental::RBDT& bdt) {
  BDTResult result;
  result.is_valid = false;
  result.scores[0] = 0.0;
  result.scores[1] = 0.0;
  result.scores[2] = 0.0;
  result.max_score = 0.0;
  result.mean_scores = 0.0;
  result.score = 0.0;
  
  int pair_indices[3][2] = {{0,1}, {2,0}, {1,2}};
  double scores[3] = {0.0, 0.0, 0.0};
  
  for (int p = 0; p < 3; ++p) {
    int i_idx = pair_indices[p][0];
    int j_idx = pair_indices[p][1];
    int unpaired_idx = -1;
    for (int k = 0; k < 3; ++k) {
      if (k != i_idx && k != j_idx) { unpaired_idx = k; break; }
    }
    if (unpaired_idx == -1) continue;
    
    std::vector<float> features = extract_features(i_idx, j_idx, unpaired_idx,
						   event.photons, ENERGY_THRESHOLD);
    TMVA::Experimental::RTensor<float> input_tensor(features.data(), {1, features.size()});
    auto bdt_result = bdt.Compute(input_tensor);
    scores[p] = bdt_result(0,0);
  }
  
  result.scores[0] = scores[0];
  result.scores[1] = scores[1];
  result.scores[2] = scores[2];
  
  int best_pair = 0;
  double max_score = scores[0];
  for (int p = 1; p < 3; ++p) {
    if (scores[p] > max_score) {
      max_score = scores[p];
      best_pair = p;
    }
  }
  
  result.max_score = scores[best_pair];
  result.score = result.max_score;
  result.mean_scores = (scores[0] + scores[1] + scores[2]) / 3.0;
  
  result.best_pair_index = best_pair;
  result.pi0_indices[0] = pair_indices[best_pair][0];
  result.pi0_indices[1] = pair_indices[best_pair][1];
  
  result.prompt_index = -1;
  for (int k = 0; k < 3; ++k) {
    if (k != result.pi0_indices[0] && k != result.pi0_indices[1]) {
      result.prompt_index = k; break;
    }
  }
  result.is_valid = (result.prompt_index != -1);
  return result;
}
