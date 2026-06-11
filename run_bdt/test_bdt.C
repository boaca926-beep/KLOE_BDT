#include <TMVA/RBDT.hxx>
#include <TMVA/RTensor.hxx>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TMath.h>
#include <vector>
#include <algorithm>

#include "../header_bdt/helper.h"
#include "../header_bdt/cut_para.h"

using namespace TMVA::Experimental;

// Constants (unchanged)
constexpr double ENERGY_THRESHOLD = 5.0;        // MeV
constexpr int N_BINS_ENERGY = 200;
constexpr int N_BINS_MASS = 200;
constexpr int N_BINS_PULL = 150;
constexpr int N_BINS_CHI2 = 100; 
constexpr int N_BINS_ANGLE = 180;
constexpr int N_BINS_BETA  = 150;
constexpr int N_BINS_BDT_SCORE  = 150;
constexpr double BETA_RANGE_MIN = 0.3;
constexpr double BETA_RANGE_MAX = 1.0;
constexpr double ENERGY_RANGE_MAX = 500.0;      // MeV
constexpr double ENERGY_PHO3_RANGE_MIN = 50.0;  // MeV
constexpr double ENERGY_PHO3_RANGE_MAX = 500.0; // MeV
constexpr double MASS_GG_RANGE_MAX = 200.0;     // MeV/c²
constexpr double MASS_GG_RANGE_MIN = 50.0;      // MeV/c²
constexpr double MASS_3PI_RANGE_MAX = 1000.0;   // MeV/c²
constexpr double MASS_3PI_RANGE_MIN = 400.0;    // MeV/c²
constexpr double MASS_2PI_RANGE_MAX = 700.0;
constexpr double MASS_2PI_RANGE_MIN = 200.0;
constexpr double PULL_RANGE_MIN = -30;          // MeV/c² / MeV
constexpr double PULL_RANGE_MAX = 30;
constexpr double CHI2_RANGE_MAX = 50.0;
constexpr double ANGLE_RANGE_MAX = 180.0;       // deg
constexpr double BDT_SCORE_MAX = 1.0;       

// Event data structure (unchanged)
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
    double score;
    int best_pair_index;
    int pi0_indices[2];
    int prompt_index;
    bool is_valid;
};

// Helper function prototypes (unchanged)
double compute_invariant_mass(int i, int j, const double photons[3][4]);
double compute_3pi_mass(int pi0_idx1, int pi0_idx2, const double photons[3][4], const double tracks[2][4]);
double compute_dipion_mass(const double tracks[2][4]);
double compute_cos_theta(int i, int j, const double photons[3][4]);
std::vector<float> extract_features(int i_idx, int j_idx, int unpaired_idx,
                                    const double photons[3][4], double energy_threshold);
BDTResult find_best_pion_pair(const EventData& event, TMVA::Experimental::RBDT& bdt);

// Histogram manager (unchanged)
class HistogramManager {
private:
    std::vector<TH1D*> histograms;
public:
    ~HistogramManager() { for (auto h : histograms) delete h; }
    TH1D* create(const char* name, const char* title, int nbins, double xmin, double xmax) {
        TH1D* h = new TH1D(name, title, nbins, xmin, xmax);
        h->SetDirectory(0);
        histograms.push_back(h);
        return h;
    }

  
};

void test_bdt() {
    gErrorIgnoreLevel = kError;
    TGaxis::SetMaxDigits(4);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetFitFormat("6.4g");

    // Configuration
    const char* model_filename = "/home/bo/Desktop/KLOE_BDT/models/bdt_pi0_TCOMB.root";
    const char* data_filename = "/home/bo/Desktop/input_bdt_TDATA_chain/cut/tree_pre_bdt.root";
    const char* tree_name = "TISR3PI_SIG";   // or "TETAGAM"
    //const char* tree_name = "TETAGAM";
    //const char* tree_name = "TDATA";
    
    gSystem->mkdir("../plots_test/", kTRUE);

    // Load BDT model
    if (gSystem->AccessPathName(model_filename)) {
        std::cerr << "ERROR: Model file not found: " << model_filename << std::endl;
        return;
    }
    RBDT bdt("BDT_pi0", model_filename);
    std::cout << "✓ BDT model loaded from " << model_filename << std::endl;

    // Open file and tree
    TFile* file = TFile::Open(data_filename);
    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: Cannot open file " << data_filename << std::endl;
        return;
    }
    TTree* tree = (TTree*)file->Get(tree_name);
    if (!tree) {
        std::cerr << "ERROR: Cannot find tree '" << tree_name << "'" << std::endl;
        return;
    }
    std::cout << "Tree " << tree_name << " has " << tree->GetEntries() << " entries." << std::endl;

    HistogramManager hists;
    HistogramManager hists2d;

    TH1D* hM2pi = hists.create("hM2pi", "", N_BINS_MASS, MASS_2PI_RANGE_MIN, MASS_2PI_RANGE_MAX);
    TH1D* hChi2 = hists.create("hChi2", "", N_BINS_CHI2, 0, CHI2_RANGE_MAX);
    
    // --- Histograms for fixed pair (photons 0,1) ---
    TH1D* hE1_fixed = hists.create("hE1_fixed", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hE2_fixed = hists.create("hE2_fixed", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hE3_fixed = hists.create("hE3_fixed", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hMgg_fixed = hists.create("hMgg_fixed", "", N_BINS_MASS, MASS_GG_RANGE_MIN, MASS_GG_RANGE_MAX);
    TH1D* hM3pi_fixed = hists.create("hM3pi_fixed", "", N_BINS_MASS, MASS_3PI_RANGE_MIN, MASS_3PI_RANGE_MAX);
    TH1D* hAngle_fixed = hists.create("hAngle_fixed", "", N_BINS_ANGLE, 0, ANGLE_RANGE_MAX);
    
    // --- Histograms for BDT-selected pair ---
    TH1D* hE1_bdt = hists.create("hE1_bdt", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hE2_bdt = hists.create("hE2_bdt", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hE3_bdt = hists.create("hE3_bdt", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);

    TH1D* hMgg_bdt = hists.create("hMgg_bdt", "", N_BINS_MASS, MASS_GG_RANGE_MIN, MASS_GG_RANGE_MAX);

    TH1D* hM3pi_bdt = hists.create("hM3pi_bdt", "", N_BINS_MASS, MASS_3PI_RANGE_MIN, MASS_3PI_RANGE_MAX);

    //
    TH1D* h1dM3pi_bdt_corr_peak = hists.create("h1dM3pi_bdt_corr_peak", "", N_BINS_MASS, MASS_3PI_RANGE_MIN, MASS_3PI_RANGE_MAX);
    TH1D* h1dM3pi_bdt_corr_non_reson = hists.create("h1dM3pi_bdt_corr_non_reson", "", N_BINS_MASS, MASS_3PI_RANGE_MIN, MASS_3PI_RANGE_MAX);

    //
    TH1D* h1dpullE1_peak = hists.create("h1dpullE1_peak", "", 100, -10, 10);
    TH1D* h1dpullE1_non_reson = hists.create("h1dpullE1_non_reson", "", 100, -10, 10);

    //
    TH1D* h1dbdtscore_peak = hists.create("h1dbdtscore_peak", "", N_BINS_BDT_SCORE, 0, BDT_SCORE_MAX);
    TH1D* h1dbdtscore_non_reson = hists.create("h1dbdtscore_non_reson", "", N_BINS_BDT_SCORE, 0, BDT_SCORE_MAX);
    
    TH1D* hAngle_bdt = hists.create("hAngle_bdt", "", N_BINS_ANGLE, 0, ANGLE_RANGE_MAX);

    // Ediff vs Eisr
    TH2D* hDeltaE_eisr_bdt_peak = new TH2D("hDeltaE_eisr_bdt_peak", "DeltaE vs. E_{isr} (omega peak)",
					   200, -700, -150,
					   N_BINS_ENERGY, ENERGY_PHO3_RANGE_MIN, ENERGY_PHO3_RANGE_MAX);
    
    TH2D* hDeltaE_eisr_bdt_non_reson = new TH2D("hDeltaE_eisr_bdt_non_reson", "DeltaE vs. E_{isr} (non_resonant)",
						200, -700, -150,
						N_BINS_ENERGY, ENERGY_PHO3_RANGE_MIN, ENERGY_PHO3_RANGE_MAX);

    
    // beta0 vs Eisr
    TH2D* hBeta0_eisr_bdt_peak = new TH2D("hBeta0_eisr_bdt_peak", "Beta0 vs. E_{isr} (omega peak)",
					  N_BINS_BETA, 0, BETA_RANGE_MAX,
					  N_BINS_ENERGY, ENERGY_PHO3_RANGE_MIN, ENERGY_PHO3_RANGE_MAX);
    
    TH2D* hBeta0_eisr_bdt_non_reson = new TH2D("hBeta0_eisr_bdt_non_reson", "Beta0 vs. E_{isr} (non-resonant)",
					  N_BINS_BETA, BETA_RANGE_MIN, BETA_RANGE_MAX,
					  N_BINS_ENERGY, ENERGY_PHO3_RANGE_MIN, ENERGY_PHO3_RANGE_MAX);
    
    // angle vs Eisr
    TH2D* hAngle_eisr_bdt_peak = new TH2D("hAngle_eisr_bdt_peak", "Opening Angle vs. E_{isr} (omega peak)",
					  N_BINS_ANGLE, 0, ANGLE_RANGE_MAX,
					  N_BINS_ENERGY, ENERGY_PHO3_RANGE_MIN, ENERGY_PHO3_RANGE_MAX);
    TH2D* hAngle_eisr_bdt_non_reson = new TH2D("hAngle_eisr_bdt_non_reson", "Opening Angle vs. E_{isr} (non-resonant)",
					  N_BINS_ANGLE, 0, ANGLE_RANGE_MAX,
					  N_BINS_ENERGY, ENERGY_PHO3_RANGE_MIN, ENERGY_PHO3_RANGE_MAX);
    

    // m2pi vs Eisr
    TH2D* hM2pi_eisr_bdt_peak = new TH2D("hM2pi_eisr_bdt_peak", "M_{2#pi} vs. E_{isr} (omega peak)",
					 N_BINS_MASS, MASS_2PI_RANGE_MIN, MASS_2PI_RANGE_MAX,
					 N_BINS_ENERGY, ENERGY_PHO3_RANGE_MIN, ENERGY_PHO3_RANGE_MAX);

    TH2D* hM2pi_eisr_bdt_non_reson = new TH2D("hM2pi_eisr_bdt_non_reson", "M_{2#pi} vs. E_{isr} (non-resonant)",
					 N_BINS_MASS, MASS_2PI_RANGE_MIN, MASS_2PI_RANGE_MAX,
					 N_BINS_ENERGY, ENERGY_PHO3_RANGE_MIN, ENERGY_PHO3_RANGE_MAX);
    
    // m3pi_bdt corr
    TH2D* hM3pi_bdt_corr = new TH2D("hM3pi_bdt_corr", "M_{3#pi} Correlation",
                                    N_BINS_MASS, MASS_3PI_RANGE_MIN, MASS_3PI_RANGE_MAX,
                                    N_BINS_MASS, MASS_3PI_RANGE_MIN, MASS_3PI_RANGE_MAX);
    TH2D* hM3pi_bdt_corr_peak = new TH2D("hM3pi_bdt_corr_peak", "M_{3#pi} Correlation (#omega peak)",
                                    N_BINS_MASS, MASS_3PI_RANGE_MIN, MASS_3PI_RANGE_MAX,
                                    N_BINS_MASS, MASS_3PI_RANGE_MIN, MASS_3PI_RANGE_MAX);
    TH2D* hM3pi_bdt_corr_non_reson = new TH2D("hM3pi_bdt_corr_reson", "M_{3#pi} Correlation (Non-resonant)",
                                    N_BINS_MASS, MASS_3PI_RANGE_MIN, MASS_3PI_RANGE_MAX,
                                    N_BINS_MASS, MASS_3PI_RANGE_MIN, MASS_3PI_RANGE_MAX);

    // m3pi true vs generated true
    TH2D* hM3pi_bdt_corr_peak_true = new TH2D("hM3pi_bdt_corr_peak_true", "M_{3#pi} Correlation (#omega peak)",
					      N_BINS_MASS, MASS_3PI_RANGE_MIN, MASS_3PI_RANGE_MAX,
					      N_BINS_MASS, MASS_3PI_RANGE_MIN, MASS_3PI_RANGE_MAX);
    TH2D* hM3pi_bdt_corr_non_reson_true = new TH2D("hM3pi_bdt_corr_reson_true", "M_{3#pi} Correlation (Non-resonant)",
						   N_BINS_MASS, MASS_3PI_RANGE_MIN, MASS_3PI_RANGE_MAX,
						   N_BINS_MASS, MASS_3PI_RANGE_MIN, MASS_3PI_RANGE_MAX);

    // --- Pull histograms (if true branches exist) ---
    TH1D* hE1_pull = hists.create("hE1_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hE2_pull = hists.create("hE2_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hE3_pull = hists.create("hE3_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hMgg_pull = hists.create("hMgg_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hM3pi_pull = hists.create("hM3pi_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hM2pi_pull = hists.create("hM2pi_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);

    // Create output ROOT file
    TFile* outfile = new TFile("../plots_test/test_bdt.root", "RECREATE");
 
    // Branch addresses (unchanged)
    double E1, px1, py1, pz1;
    double E2, px2, py2, pz2;
    double E3, px3, py3, pz3;
    double ppl_E, ppl_px, ppl_py, ppl_pz;
    double pmi_E, pmi_px, pmi_py, pmi_pz;
    double lagvalue_min_7C, deltaE, betapi0, angle_pi0gam12, ppIM;
    double betapi0_bdt;
    double angle;
    double E1_true, px1_true, py1_true, pz1_true;
    double E2_true, px2_true, py2_true, pz2_true;
    double E3_true, px3_true, py3_true, pz3_true;
    double ppl_E_true, ppl_px_true, ppl_py_true, ppl_pz_true;
    double pmi_E_true, pmi_px_true, pmi_py_true, pmi_pz_true;
    double bdt_score = 0;
    double true_m3pi = 0;
    double pull_E3 = 0;
    
    int recon_indx_bdt = 0, recon_indx = 0;
    int bkg_indx = 0;

    tree->SetBranchAddress("Br_recon_indx_bdt", &recon_indx_bdt);
    tree->SetBranchAddress("Br_bkg_indx", &bkg_indx);
    
    tree->SetBranchAddress("Br_E1", &E1);
    tree->SetBranchAddress("Br_px1", &px1);
    tree->SetBranchAddress("Br_py1", &py1);
    tree->SetBranchAddress("Br_pz1", &pz1);
    tree->SetBranchAddress("Br_E2", &E2);
    tree->SetBranchAddress("Br_px2", &px2);
    tree->SetBranchAddress("Br_py2", &py2);
    tree->SetBranchAddress("Br_pz2", &pz2);
    tree->SetBranchAddress("Br_E3", &E3);
    tree->SetBranchAddress("Br_px3", &px3);
    tree->SetBranchAddress("Br_py3", &py3);
    tree->SetBranchAddress("Br_pz3", &pz3);
    tree->SetBranchAddress("Br_ppl_E", &ppl_E);
    tree->SetBranchAddress("Br_ppl_px", &ppl_px);
    tree->SetBranchAddress("Br_ppl_py", &ppl_py);
    tree->SetBranchAddress("Br_ppl_pz", &ppl_pz);
    tree->SetBranchAddress("Br_pmi_E", &pmi_E);
    tree->SetBranchAddress("Br_pmi_px", &pmi_px);
    tree->SetBranchAddress("Br_pmi_py", &pmi_py);
    tree->SetBranchAddress("Br_pmi_pz", &pmi_pz);
    tree->SetBranchAddress("Br_lagvalue_min_7C", &lagvalue_min_7C);
    tree->SetBranchAddress("Br_deltaE", &deltaE);
    tree->SetBranchAddress("Br_betapi0", &betapi0);
    tree->SetBranchAddress("Br_betapi0_bdt", &betapi0_bdt);
    tree->SetBranchAddress("Br_angle_pi0gam12", &angle_pi0gam12);
    tree->SetBranchAddress("Br_angle_pi0gam12_bdt", &angle);
    tree->SetBranchAddress("Br_ppIM", &ppIM);
    tree->SetBranchAddress("Br_bdt_score", &bdt_score);
    tree->SetBranchAddress("Br_true_m3pi", &true_m3pi);
    tree->SetBranchAddress("Br_pull_E3", &pull_E3);
 
    bool hasTrue = (tree->GetBranch("Br_E1_true") != nullptr);
    if (hasTrue) {
        tree->SetBranchAddress("Br_E1_true", &E1_true);
        tree->SetBranchAddress("Br_px1_true", &px1_true);
        tree->SetBranchAddress("Br_py1_true", &py1_true);
        tree->SetBranchAddress("Br_pz1_true", &pz1_true);
        tree->SetBranchAddress("Br_E2_true", &E2_true);
        tree->SetBranchAddress("Br_px2_true", &px2_true);
        tree->SetBranchAddress("Br_py2_true", &py2_true);
        tree->SetBranchAddress("Br_pz2_true", &pz2_true);
        tree->SetBranchAddress("Br_E3_true", &E3_true);
        tree->SetBranchAddress("Br_px3_true", &px3_true);
        tree->SetBranchAddress("Br_py3_true", &py3_true);
        tree->SetBranchAddress("Br_pz3_true", &pz3_true);
        tree->SetBranchAddress("Br_ppl_E_true", &ppl_E_true);
        tree->SetBranchAddress("Br_ppl_px_true", &ppl_px_true);
        tree->SetBranchAddress("Br_ppl_py_true", &ppl_py_true);
        tree->SetBranchAddress("Br_ppl_pz_true", &ppl_pz_true);
        tree->SetBranchAddress("Br_pmi_E_true", &pmi_E_true);
        tree->SetBranchAddress("Br_pmi_px_true", &pmi_px_true);
        tree->SetBranchAddress("Br_pmi_py_true", &pmi_py_true);
        tree->SetBranchAddress("Br_pmi_pz_true", &pmi_pz_true);
    }

    // Event loop (unchanged)
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        EventData event;
        event.photons[0][0] = E1; event.photons[0][1] = px1; event.photons[0][2] = py1; event.photons[0][3] = pz1;
        event.photons[1][0] = E2; event.photons[1][1] = px2; event.photons[1][2] = py2; event.photons[1][3] = pz2;
        event.photons[2][0] = E3; event.photons[2][1] = px3; event.photons[2][2] = py3; event.photons[2][3] = pz3;
        event.tracks[0][0] = ppl_E; event.tracks[0][1] = ppl_px; event.tracks[0][2] = ppl_py; event.tracks[0][3] = ppl_pz;
        event.tracks[1][0] = pmi_E; event.tracks[1][1] = pmi_px; event.tracks[1][2] = pmi_py; event.tracks[1][3] = pmi_pz;
        event.lagvalue_min_7C = lagvalue_min_7C;
        event.deltaE = deltaE;
        event.betapi0 = betapi0;
        event.angle_pi0gam12 = angle_pi0gam12;
        event.ppIM = ppIM;

	double m2pi = compute_dipion_mass(event.tracks);
        hM2pi->Fill(m2pi);
	hChi2->Fill(lagvalue_min_7C);
	
	// Fixed pair
        double m_gg_fixed = compute_invariant_mass(0, 1, event.photons);
        double m3pi_fixed = compute_3pi_mass(0, 1, event.photons, event.tracks);
        hE1_fixed->Fill(event.photons[0][0]);
        hE2_fixed->Fill(event.photons[1][0]);
        hE3_fixed->Fill(event.photons[2][0]);
        hMgg_fixed->Fill(m_gg_fixed);
        hM3pi_fixed->Fill(m3pi_fixed);
        hAngle_fixed->Fill(angle_pi0gam12);
	
        // BDT selection
        BDTResult result = find_best_pion_pair(event, bdt);
        if (!result.is_valid) continue;

        double e1 = event.photons[result.pi0_indices[0]][0];
        double e2 = event.photons[result.pi0_indices[1]][0];
        double e3 = event.photons[result.prompt_index][0];
        double m_gg = compute_invariant_mass(result.pi0_indices[0], result.pi0_indices[1], event.photons);
        double m3pi = compute_3pi_mass(result.pi0_indices[0], result.pi0_indices[1], event.photons, event.tracks);
	
        hE1_bdt->Fill(e1);
        hE2_bdt->Fill(e2);
        hE3_bdt->Fill(e3);
        hMgg_bdt->Fill(m_gg);
        hM3pi_bdt->Fill(m3pi);
        hAngle_bdt->Fill(angle);
	
        if (hasTrue) {
            EventData event_true;
            event_true.photons[0][0] = E1_true; event_true.photons[0][1] = px1_true; event_true.photons[0][2] = py1_true; event_true.photons[0][3] = pz1_true;
            event_true.photons[1][0] = E2_true; event_true.photons[1][1] = px2_true; event_true.photons[1][2] = py2_true; event_true.photons[1][3] = pz2_true;
            event_true.photons[2][0] = E3_true; event_true.photons[2][1] = px3_true; event_true.photons[2][2] = py3_true; event_true.photons[2][3] = pz3_true;
            event_true.tracks[0][0] = ppl_E_true; event_true.tracks[0][1] = ppl_px_true; event_true.tracks[0][2] = ppl_py_true; event_true.tracks[0][3] = ppl_pz_true;
            event_true.tracks[1][0] = pmi_E_true; event_true.tracks[1][1] = pmi_px_true; event_true.tracks[1][2] = pmi_py_true; event_true.tracks[1][3] = pmi_pz_true;

            double e1_true = event_true.photons[result.pi0_indices[0]][0];
            double e2_true = event_true.photons[result.pi0_indices[1]][0];
            double e3_true = event_true.photons[result.prompt_index][0];
            double m_gg_true = compute_invariant_mass(result.pi0_indices[0], result.pi0_indices[1], event_true.photons);
            double m3pi_true = compute_3pi_mass(result.pi0_indices[0], result.pi0_indices[1], event_true.photons, event_true.tracks);
            double m2pi_true = compute_dipion_mass(event_true.tracks);

	    //cout << m3pi << ", " << m3pi_true << endl;
	    //cout << recon_indx_bdt << ", " << bkg_indx << endl;
	    //cout << m2pi << ", " << angle << endl;
	    //cout << betapi0_bdt << endl;
	    //cout << true_m3pi << endl;
	    
	    if (recon_indx_bdt == 2 && bkg_indx == 1) {
	      hM3pi_bdt_corr_peak->Fill(m3pi_true, m3pi);
	      hM3pi_bdt_corr_peak_true->Fill(m3pi_true, true_m3pi);
	      h1dM3pi_bdt_corr_peak->Fill(m3pi);
	      hM2pi_eisr_bdt_peak->Fill(m2pi, e3);
	      hAngle_eisr_bdt_peak->Fill(angle, e3);
	      hBeta0_eisr_bdt_peak->Fill(betapi0_bdt, e3);
	      hDeltaE_eisr_bdt_peak->Fill(deltaE, e3);
	      h1dbdtscore_peak->Fill(bdt_score);
	      h1dpullE1_peak->Fill(pull_E3);
	    }
	    else {
	      hM3pi_bdt_corr_non_reson->Fill(m3pi_true, m3pi);
	      hM3pi_bdt_corr_non_reson_true->Fill(m3pi_true, true_m3pi);
	      h1dM3pi_bdt_corr_non_reson->Fill(m3pi);
	      hM2pi_eisr_bdt_non_reson->Fill(m2pi, e3);
	      hAngle_eisr_bdt_non_reson->Fill(angle, e3);
	      hBeta0_eisr_bdt_non_reson->Fill(betapi0_bdt, e3);
	      hDeltaE_eisr_bdt_non_reson->Fill(deltaE, e3);
	      h1dbdtscore_non_reson->Fill(bdt_score);
	      h1dpullE1_non_reson->Fill(pull_E3);
	    }
	    
            hE1_pull->Fill(e1 - e1_true);
            hE2_pull->Fill(e2 - e2_true);
            hE3_pull->Fill(e3 - e3_true);
            hMgg_pull->Fill(m_gg - m_gg_true);
            hM3pi_pull->Fill(m3pi - m3pi_true);
            hM2pi_pull->Fill(m2pi - m2pi_true);
        }

    }

    // After event loop, before drawing
    if (hasTrue) {
        // Fit the 2D correlation histogram
        TFitResultPtr fitRes = hM3pi_bdt_corr_peak->Fit("pol1", "S");
        if (fitRes.Get() && fitRes->IsValid()) {
            double p0 = fitRes->Parameter(0);
            double p1 = fitRes->Parameter(1);
            double err0 = fitRes->ParError(0);
            double err1 = fitRes->ParError(1);
            double corr = hM3pi_bdt_corr_peak->GetCorrelationFactor();
            std::cout << "\n=== M3pi Correlation (BDT-selected) ===" << std::endl;
            std::cout << "Linear fit: M_reco = p0 + p1 * M_true" << std::endl;
            std::cout << "p0 = " << p0 << " +/- " << err0 << " MeV/c²" << std::endl;
            std::cout << "p1 = " << p1 << " +/- " << err1 << std::endl;
            std::cout << "Pearson correlation coefficient r = " << corr << std::endl;
            std::cout << "Pull RMS = " << hM3pi_pull->GetRMS() << " MeV/c²" << std::endl;
            std::cout << "=======================================\n" << std::endl;
        } else {
            std::cout << "Fit failed!" << std::endl;
        }
    }
    
    // Normalise pull histograms
    auto safeNormalize = [](TH1D* h) { if (h && h->Integral() > 0) h->Scale(1.0 / h->Integral()); };
    if (hasTrue) {
        safeNormalize(hE1_pull); safeNormalize(hE2_pull); safeNormalize(hE3_pull);
        safeNormalize(hMgg_pull); safeNormalize(hM3pi_pull); safeNormalize(hM2pi_pull);
    }

    // --- Drawing functions ---
    auto setHistStyle = [](TH1D* h, Color_t color, const char* xTitle, const char* yTitle) {
        h->SetLineColor(color);
        h->SetLineWidth(1);
        h->SetFillColor(color);
        h->SetFillStyle(3001);
        h->GetXaxis()->SetTitle(xTitle);
        h->GetYaxis()->SetTitle(yTitle);
        h->GetXaxis()->SetTitleSize(0.05);
        h->GetYaxis()->SetTitleSize(0.05);
        h->GetXaxis()->SetLabelSize(0.04);
        h->GetYaxis()->SetLabelSize(0.04);
        h->GetXaxis()->SetTitleOffset(1.3);
	h->GetXaxis()->SetNdivisions(6, kTRUE);
        h->GetXaxis()->CenterTitle();
        h->GetYaxis()->CenterTitle();
    };

    //
    auto drawDoubleCompr = [&](const char* name, const char* title,
			       TH1D* h1, TH1D* h2,
			       const char* xTitle, const char* yTitle,
			       bool logy = false) {
      TCanvas* c = new TCanvas(name, title, 600, 600);
      c->cd();
      gPad->SetBottomMargin(0.14);
      gPad->SetLeftMargin(0.16);
      
      setHistStyle(h1, kRed, xTitle, yTitle);
      setHistStyle(h2, kBlue, xTitle, yTitle);
      double max = std::max(h1->GetMaximum(), h2->GetMaximum());
      if (logy) {
        h1->GetYaxis()->SetRangeUser(0.5, max * 1.5);
        h2->GetYaxis()->SetRangeUser(0.5, max * 1.5);
      } else {
        h1->GetYaxis()->SetRangeUser(0, max * 1.2);
        h2->GetYaxis()->SetRangeUser(0, max * 1.2);
      }
      h1->Draw("HIST");
      h2->Draw("HIST SAME");
      if (logy) gPad->SetLogy();
      TLegend* leg = new TLegend(0.2, 0.85, 0.6, 0.9);
      leg->SetNColumns(2);
      //leg->AddEntry(h2, "#eta peak", "f");
      leg->AddEntry(h2, "#omega peak", "f");
      leg->AddEntry(h1, "Combinatorial", "f");
      leg->Draw();
      c->SaveAs(Form("../plots_test/%s.pdf", name));
      c->Write();
      delete c;
    };
    
    // Modified drawTripleOverlay with logy parameter
    auto drawTripleOverlay = [&](const char* name, const char* title,
				 TH1D* h1_fixed, TH1D* h1_bdt,
				 TH1D* h2_fixed, TH1D* h2_bdt,
				 TH1D* h3_fixed, TH1D* h3_bdt,
				 const char* xTitle1, const char* xTitle2, const char* xTitle3,
				 const char* yTitle, bool logy = false) {
      TCanvas* c = new TCanvas(name, title, 1800, 600);
      c->Divide(3,1);
      for (int i=1; i<=3; ++i) {
        TPad* pad = (TPad*)c->GetPad(i);
        pad->SetBottomMargin(0.14);
        pad->SetLeftMargin(0.16);
      }
      auto drawPad = [&](int iPad, TH1D* h_fixed, TH1D* h_bdt, const char* xtitle) {
        c->cd(iPad);
        setHistStyle(h_fixed, kBlue, xtitle, yTitle);
        setHistStyle(h_bdt, kRed, xtitle, yTitle);
        double max = std::max(h_fixed->GetMaximum(), h_bdt->GetMaximum());
        if (logy) {
	  // For log scale, avoid zero: set lower bound to a small positive number
	  h_fixed->GetYaxis()->SetRangeUser(0.5, max * 1.5);
	  h_bdt->GetYaxis()->SetRangeUser(0.5, max * 1.5);
        } else {
	  h_fixed->GetYaxis()->SetRangeUser(0, max * 1.2);
	  h_bdt->GetYaxis()->SetRangeUser(0, max * 1.2);
        }
        h_fixed->Draw("HIST");
        h_bdt->Draw("HIST SAME");
        if (logy) gPad->SetLogy();
        TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
        leg->AddEntry(h_fixed, "#chi^{2}-selected pair", "f");
        leg->AddEntry(h_bdt, "BDT-selected pair", "f");
        leg->Draw();
      };
      drawPad(1, h1_fixed, h1_bdt, xTitle1);
      drawPad(2, h2_fixed, h2_bdt, xTitle2);
      drawPad(3, h3_fixed, h3_bdt, xTitle3);
      c->SaveAs(Form("../plots_test/%s.pdf", name));
      c->Write();
      delete c;
    };
    
    // drawTriple unchanged (for pulls)
    auto drawTriple = [&](const char* name, const char* title,
                          TH1D* h1, TH1D* h2, TH1D* h3,
                          const char* xTitle1, const char* xTitle2, const char* xTitle3,
                          const char* yTitle, Color_t color = kRed) {
        TCanvas* c = new TCanvas(name, title, 1800, 600);
        c->Divide(3,1);
        for (int i=1; i<=3; ++i) {
            TPad* pad = (TPad*)c->GetPad(i);
            pad->SetBottomMargin(0.14);
            pad->SetLeftMargin(0.16);
        }
        auto drawPad = [&](int iPad, TH1D* h, const char* xtitle) {
            c->cd(iPad);
            setHistStyle(h, color, xtitle, yTitle);
            double max = h->GetMaximum();
            if (max > 0) h->GetYaxis()->SetRangeUser(0, max * 1.2);
            h->Draw("HIST");
            TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
            leg->AddEntry(h, "BDT-selected", "f");
            leg->Draw();
        };
        drawPad(1, h1, xTitle1);
        drawPad(2, h2, xTitle2);
        drawPad(3, h3, xTitle3);
        c->SaveAs(Form("../plots_test/%s.pdf", name));
	c->Write();
        delete c;
    };

    // Drawing function for 2D histogram
    auto draw2D = [&](const char* name, const char* title,
		      TH2D* h2, const char* xTitle, const char* yTitle,
		      bool logz = false,
		      std::vector<double> hlines = {},
		      std::vector<double> vlines = {}) {
      TCanvas* c = new TCanvas(name, title, 800, 800);
      c->SetLeftMargin(0.15);
      c->SetRightMargin(0.15);
      c->SetBottomMargin(0.15);
      h2->GetXaxis()->SetTitle(xTitle);
      h2->GetYaxis()->SetTitle(yTitle);
      h2->GetXaxis()->SetTitleSize(0.05);
      h2->GetYaxis()->SetTitleSize(0.05);
      h2->GetYaxis()->SetTitleOffset(1.3);
      h2->GetXaxis()->SetLabelSize(0.04);
      h2->GetYaxis()->SetLabelSize(0.04);
      h2->GetXaxis()->CenterTitle();
      h2->GetYaxis()->CenterTitle();
      h2->Draw("COLZ");
      if (logz) c->SetLogz();
      c->Update();
      
      // Draw the diagonal fit line (pol1) if it exists
      TF1 *fit = h2->GetFunction("pol1");
      if (fit) {
        fit->SetLineColor(kRed);
        fit->SetLineWidth(2);
        fit->SetLineStyle(1);
        fit->Draw("same");
      }
      
      // Draw horizontal lines
      if (!hlines.empty()) {
        double xmin = gPad->GetUxmin();
        double xmax = gPad->GetUxmax();
	for (double y : hlines) {
	  TLine *line = new TLine(xmin, y, xmax, y);
	  line->SetLineColor(kBlack);
	  line->SetLineWidth(3);
	  //line->SetLineStyle(2);   // dashed
	  line->Draw("same");      // <-- lower‑case "same"
        }
	
      }

      // Draw vertical lines
      if (!vlines.empty()) {
       	double ymin = gPad->GetUymin();
	double ymax = gPad->GetUymax();
	for (double x : vlines) {
	  TLine *line1 = new TLine(x, ymin, x, ymax);
	  line1->SetLineColor(kBlack);
	  line1->SetLineWidth(3);
	  //line1->SetLineStyle(2);   // dashed
	  line1->Draw("same");      // <-- lower‑case "same"
	}
	
      }
      c->SaveAs(Form("../plots_test/%s.pdf", name));
      c->Write();
      delete c;
    };
    
    // Produce plots
    // m3pi_bdt_peak and m3pi_bdt_non_reson comparison

    drawDoubleCompr(Form("m3pi_bdt_compr_%s", tree_name), "3pi Invaraint Mass Comparsion",
		    h1dM3pi_bdt_corr_non_reson, h1dM3pi_bdt_corr_peak,
		    "M_{3#pi} [MeV/c^{2}]", "Entries", false);

    drawDoubleCompr(Form("pullE1_compr_%s", tree_name), "Pull E1 Mass Comparsion",
		    h1dpullE1_non_reson, h1dpullE1_peak,
		    "Pull E_{3}", "Entries", false);
    
    drawDoubleCompr(Form("bdtscore_compr_%s", tree_name), "BDT score Comparsion",
		    h1dbdtscore_non_reson, h1dbdtscore_peak,
		    "BDT value", "Entries", true);
    

    // Photon energies: linear y-axis
    drawTripleOverlay("photon_energies", "Photon Energies",
                      hE1_fixed, hE1_bdt, hE2_fixed, hE2_bdt, hE3_fixed, hE3_bdt,
                      "E_{1} [MeV]", "E_{2} [MeV]", "E_{3} [MeV]", "Entries", true);

    // Kinematic variables: logarithmic y-axis
    drawTripleOverlay("kine_vars", "Kinematic Variables",
                      hMgg_fixed, hMgg_bdt, hM3pi_fixed, hM3pi_bdt, hAngle_fixed, hAngle_bdt,
                      "M_{#gamma#gamma} [MeV/c^{2}]", "M_{3#pi} [MeV/c^{2}]", "#angle_{#gamma#gamma} [#circ]", "Entries", true);

    // M_{2#pi} vs. E_{#gamma_{3}}
    draw2D(Form("angle_eisr_bdt_peak_%s", tree_name), "#angle_{#gamma#gamma} vs. E_{#gamma_{3}} (omega-peak)",
	   hAngle_eisr_bdt_peak, "#angle_{#gamma#gamma} [#circ]", "E_{#gamma_{3}} [MeV]", true, {}, {66.,180.});
    
    draw2D(Form("angle_eisr_bdt_non_reson_%s", tree_name), "#angle_{#gamma#gamma} vs. E_{#gamma_{3}} (Non-resonant)",
	   hAngle_eisr_bdt_non_reson, "#angle_{#gamma#gamma} [#circ]", "E_{#gamma_{3}} [MeV]", true, {}, {66.,180.});

    if (hasTrue) {
        drawTriple("energy_pulls", "Energy Pulls (BDT-selected)",
                   hE1_pull, hE2_pull, hE3_pull,
                   "E_{1} pull [MeV]", "E_{2} pull [MeV]", "E_{3} pull [MeV]",
                   "Normalized entries", kRed);
        drawTriple("mass_pulls", "Mass Pulls (BDT-selected)",
                   hMgg_pull, hM3pi_pull, hM2pi_pull,
                   "M_{#gamma#gamma} pull [MeV/c^{2}]", "M_{3#pi} pull [MeV/c^{2}]", "M_{2#pi} pull [MeV/c^{2}]",
                   "Normalized entries", kRed);
	
	draw2D(Form("m3pi_correlation_peak_%s_true", tree_name), "M_{3#pi} Correlation (BDT-selected)",
               hM3pi_bdt_corr_peak_true, "M^{true, bdt}_{3#pi} [MeV/c^{2}]", "M^{true, gen}_{3#pi} [MeV/c^{2}]", true, {760, 800}, {760, 800});
	draw2D(Form("m3pi_correlation_peak_%s", tree_name), "M_{3#pi} Correlation (BDT-selected)",
               hM3pi_bdt_corr_peak, "M^{true}_{3#pi} [MeV/c^{2}]", "M^{rec}_{3#pi} [MeV/c^{2}]", true, {760, 800}, {760, 800});

	draw2D(Form("m3pi_correlation_non_reso_%s", tree_name), "M_{3#pi} Correlation Non-resonance (BDT-selected)",
               hM3pi_bdt_corr_non_reson, "M^{true}_{3#pi} [MeV/c^{2}]", "M^{rec}_{3#pi} [MeV/c^{2}]", true, {760, 800}, {760, 800});
	draw2D(Form("m3pi_correlation_non_reso_%s_true", tree_name), "M_{3#pi} Correlation Non-resonance (BDT-selected)",
               hM3pi_bdt_corr_non_reson_true, "M^{true, bdt}_{3#pi} [MeV/c^{2}]", "M^{true, gen}_{3#pi} [MeV/c^{2}]", true, {760, 800}, {760, 800});
    }

    // Write histograms to output ROOT file
    outfile->cd();
    hM2pi->Write(); hChi2->Write();
    hE1_fixed->Write(); hE2_fixed->Write(); hE3_fixed->Write();
    hMgg_fixed->Write(); hM3pi_fixed->Write();
    hE1_bdt->Write(); hE2_bdt->Write(); hE3_bdt->Write();
    hMgg_bdt->Write(); hM3pi_bdt->Write(); 
    hM2pi_eisr_bdt_peak->Write(); hM2pi_eisr_bdt_non_reson->Write();
    hAngle_eisr_bdt_peak->Write(); hAngle_eisr_bdt_non_reson->Write();
    hBeta0_eisr_bdt_peak->Write(); hBeta0_eisr_bdt_non_reson->Write();
    hDeltaE_eisr_bdt_peak->Write(); hDeltaE_eisr_bdt_non_reson->Write();
	   
    if (hasTrue) {
      hE1_pull->Write(); hE2_pull->Write(); hE3_pull->Write();
      hMgg_pull->Write(); hM3pi_pull->Write(); hM2pi_pull->Write();
      hM3pi_bdt_corr_peak->Write();
    }
    outfile->Close();
 
    file->Close();
    std::cout << "\nAll plots saved to ../plots_test/" << std::endl;

    gROOT->GetListOfCanvases()->Delete();
    
}

// ---------- Helper function implementations (unchanged) ----------
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
    int best_pair = 0;
    for (int p = 1; p < 3; ++p) if (scores[p] > scores[best_pair]) best_pair = p;
    result.score = scores[best_pair];
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
