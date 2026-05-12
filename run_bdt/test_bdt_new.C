#include <TMVA/RBDT.hxx>
#include <TMVA/RTensor.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TMath.h>
#include <memory>
#include <array>
#include <vector>
#include <string>
#include <algorithm>

#include "../header_bdt/helper.h"
#include "../header_bdt/cut_para.h"

using namespace TMVA::Experimental;

// Constants
constexpr double ENERGY_THRESHOLD = 5.0;        // MeV
constexpr double BDT_CUT_VALUE = 0.4;           // BDT score threshold
constexpr int N_BINS_ENERGY = 200;
constexpr int N_BINS_MASS = 200;
constexpr int N_BINS_CHI2 = 100;
constexpr int N_BINS_PULL = 150;
constexpr double ENERGY_RANGE_MAX = 500.0;      // MeV
constexpr double MASS_RANGE_MAX = 1000.0;       // MeV/c²
constexpr double MASS_GG_RANGE_MAX = 200.0;     // MeV/c²
constexpr double MASS_GG_RANGE_MIN = 50.0;      // MeV/c²
constexpr double CHI2_RANGE_MAX = 0.0;      
constexpr double COS_THETA_RANGE_MIN = -1.0;
constexpr double COS_THETA_RANGE_MAX = 1.0;
constexpr double PULL_RANGE_MIN = -30;          // MeV/c²
constexpr double PULL_RANGE_MAX = 30;           // MeV/c²

constexpr double PI = TMath::Pi();

// Event data structure for better organization
struct EventData {
    double photons[3][4];  // [photon_index][E, px, py, pz]
    double tracks[2][4];   // [track_index][E, px, py, pz]
    double lagvalue_min_7C;
    double deltaE;
    double betapi0;
    double angle_pi0gam12;
    double ppIM;
    int bkg_indx;
    int recon_indx;
};

// BDT result structure
struct BDTResult {
    double score;
    int best_pair_index;
    int pi0_indices[2];
    int prompt_index;
    bool is_valid;
};

// Histogram manager to ensure proper cleanup
class HistogramManager {
private:
    std::vector<TH1D*> histograms;
    
public:
    HistogramManager() = default;  // Add default constructor
    
    ~HistogramManager() {
        for (auto hist : histograms) {
            delete hist;
        }
    }
    
    TH1D* create(const char* name, const char* title, int nbins, double xmin, double xmax) {
        TH1D* hist = new TH1D(name, title, nbins, xmin, xmax);
        histograms.push_back(hist);
        return hist;
    }
    
    // Prevent copying
    HistogramManager(const HistogramManager&) = delete;
    HistogramManager& operator=(const HistogramManager&) = delete;
};

// Helper function prototypes
double compute_invariant_mass(int i, int j, const double photons[3][4]);
double compute_3pi_mass(int pi0_idx1, int pi0_idx2, const double photons[3][4], const double tracks[2][4]);
double compute_cos_theta(int i, int j, const double photons[3][4]);
double get_fbeta(double a1, double b1, double c1, double m2pi);
double compute_dipion_mass(const double tracks[2][4]);
std::vector<float> extract_features(int i_idx, int j_idx, int unpaired_idx, 
                                    const double photons[3][4], double energy_threshold);
BDTResult find_best_pion_pair(const EventData& event, TMVA::Experimental::RBDT& bdt);

/**
 * Apply BDT model to KLOE detector data for η→3π⁰ decay analysis
 * @param model_filename Path to the BDT model file
 * @param data_filename Path to the input ROOT data file
 * 
 * The function selects the best γγ pair among the three photons
 * using a BDT trained on kinematic variables, then applies both
 * KLOE χ² selection and BDT selection for comparison.
 */
void test_bdt_new() {

    // ROOT Style settings
    gErrorIgnoreLevel = kError;
    TGaxis::SetMaxDigits(4);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetFitFormat("6.4g");

    const char* model_filename = "/home/kloe/Desktop/KLOE_BDT/models/bdt_pi0_TCOMB.root";
    //const char* data_filename = "/home/kloe/Desktop/KLOE_BDT/dataset/kloe_bdt_norm.root";
    const char* data_filename = "/home/kloe/Desktop/input_bdt_TDATA_chain/cut/tree_pre.root";
  
    
    // Manually load libraries
    gSystem->Load("libTMVA");
    gSystem->Load("libTMVAUtils");
    
    std::cout << "Testing model_KLOE ..." << std::endl;

    // 1. Validate model file exists
    if (gSystem->AccessPathName(model_filename)) {
        std::cerr << "ERROR: Model file " << model_filename << " does not exist!" << std::endl;
        return;
    }

    // 2. Load the BDT model
    std::cout << "Loading model from " << model_filename << std::endl;
    TMVA::Experimental::RBDT bdt("BDT_pi0", model_filename);
    
    std::cout << "✓ Model loaded successfully!" << std::endl;

    //const TString phys_ch[2] = {"TDATA", "Data"};
    //const TString phys_ch[2] = {"TETAGAM", "Etagam"};
    const TString phys_ch[2] = {"TISR3PI_SIG", "Signal"};
    
    const TString ch_nm = phys_ch[0];
    const TString ch_type = phys_ch[1];

    // 3. Initialize histograms using manager
    HistogramManager hists;
    
    // Energy histograms (for checking all 3 photons)
    TH1D* he1 = hists.create("he1", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* he2 = hists.create("he2", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* he3 = hists.create("he3", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    
    // Kinematic histograms (for debugging)
    TH1D* hcos_theta = hists.create("hcos_theta", "", N_BINS_ENERGY, COS_THETA_RANGE_MIN, COS_THETA_RANGE_MAX);
    TH1D* hopen_angle = hists.create("hopen_angle", "", N_BINS_ENERGY, 0, PI);
    TH1D* hE_asym = hists.create("hE_asym", "", N_BINS_ENERGY, 0, 1);
    TH1D* he_min_x_angle = hists.create("he_min_x_angle", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hE_diff = hists.create("hE_diff", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hasym_x_angle = hists.create("hasym_x_angle", "", N_BINS_ENERGY, 0, PI);
    TH1D* hM2pi_pull = hists.create("hM2pi_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    
    // Histograms for BDT selection
    TH1D* hchi2_BDT_good = hists.create("hchi2_BDT_good", "", N_BINS_CHI2, 0, CHI2_RANGE_MAX);
    TH1D* hchi2_BDT_bad = hists.create("hchi2_BDT_bad", "", N_BINS_CHI2, 0, CHI2_RANGE_MAX);
    TH1D* hchi2_BDT = hists.create("hchi2_BDT", "", N_BINS_CHI2, 0, CHI2_RANGE_MAX);
    
    TH1D* hM3pi_BDT_good = hists.create("hM3pi_BDT_good", "", N_BINS_MASS, 400, MASS_RANGE_MAX);
    TH1D* hM3pi_BDT_bad = hists.create("hM3pi_BDT_bad", "", N_BINS_MASS, 400, MASS_RANGE_MAX);
    TH1D* hM3pi_BDT = hists.create("hM3pi_BDT", "", N_BINS_MASS, 400, MASS_RANGE_MAX);

    TH1D* hM_gg_BDT_good = hists.create("hM_gg_BDT_good", "", N_BINS_ENERGY, MASS_GG_RANGE_MIN, MASS_GG_RANGE_MAX);
    TH1D* hM_gg_BDT_bad = hists.create("hM_gg_BDT_bad", "", N_BINS_ENERGY, MASS_GG_RANGE_MIN, MASS_GG_RANGE_MAX);
    TH1D* hM_gg_BDT = hists.create("hM_gg_BDT", "", N_BINS_ENERGY, MASS_GG_RANGE_MIN, MASS_GG_RANGE_MAX);
    
    TH1D* hE1_BDT_good = hists.create("hE1_BDT_good", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hE1_BDT_bad = hists.create("hE1_BDT_bad", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hE1_BDT = hists.create("hE1_BDT", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    
    TH1D* hE2_BDT_good = hists.create("hE2_BDT_good", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hE2_BDT_bad = hists.create("hE2_BDT_bad", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hE2_BDT = hists.create("hE2_BDT", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    
    TH1D* hE3_BDT_good = hists.create("hE3_BDT_good", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hE3_BDT_bad = hists.create("hE3_BDT_bad", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hE3_BDT = hists.create("hE3_BDT", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);

    TH1D* hE1_pull_good = hists.create("hE1_pull_good", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hE1_pull_bad = hists.create("hE1_pull_bad", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hE1_pull = hists.create("hE1_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);

    TH1D* hE2_pull_good = hists.create("hE2_pull_good", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hE2_pull_bad = hists.create("hE2_pull_bad", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hE2_pull = hists.create("hE2_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);

    TH1D* hE3_pull_good = hists.create("hE3_pull_good", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hE3_pull_bad = hists.create("hE3_pull_bad", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hE3_pull = hists.create("hE3_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);

    // px pulls
    TH1D* hPx1_pull_good = hists.create("hPx1_pull_good", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPx1_pull_bad  = hists.create("hPx1_pull_bad", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPx1_pull  = hists.create("hPx1_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    
    TH1D* hPx2_pull_good = hists.create("hPx2_pull_good", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPx2_pull_bad  = hists.create("hPx2_pull_bad", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPx2_pull  = hists.create("hPx2_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    
    TH1D* hPx3_pull_good = hists.create("hPx3_pull_good", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPx3_pull_bad  = hists.create("hPx3_pull_bad", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPx3_pull  = hists.create("hPx3_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    
    // py pulls
    TH1D* hPy1_pull_good = hists.create("hPy1_pull_good", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPy1_pull_bad  = hists.create("hPy1_pull_bad", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPy1_pull  = hists.create("hPy1_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    
    TH1D* hPy2_pull_good = hists.create("hPy2_pull_good", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPy2_pull_bad  = hists.create("hPy2_pull_bad", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPy2_pull  = hists.create("hPy2_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    
    TH1D* hPy3_pull_good = hists.create("hPy3_pull_good", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPy3_pull_bad  = hists.create("hPy3_pull_bad", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPy3_pull  = hists.create("hPy3_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    
    // pz pulls
    TH1D* hPz1_pull_good = hists.create("hPz1_pull_good", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPz1_pull_bad  = hists.create("hPz1_pull_bad", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPz1_pull  = hists.create("hPz1_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    
    TH1D* hPz2_pull_good = hists.create("hPz2_pull_good", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPz2_pull_bad  = hists.create("hPz2_pull_bad", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPz2_pull  = hists.create("hPz2_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    
    TH1D* hPz3_pull_good = hists.create("hPz3_pull_good", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPz3_pull_bad  = hists.create("hPz3_pull_bad", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPz3_pull  = hists.create("hPz3_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);

    //
    TH1D* hM3pi_pull_good = hists.create("hM3pi_pull_good", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hM3pi_pull_bad = hists.create("hM3pi_pull_bad", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hM3pi_pull = hists.create("hM3pi_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    
    TH1D* hM_gg_pull_good = hists.create("hM_gg_pull_good", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hM_gg_pull_bad = hists.create("hM_gg_pull_bad", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hM_gg_pull = hists.create("hM_gg_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    
    // Mass histograms for KLOE selection (using fixed pair 0,1)
    TH1D* hM_gg = hists.create("hM_gg", "", N_BINS_ENERGY, MASS_GG_RANGE_MIN, MASS_GG_RANGE_MAX);
    TH1D* hM_gg_good = hists.create("hM_gg_good", "", N_BINS_ENERGY, MASS_GG_RANGE_MIN, MASS_GG_RANGE_MAX);
    TH1D* hM_gg_bad = hists.create("hM_gg_bad", "", N_BINS_ENERGY, MASS_GG_RANGE_MIN, MASS_GG_RANGE_MAX);
    
    TH1D* hM3pi = hists.create("hM3pi", "", N_BINS_MASS, 400, MASS_RANGE_MAX);
    TH1D* hM3pi_good = hists.create("hM3pi_good", "", N_BINS_MASS, 400, MASS_RANGE_MAX);
    TH1D* hM3pi_bad = hists.create("hM3pi_bad", "", N_BINS_MASS, 400, MASS_RANGE_MAX);
    
    TH1D* hE1 = hists.create("hE1", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hE1_good = hists.create("hE1_good", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hE1_bad = hists.create("hE1_bad", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    
    TH1D* hE2 = hists.create("hE2", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hE2_good = hists.create("hE2_good", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hE2_bad = hists.create("hE2_bad", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    
    // Event counters
    int evnt_KLOE = 0;
    int evnt_good = 0;
    int evnt_bad = 0;
    int n_found = 0;
    int n_discarded = 0;

    // 4. Process data file if it exists
    if (!gSystem->AccessPathName(data_filename)) {
        std::cout << "\nProcessing data file: " << data_filename << std::endl;

        // Open ROOT file
        TFile* file = TFile::Open(data_filename);
        if (!file || file->IsZombie()) {
            std::cerr << "ERROR: Cannot open file " << data_filename << std::endl;
            return;
        }

        // Get the tree
        std::cout << "Reading tree: " << ch_nm << std::endl;
        TTree* tree = (TTree*)file->Get(ch_nm);
        
        if (!tree) {
            std::cerr << "ERROR: Cannot find tree '" << ch_nm << "' in file" << std::endl;
            file->Close();
            return;
        }

        int nentries = tree->GetEntries();
        std::cout << "Tree has " << nentries << " entries" << std::endl;

        // Branch variables
        double lagvalue_min_7C = 0., deltaE = 0., betapi0 = 0., angle_pi0gam12 = 0.;
        double ppIM = 0.;
        double E1 = 0., px1 = 0., py1 = 0., pz1 = 0.;
        double E2 = 0., px2 = 0., py2 = 0., pz2 = 0.;
        double E3 = 0., px3 = 0., py3 = 0., pz3 = 0.;
	double E1_true = 0., px1_true = 0., py1_true = 0., pz1_true = 0.;
	double E2_true = 0., px2_true = 0., py2_true = 0., pz2_true = 0.;
  	double E3_true = 0., px3_true = 0., py3_true = 0., pz3_true = 0.;
  	double ppl_E = 0., ppl_px = 0., ppl_py = 0., ppl_pz = 0.;
	double pmi_E = 0., pmi_px = 0., pmi_py = 0., pmi_pz = 0.;
	double ppl_E_true = 0., ppl_px_true = 0., ppl_py_true = 0., ppl_pz_true = 0.;
	double pmi_E_true = 0., pmi_px_true = 0., pmi_py_true = 0., pmi_pz_true = 0.;
        int bkg_indx = 0, recon_indx = 0;

        // Set branch addresses
        tree->SetBranchAddress("Br_deltaE", &deltaE);
        tree->SetBranchAddress("Br_angle_pi0gam12", &angle_pi0gam12);
        tree->SetBranchAddress("Br_betapi0", &betapi0);
        tree->SetBranchAddress("Br_lagvalue_min_7C", &lagvalue_min_7C);
        tree->SetBranchAddress("Br_ppIM", &ppIM);
        tree->SetBranchAddress("Br_bkg_indx", &bkg_indx);
        tree->SetBranchAddress("Br_recon_indx", &recon_indx);
        tree->SetBranchAddress("Br_ppl_E", &ppl_E);
        tree->SetBranchAddress("Br_ppl_px", &ppl_px);
        tree->SetBranchAddress("Br_ppl_py", &ppl_py);
        tree->SetBranchAddress("Br_ppl_pz", &ppl_pz);
        tree->SetBranchAddress("Br_pmi_E", &pmi_E);
        tree->SetBranchAddress("Br_pmi_px", &pmi_px);
        tree->SetBranchAddress("Br_pmi_py", &pmi_py);
        tree->SetBranchAddress("Br_pmi_pz", &pmi_pz);
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
        
        // Create output file and tree
        TFile* outfile = TFile::Open("output_with_bdt.root", "RECREATE");
        TTree* outtree = new TTree("new_tree", "Tree with BDT response");
        
        double bdt_score = 0.0;
        double m_gg_bdt = 0.0, m3pi_bdt = 0.0;
	double m_gg_bdt_true = 0.0, m3pi_bdt_true = 0.0;
	double m_gg_pull = 0.0, m3pi_pull = 0.0;
	double m2pi = 0.0, m2pi_true = 0., m2pi_pull = 0.0;
	double e1_bdt = 0.0, e2_bdt = 0.0, e3_bdt = 0.0;
	double px1_bdt = 0.0, px2_bdt = 0.0, px3_bdt = 0.0;
	double py1_bdt = 0.0, py2_bdt = 0.0, py3_bdt = 0.0;
	double pz1_bdt = 0.0, pz2_bdt = 0.0, pz3_bdt = 0.0;
	double e1_bdt_true = 0.0, e2_bdt_true = 0.0, e3_bdt_true = 0.0;
	double px1_bdt_true = 0.0, px2_bdt_true = 0.0, px3_bdt_true = 0.0;
	double py1_bdt_true = 0.0, py2_bdt_true = 0.0, py3_bdt_true = 0.0;
	double pz1_bdt_true = 0.0, pz2_bdt_true = 0.0, pz3_bdt_true = 0.0;
	double e1_pull = 0.0, e2_pull = 0.0, e3_pull = 0.0;
	double px1_pull = 0.0, px2_pull = 0.0, px3_pull = 0.0;
	double py1_pull = 0.0, py2_pull = 0.0, py3_pull = 0.0;
	double pz1_pull = 0.0, pz2_pull = 0.0, pz3_pull = 0.0;
	
        int event_id = 0;
        
        outtree->Branch("event_id", &event_id);
        outtree->Branch("bdt_score", &bdt_score);
        outtree->Branch("m_gg_bdt", &m_gg_bdt);
	outtree->Branch("m_gg_bdt_true", &m_gg_bdt_true);
	outtree->Branch("m_gg_pull", &m_gg_pull);
        outtree->Branch("m3pi_bdt", &m3pi_bdt);
	outtree->Branch("m3pi_bdt_true", &m3pi_bdt_true);
	outtree->Branch("m3pi_pull", &m3pi_pull);
	outtree->Branch("m2pi", &m2pi);
	outtree->Branch("m2pi_true", &m2pi_true);
	outtree->Branch("m2pi_pull", &m2pi_pull);
        outtree->Branch("e1_bdt", &e1_bdt);
        outtree->Branch("e2_bdt", &e2_bdt);
        outtree->Branch("e3_bdt", &e3_bdt);
	outtree->Branch("e1_bdt_true", &e1_bdt_true);
	outtree->Branch("e2_bdt_true", &e2_bdt_true);
        outtree->Branch("e3_bdt_true", &e3_bdt_true);
        outtree->Branch("e1_pull", &e1_pull);
	outtree->Branch("e2_pull", &e2_pull);
	outtree->Branch("e3_pull", &e3_pull);
	outtree->Branch("px1_pull", &px1_pull);
	outtree->Branch("px2_pull", &px2_pull);
	outtree->Branch("px3_pull", &px3_pull);
	outtree->Branch("py1_pull", &py1_pull);
	outtree->Branch("py2_pull", &py2_pull);
	outtree->Branch("py3_pull", &py3_pull);
	outtree->Branch("pz1_pull", &pz1_pull);
	outtree->Branch("pz2_pull", &pz2_pull);
	outtree->Branch("pz3_pull", &pz3_pull);

        // Main event loop
        for (int i = 0; i < nentries; i++) {
            tree->GetEntry(i);

	    // Test variables
	    //cout << E1_true << ", " << px1_true << ", " << py1_true << ", " << pz1_true << endl;
	    //cout << lagvalue_min_7C << endl;
	    //cout << ppl_px_true << ", " << ppl_py_true << ", " << ppl_pz_true << endl;
	    //cout << pmi_px_true << ", " << pmi_py_true << ", " << pmi_pz_true << endl;
	    
            // Apply kinematic cuts
            if (lagvalue_min_7C > chi2_cut) continue;
            if (deltaE > deltaE_cut) continue;
            if (angle_pi0gam12 > angle_cut) continue;
            if (betapi0 > get_fbeta(beta_cut, c0, c1, ppIM)) continue;

            // Clean data - remove NaN values
            if (TMath::IsNaN(E1) || TMath::IsNaN(E2) || TMath::IsNaN(E3)) continue;
            if (TMath::IsNaN(ppl_E) || TMath::IsNaN(pmi_E)) continue;

            // Fill energy histograms for all three photons (for debugging)
            he1->Fill(E1);
            he2->Fill(E2);
            he3->Fill(E3);

            // Prepare event data structure
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
            event.bkg_indx = bkg_indx;
            event.recon_indx = recon_indx;

	    // Prepare event photon true data structure
	    EventData event_true;
            event_true.photons[0][0] = E1_true; event_true.photons[0][1] = px1_true; event_true.photons[0][2] = py1_true; event_true.photons[0][3] = pz1_true;
            event_true.photons[1][0] = E2_true; event_true.photons[1][1] = px2_true; event_true.photons[1][2] = py2_true; event_true.photons[1][3] = pz2_true;
            event_true.photons[2][0] = E3_true; event_true.photons[2][1] = px3_true; event_true.photons[2][2] = py3_true; event_true.photons[2][3] = pz3_true;
            event_true.tracks[0][0] = ppl_E_true; event_true.tracks[0][1] = ppl_px_true; event_true.tracks[0][2] = ppl_py_true; event_true.tracks[0][3] = ppl_pz_true;
            event_true.tracks[1][0] = pmi_E_true; event_true.tracks[1][1] = pmi_px_true; event_true.tracks[1][2] = pmi_py_true; event_true.tracks[1][3] = pmi_pz_true;
            
	    // Compute using fixed pair (0,1) for KLOE comparison
            double m_gg_fixed = compute_invariant_mass(0, 1, event.photons);
            double m3pi_fixed = compute_3pi_mass(0, 1, event.photons, event.tracks);
            m2pi = compute_dipion_mass(event.tracks);
	    m2pi_true = compute_dipion_mass(event_true.tracks);
	    m2pi_pull = m2pi - m2pi_true;
	    
            // Fill histograms for fixed pair
            hM_gg->Fill(m_gg_fixed);
            hM3pi->Fill(m3pi_fixed);
            hE1->Fill(event.photons[0][0]);
            hE2->Fill(event.photons[1][0]);
	    hM2pi_pull->Fill(m2pi_pull);
	      
            // Find best pair using BDT
            BDTResult result = find_best_pion_pair(event, bdt);
            
            if (!result.is_valid) continue;

            // Extract best pair information
            e1_bdt = event.photons[result.pi0_indices[0]][0];
            e2_bdt = event.photons[result.pi0_indices[1]][0];
            e3_bdt = event.photons[result.prompt_index][0];

	    px1_bdt = event.photons[result.pi0_indices[0]][1];
	    py1_bdt = event.photons[result.pi0_indices[0]][2];
	    pz1_bdt = event.photons[result.pi0_indices[0]][3];

	    px2_bdt = event.photons[result.pi0_indices[1]][1];
	    py2_bdt = event.photons[result.pi0_indices[1]][2];
	    pz2_bdt = event.photons[result.pi0_indices[1]][3];

	    px3_bdt = event.photons[result.prompt_index][1];
	    py3_bdt = event.photons[result.prompt_index][2];
	    pz3_bdt = event.photons[result.prompt_index][3];
 
	    e1_bdt_true = event_true.photons[result.pi0_indices[0]][0];
            e2_bdt_true = event_true.photons[result.pi0_indices[1]][0];
            e3_bdt_true = event_true.photons[result.prompt_index][0];

	    px1_bdt_true = event_true.photons[result.pi0_indices[0]][1];
	    py1_bdt_true = event_true.photons[result.pi0_indices[0]][2];
	    pz1_bdt_true = event_true.photons[result.pi0_indices[0]][3];

	    px2_bdt_true = event_true.photons[result.pi0_indices[1]][1];
	    py2_bdt_true = event_true.photons[result.pi0_indices[1]][2];
	    pz2_bdt_true = event_true.photons[result.pi0_indices[1]][3];
	    
	    px3_bdt_true = event_true.photons[result.prompt_index][1];
	    py3_bdt_true = event_true.photons[result.prompt_index][2];
	    pz3_bdt_true = event_true.photons[result.prompt_index][3];

	    e1_pull = e1_bdt - e1_bdt_true;
	    e2_pull = e2_bdt - e2_bdt_true;
	    e3_pull = e3_bdt - e3_bdt_true;

	    px1_pull = px1_bdt - px1_bdt_true;
	    py1_pull = py1_bdt - py1_bdt_true;
	    pz1_pull = pz1_bdt - pz1_bdt_true;

	    px2_pull = px2_bdt - px2_bdt_true;
	    py2_pull = py2_bdt - py2_bdt_true;
	    pz2_pull = pz2_bdt - pz2_bdt_true;

	    px3_pull = px3_bdt - px3_bdt_true;
	    py3_pull = py3_bdt - py3_bdt_true;
	    pz3_pull = pz3_bdt - pz3_bdt_true;

	    m_gg_bdt = compute_invariant_mass(result.pi0_indices[0], result.pi0_indices[1], event.photons);
	    m_gg_bdt_true = compute_invariant_mass(result.pi0_indices[0], result.pi0_indices[1], event_true.photons);

            m3pi_bdt = compute_3pi_mass(result.pi0_indices[0], result.pi0_indices[1], event.photons, event.tracks);
	    m_gg_pull = m_gg_bdt - m_gg_bdt_true;
	    
	    m3pi_bdt_true = compute_3pi_mass(result.pi0_indices[0], result.pi0_indices[1], event_true.photons, event_true.tracks);
	    m3pi_pull = m3pi_bdt - m3pi_bdt_true;
            bdt_score = result.score;
            event_id = i;

            // Fill output tree
            outtree->Fill();

            // BDT selection
            if (bdt_score > BDT_CUT_VALUE) {
                n_found++;
		hchi2_BDT_good->Fill(lagvalue_min_7C);
                hE1_BDT_good->Fill(e1_bdt);
                hE2_BDT_good->Fill(e2_bdt);
                hE3_BDT_good->Fill(e3_bdt);
                hM_gg_BDT_good->Fill(m_gg_bdt);
                hM3pi_BDT_good->Fill(m3pi_bdt);

		hE1_pull_good->Fill(e1_pull);
		hE2_pull_good->Fill(e2_pull);
		hE3_pull_good->Fill(e3_pull);
		
		hPx1_pull_good->Fill(px1_pull);
		hPx2_pull_good->Fill(px2_pull);
		hPx3_pull_good->Fill(px3_pull);

		hPy1_pull_good->Fill(py1_pull);
		hPy2_pull_good->Fill(py2_pull);
		hPy3_pull_good->Fill(py3_pull);

		hPz1_pull_good->Fill(pz1_pull);
		hPz2_pull_good->Fill(pz2_pull);
		hPz3_pull_good->Fill(pz3_pull);

		hM3pi_pull_good->Fill(m3pi_pull);
		hM_gg_pull_good->Fill(m_gg_pull);
            } else {
                n_discarded++;
		hchi2_BDT_bad->Fill(lagvalue_min_7C);
                hE1_BDT_bad->Fill(e1_bdt);
                hE2_BDT_bad->Fill(e2_bdt);
                hE3_BDT_bad->Fill(e3_bdt);
                hM_gg_BDT_bad->Fill(m_gg_bdt);
                hM3pi_BDT_bad->Fill(m3pi_bdt);
            }

            // Fill total histograms
	    hchi2_BDT->Fill(lagvalue_min_7C);
            hE1_BDT->Fill(e1_bdt);
            hE2_BDT->Fill(e2_bdt);
            hE3_BDT->Fill(e3_bdt);
            hM_gg_BDT->Fill(m_gg_bdt);
	    hM3pi_BDT->Fill(m3pi_bdt);
		
            evnt_KLOE++;
        }

        // Write histograms and close files
        outfile->cd();
        //outtree->Write();
        
        //hM3pi->Write();
        //hM3pi_good->Write();
        //hM3pi_bad->Write();
        hM3pi_BDT_good->Write();
        hM3pi_BDT_bad->Write();
        hM3pi_BDT->Write();
        
        outfile->Close();
        file->Close();

        std::cout << "\n✓ pi0 finding complete!" << std::endl;
        std::cout << "  Total events after cuts: " << evnt_KLOE << std::endl;
        std::cout << "  KLOE selected: " << evnt_good << std::endl;
        std::cout << "  KLOE discarded: " << evnt_bad << std::endl;
        std::cout << "  BDT selected (score > " << BDT_CUT_VALUE << "): " << n_found << std::endl;
        std::cout << "  BDT discarded (score < " << BDT_CUT_VALUE << "): " << n_discarded << std::endl;
    }

    // 5. Make comparison plots
    if (evnt_KLOE > 0) {
        // Original raw counts canvas
        TCanvas* cv0 = new TCanvas("c1", "BDT Selection (" + ch_nm + ")", 2000, 1200);
        //cv0->SetLeftMargin(0.1);
        //cv0->SetBottomMargin(0.1);
        cv0->Divide(3, 2);
        
        // Energy of gamma 1: E1
        cv0->cd(1);
        double ymax_e1 = hE1_BDT->GetBinContent(hE1_BDT->GetMaximumBin());
        double ymax_e2 = hE2_BDT->GetBinContent(hE2_BDT->GetMaximumBin());
        double ymax_e3 = hE3_BDT->GetBinContent(hE3_BDT->GetMaximumBin());
        
        TPaveText* pt1 = new TPaveText(0.11, 0.87, 0.80, 0.89, "NDC");
        PteAttr(pt1);
        pt1->SetTextSize(0.03);
        pt1->SetTextColor(kBlack);
        pt1->AddText(Form("Events=%d, BDT Cut Value=%.1f, BDT Selected=%d, Discarded=%d", evnt_KLOE, BDT_CUT_VALUE, n_found, evnt_KLOE - n_found));
        
        format_h(hE1, 1, 2);
	format_h(hE1_BDT, 1, 1);
        //format_h(hE1_good, 4, 1);
        //format_h(hE1_bad, 2, 1);
        formatfill_h(hE1_BDT_good, 3, 3001);
        formatfill_h(hE1_BDT_bad, 2, 3001);

	hE1_BDT_good->SetLineColor(kGreen+2);
        hE1_BDT_good->SetFillColorAlpha(kGreen+1, 0.4);
        hE1_BDT_bad->SetLineColor(kRed+2);
        hE1_BDT_bad->SetFillColorAlpha(kRed+1, 0.4);
        
        hE1_BDT->GetYaxis()->SetNdivisions(505);
        hE1_BDT->GetYaxis()->SetRangeUser(0.1, ymax_e1 * 1.2);
	hE1_BDT->GetYaxis()->SetTitle("Entries");
	hE1_BDT->GetYaxis()->CenterTitle();
	hE1_BDT->GetYaxis()->SetTitleSize(0.03);
        
	hE1_BDT->GetXaxis()->SetTitle("E_{1} [MeV]");
        hE1_BDT->GetXaxis()->CenterTitle();
        hE1_BDT->GetXaxis()->SetTitleSize(0.04);
        
        hE1_BDT->Draw();
        hE1_BDT_good->Draw("Same");
        hE1_BDT_bad->Draw("Same");
        pt1->Draw("Same");
        
        TLegend* legd_cv = new TLegend(0.6, 0.6, 0.9, 0.85);
        legd_cv->SetTextFont(132);
        legd_cv->SetFillStyle(0);
        legd_cv->SetBorderSize(0);
        legd_cv->SetNColumns(1);
        legd_cv->AddEntry(hE1_BDT, "BDT Sum", "l");
        legd_cv->AddEntry(hE1_BDT_good, "BDT Selected", "f");
        legd_cv->AddEntry(hE1_BDT_bad, "BDT Discarded", "f");
        legd_cv->Draw("Same");
        legtextsize(legd_cv, 0.04);
        
        // Energy of gamma 2: E2
        cv0->cd(2);
        format_h(hE2, 1, 2);
        format_h(hE2_BDT, 1, 1);
        formatfill_h(hE2_BDT_good, 3, 3001);
        formatfill_h(hE2_BDT_bad, 2, 3001);

	hE2_BDT_good->SetLineColor(kGreen+2);
        hE2_BDT_good->SetFillColorAlpha(kGreen+1, 0.4);
        hE2_BDT_bad->SetLineColor(kRed+2);
        hE2_BDT_bad->SetFillColorAlpha(kRed+1, 0.4);
        
        hE2_BDT->GetYaxis()->SetNdivisions(505);
        hE2_BDT->GetYaxis()->SetRangeUser(0.1, ymax_e2 * 1.2);
	hE2_BDT->GetYaxis()->SetTitle("Entries");
	hE2_BDT->GetYaxis()->SetTitleSize(0.03);
        hE2_BDT->GetYaxis()->CenterTitle();
        hE2_BDT->GetXaxis()->SetTitle("E_{2} [MeV]");
        hE2_BDT->GetXaxis()->CenterTitle();
        hE2_BDT->GetXaxis()->SetTitleSize(0.04);
        
        hE2_BDT->Draw();
        hE2_BDT_good->Draw("Same");
        hE2_BDT_bad->Draw("Same");
        
        // Energy of gamma3: E3
        cv0->cd(3);
        format_h(hE3_BDT, 1, 1);
        formatfill_h(hE3_BDT_good, 3, 3001);
        formatfill_h(hE3_BDT_bad, 2, 3001);

	hE3_BDT_good->SetLineColor(kGreen+2);
        hE3_BDT_good->SetFillColorAlpha(kGreen+1, 0.4);
        hE3_BDT_bad->SetLineColor(kRed+2);
        hE3_BDT_bad->SetFillColorAlpha(kRed+1, 0.4);
        
        hE3_BDT->GetYaxis()->SetNdivisions(505);
        hE3_BDT->GetYaxis()->SetRangeUser(0.1, ymax_e3 * 1.2);
	hE3_BDT->GetYaxis()->SetTitle("Entries");
	hE3_BDT->GetYaxis()->CenterTitle();
	hE3_BDT->GetYaxis()->SetTitleSize(0.03);
        hE3_BDT->GetXaxis()->SetTitle("E_{3} [MeV]");
        hE3_BDT->GetXaxis()->CenterTitle();
        hE3_BDT->GetXaxis()->SetTitleSize(0.04);
        
        hE3_BDT->Draw();
        hE3_BDT_good->Draw("Same");
        hE3_BDT_bad->Draw("Same");
        
        // Invariant mass of gamma gamma: M_gg
        cv0->cd(4);
        double ymax_m_gg = hM_gg_BDT->GetBinContent(hM_gg_BDT->GetMaximumBin());
        
        format_h(hM_gg_BDT, 1, 1);
        formatfill_h(hM_gg_BDT_good, 3, 3001);
        formatfill_h(hM_gg_BDT_bad, 2, 3001);

	hM_gg_BDT_good->SetLineColor(kGreen+2);
        hM_gg_BDT_good->SetFillColorAlpha(kGreen+1, 0.4);
        hM_gg_BDT_bad->SetLineColor(kRed+2);
        hM_gg_BDT_bad->SetFillColorAlpha(kRed+1, 0.4);
        
        hM_gg_BDT->GetYaxis()->SetNdivisions(505);
        hM_gg_BDT->GetYaxis()->SetRangeUser(0.1, ymax_m_gg * 1.2);
	hM_gg_BDT->GetYaxis()->SetTitle("Entries");
	hM_gg_BDT->GetYaxis()->CenterTitle();
	hM_gg_BDT->GetXaxis()->SetTitleSize(0.03);
        hM_gg_BDT->GetXaxis()->SetTitle("M(#gamma#gamma) [MeV/c^{2}]");
        hM_gg_BDT->GetXaxis()->CenterTitle();
        hM_gg_BDT->GetXaxis()->SetTitleSize(0.04);
        
        hM_gg_BDT->Draw();
        hM_gg_BDT_good->Draw("Same");
        hM_gg_BDT_bad->Draw("Same");
	gPad->SetLogy(1);
        
        // Invariant mass of 3pi:  M_3pi
        cv0->cd(5);
        double ymax_m3pi = hM3pi_BDT->GetBinContent(hM3pi_BDT->GetMaximumBin());
        
        format_h(hM3pi, 1, 2);
	format_h(hM3pi_BDT, 1, 1);
        formatfill_h(hM3pi_BDT_good, 3, 3001);
        formatfill_h(hM3pi_BDT_bad, 2, 3001);

	hM3pi_BDT_good->SetLineColor(kGreen+2);
        hM3pi_BDT_good->SetFillColorAlpha(kGreen+1, 0.4);
        hM3pi_BDT_bad->SetLineColor(kRed+2);
        hM3pi_BDT_bad->SetFillColorAlpha(kRed+1, 0.4);
        
        hM3pi_BDT->GetYaxis()->SetNdivisions(505);
        hM3pi_BDT->GetYaxis()->SetRangeUser(0.1, ymax_m3pi * 1.5);
	hM3pi_BDT->GetYaxis()->SetTitle("Entries");
	hM3pi_BDT->GetYaxis()->CenterTitle();
	hM3pi_BDT->GetYaxis()->SetTitleSize(0.03);
        hM3pi_BDT->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
        hM3pi_BDT->GetXaxis()->CenterTitle();
        hM3pi_BDT->GetXaxis()->SetTitleSize(0.04);
        
        hM3pi_BDT->Draw();
        hM3pi_BDT_good->Draw("Same");
        hM3pi_BDT_bad->Draw("Same");
        gPad->SetLogy(1);

	// Chi2
        cv0->cd(6);
        double ymax_chi2 = hchi2_BDT->GetBinContent(hchi2_BDT->GetMaximumBin());
        
        format_h(hchi2_BDT, 1, 1);
        formatfill_h(hchi2_BDT_good, 3, 3001);
        formatfill_h(hchi2_BDT_bad, 2, 3001);

	hchi2_BDT_good->SetLineColor(kGreen+2);
        hchi2_BDT_good->SetFillColorAlpha(kGreen+1, 0.4);
        hchi2_BDT_bad->SetLineColor(kRed+2);
        hchi2_BDT_bad->SetFillColorAlpha(kRed+1, 0.4);
        
        hchi2_BDT->GetYaxis()->SetNdivisions(505);
        hchi2_BDT->GetYaxis()->SetRangeUser(0.1, ymax_chi2 * 1.5);
	hchi2_BDT->GetYaxis()->SetTitle("Entries");
	hchi2_BDT->GetYaxis()->CenterTitle();
	hchi2_BDT->GetYaxis()->SetTitleSize(0.03);
        hchi2_BDT->GetXaxis()->SetTitle("#chi^{2}");
        hchi2_BDT->GetXaxis()->CenterTitle();
        hchi2_BDT->GetXaxis()->SetTitleSize(0.04);
        
        hchi2_BDT->Draw();
        hchi2_BDT_good->Draw("Same");
        hchi2_BDT_bad->Draw("Same");
        //gPad->SetLogy(1);
        
        cv0->SaveAs(Form("../plots_bdt/bdt_spectra_%s.pdf", ch_nm.Data()));
        delete cv0;
        delete pt1;
        delete legd_cv;

	// ==================================================================
	// Normalized pull histograms
	// ==================================================================
	auto safeNormalize = [](TH1D* h) {
	  if (h && h->Integral() > 0) h->Scale(1.0 / h->Integral());
	};

	TH1D* hE1_pull_good_norm = (TH1D*)hE1_pull_good->Clone("hE1_pull_good_norm");
	TH1D* hE2_pull_good_norm = (TH1D*)hE2_pull_good->Clone("hE2_pull_good_norm");
	TH1D* hE3_pull_good_norm = (TH1D*)hE3_pull_good->Clone("hE3_pull_good_norm");
	TH1D* hPx1_pull_good_norm = (TH1D*)hPx1_pull_good->Clone("hPx1_pull_good_norm");
	TH1D* hPx2_pull_good_norm = (TH1D*)hPx2_pull_good->Clone("hPx2_pull_good_norm");
	TH1D* hPx3_pull_good_norm = (TH1D*)hPx3_pull_good->Clone("hPx3_pull_good_norm");
	TH1D* hPy1_pull_good_norm = (TH1D*)hPy1_pull_good->Clone("hPy1_pull_good_norm");
	TH1D* hPy2_pull_good_norm = (TH1D*)hPy2_pull_good->Clone("hPy2_pull_good_norm");
	TH1D* hPy3_pull_good_norm = (TH1D*)hPy3_pull_good->Clone("hPy3_pull_good_norm");
	TH1D* hPz1_pull_good_norm = (TH1D*)hPz1_pull_good->Clone("hPz1_pull_good_norm");
	TH1D* hPz2_pull_good_norm = (TH1D*)hPz2_pull_good->Clone("hPz2_pull_good_norm");
	TH1D* hPz3_pull_good_norm = (TH1D*)hPz3_pull_good->Clone("hPz3_pull_good_norm");
	TH1D* hM3pi_pull_good_norm = (TH1D*)hM3pi_pull_good->Clone("hM3pi_pull_good_norm");
	TH1D* hM_gg_pull_good_norm = (TH1D*)hM_gg_pull_good->Clone("hM_gg_pull_good_norm");
 
	safeNormalize(hE1_pull_good_norm);
	safeNormalize(hE2_pull_good_norm);
	safeNormalize(hE3_pull_good_norm);
	safeNormalize(hPx1_pull_good_norm);
	safeNormalize(hPx2_pull_good_norm);
	safeNormalize(hPx3_pull_good_norm);
	safeNormalize(hPy1_pull_good_norm);
	safeNormalize(hPy2_pull_good_norm);
	safeNormalize(hPy3_pull_good_norm);
	safeNormalize(hPz1_pull_good_norm);
	safeNormalize(hPz2_pull_good_norm);
	safeNormalize(hPz3_pull_good_norm);
	safeNormalize(hM3pi_pull_good_norm);
	safeNormalize(hM_gg_pull_good_norm);

	// ==================================================================
	// Improved layout: separate canvases for each variable type
	// ==================================================================

	// Common style settings
	auto setHistStyle = [](TH1D* h, Color_t color, const char* xTitle, const char* yTitle) {
	  h->SetLineColor(color);
	  h->SetLineWidth(1);
	  h->SetFillColorAlpha(color, 0.3);
	  h->GetYaxis()->SetTitle(yTitle);
	  h->GetXaxis()->SetTitle(xTitle);
	  h->GetYaxis()->SetTitleSize(0.05);
	  h->GetXaxis()->SetTitleSize(0.05);
	  h->GetYaxis()->SetLabelSize(0.045);
	  h->GetXaxis()->SetLabelSize(0.045);
	  h->GetXaxis()->SetTitleOffset(1.3);
	  h->GetXaxis()->CenterTitle();
	  h->GetYaxis()->CenterTitle();
	};
	
	auto drawSet = [&](const char* canvasName, const char* canvasTitle, const char* lgdTitle,
			   TH1D* h1, TH1D* h2, TH1D* h3,
			   const char* xTitle1, const char* xTitle2, const char* xTitle3,
			   const char* yTitle1, const char* yTitle2, const char* yTitle3,
			   double yMaxFactor = 1.2) {
	  TCanvas* c = new TCanvas(canvasName, canvasTitle, 1800, 600);
	  c->Divide(3,1);
	  
	  double max1 = h1->GetMaximum();
	  double max2 = h2->GetMaximum();
	  double max3 = h3->GetMaximum();
	  //double yMax = TMath::Max(max1, TMath::Max(max2, max3)) * yMaxFactor;
	  
	  // Set margins for each pad individually
	  for (int iPad = 1; iPad <= 3; ++iPad) {
	    TPad* pad = (TPad*)c->GetPad(iPad);
	    pad->SetBottomMargin(0.22);   // more space for x‑title
	    pad->SetLeftMargin(0.15);     // more space for y‑title
	  }

	  TPaveText* pt2 = new TPaveText(0.2, 0.80, 0.80, 0.89, "NDC");
	  PteAttr(pt2);
	  pt2->SetTextSize(0.05);
	  pt2->SetTextColor(kBlack);
	  pt2->AddText("Signal");
        
	  c->cd(1);
	  setHistStyle(h1, kGreen+2, xTitle1, yTitle1);
	  h1->GetYaxis()->SetRangeUser(0, h1->GetBinContent(h1->GetMaximumBin()) * yMaxFactor);
	  h1->Draw("HIST");
	  pt2->Draw("Same");
	  TLegend* leg1 = new TLegend(0.65, 0.8, 0.9, 0.9);
	  leg1->SetBorderSize(0);
	  leg1->SetFillStyle(0);
	  leg1->AddEntry(h1, lgdTitle, "f");
	  leg1->Draw();
	  
	  c->cd(2);
	  setHistStyle(h2, kGreen+2, xTitle2, yTitle2);
	  h2->GetYaxis()->SetRangeUser(0, h2->GetBinContent(h2->GetMaximumBin()) * yMaxFactor);
	  h2->Draw("HIST");
	  //TLegend* leg2 = new TLegend(0.65, 0.75, 0.9, 0.9);
	  //leg2->SetBorderSize(0);
	  //leg2->SetFillStyle(0);
	  //leg2->AddEntry(h2, "BDT Selected", "f");
	  //leg2->Draw();
	  
	  c->cd(3);
	  setHistStyle(h3, kGreen+2, xTitle3, yTitle3);
	  h3->GetYaxis()->SetRangeUser(0, h3->GetBinContent(h3->GetMaximumBin()) * yMaxFactor);
	  h3->Draw("HIST");
	  //TLegend* leg3 = new TLegend(0.65, 0.75, 0.9, 0.9);
	  //leg3->SetBorderSize(0);
	  //leg3->SetFillStyle(0);
	  //leg3->AddEntry(h3, "BDT Selected", "f");
	  //leg3->Draw();
	  
	  c->SaveAs(Form("../plots_bdt/%s.pdf", canvasName));
	  delete c;
	};

	// Draw Energy pulls
	drawSet("energy_good_pulls", "Energy Pulls (BDT Selected)", "BDT Selected",
		hE1_pull_good_norm, hE2_pull_good_norm, hE3_pull_good_norm,
		"E_{1} pull [MeV]", "E_{2} pull [MeV]", "E_{3} pull [MeV]",
		"Normalized entries", "Normalized entries", "Normalized entries");
	
	// Draw Px pulls
	drawSet("px_good_pulls", "P_{x} Pulls (BDT Selected)", "BDT Selected",
		hPx1_pull_good_norm, hPx2_pull_good_norm, hPx3_pull_good_norm,
		"p_{x,1} pull [MeV/c]", "p_{x,2} pull [MeV/c]", "p_{x,3} pull [MeV/c]",
		"Normalized entries", "Normalized entries", "Normalized entries");
	
	// Draw Py pulls
	drawSet("py_good_pulls", "P_{y} Pulls (BDT Selected)", "BDT Selected",
		hPy1_pull_good_norm, hPy2_pull_good_norm, hPy3_pull_good_norm,
		"p_{y,1} pull [MeV/c]", "p_{y,2} pull [MeV/c]", "p_{y,3} pull [MeV/c]",
		"Normalized entries", "Normalized entries", "Normalized entries");
	
	// Draw Pz pulls
	drawSet("pz_good_pulls", "P_{z} Pulls (BDT Selected)", "BDT Selected",
		hPz1_pull_good_norm, hPz2_pull_good_norm, hPz3_pull_good_norm,
		"p_{z,1} pull [MeV/c]", "p_{z,2} pull [MeV/c]", "p_{z,3} pull [MeV/c]",
		"Normalized entries", "Normalized entries", "Normalized entries");

	// Draw m3pi pulls
	drawSet("mass_good_pulls", "M_{3#pi} Pulls (BDT Selected)", "BDT Selected",
		hM3pi_pull_good_norm, hM_gg_pull_good_norm, hM2pi_pull,
		"M_{3#pi} pull [MeV/c^{2}]", "M_{#gamma#gamma} pull [MeV/c^{2}]", "M_{2#pi} pull [MeV/c^{2}]",
		"Normalized Entries", "Normalized Entries", "Entries");
    }
}

// Helper function implementations

double get_fbeta(double a1, double b1, double c1, double m2pi) {
    m2pi = m2pi / 1000.;
    double fbeta = a1 + 1. / (exp((m2pi - c1) / b1) - 1.);
    return fbeta;
}

double compute_dipion_mass(const double tracks[2][4]) {
    double E_sum = tracks[0][0] + tracks[1][0];
    double px_sum = tracks[0][1] + tracks[1][1];
    double py_sum = tracks[0][2] + tracks[1][2];
    double pz_sum = tracks[0][3] + tracks[1][3];
    double mass2 = E_sum*E_sum - (px_sum*px_sum + py_sum*py_sum + pz_sum*pz_sum);
    return (mass2 > 0) ? sqrt(mass2) : 0.0;
}

double compute_invariant_mass(int i, int j, const double photons[3][4]) {
    double E_sum = photons[i][0] + photons[j][0];
    double px_sum = photons[i][1] + photons[j][1];
    double py_sum = photons[i][2] + photons[j][2];
    double pz_sum = photons[i][3] + photons[j][3];
    
    double mass2 = E_sum * E_sum - (px_sum * px_sum + py_sum * py_sum + pz_sum * pz_sum);
    return (mass2 > 0) ? TMath::Sqrt(mass2) : 0.0;
}

double compute_3pi_mass(int pi0_idx1, int pi0_idx2, const double photons[3][4], const double tracks[2][4]) {
    // Sum pi0 and track momenta
    double E_sum = photons[pi0_idx1][0] + photons[pi0_idx2][0] + tracks[0][0] + tracks[1][0];
    double px_sum = photons[pi0_idx1][1] + photons[pi0_idx2][1] + tracks[0][1] + tracks[1][1];
    double py_sum = photons[pi0_idx1][2] + photons[pi0_idx2][2] + tracks[0][2] + tracks[1][2];
    double pz_sum = photons[pi0_idx1][3] + photons[pi0_idx2][3] + tracks[0][3] + tracks[1][3];
    
    double mass2 = E_sum * E_sum - (px_sum * px_sum + py_sum * py_sum + pz_sum * pz_sum);
    return (mass2 > 0) ? TMath::Sqrt(mass2) : 0.0;
}

double compute_cos_theta(int i, int j, const double photons[3][4]) {
    double px_sum = photons[i][1] + photons[j][1];
    double py_sum = photons[i][2] + photons[j][2];
    double pz_sum = photons[i][3] + photons[j][3];
    double p_mag = TMath::Sqrt(px_sum*px_sum + py_sum*py_sum + pz_sum*pz_sum);
    
    if (p_mag < 1e-10) return 0.0;
    return pz_sum / p_mag;
}

std::vector<float> extract_features(int i_idx, int j_idx, int unpaired_idx,
                                    const double photons[3][4], double energy_threshold) {
    std::vector<float> features(10, 0.0f);
    
    double e1 = photons[i_idx][0];
    double e2 = photons[j_idx][0];
    double e3 = photons[unpaired_idx][0];
    
    // Default values (will be overwritten if energy threshold is met)
    features[5] = (float)e1;
    features[6] = (float)e2;
    features[7] = (float)e3;
    
    if (e1 >= energy_threshold && e2 >= energy_threshold) {
        double m_gg = compute_invariant_mass(i_idx, j_idx, photons);
        double cos_theta = compute_cos_theta(i_idx, j_idx, photons);
        double opening_angle = TMath::ACos(cos_theta);
        double denominator = e1 + e2;
        double E_asym = (denominator > 1e-10) ? TMath::Abs(e1 - e2) / denominator : 0.0;
        E_asym = TMath::Max(0.0, TMath::Min(1.0, E_asym));
        double e_min_x_angle = TMath::Min(e1, e2) * opening_angle;
        double E_diff = TMath::Abs(e1 - e2);
        double asym_x_angle = E_asym * opening_angle;
        
        features[0] = (float)m_gg;
        features[1] = (float)opening_angle;
        features[2] = (float)cos_theta;
        features[3] = (float)E_asym;
        features[4] = (float)e_min_x_angle;
        features[5] = (float)e1;          // E1
        features[6] = (float)e2;          // E2
        features[7] = (float)e3;          // E3
        features[8] = (float)asym_x_angle;
        features[9] = (float)E_diff;
    }
    
    return features;
}

BDTResult find_best_pion_pair(const EventData& event, TMVA::Experimental::RBDT& bdt) {
    BDTResult result;
    result.is_valid = false;
    
    // All 3 possible pairs of photons
    int pair_indices[3][2] = {{0, 1}, {2, 0}, {1, 2}};
    double scores[3] = {0.0, 0.0, 0.0};
    
    for (int p = 0; p < 3; p++) {
        int i_idx = pair_indices[p][0];
        int j_idx = pair_indices[p][1];
        
        // Find the unpaired photon
        int unpaired_idx = -1;
        for (int k = 0; k < 3; k++) {
            if (k != i_idx && k != j_idx) {
                unpaired_idx = k;
                break;
            }
        }
        
        if (unpaired_idx == -1) continue;
        
        // Extract features and compute BDT score
        std::vector<float> features = extract_features(i_idx, j_idx, unpaired_idx, 
                                                       event.photons, ENERGY_THRESHOLD);
        
        TMVA::Experimental::RTensor<float> input_tensor(features.data(), {1, features.size()});
        auto bdt_result = bdt.Compute(input_tensor);
        scores[p] = bdt_result(0, 0);
    }
    
    // Find best pair (highest score)
    int best_pair = 0;
    for (int p = 1; p < 3; p++) {
        if (scores[p] > scores[best_pair]) best_pair = p;
    }
    
    result.score = scores[best_pair];
    result.best_pair_index = best_pair;
    result.pi0_indices[0] = pair_indices[best_pair][0];
    result.pi0_indices[1] = pair_indices[best_pair][1];
    
    // Find prompt photon (the one not in the best pair)
    result.prompt_index = -1;
    for (int k = 0; k < 3; k++) {
        if (k != result.pi0_indices[0] && k != result.pi0_indices[1]) {
            result.prompt_index = k;
            break;
        }
    }
    
    result.is_valid = (result.prompt_index != -1);
    return result;
}
