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
constexpr double ENERGY_RANGE_MAX = 500.0;      // MeV
constexpr double MASS_RANGE_MAX = 1000.0;       // MeV/c²
constexpr double MASS_GG_RANGE_MAX = 200.0;     // MeV/c²
constexpr double MASS_GG_RANGE_MIN = 50.0;      // MeV/c²
constexpr double CHI2_RANGE_MAX = 0.0;      
constexpr double COS_THETA_RANGE_MIN = -1.0;
constexpr double COS_THETA_RANGE_MAX = 1.0;
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

    const TString phys_ch[2] = {"TDATA", "Data"};
    //const TString phys_ch[2] = {"TETAGAM", "Etagam"};
    //const TString phys_ch[2] = {"TISR3PI_SIG", "Signal"};
    
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
        double ppl_E = 0., ppl_px = 0., ppl_py = 0., ppl_pz = 0.;
        double pmi_E = 0., pmi_px = 0., pmi_py = 0., pmi_pz = 0.;
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

        // Create output file and tree
        TFile* outfile = TFile::Open("output_with_bdt.root", "RECREATE");
        TTree* outtree = new TTree("new_tree", "Tree with BDT response");
        
        double bdt_score = 0.0;
        double m_gg_bdt = 0.0, m3pi_bdt = 0.0;
        double e1_bdt = 0.0, e2_bdt = 0.0, e3_bdt = 0.0;
        int event_id = 0;
        
        outtree->Branch("event_id", &event_id);
        outtree->Branch("bdt_score", &bdt_score);
        outtree->Branch("m_gg_bdt", &m_gg_bdt);
        outtree->Branch("m3pi_bdt", &m3pi_bdt);
        outtree->Branch("e1_bdt", &e1_bdt);
        outtree->Branch("e2_bdt", &e2_bdt);
        outtree->Branch("e3_bdt", &e3_bdt);

        // Main event loop
        for (int i = 0; i < nentries; i++) {
            tree->GetEntry(i);

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

            // Compute using fixed pair (0,1) for KLOE comparison
            double m_gg_fixed = compute_invariant_mass(0, 1, event.photons);
            double m3pi_fixed = compute_3pi_mass(0, 1, event.photons, event.tracks);
            
            // Fill histograms for fixed pair
            hM_gg->Fill(m_gg_fixed);
            hM3pi->Fill(m3pi_fixed);
            hE1->Fill(event.photons[0][0]);
            hE2->Fill(event.photons[1][0]);

            // Find best pair using BDT
            BDTResult result = find_best_pion_pair(event, bdt);
            
            if (!result.is_valid) continue;

            // Extract best pair information
            e1_bdt = event.photons[result.pi0_indices[0]][0];
            e2_bdt = event.photons[result.pi0_indices[1]][0];
            e3_bdt = event.photons[result.prompt_index][0];
            
            m_gg_bdt = compute_invariant_mass(result.pi0_indices[0], result.pi0_indices[1], event.photons);
            m3pi_bdt = compute_3pi_mass(result.pi0_indices[0], result.pi0_indices[1], event.photons, event.tracks);
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
        outtree->Write();
        
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

	/*
        // ==================== NEW: Normalized (differential) shapes ====================
        // Create normalized clones
        TH1D* hE1_good_norm = (TH1D*)hE1_BDT_good->Clone("hE1_good_norm");
        TH1D* hE1_bad_norm  = (TH1D*)hE1_BDT_bad->Clone("hE1_bad_norm");
        hE1_good_norm->Scale(1.0 / hE1_good_norm->Integral());
        hE1_bad_norm->Scale(1.0 / hE1_bad_norm->Integral());

        TH1D* hE2_good_norm = (TH1D*)hE2_BDT_good->Clone("hE2_good_norm");
        TH1D* hE2_bad_norm  = (TH1D*)hE2_BDT_bad->Clone("hE2_bad_norm");
        hE2_good_norm->Scale(1.0 / hE2_good_norm->Integral());
        hE2_bad_norm->Scale(1.0 / hE2_bad_norm->Integral());

        TH1D* hE3_good_norm = (TH1D*)hE3_BDT_good->Clone("hE3_good_norm");
        TH1D* hE3_bad_norm  = (TH1D*)hE3_BDT_bad->Clone("hE3_bad_norm");
        hE3_good_norm->Scale(1.0 / hE3_good_norm->Integral());
        hE3_bad_norm->Scale(1.0 / hE3_bad_norm->Integral());

        TH1D* hMgg_good_norm = (TH1D*)hM_gg_BDT_good->Clone("hMgg_good_norm");
        TH1D* hMgg_bad_norm  = (TH1D*)hM_gg_BDT_bad->Clone("hMgg_bad_norm");
        hMgg_good_norm->Scale(1.0 / hMgg_good_norm->Integral());
        hMgg_bad_norm->Scale(1.0 / hMgg_bad_norm->Integral());

        TH1D* hM3pi_good_norm = (TH1D*)hM3pi_BDT_good->Clone("hM3pi_good_norm");
        TH1D* hM3pi_bad_norm  = (TH1D*)hM3pi_BDT_bad->Clone("hM3pi_bad_norm");
        hM3pi_good_norm->Scale(1.0 / hM3pi_good_norm->Integral());
        hM3pi_bad_norm->Scale(1.0 / hM3pi_bad_norm->Integral());

        // New canvas for normalized comparison
        TCanvas* cv_norm = new TCanvas("cv_norm", "Normalized BDT Comparison", 1800, 1200);
        cv_norm->Divide(2,3);
        
        // E1
        cv_norm->cd(1);
        hE1_good_norm->SetLineColor(kGreen+2);
        hE1_good_norm->SetFillColorAlpha(kGreen+1, 0.4);
        hE1_bad_norm->SetLineColor(kRed+2);
        hE1_bad_norm->SetFillColorAlpha(kRed+1, 0.4);
        hE1_good_norm->GetYaxis()->SetTitle("Normalized entries");
        hE1_good_norm->GetXaxis()->SetTitle("E_{1} [MeV]");
        hE1_good_norm->Draw("HIST");
        hE1_bad_norm->Draw("HIST SAME");
        TLegend* leg1 = new TLegend(0.7,0.7,0.9,0.9);
        leg1->AddEntry(hE1_good_norm, "BDT Selected", "f");
        leg1->AddEntry(hE1_bad_norm, "BDT Discarded", "f");
        leg1->Draw();
        
        // E2
        cv_norm->cd(2);
        hE2_good_norm->GetYaxis()->SetTitle("Normalized entries");
        hE2_good_norm->GetXaxis()->SetTitle("E_{2} [MeV]");
        hE2_good_norm->Draw("HIST");
        hE2_bad_norm->Draw("HIST SAME");
        TLegend* leg2 = new TLegend(0.7,0.7,0.9,0.9);
        leg2->AddEntry(hE2_good_norm, "BDT Selected", "f");
        leg2->AddEntry(hE2_bad_norm, "BDT Discarded", "f");
        leg2->Draw();
        
        // E3
        cv_norm->cd(3);
        hE3_good_norm->GetYaxis()->SetTitle("Normalized entries");
        hE3_good_norm->GetXaxis()->SetTitle("E_{3} [MeV]");
        hE3_good_norm->Draw("HIST");
        hE3_bad_norm->Draw("HIST SAME");
        TLegend* leg3 = new TLegend(0.7,0.7,0.9,0.9);
        leg3->AddEntry(hE3_good_norm, "BDT Selected", "f");
        leg3->AddEntry(hE3_bad_norm, "BDT Discarded", "f");
        leg3->Draw();
        
        // M_gg
        cv_norm->cd(4);
        hMgg_good_norm->GetYaxis()->SetTitle("Normalized entries");
        hMgg_good_norm->GetXaxis()->SetTitle("M(#gamma#gamma) [MeV/c^{2}]");
        hMgg_good_norm->Draw("HIST");
        hMgg_bad_norm->Draw("HIST SAME");
        TLegend* leg4 = new TLegend(0.7,0.7,0.9,0.9);
        leg4->AddEntry(hMgg_good_norm, "BDT Selected", "f");
        leg4->AddEntry(hMgg_bad_norm, "BDT Discarded", "f");
        leg4->Draw();
        
        // M_3pi
        cv_norm->cd(5);
        hM3pi_good_norm->GetYaxis()->SetTitle("Normalized entries");
        hM3pi_good_norm->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
        hM3pi_good_norm->Draw("HIST");
        hM3pi_bad_norm->Draw("HIST SAME");
        TLegend* leg5 = new TLegend(0.7,0.7,0.9,0.9);
        leg5->AddEntry(hM3pi_good_norm, "BDT Selected", "f");
        leg5->AddEntry(hM3pi_bad_norm, "BDT Discarded", "f");
        leg5->Draw();
        
        cv_norm->SaveAs(Form("../plots_bdt/bdt_norm_%s.pdf", ch_nm.Data()));
        delete cv_norm;
        delete leg1; delete leg2; delete leg3; delete leg4; delete leg5;
	*/
    }
}

// Helper function implementations

double get_fbeta(double a1, double b1, double c1, double m2pi) {
    m2pi = m2pi / 1000.;
    double fbeta = a1 + 1. / (exp((m2pi - c1) / b1) - 1.);
    return fbeta;
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
