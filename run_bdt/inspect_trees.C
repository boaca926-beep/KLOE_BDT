#include "../header_bdt/helper.h"

// Constants
constexpr double BDT_CUT_VALUE = 0.4;           // BDT score threshold
constexpr int N_BINS_ENERGY = 200;
constexpr double ENERGY_RANGE_MAX = 500.0;      // MeV

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

// ----------------------------------------------------------------------
// Configuration
// ----------------------------------------------------------------------
//#define INPUT_FILE_PATH "/home/kloe/Desktop/input_bdt_TDATA_chain/cut/tree_pre_bdt.root"
#define INPUT_FILE_PATH "/home/kloe/Desktop/input_tmp_TDATA_chain/cut/tree_pre_bdt.root"
#define OUTPUT_FILE_PATH "./output_with_bdt.root"

struct DataAttr {
  const char* tree_name;      // name of the TTree to read
  const char* label;          // legend label (e.g., "Selected")
  const char* title_addon;    // additional title (e.g., "#eta#gamma")
  Color_t     color;          // colour for fill/line
  TString     cut_label;      // optional cut description (e.g., "BDT > 0.4")
  TString     entry_label;
};

// Example for TETAGAM (BDT‑selected)
/*
DataAttr myData = {
    "TETAGAM",          // tree_name
    "Selected",         // label
    "#eta#gamma",       // title_addon
    kGreen,             // colour
    Form("BDT > %.1f", BDT_CUT_VALUE)
};
*/

DataAttr myData = {
    "TETAGAM_COMB",          // tree_name
    "Discarded",             // label
    "#eta#gamma",            // title_addon
    kRed,                    // colour
    Form("BDT < %.1f", BDT_CUT_VALUE)
};

/**
 * Inspect trees in tree_pre_bdt.root
 * Apply BDT model to KLOE detector data for e+e-→3π decay analysis
 * @param INPUT_FILE_PATH Path to the input ROOT data file
 * @param OUTPUT_FILE_PATH Path to the output ROOT data file
 * @param DATA_TYPE Type of data to inspect
 *
 * The function inspect kinematic variables used in the analysis for comparison.
 */

int inspect_trees(){
  TStopwatch timer;
  timer.Start();
            
  // ROOT Style settings
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFitFormat("6.4g");

  // Initialize histograms
  HistogramManager hists;
    
  TH1D* hE1_BDT = hists.create("hE1_BDT", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
  TH1D* hE2_BDT = hists.create("hE2_BDT", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
  TH1D* hE3_BDT = hists.create("hE3_BDT", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    
  // 1. Process data file if it exists
  if (!gSystem->AccessPathName(INPUT_FILE_PATH)) {
    std::cout << "\nProcessing data file: " << INPUT_FILE_PATH << std::endl;

    TFile *f_input = new TFile(INPUT_FILE_PATH);
    TTree *DATA_TREE = (TTree*)f_input->Get(myData.tree_name);
    if (!DATA_TREE) {
      cerr << "ERROR: Cannot find tree " << myData.tree_name << endl;
      return 1;
    }
    else {
      cout << myData.tree_name << " exists!" << endl;
    }

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
    double betapi0_bdt = 0.0;
    double angle_pi0gam12_bdt = 0.0;
	
    TLorentzVector pi0gam1_bdt(0.,0.,0.,0.);
    TLorentzVector pi0gam2_bdt(0.,0.,0.,0.);
    TLorentzVector pi0_bdt;
	
    int event_id = 0;
    
    // Set branch addresses
    DATA_TREE->SetBranchAddress("Br_deltaE", &deltaE);
    DATA_TREE->SetBranchAddress("Br_angle_pi0gam12", &angle_pi0gam12);
    DATA_TREE->SetBranchAddress("Br_betapi0", &betapi0);
    DATA_TREE->SetBranchAddress("Br_lagvalue_min_7C", &lagvalue_min_7C);
    DATA_TREE->SetBranchAddress("Br_ppIM", &ppIM);
    DATA_TREE->SetBranchAddress("Br_bkg_indx", &bkg_indx);
    DATA_TREE->SetBranchAddress("Br_recon_indx", &recon_indx);
    DATA_TREE->SetBranchAddress("Br_ppl_E", &ppl_E);
    DATA_TREE->SetBranchAddress("Br_ppl_px", &ppl_px);
    DATA_TREE->SetBranchAddress("Br_ppl_py", &ppl_py);
    DATA_TREE->SetBranchAddress("Br_ppl_pz", &ppl_pz);
    DATA_TREE->SetBranchAddress("Br_pmi_E", &pmi_E);
    DATA_TREE->SetBranchAddress("Br_pmi_px", &pmi_px);
    DATA_TREE->SetBranchAddress("Br_pmi_py", &pmi_py);
    DATA_TREE->SetBranchAddress("Br_pmi_pz", &pmi_pz);
    DATA_TREE->SetBranchAddress("Br_E1", &E1);
    DATA_TREE->SetBranchAddress("Br_px1", &px1);
    DATA_TREE->SetBranchAddress("Br_py1", &py1);
    DATA_TREE->SetBranchAddress("Br_pz1", &pz1);
    DATA_TREE->SetBranchAddress("Br_E2", &E2);
    DATA_TREE->SetBranchAddress("Br_px2", &px2);
    DATA_TREE->SetBranchAddress("Br_py2", &py2);
    DATA_TREE->SetBranchAddress("Br_pz2", &pz2);
    DATA_TREE->SetBranchAddress("Br_E3", &E3);
    DATA_TREE->SetBranchAddress("Br_px3", &px3);
    DATA_TREE->SetBranchAddress("Br_py3", &py3);
    DATA_TREE->SetBranchAddress("Br_pz3", &pz3);
    DATA_TREE->SetBranchAddress("Br_E1_true", &E1_true);
    DATA_TREE->SetBranchAddress("Br_px1_true", &px1_true);
    DATA_TREE->SetBranchAddress("Br_py1_true", &py1_true);
    DATA_TREE->SetBranchAddress("Br_pz1_true", &pz1_true);
    DATA_TREE->SetBranchAddress("Br_E2_true", &E2_true);
    DATA_TREE->SetBranchAddress("Br_px2_true", &px2_true);
    DATA_TREE->SetBranchAddress("Br_py2_true", &py2_true);
    DATA_TREE->SetBranchAddress("Br_pz2_true", &pz2_true);
    DATA_TREE->SetBranchAddress("Br_E3_true", &E3_true);
    DATA_TREE->SetBranchAddress("Br_px3_true", &px3_true);
    DATA_TREE->SetBranchAddress("Br_py3_true", &py3_true);
    DATA_TREE->SetBranchAddress("Br_pz3_true", &pz3_true);
    DATA_TREE->SetBranchAddress("Br_ppl_E_true", &ppl_E_true);
    DATA_TREE->SetBranchAddress("Br_ppl_px_true", &ppl_px_true);
    DATA_TREE->SetBranchAddress("Br_ppl_py_true", &ppl_py_true);
    DATA_TREE->SetBranchAddress("Br_ppl_pz_true", &ppl_pz_true);
    DATA_TREE->SetBranchAddress("Br_pmi_E_true", &pmi_E_true);
    DATA_TREE->SetBranchAddress("Br_pmi_px_true", &pmi_px_true);
    DATA_TREE->SetBranchAddress("Br_pmi_py_true", &pmi_py_true);
    DATA_TREE->SetBranchAddress("Br_pmi_pz_true", &pmi_pz_true);
    DATA_TREE->SetBranchAddress("Br_e1_bdt", &e1_bdt);
    DATA_TREE->SetBranchAddress("Br_e2_bdt", &e2_bdt);
    DATA_TREE->SetBranchAddress("Br_e3_bdt", &e3_bdt);

    TFile* outfile = TFile::Open(OUTPUT_FILE_PATH, "RECREATE");
    /*
    TTree* outtree = new TTree("new_tree", "Tree with BDT response");
        
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
    outtree->Branch("betapi0_bdt", &betapi0_bdt);
    outtree->Branch("betapi0", &betapi0);
    outtree->Branch("angle_pi0gam12_bdt", &angle_pi0gam12_bdt);
    outtree->Branch("angle_pi0gam12", &angle_pi0gam12);
    */
    
    // Main event loop
    int nentries = DATA_TREE->GetEntries();
        
    for (int i = 0; i < nentries; i++) {
      DATA_TREE->GetEntry(i);
      
      // Test variables
      //cout << E1_true << ", " << px1_true << ", " << py1_true << ", " << pz1_true << endl;
      //cout << lagvalue_min_7C << endl;
      //cout << ppl_px_true << ", " << ppl_py_true << ", " << ppl_pz_true << endl;
      //cout << pmi_px_true << ", " << pmi_py_true << ", " << pmi_pz_true << endl;
      //cout << e1_bdt << endl;
      
      hE1_BDT->Fill(e1_bdt);
      hE2_BDT->Fill(e2_bdt);
      hE3_BDT->Fill(e3_bdt);
                
    }

    
    // Write histograms and close files
    outfile->cd();
    hE1_BDT->Write();
    outfile->Close();
    f_input->Close();

    // Common style settings
    auto setHistStyle = [](TH1D* h, Color_t color, const char* xTitle, const char* yTitle) {
      h->SetLineColor(color);
      h->SetLineWidth(1);
      h->SetFillColor(color);
      h->SetFillStyle(3001);
      //h->SetFillColorAlpha(color, 0.3);
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
		       double yMaxFactor = 1.2, Color_t color = kGreen) {
      TCanvas* c = new TCanvas(canvasName, canvasTitle, 1800, 600);
      c->Divide(3,1);
      
      double max1 = h1->GetMaximum();
      double max2 = h2->GetMaximum();
      double max3 = h3->GetMaximum();
      //double yMax = TMath::Max(max1, TMath::Max(max2, max3)) * yMaxFactor;
      
      // Set margins for each pad individually
      for (int iPad = 1; iPad <= 3; ++iPad) {
	TPad* pad = (TPad*)c->GetPad(iPad);
	pad->SetBottomMargin(0.13);   // more space for x‑title
	pad->SetLeftMargin(0.16);     // more space for y‑title
      }
      
      TPaveText* pt1 = new TPaveText(0.2, 0.80, 0.80, 0.89, "NDC");
      PteAttr(pt1);
      pt1->SetTextSize(0.05);
      pt1->SetTextColor(kBlack);
      pt1->AddText(myData.title_addon);

      c->cd(1);
      setHistStyle(h1, color+2, xTitle1, yTitle1);
      h1->GetYaxis()->SetRangeUser(0, h1->GetBinContent(h1->GetMaximumBin()) * yMaxFactor);
      h1->Draw("HIST F");
      pt1->Draw("Same");
      TLegend* leg1 = new TLegend(0.65, 0.8, 0.9, 0.9);
      leg1->SetBorderSize(0);
      leg1->SetFillStyle(0);
      leg1->AddEntry(h1, lgdTitle, "f");
      leg1->Draw();

      TPaveText* pt2 = new TPaveText(0.2, 0.80, 0.80, 0.89, "NDC");
      PteAttr(pt2);
      pt2->SetTextSize(0.05);
      pt2->SetTextColor(kBlack);
      pt2->AddText(myData.cut_label);
        
      c->cd(2);
      setHistStyle(h2, color+2, xTitle2, yTitle2);
      h2->GetYaxis()->SetRangeUser(0, h2->GetBinContent(h2->GetMaximumBin()) * yMaxFactor);
      h2->Draw("HIST F");
      pt2->Draw("Same");
      //TLegend* leg2 = new TLegend(0.65, 0.75, 0.9, 0.9);
      //leg2->SetBorderSize(0);
      //leg2->SetFillStyle(0);
      //leg2->AddEntry(h2, "BDT Selected", "f");
      //leg2->Draw();

      double entries = hE1_BDT->GetEntries();
      myData.entry_label = Form("Events=%.0f", entries);
      cout << myData.entry_label << endl;

      TPaveText* pt3 = new TPaveText(0.2, 0.80, 0.80, 0.89, "NDC");
      PteAttr(pt3);
      pt3->SetTextSize(0.05);
      pt3->SetTextColor(kBlack);
      pt3->AddText(myData.entry_label);
      
      c->cd(3);
      setHistStyle(h3, color+2, xTitle3, yTitle3);
      h3->GetYaxis()->SetRangeUser(0, h3->GetBinContent(h3->GetMaximumBin()) * yMaxFactor);
      h3->Draw("HIST F");
      pt3->Draw("Same");
      //TLegend* leg3 = new TLegend(0.65, 0.75, 0.9, 0.9);
      //leg3->SetBorderSize(0);
      //leg3->SetFillStyle(0);
      //leg3->AddEntry(h3, "BDT Selected", "f");
      //leg3->Draw();
      
      c->SaveAs(Form("../plots_bdt/%s.pdf", canvasName));
      delete c;
    };


    // Draw Energy pulls
    TString phoE_BDT_cv_nm = TString("phoE_BDT_") + myData.label + "_" + myData.tree_name;
    TString phoE_BDT_title = TString("Photon Energy ") + myData.tree_name;
    drawSet(phoE_BDT_cv_nm.Data(), phoE_BDT_title.Data(), myData.label,
	    hE1_BDT, hE2_BDT, hE3_BDT,
	    "E_{1} [MeV]", "E_{2} [MeV]", "E_{3} [MeV]",
	    "Entries", "Entries", "Entries",
	    1.2,
	    myData.color);
	
    // ==================================================================
    // Normalized pull histograms
    // ==================================================================
    //auto safeNormalize = [](TH1D* h) {
    //  if (h && h->Integral() > 0) h->Scale(1.0 / h->Integral());
    //};

  }

  /*
  // ---------- Event loop ----------
  Long64_t nentries = ALLCHAIN_CUT->GetEntries();
  cout << "Processing " << nentries << " events" << endl;

  for (Long64_t irow = 0; irow < nentries; irow++) {
    ALLCHAIN_CUT->GetEntry(irow);
    if (irow % 100000 == 0) cout << "Event " << irow << endl;
  }
  */

  timer.Stop();
  timer.Print();
  return 0;
}
