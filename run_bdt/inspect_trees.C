// Constants
constexpr double BDT_CUT_VALUE = 0.4;           // BDT score threshold (not used)
constexpr int N_BINS_ENERGY = 200;
constexpr int N_BINS_MASS = 200;
constexpr int N_BINS_PULL = 150;
constexpr int N_BINS_CHI2 = 100;
constexpr int N_BINS_ANGLE = 180;
constexpr double ENERGY_RANGE_MAX = 500.0;      // MeV
constexpr double MASS_GG_RANGE_MAX = 200.0;     // MeV/c²
constexpr double MASS_GG_RANGE_MIN = 50.0;      // MeV/c²
constexpr double MASS_3PI_RANGE_MAX = 1000.0;   // MeV/c²
constexpr double MASS_3PI_RANGE_MIN = 400.0;    // MeV/c²
constexpr double MASS_2PI_RANGE_MAX = 800.0;    // MeV/c²
constexpr double MASS_2PI_RANGE_MIN = 200.0;    // MeV/c²
constexpr double ANGLE_RANGE_MAX = 180.0;       // deg
constexpr double CHI2_RANGE_MAX = 50.0;      
constexpr double PULL_RANGE_MIN = -30;          // MeV/c² / MeV
constexpr double PULL_RANGE_MAX = 30;           // MeV/c² / MeV

// Histogram manager (with SetDirectory(0) to avoid ROOT ownership)
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

// ----------------------------------------------------------------------
// Configuration – single tree, no BDT score cut
// ----------------------------------------------------------------------
#define INPUT_FILE_PATH "/home/kloe/Desktop/input_bdt_TDATA_chain/cut/tree_pre_bdt.root"
#define OUTPUT_FILE_PATH "./output_with_bdt.root"

struct DataAttr {
  const char* tree_name;      // single tree containing BDT-selected branches
  const char* label;          // label for plots
  Color_t     color;          // colour for histograms
  TString     cut_label;      // optional description (e.g., "BDT selection")
};

// Example for signal 3π channel
DataAttr myData = {
  "TISR3PI_SIG",              // tree_name
  "BDT-selected pair",        // label
  kGreen,                     // colour
  "BDT selection (no score cut)"
};

// ----------------------------------------------------------------------
int inspect_trees() {
    TStopwatch timer;
    timer.Start();

    // ROOT Style settings
    gErrorIgnoreLevel = kError;
    TGaxis::SetMaxDigits(4);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetFitFormat("6.4g");

    gSystem->mkdir("../plots_trees/", kTRUE);

    if (gSystem->AccessPathName(INPUT_FILE_PATH)) {
        std::cerr << "ERROR: Input file not found: " << INPUT_FILE_PATH << std::endl;
        return 1;
    }

    TFile *f_input = new TFile(INPUT_FILE_PATH);
    TTree *tree = (TTree*)f_input->Get(myData.tree_name);
    if (!tree) {
        std::cerr << "ERROR: Cannot find tree " << myData.tree_name << std::endl;
        return 1;
    }
    std::cout << "Tree " << myData.tree_name << " has " << tree->GetEntries() << " entries." << std::endl;

    HistogramManager hists;

    // ---------- Histograms (one set for all events) ----------
    // Energies
    TH1D* hE1 = hists.create("hE1", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hE2 = hists.create("hE2", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);
    TH1D* hE3 = hists.create("hE3", "", N_BINS_ENERGY, 0, ENERGY_RANGE_MAX);

    // Kinematic variables
    TH1D* hMgg = hists.create("hMgg", "", N_BINS_MASS, MASS_GG_RANGE_MIN, MASS_GG_RANGE_MAX);
    TH1D* hM3pi = hists.create("hM3pi", "", N_BINS_MASS, MASS_3PI_RANGE_MIN, MASS_3PI_RANGE_MAX);
    TH1D* hM2pi = hists.create("hM2pi", "", N_BINS_MASS, MASS_2PI_RANGE_MIN, MASS_2PI_RANGE_MAX);
    TH1D* hAngle = hists.create("hAngle", "", N_BINS_ANGLE, 0, ANGLE_RANGE_MAX);
    TH1D* hChi2 = hists.create("hChi2", "", N_BINS_CHI2, 0, CHI2_RANGE_MAX);
    
    // Pulls
    TH1D* hE1_pull = hists.create("hE1_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hE2_pull = hists.create("hE2_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hE3_pull = hists.create("hE3_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPx1_pull = hists.create("hPx1_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPx2_pull = hists.create("hPx2_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPx3_pull = hists.create("hPx3_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPy1_pull = hists.create("hPy1_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPy2_pull = hists.create("hPy2_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPy3_pull = hists.create("hPy3_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPz1_pull = hists.create("hPz1_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPz2_pull = hists.create("hPz2_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hPz3_pull = hists.create("hPz3_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hMgg_pull = hists.create("hMgg_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hM3pi_pull = hists.create("hM3pi_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);
    TH1D* hM2pi_pull = hists.create("hM2pi_pull", "", N_BINS_PULL, PULL_RANGE_MIN, PULL_RANGE_MAX);

    // ---------- Branch addresses ----------
    double e1, e2, e3, mgg, m3pi, m2pi, angle, chi2;
    double e1p, e2p, e3p;
    double px1p, px2p, px3p, py1p, py2p, py3p, pz1p, pz2p, pz3p;
    double mggp, m3pip, m2pip;

    tree->SetBranchAddress("Br_e1_bdt", &e1);
    tree->SetBranchAddress("Br_e2_bdt", &e2);
    tree->SetBranchAddress("Br_e3_bdt", &e3);
    tree->SetBranchAddress("Br_m_gg_bdt", &mgg);
    tree->SetBranchAddress("Br_m3pi_bdt", &m3pi);
    tree->SetBranchAddress("Br_ppIM", &m2pi);          // dipion mass
    tree->SetBranchAddress("Br_angle_pi0gam12_bdt", &angle);
    tree->SetBranchAddress("Br_lagvalue_min_7C", &chi2);
    tree->SetBranchAddress("Br_e1_pull", &e1p);
    tree->SetBranchAddress("Br_e2_pull", &e2p);
    tree->SetBranchAddress("Br_e3_pull", &e3p);
    tree->SetBranchAddress("Br_px1_pull", &px1p);
    tree->SetBranchAddress("Br_px2_pull", &px2p);
    tree->SetBranchAddress("Br_px3_pull", &px3p);
    tree->SetBranchAddress("Br_py1_pull", &py1p);
    tree->SetBranchAddress("Br_py2_pull", &py2p);
    tree->SetBranchAddress("Br_py3_pull", &py3p);
    tree->SetBranchAddress("Br_pz1_pull", &pz1p);
    tree->SetBranchAddress("Br_pz2_pull", &pz2p);
    tree->SetBranchAddress("Br_pz3_pull", &pz3p);
    tree->SetBranchAddress("Br_m_gg_pull", &mggp);
    tree->SetBranchAddress("Br_m3pi_pull", &m3pip);
    tree->SetBranchAddress("Br_m2pi_pull", &m2pip);

    // ---------- Event loop (no BDT score cut) ----------
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

	//cout << chi2 << endl;
	
        // Fill all histograms (no event rejection)
        hE1->Fill(e1); hE2->Fill(e2); hE3->Fill(e3);
        hMgg->Fill(mgg); hM3pi->Fill(m3pi); hM2pi->Fill(m2pi);
        hAngle->Fill(angle); hChi2->Fill(chi2);
        hE1_pull->Fill(e1p); hE2_pull->Fill(e2p); hE3_pull->Fill(e3p);
        hPx1_pull->Fill(px1p); hPx2_pull->Fill(px2p); hPx3_pull->Fill(px3p);
        hPy1_pull->Fill(py1p); hPy2_pull->Fill(py2p); hPy3_pull->Fill(py3p);
        hPz1_pull->Fill(pz1p); hPz2_pull->Fill(pz2p); hPz3_pull->Fill(pz3p);
        hMgg_pull->Fill(mggp); hM3pi_pull->Fill(m3pip); hM2pi_pull->Fill(m2pip);
    }

    // ---------- Normalisation of pull histograms (shape comparison) ----------
    auto safeNormalize = [](TH1D* h) { if (h && h->Integral()>0) h->Scale(1.0/h->Integral()); };
    std::vector<TH1D*> pullHists = {
        hE1_pull, hE2_pull, hE3_pull,
        hPx1_pull, hPx2_pull, hPx3_pull,
        hPy1_pull, hPy2_pull, hPy3_pull,
        hPz1_pull, hPz2_pull, hPz3_pull,
        hMgg_pull, hM3pi_pull
    };
    for (auto h : pullHists) safeNormalize(h);
    // hM2pi_pull left as raw counts (optional)

    // ---------- Drawing functions (single histogram per pad) ----------
    auto setHistStyle = [](TH1D* h, Color_t color, const char* xTitle, const char* yTitle) {
        h->SetLineColor(color);
        h->SetLineWidth(1);
	h->SetFillColor(color);
        //h->SetFillColorAlpha(color, 0.2);
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

    auto drawTriple = [&](const char* name, const char* title,
                          TH1D* h1, TH1D* h2, TH1D* h3,
                          const char* xTitle1, const char* xTitle2, const char* xTitle3,
                          const char* yTitle, bool logy = false) {
        TCanvas* c = new TCanvas(name, title, 1800, 600);
        c->Divide(3,1);
        for (int i=1; i<=3; ++i) {
            TPad* pad = (TPad*)c->GetPad(i);
            pad->SetBottomMargin(0.14);
            pad->SetLeftMargin(0.16);
        }
        auto drawPad = [&](int iPad, TH1D* h, const char* xtitle) {
            c->cd(iPad);
            setHistStyle(h, myData.color, xtitle, yTitle);
            double max = h->GetMaximum();
            if (max > 0) h->GetYaxis()->SetRangeUser(0, max * 1.2);
            h->Draw("HIST");
            //if (logy) c->SetLogy();
            TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
            leg->AddEntry(h, myData.label, "f");
            leg->Draw();
        };
        drawPad(1, h1, xTitle1);
        drawPad(2, h2, xTitle2);
        drawPad(3, h3, xTitle3);
        c->SaveAs(Form("../plots_trees/%s.pdf", name));
        delete c;
    };

    auto drawSingle = [&](const char* name, const char* title,
			  TH1D* h, const char* xTitle, const char* yTitle,
			  bool logy = false) {
      TCanvas* c = new TCanvas(name, title, 1000, 800); // larger canvas
      c->SetBottomMargin(0.15);
      c->SetLeftMargin(0.12);
      setHistStyle(h, myData.color, xTitle, yTitle);
      h->GetXaxis()->SetNdivisions(5, kTRUE);
      double max = h->GetMaximum();
      if (max > 0) {
        if (logy) h->GetYaxis()->SetRangeUser(0.5, max * 1.5);
        else h->GetYaxis()->SetRangeUser(0, max * 1.2);
      }
      h->Draw("HIST");
      //if (logy) c->SetLogy();
      TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
      leg->AddEntry(h, myData.label, "f");
      leg->Draw();
      c->SaveAs(Form("../plots_trees/%s.pdf", name));
      delete c;
    };
    
    // ---------- Produce plots ----------
    // 1. Photon energies
    drawTriple("photon_energies", "Photon Energies (BDT-selected)",
               hE1, hE2, hE3,
               "E_{1} [MeV]", "E_{2} [MeV]", "E_{3} [MeV]", "Entries");

    // 2. Energy pulls (normalised)
    drawTriple("E_pulls", "Energy Pulls (BDT-selected)",
               hE1_pull, hE2_pull, hE3_pull,
               "E_{1} pull [MeV]", "E_{2} pull [MeV]", "E_{3} pull [MeV]", "Normalized entries");

    // 3. Momentum pulls
    drawTriple("px_pulls", "Px Pulls (BDT-selected)",
               hPx1_pull, hPx2_pull, hPx3_pull,
               "p_{x,1} pull [MeV/c]", "p_{x,2} pull [MeV/c]", "p_{x,3} pull [MeV/c]", "Normalized entries");
    drawTriple("py_pulls", "Py Pulls (BDT-selected)",
               hPy1_pull, hPy2_pull, hPy3_pull,
               "p_{y,1} pull [MeV/c]", "p_{y,2} pull [MeV/c]", "p_{y,3} pull [MeV/c]", "Normalized entries");
    drawTriple("pz_pulls", "Pz Pulls (BDT-selected)",
               hPz1_pull, hPz2_pull, hPz3_pull,
               "p_{z,1} pull [MeV/c]", "p_{z,2} pull [MeV/c]", "p_{z,3} pull [MeV/c]", "Normalized entries");

    // 4. Mass pulls (Mgg, M3pi normalised; M2pi raw counts)
    drawTriple("mass_pulls", "Mass Pulls (BDT-selected)",
               hMgg_pull, hM3pi_pull, hM2pi_pull,
               "M_{#gamma#gamma} pull [MeV/c^{2}]", "M_{3#pi} pull [MeV/c^{2}]", "M_{2#pi} pull [MeV/c^{2}]",
               "Normalized entries");

    // 5. Kinematic variables (raw counts)
    drawTriple("kine_vars", "Kinematic Variables (BDT-selected)",
               hMgg, hM3pi, hAngle,
               "M_{#gamma#gamma} [MeV/c^{2}]", "M_{3#pi} [MeV/c^{2}]", "#angle_{#gamma#gamma} [#circ]", "Entries");

    // 6. Single plots: M2pi (dipion mass) and chi2
    drawSingle("m2pi_distribution", "Dipion Mass (BDT-selected)",
               hM2pi, "M_{2#pi} [MeV/c^{2}]", "Entries");
    drawSingle("chi2_distribution", "#chi^{2} Distribution (BDT-selected)",
               hChi2, "#chi^{2}", "Entries", true); // log y

    cout << hChi2->GetEntries() << endl;

    std::cout << "Chi2 max = " << hChi2->GetBinContent(hChi2->GetMaximumBin()) << std::endl;
    std::cout << "Chi2 entries = " << hChi2->GetEntries() << std::endl;
    std::cout << "Chi2 overflow = " << hChi2->GetBinContent(hChi2->GetNbinsX() + 1) << std::endl;
    
    f_input->Close();
    timer.Stop();
    timer.Print();
    return 0;
}
