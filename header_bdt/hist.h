TFile *f_cut = new TFile(outputCut + "tree_pre_bdt.root");
TFile *f_gen = new TFile(outputGen + "tree_gen.root");

//cout << f_gen -> GetName() << endl;

const double eeg_lsf = 2.;
const double mass_sigma_nb = 1;
const double sfw2d_sigma_nb = 1;

//IM3pi (analysis):
const double IM3pi_min = 380;
const double IM3pi_max = 1020;
const double IM3pi_sigma = 2.65;
const int IM3pi_bin = TMath::Nint((IM3pi_max - IM3pi_min) / mass_sigma_nb / IM3pi_sigma);

//Angle_gamm12
const double angle_min = 20;
const double angle_max = 140;
const double angle_sigma = 0.53;
const int angle_bin = TMath::Nint((angle_max - angle_min) / sfw2d_sigma_nb / angle_sigma);

//mgg
//mgg_bin, mgg_min, mgg_max
const double mgg_min = 100;
const double mgg_max = 180;
const double mgg_sigma = 2.09;
const int mgg_bin = TMath::Nint((mgg_max - mgg_min) / sfw2d_sigma_nb / mgg_sigma);

//Eisr:
const double Eisr_min = 50;
const double Eisr_max = 500;
const double Eisr_sigma = 2.48;
const int Eisr_bin = TMath::Nint((Eisr_max - Eisr_min) / sfw2d_sigma_nb / Eisr_sigma);

//bdt_score
//bdt_score_bin, bdt_score_min, bdt_score_max
const double bdt_score_min = 0;
const double bdt_score_max = 1;
const int bdt_score_bin = 150;
  
//ppIM:
const double ppIM_min = 200; 
const double ppIM_max = 700; 
const double ppIM_sigma = 2.30;
const int ppIM_bin = TMath::Nint((ppIM_max - ppIM_min) / sfw2d_sigma_nb / ppIM_sigma);

// betapi0:
const double betapi0_min = 0.;
const double betapi0_max = 1.0;
const int betapi0_bin = 50;

//sfw1d
const double xmin = 770;
const double xmax = 800;
const int xbins = 60; 

//crx3pi
const double hmin = 740; //770, 700
const double hmax = 820; //800, 820
const int hbins = TMath::Nint((hmax - hmin) / 0.25 / IM3pi_sigma);

TList *HIM3pi_fit = new TList(); // IM3pi distr. for fit omega parameters
TList *HPeakNonReson = new TList(); // IM3pi distr. for peak and non-resonant
TList *HSFW2D = new TList(); // Eisr vs. ppIM distr. for MC normalization
TList *HSFW1D = new TList(); // IM3pi distr. for signal MC tuning
TList *HSIG = new TList(); // IM3pi distr. signal true and generated
TList *HIM3pi_crx = new TList(); // IM3pi distr. for crx3pi obs.
TList *HppIM_vs_betapi0 = new TList();   // 2D for ppIM vs betapi0
TList *Heisr_vs_angle = new TList();   // 2D for Eisr vs angle

//TRandom *rnd=0;

// methods
void fillHist() {

  // data and MC background
  TIter next_tree(f_cut -> GetListOfKeys());
  TKey *key;

  while ( (key = (TKey *) next_tree() ) ) {// start tree while loop

    TString objnm_tree = key -> GetName();
    TString classnm_tree = key -> GetClassName();
    key -> GetSeekKey();

    cout << "classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;

    TTree *tree_tmp = (TTree*)f_cut -> Get(objnm_tree);
    //cout << tree_tmp -> GetName() << endl;
    if (!tree_tmp) continue;

    // ---- Set branch addresses for this tree ----
    double m3pi_bdt = 0., m3pi_true_bdt = 0., Eisr = 0., ppIM = 0., betapi0_bdt = 0.;
    double angle_bdt = 0., m_gg = 0., bdt_score = 0.;
    int recon_indx_bdt = 0, bkg_indx = 0;
    
    tree_tmp->SetBranchAddress("Br_m3pi_bdt", &m3pi_bdt);
    tree_tmp->SetBranchAddress("Br_m3pi_true_bdt", &m3pi_true_bdt);
    tree_tmp->SetBranchAddress("Br_e3_bdt", &Eisr);
    tree_tmp->SetBranchAddress("Br_ppIM", &ppIM);
    tree_tmp->SetBranchAddress("Br_m_gg_bdt", &m_gg);
    tree_tmp->SetBranchAddress("Br_bdt_score", &bdt_score);
    tree_tmp->SetBranchAddress("Br_betapi0_bdt", &betapi0_bdt);
    tree_tmp->SetBranchAddress("Br_angle_pi0gam12_bdt", &angle_bdt);
    tree_tmp->SetBranchAddress("Br_recon_indx_bdt", &recon_indx_bdt);
    tree_tmp->SetBranchAddress("Br_bkg_indx", &bkg_indx);
    	
    // fill histos
    TH1D * h1d_tmp_crx = new TH1D("h1d_IM3pi_" + objnm_tree + "_CRX", "", hbins, hmin, hmax);
    h1d_tmp_crx -> Sumw2();

    TH1D * h1d_tmp_sfw = new TH1D("h1d_IM3pi_sfw_" + objnm_tree, "", xbins, xmin, xmax);
    h1d_tmp_sfw -> Sumw2();

    TH1D * h1d_tmp = new TH1D("h1d_IM3pi_" + objnm_tree, "", IM3pi_bin, IM3pi_min, IM3pi_max);
    h1d_tmp -> Sumw2();

    TH1D * h1d_tmp_true = new TH1D("h1d_IM3pi_" + objnm_tree + "_TRUE", "", IM3pi_bin, IM3pi_min, IM3pi_max);
    h1d_tmp_true -> Sumw2();

    //
    TH1D * h1d_tmp_peak = new TH1D("h1d_IM3pi_" + objnm_tree + "_PEAK", "", IM3pi_bin, IM3pi_min, IM3pi_max);
    h1d_tmp_peak -> Sumw2();

    TH1D * h1d_tmp_peak_true = new TH1D("h1d_IM3pi_" + objnm_tree + "_PEAK_TRUE", "", IM3pi_bin, IM3pi_min, IM3pi_max);
    h1d_tmp_peak_true -> Sumw2();

    //
    TH1D * h1d_tmp_non_reson = new TH1D("h1d_IM3pi_" + objnm_tree + "_NON_RESON", "", IM3pi_bin, IM3pi_min, IM3pi_max);
    h1d_tmp_non_reson -> Sumw2();

    TH1D * h1d_tmp_non_reson_true = new TH1D("h1d_IM3pi_" + objnm_tree + "_NON_RESON_TRUE", "", IM3pi_bin, IM3pi_min, IM3pi_max);
    h1d_tmp_non_reson_true -> Sumw2();

    //
    TH2D * h2d_tmp = new TH2D("h2d_sfw_" + objnm_tree, "", ppIM_bin, ppIM_min, ppIM_max, Eisr_bin, Eisr_min, Eisr_max);
    h2d_tmp -> Sumw2();

    TH2D * h2d_tmp_peak = new TH2D("h2d_sfw_" + objnm_tree + "_peak", "", ppIM_bin, ppIM_min, ppIM_max, Eisr_bin, Eisr_min, Eisr_max);
    h2d_tmp_peak -> Sumw2();

    TH2D * h2d_tmp_non_reson = new TH2D("h2d_sfw_" + objnm_tree + "_non_reson", "", ppIM_bin, ppIM_min, ppIM_max, Eisr_bin, Eisr_min, Eisr_max);
    h2d_tmp_non_reson -> Sumw2();
    
    /*
    TH2D * h2d_tmp = new TH2D("h2d_sfw_" + objnm_tree, "", angle_bin, angle_min, angle_max, mgg_bin, mgg_min, mgg_max);
    h2d_tmp -> Sumw2();

    TH2D * h2d_tmp_peak = new TH2D("h2d_sfw_" + objnm_tree + "_peak", "", angle_bin, angle_min, angle_max, mgg_bin, mgg_min, mgg_max);
    h2d_tmp_peak -> Sumw2();

    TH2D * h2d_tmp_non_reson = new TH2D("h2d_sfw_" + objnm_tree + "_non_reson", "", angle_bin, angle_min, angle_max, mgg_bin, mgg_min, mgg_max);
    h2d_tmp_non_reson -> Sumw2();
    */

    /*
    TH2D * h2d_tmp = new TH2D("h2d_sfw_" + objnm_tree, "", bdt_score_bin, bdt_score_min, bdt_score_max, Eisr_bin, Eisr_min, Eisr_max);
    h2d_tmp -> Sumw2();

    TH2D * h2d_tmp_peak = new TH2D("h2d_sfw_" + objnm_tree + "_peak", "", bdt_score_bin, bdt_score_min, bdt_score_max, Eisr_bin, Eisr_min, Eisr_max);
    h2d_tmp_peak -> Sumw2();

    TH2D * h2d_tmp_non_reson = new TH2D("h2d_sfw_" + objnm_tree + "_non_reson", "", bdt_score_bin, bdt_score_min, bdt_score_max, Eisr_bin, Eisr_min, Eisr_max);
    h2d_tmp_non_reson -> Sumw2();
    */
    //
    TH2D *h2d_ppIM_vs_beta = new TH2D("h2d_ppIM_vs_betapi0_" + objnm_tree, "", 200, 0.25, 0.65, 200, 0.3, 1.);  // range in GeV  
    h2d_ppIM_vs_beta -> Sumw2();

    TH2D *h2d_ppIM_vs_beta_peak = new TH2D("h2d_ppIM_vs_betapi0_" + objnm_tree + "_peak", "", 200, 0.25, 0.65, 200, 0.3, 1.);  // range in GeV  
    h2d_ppIM_vs_beta_peak -> Sumw2();

    TH2D *h2d_ppIM_vs_beta_non_reson = new TH2D("h2d_ppIM_vs_betapi0_" + objnm_tree + "_non_reson", "", 200, 0.25, 0.65, 200, 0.3, 1.);  // range in GeV  
    h2d_ppIM_vs_beta_non_reson -> Sumw2();

    //
    TH2D *h2d_angle_vs_eisr = new TH2D("h2d_angle_vs_eisr_" + objnm_tree, "", angle_bin, angle_min, angle_max, Eisr_bin, Eisr_min, Eisr_max);    
    h2d_angle_vs_eisr -> Sumw2();

    for (Int_t irow = 0; irow < tree_tmp -> GetEntries(); irow++) {// loop chain

      tree_tmp -> GetEntry(irow);
      
      // filling
      h1d_tmp -> Fill(m3pi_bdt);
      h1d_tmp_sfw -> Fill(m3pi_bdt);
      h1d_tmp_crx -> Fill(m3pi_bdt);
      h1d_tmp_true -> Fill(m3pi_true_bdt);
      h2d_tmp -> Fill(ppIM, Eisr);
      //h2d_tmp -> Fill(angle_bdt, m_gg);
      //h2d_tmp -> Fill(bdt_score, Eisr);
      //cout << "ppIM = " << ppIM << ", Eisr = " << Eisr << endl;
      //cout << "angle_bdt = " << angle_bdt << ", m_gg = " << m_gg << endl;
      //cout << "bdt_score = " << bdt_score << ", m_gg = " << m_gg << endl;
      
      h2d_ppIM_vs_beta -> Fill(ppIM * 1e-3, betapi0_bdt);
      h2d_angle_vs_eisr->Fill(angle_bdt, Eisr);

      if (recon_indx_bdt == 2 && bkg_indx == 1) {
        h1d_tmp_peak->Fill(m3pi_bdt);
        h1d_tmp_peak_true->Fill(m3pi_true_bdt);
        h2d_tmp_peak->Fill(ppIM, Eisr);
	//h2d_tmp_peak->Fill(angle_bdt, m_gg);
	//h2d_tmp_peak -> Fill(bdt_score, Eisr);
        h2d_ppIM_vs_beta_peak -> Fill(ppIM * 1e-3, betapi0_bdt);
      }
      else {
        h1d_tmp_non_reson->Fill(m3pi_bdt);
        h1d_tmp_non_reson_true->Fill(m3pi_true_bdt);
        h2d_tmp_non_reson->Fill(ppIM, Eisr);
	//h2d_tmp_non_reson->Fill(angle_bdt, m_gg);
	//h2d_tmp_non_reson->Fill(bdt_score, Eisr);
        h2d_ppIM_vs_beta_non_reson -> Fill(ppIM * 1e-3, betapi0_bdt);
      }
      
    }

    // scale EEG by eeg_lsf
    if (objnm_tree == "TEEG") {
       h1d_tmp -> Scale(eeg_lsf);
       h1d_tmp_peak->Scale(eeg_lsf);
       h1d_tmp_non_reson->Scale(eeg_lsf);
       
       h2d_tmp -> Scale(eeg_lsf);
       h2d_tmp_peak -> Scale(eeg_lsf);
       h2d_tmp_non_reson -> Scale(eeg_lsf);

       h2d_ppIM_vs_beta -> Scale(eeg_lsf);
       h2d_ppIM_vs_beta_peak -> Scale(eeg_lsf);
       h2d_ppIM_vs_beta_non_reson -> Scale(eeg_lsf);
       
       h2d_angle_vs_eisr->Scale(eeg_lsf);
    }

    HIM3pi_fit -> Add(h1d_tmp);
    HIM3pi_fit -> Add(h1d_tmp_true);
    
    HPeakNonReson -> Add(h1d_tmp_peak);
    HPeakNonReson -> Add(h1d_tmp_peak_true);
    HPeakNonReson -> Add(h1d_tmp_non_reson);
    HPeakNonReson -> Add(h1d_tmp_non_reson_true);
    
    HSFW1D -> Add(h1d_tmp_sfw);

    HSFW2D -> Add(h2d_tmp);
    HSFW2D -> Add(h2d_tmp_peak);
    HSFW2D -> Add(h2d_tmp_non_reson);
    
    HppIM_vs_betapi0 -> Add(h2d_ppIM_vs_beta);
    HppIM_vs_betapi0 -> Add(h2d_ppIM_vs_beta_peak);
    HppIM_vs_betapi0 -> Add(h2d_ppIM_vs_beta_non_reson);
    
    Heisr_vs_angle->Add(h2d_angle_vs_eisr);
    HIM3pi_crx -> Add(h1d_tmp_crx);
    
  }  // Close while loop

  // ===== SIGNAL TRUE (from f_cut) =====
  TH1D * hsig_true = new TH1D("hsig_true", "", IM3pi_bin, IM3pi_min, IM3pi_max);
  hsig_true -> Sumw2();

  TH1D * hsig_true_crx = new TH1D("h1d_IM3pi_TISR3PI_SIG_TRUE_CRX", "", hbins, hmin, hmax);
  hsig_true_crx -> Sumw2();

  TTree *TISR3PI_SIG = (TTree*)f_cut -> Get("TISR3PI_SIG");

  if (TISR3PI_SIG) {
    double m3pi_true = 0.;
    TISR3PI_SIG->SetBranchAddress("Br_m3pi_true_bdt", &m3pi_true);
    
    for (Int_t irow = 0; irow < TISR3PI_SIG -> GetEntries(); irow++) {
      TISR3PI_SIG -> GetEntry(irow);
      hsig_true -> Fill(m3pi_true);
      hsig_true_crx -> Fill(m3pi_true);
    }
  }

  HSIG -> Add(hsig_true);
  HIM3pi_crx -> Add(hsig_true_crx);

  // ===== SIGNAL GENERATED (from f_gen) =====
  TTree *TISR3PI_SIG_GEN = (TTree*)f_gen -> Get("TISR3PI_SIG_GEN");

  if (TISR3PI_SIG_GEN) {
    double m3pi_gen = 0.;

    TH1D * hsig_gen = new TH1D("hsig_gen", "", IM3pi_bin, IM3pi_min, IM3pi_max);
    hsig_gen -> Sumw2();

    TH1D * hsig_gen_crx = new TH1D("h1d_IM3pi_TISR3PI_SIG_GEN_CRX", "", hbins, hmin, hmax);
    hsig_gen_crx -> Sumw2();

    TISR3PI_SIG_GEN->SetBranchAddress("Br_IM3pi_gen", &m3pi_gen);

    for (Int_t irow = 0; irow < TISR3PI_SIG_GEN -> GetEntries(); irow++) {
      TISR3PI_SIG_GEN -> GetEntry(irow);
      hsig_gen -> Fill(m3pi_gen);
      hsig_gen_crx -> Fill(m3pi_gen);
    }

    HSIG -> Add(hsig_gen);
    HIM3pi_crx -> Add(hsig_gen_crx);
  }

}  // Close fillHist() function
