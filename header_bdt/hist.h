TFile *f_cut = new TFile(outputCut + "tree_pre.root");
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

//Eisr:
const double Eisr_min = 50;
const double Eisr_max = 500;
const double Eisr_sigma = 2.48;
const int Eisr_bin = TMath::Nint((Eisr_max - Eisr_min) / sfw2d_sigma_nb / Eisr_sigma);

//ppIM:
const double ppIM_min = 200;
const double ppIM_max = 700;
const double ppIM_sigma = 2.30;
const int ppIM_bin = TMath::Nint((ppIM_max - ppIM_min) / sfw2d_sigma_nb / ppIM_sigma);

//sfw1d
const double xmin = 770;
const double xmax = 800;
const int xbins = 60; 

//crx3pi
const double hmin = 740; //770, 700
const double hmax = 820; //800, 820
const int hbins = TMath::Nint((hmax - hmin) / 0.25 / IM3pi_sigma);


TList *HIM3pi_fit = new TList(); // IM3pi distr. for fit omega parameters
TList *HSFW2D = new TList(); // Eisr vs. ppIM distr. for MC normalization
TList *HSFW1D = new TList(); // IM3pi distr. for signal MC tuning
TList *HSIG = new TList(); // IM3pi distr. signal true and generated
TList *HIM3pi_crx = new TList(); // IM3pi distr. for crx3pi obs.
TRandom *rnd=0;


// methods
void fillHist() {

  // data and MC background
  TIter next_tree(f_cut -> GetListOfKeys());

  TString objnm_tree, classnm_tree;

  double m3pi = 0., m3pi_true = 0., m3pi_corred = 0.;
  double Eisr = 0.;
  double ppIM = 0.;

  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {// start tree while loop

    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();

    cout << "classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;

    TTree *tree_tmp = (TTree*)f_cut -> Get(objnm_tree);
    //cout << tree_tmp -> GetName() << endl;

    // fill histos
    TH1D * h1d_tmp_crx = new TH1D("h1d_IM3pi_" + objnm_tree + "_CRX", "", hbins, hmin, hmax);
    h1d_tmp_crx -> Sumw2();

    TH1D * h1d_tmp_sfw = new TH1D("h1d_IM3pi_sfw_" + objnm_tree, "", xbins, xmin, xmax);
    h1d_tmp_sfw -> Sumw2();

    TH1D * h1d_tmp = new TH1D("h1d_IM3pi_" + objnm_tree, "", IM3pi_bin, IM3pi_min, IM3pi_max);
    h1d_tmp -> Sumw2();

    TH1D * h1d_tmp_true = new TH1D("h1d_IM3pi_" + objnm_tree + "_TRUE", "", IM3pi_bin, IM3pi_min, IM3pi_max);
    h1d_tmp_true -> Sumw2();

    TH2D * h2d_tmp = new TH2D("h2d_sfw_" + objnm_tree, "", ppIM_bin, ppIM_min, ppIM_max, Eisr_bin, Eisr_min, Eisr_max);
    h2d_tmp -> Sumw2();

    for (Int_t irow = 0; irow < tree_tmp -> GetEntries(); irow++) {// loop chain

      tree_tmp -> GetEntry(irow);
      
      m3pi = tree_tmp -> GetLeaf("Br_IM3pi_7C") -> GetValue(0);
      m3pi_true = tree_tmp -> GetLeaf("Br_IM3pi_true") -> GetValue(0);

      Eisr = tree_tmp -> GetLeaf("Br_Eisr") -> GetValue(0);
      ppIM = tree_tmp -> GetLeaf("Br_ppIM") -> GetValue(0);

      // filling
      h1d_tmp -> Fill(m3pi);

      h1d_tmp_sfw -> Fill(m3pi);

      h1d_tmp_crx -> Fill(m3pi);

      h1d_tmp_true -> Fill(m3pi_true);

      h2d_tmp -> Fill(ppIM, Eisr);
      
    }

    // scale EEG by eeg_lsf
    if (objnm_tree == "TEEG") {
       h1d_tmp -> Scale(eeg_lsf);
       h2d_tmp -> Scale(eeg_lsf);
    }

    HIM3pi_fit -> Add(h1d_tmp);
    HIM3pi_fit -> Add(h1d_tmp_true);
    
    HSFW1D -> Add(h1d_tmp_sfw);
    HSFW2D -> Add(h2d_tmp);

    HIM3pi_crx -> Add(h1d_tmp_crx);
    
  }

  // hsig_true

  m3pi_true = 0.;

  TH1D * hsig_true = new TH1D("hsig_true", "", IM3pi_bin, IM3pi_min, IM3pi_max);
  hsig_true -> Sumw2();

  TH1D * hsig_true_crx = new TH1D("h1d_IM3pi_TISR3PI_SIG_TRUE_CRX", "", hbins, hmin, hmax);
  hsig_true_crx -> Sumw2();

  TTree *TISR3PI_SIG = (TTree*)f_cut -> Get("TISR3PI_SIG");

  for (Int_t irow = 0; irow < TISR3PI_SIG -> GetEntries(); irow++) {

    TISR3PI_SIG -> GetEntry(irow);

    m3pi_true = TISR3PI_SIG -> GetLeaf("Br_IM3pi_true") -> GetValue(0);

    //cout << m3pi_true << endl;

    hsig_true -> Fill(m3pi_true);
    hsig_true_crx -> Fill(m3pi_true);
    
  }

  HSIG -> Add(hsig_true);
  HIM3pi_crx -> Add(hsig_true_crx);
  
  // hsig_gen
  TTree *TISR3PI_SIG_GEN = (TTree*)f_gen -> Get("TISR3PI_SIG_GEN");

  double m3pi_gen = 0.;
  
  TH1D * hsig_gen = new TH1D("hsig_gen", "", IM3pi_bin, IM3pi_min, IM3pi_max);
  hsig_gen -> Sumw2();

  TH1D * hsig_gen_crx = new TH1D("h1d_IM3pi_TISR3PI_SIG_GEN_CRX", "", hbins, hmin, hmax);
  hsig_gen_crx -> Sumw2();

  for (Int_t irow = 0; irow < TISR3PI_SIG_GEN -> GetEntries(); irow++) {
    
    TISR3PI_SIG_GEN -> GetEntry(irow);
    
    m3pi_gen = TISR3PI_SIG_GEN -> GetLeaf("Br_IM3pi_gen") -> GetValue(0);
    
    hsig_gen -> Fill(m3pi_gen);
    hsig_gen_crx -> Fill(m3pi_gen);
    
  }

  HSIG -> Add(hsig_gen);
  HIM3pi_crx -> Add(hsig_gen_crx);
  
}

