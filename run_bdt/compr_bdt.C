#include "../header_method/method.h"
#include "../header_bdt/compr.h"   // defines var_nm, binsize, var_min, var_max
#include "../header_plot/plot.h"

TRandom3 *rnd = new TRandom3();

int compr_bdt() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);

  TH1::SetDefaultSumw2();

  const TString tree_file_nm = "/home/kloe/Desktop/input_bdt_TDATA_norm/cut/tree_pre_bdt.root";
  
  TFile* tree_file = new TFile(tree_file_nm);
  if (!tree_file || tree_file->IsZombie()) {
    cerr << "ERROR: Cannot open " << tree_file_nm << endl;   // fixed
    return 1;
  }

  getObj(tree_file);   // ensure defined in method.h

  const int TLSize = 9;
  TTree *TDATA      = static_cast<TTree*>(tree_file->Get("TDATA"));
  TTree *TEEG       = static_cast<TTree*>(tree_file->Get("TEEG"));
  TTree *TOMEGAPI   = static_cast<TTree*>(tree_file->Get("TOMEGAPI"));
  TTree *TKSL       = static_cast<TTree*>(tree_file->Get("TKSL"));
  TTree *TKPM       = static_cast<TTree*>(tree_file->Get("TKPM"));
  TTree *TRHOPI     = static_cast<TTree*>(tree_file->Get("TRHOPI"));
  TTree *TETAGAM    = static_cast<TTree*>(tree_file->Get("TETAGAM"));
  TTree *TBKGREST   = static_cast<TTree*>(tree_file->Get("TBKGREST"));
  TTree *TISR3PI_SIG = static_cast<TTree*>(tree_file->Get("TISR3PI_SIG"));

  TTree *TrList[TLSize] = {TDATA, TEEG, TOMEGAPI, TKSL, TKPM, TRHOPI, TETAGAM, TBKGREST, TISR3PI_SIG};
  const TString TrNm[TLSize] = {"data", "eeg", "omegapi", "ksl", "kpm", "rhopi", "etagam", "bkgrest", "isr3pi"};
  int color_list[TLSize] = {1, 6, 7, 28, 46, 42, 3, 37, 4};

  TList *Hlist = new TList();   // use TList for proper FindObject

  double var_value = 0.;

  for (int i = 0; i < TLSize; ++i) {
    if (!TrList[i]) {
      cerr << "WARNING: Tree " << TrNm[i] << " not found!" << endl;
      continue;
    }

    TH1D *h = new TH1D("hist_" + TrNm[i], "", binsize, var_min, var_max);
    h->Sumw2();

    // Set branch address for efficiency
    TrList[i]->SetBranchAddress("Br_" + var_nm, &var_value);

    Long64_t nentries = TrList[i]->GetEntries();
    for (Long64_t irow = 0; irow < nentries; ++irow) {
      TrList[i]->GetEntry(irow);
      h->Fill(var_value);
    }

    format_h(h, color_list[i], 2);
    Hlist->Add(h);
    cout << h->GetName() << endl;
  }

  // Helper to retrieve histograms from the list
  auto getHist = [&](const char* name) -> TH1D* {
    TH1D* h = (TH1D*)Hlist->FindObject(name);
    if (!h) cerr << "WARNING: histogram " << name << " not found in list" << endl;
    return h;
  };

  TH1D *hist_eeg     = getHist("hist_eeg");
  TH1D *hist_omegapi = getHist("hist_omegapi");
  TH1D *hist_ksl     = getHist("hist_ksl");
  TH1D *hist_kpm     = getHist("hist_kpm");
  TH1D *hist_rhopi   = getHist("hist_rhopi");
  TH1D *hist_etagam  = getHist("hist_etagam");
  TH1D *hist_bkgrest = getHist("hist_bkgrest");
  TH1D *hist_isr3pi  = getHist("hist_isr3pi");

  // MC rest (background without signal)
  TH1D *hist_mcrest = (TH1D*)hist_bkgrest->Clone();
  hist_mcrest->Add(hist_kpm, 1.);
  hist_mcrest->Add(hist_rhopi, 1.);
  hist_mcrest->SetName("hist_mcrest");
  Hlist->Add(hist_mcrest);

  // Clone histograms for scaling (no scaling applied yet)
  TH1D *hist_eeg_sc      = (TH1D*)hist_eeg->Clone();
  hist_eeg_sc->SetName("hist_eeg_sc");
  TH1D *hist_isr3pi_sc   = (TH1D*)hist_isr3pi->Clone();
  hist_isr3pi_sc->SetName("hist_isr3pi_sc");
  TH1D *hist_omegapi_sc  = (TH1D*)hist_omegapi->Clone();
  hist_omegapi_sc->SetName("hist_omegapi_sc");
  TH1D *hist_etagam_sc   = (TH1D*)hist_etagam->Clone();
  hist_etagam_sc->SetName("hist_etagam_sc");
  TH1D *hist_ksl_sc      = (TH1D*)hist_ksl->Clone();
  hist_ksl_sc->SetName("hist_ksl_sc");
  TH1D *hist_mcrest_sc   = (TH1D*)hist_mcrest->Clone();
  hist_mcrest_sc->SetName("hist_mcrest_sc");

  // Background sum (EEG + all backgrounds)
  TH1D *hist_bkgsum_sc = (TH1D*)hist_eeg_sc->Clone();
  hist_bkgsum_sc->Add(hist_omegapi_sc, 1.);
  hist_bkgsum_sc->Add(hist_ksl_sc, 1.);
  hist_bkgsum_sc->Add(hist_etagam_sc, 1.);
  hist_bkgsum_sc->Add(hist_mcrest_sc, 1.);
  hist_bkgsum_sc->SetName("hist_bkgsum_sc");
  format_h(hist_bkgsum_sc, 6, 2);
  
  Hlist->Add(hist_eeg_sc);
  Hlist->Add(hist_isr3pi_sc);
  Hlist->Add(hist_omegapi_sc);
  Hlist->Add(hist_etagam_sc);
  Hlist->Add(hist_ksl_sc);
  Hlist->Add(hist_mcrest_sc);
  Hlist->Add(hist_bkgsum_sc);

  // Create output directory
  TString out_dir = "../output_" + TString(var_nm);
  gSystem->mkdir(out_dir, kTRUE);   // creates directory and parents if needed

  TString outfile_name = out_dir + "/hist_" + var_nm + ".root";
  TFile *f_out = new TFile(outfile_name, "recreate");
  Hlist->Write("Hlist", TObject::kSingleKey);
  f_out->Close();
  delete Hlist;   
  tree_file->Close();
  return 0;
}
