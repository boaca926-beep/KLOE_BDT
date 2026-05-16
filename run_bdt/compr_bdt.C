#include "../header_method/method.h"
#include "../header_bdt/compr.h"
#include "../header_plot/plot.h"

TRandom3 *rnd = new TRandom3();

int compr_bdt(){

  //gROOT->SetBatch(kTRUE);  
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);

  // switch on histogram errors
  TH1::SetDefaultSumw2();

  // random generator
  rnd=new TRandom3();

  // Inspect the input tree
  //cout << infile_nm << endl;
  const TString tree_file_nm = "/home/kloe/Desktop/input_bdt_TDATA_chain/cut/tree_pre_bdt.root";

  TFile* tree_file = new TFile(tree_file_nm);
  
  if (!tree_file || tree_file->IsZombie()) {
    cerr << "ERROR: Cannot open " << tree_file << endl;
    return 1;
  }

  getObj(tree_file);

  const int TLSize = 9;
  char name[TLSize], title[TLSize];

  TTree *TDATA = static_cast<TTree*>(tree_file -> Get("TDATA"));
  TTree *TEEG = static_cast<TTree*>(tree_file -> Get("TEEG"));
  TTree *TOMEGAPI = static_cast<TTree*>(tree_file -> Get("TOMEGAPI"));
  TTree *TKSL = static_cast<TTree*>(tree_file -> Get("TKSL"));
  TTree *TKPM = static_cast<TTree*>(tree_file -> Get("TKPM"));
  TTree *TRHOPI = static_cast<TTree*>(tree_file -> Get("TRHOPI"));
  TTree *TETAGAM = static_cast<TTree*>(tree_file -> Get("TETAGAM"));
  TTree *TBKGREST = static_cast<TTree*>(tree_file -> Get("TBKGREST"));
  TTree *TISR3PI_SIG = static_cast<TTree*>(tree_file -> Get("TISR3PI_SIG"));

  TTree *TrList[TLSize] = {TDATA, TEEG, TOMEGAPI, TKSL, TKPM, TRHOPI, TETAGAM, TBKGREST, TISR3PI_SIG};
  const TString TrNm[TLSize] = {"data", "eeg", "omegapi" , "ksl", "kpm", "rhopi" , "etagam" , "bkgrest", "isr3pi"};


  int color_list[TLSize] = {1, 6, 7, 28, 46, 42, 3, 37, 4};

  TObjArray *Hlist = new TObjArray();

  TH1D *h;

  double var_value = 0.;

  for (int i = 0; i < TLSize; i ++) {// start MC type loop

    h = new TH1D("hist_" + TrNm[i], "", binsize, var_min, var_max);

    for (Int_t irow = 0; irow < TrList[i] -> GetEntries(); irow++) {

      //if (irow > 1e3) break;
      
      TrList[i] -> GetEntry(irow);
      
      var_value = TrList[i] -> GetLeaf("Br_" + var_nm) -> GetValue(0);

      //cout << var_value << endl;

      h->Fill(var_value);
      
    }

    format_h(h, color_list[i], 2);
    Hlist -> Add(h);
    h -> Draw("same");
    
  }
  
  return 0;

}
