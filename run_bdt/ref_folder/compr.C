#include "compr.h"

#include "../plot.h"
#include "../hist.h"
#include "../header_method/method.h"

int compr(){

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
  const TString gen_file = "/home/bo/Desktop/analysis/chains_norm/sig_gen.root";
  const TString tree_file = "/home/bo/Desktop/analysis/crx3pi/output_norm/tree_cut0.root";

  if (!infile || gen_file->IsZombie()) {
    cerr << "ERROR: Cannot open " << gen_file << endl;
    return 1;
  }

  if (!infile || tree_file->IsZombie()) {
    cerr << "ERROR: Cannot open " << tree_file << endl;
    return 1;
  }

  getObj(gen_file);
  getObj(tree_file);
  
  // get generated signal
  TFile* intree_gen = new TFile(gen_file);
  
   TTree * ALLCHAIN_GEN = static_cast<TTree*>(intree_gen -> Get("ALLCHAIN_GEN")); 

   const double evnb_sig_gen = ALLCHAIN_GEN -> GetEntries(); //number of generated signal events

   cout << evnb_sig_gen << endl;

   // get MC recon.
   
   TFile* intree = new TFile(tree_file);
  
   TIter next_tree(intree -> GetListOfKeys());

   //TString objnm_tree, classnm_tree;

  i = 0;
  
  while ( (key = (TKey *) next_tree() ) ) {
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    //key -> GetSeekKey();
    
    //cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  
  // check branches
  //TDATA -> GetListOfLeaves() -> Print();

  //bkg
  const int TLSize = 9;
  char name[TLSize], title[TLSize];

  TTree *TDATA = static_cast<TTree*>(intree -> Get("TDATA"));
  TTree *TEEG = static_cast<TTree*>(intree -> Get("TEEG"));
  TTree *TOMEGAPI = static_cast<TTree*>(intree -> Get("TOMEGAPI"));
  TTree *TKSL = static_cast<TTree*>(intree -> Get("TKSL"));
  TTree *TKPM = static_cast<TTree*>(intree -> Get("TKPM"));
  TTree *TRHOPI = static_cast<TTree*>(intree -> Get("TRHOPI"));
  TTree *TETAGAM = static_cast<TTree*>(intree -> Get("TETAGAM"));
  TTree *TBKGREST = static_cast<TTree*>(intree -> Get("TBKGREST"));
  TTree *TISR3PI_SIG = static_cast<TTree*>(intree -> Get("TISR3PI_SIG"));

  TTree *TrList[TLSize] = {TDATA, TEEG, TOMEGAPI, TKSL, TKPM, TRHOPI, TETAGAM, TBKGREST, TISR3PI_SIG};
  const TString TrNm[TLSize] = {"data", "eeg", "omegapi" , "ksl", "kpm", "rhopi" , "etagam" , "bkgrest", "isr3pi"};


  int color_list[TLSize] = {1, 6, 7, 28, 46, 42, 3, 37, 4};
  
  // Create arrays of histograms.
  TObjArray *Hlist = new TObjArray();
  TObjArray *H2dlist = new TObjArray();

  TH1D *h;
  TH2D *h2d;
  TH2D *h2d_pchi2;

  // fixed variables
  double var_value = 0., IM3pi = 0., IM3pi_true = 0., IM3pi_det = 0.;
  double evnb_data = 0., evnb_eeg = 0., evnb_omegapi = 0., evnb_ksl = 0., evnb_kpm = 0., evnb_rhopi = 0., evnb_etagam = 0., evnb_bkgrest = 0., evnb_isr3pi = 0.;
  double betapi0 = 0., IM2pi = 0.;
  double chi2 = 0., pvalue = 0.;
  
  //if (var_nm.Contains("IM3pi_7C")) cout << var_nm << endl;

  for (int i = 0; i < TLSize; i ++) {// start MC type loop

    //sprintf(name,"h%d",i);
    //sprintf(title,"histo nr:%d",i);
    //cout << TrNm[i] << ", histo nr = " << i << ", name = " << TrList[i] -> GetName() << endl;
 
    h = new TH1D("hist_" + TrNm[i], "", binsize, var_min, var_max);
    h2d = new TH2D("h2d_discrp_" + TrNm[i], "", 200, 0.25, 0.65, 200, 0.3, 1.);
    h2d_pchi2 = new TH2D("h2d_pchi2_" + TrNm[i], "", 200, 0., 1., 200, 0., 100.);

    //cout << TLSize << endl;
    //cout << TrNm[i] << endl;
    //cout << "hist_" + TrNm[i] << ", i = " << i << endl;
    
    for (Int_t irow = 0; irow < TrList[i] -> GetEntries(); irow++) {

      //if (irow > 1e3) break;
      
      TrList[i] -> GetEntry(irow);
      
      var_value = TrList[i] -> GetLeaf("Br_" + var_nm) -> GetValue(0);

      //cout << var_value << endl;
      
      if (TrNm[i] == "data") evnb_data ++;
      else if (TrNm[i] == "eeg") evnb_eeg ++;
      else if (TrNm[i] == "omegapi") evnb_omegapi ++;
      else if (TrNm[i] == "ksl") evnb_ksl ++;
      else if (TrNm[i] == "kpm") evnb_kpm ++;
      else if (TrNm[i] == "rhopi") evnb_rhopi ++;
      else if (TrNm[i] == "etagam") evnb_etagam ++;
      else if (TrNm[i] == "bkgrest") evnb_bkgrest ++;
      else if (TrNm[i] == "isr3pi") evnb_isr3pi ++;
      
      IM3pi_true = TrList[i] -> GetLeaf("Br_IM3pi_true") -> GetValue(0);
      IM3pi = TrList[i] -> GetLeaf("Br_IM3pi_7C") -> GetValue(0);
      IM3pi_det = DetectorEvent(TMath::Abs(IM3pi_true));
      
      betapi0 = TrList[i] -> GetLeaf("Br_betapi0") -> GetValue(0);
      IM2pi = TrList[i] -> GetLeaf("Br_ppIM") -> GetValue(0) * 1e-3;
      
      chi2 = TrList[i] -> GetLeaf("Br_lagvalue_min_7C") -> GetValue(0);
      pvalue = TrList[i] -> GetLeaf("Br_pvalue") -> GetValue(0);

      h2d -> Fill(IM2pi, betapi0);
      h2d_pchi2 -> Fill(pvalue, chi2);
      h -> Fill(var_value);
	
    }
    
    format_h(h, color_list[i], 2);
    Hlist -> Add(h);
    h -> Draw();
      
    H2dlist -> Add(h2d);
    H2dlist -> Add(h2d_pchi2);
    

    
  }

  
  //double evnb_data = 0., evnb_eeg = 0., evnb_omegapi = 0., evnb_ksl = 0., evnb_kpm = 0., evnb_rhopi = 0., evnb_etagam = 0., evnb_bkgrest = 0., evnb_isr3pi = 0.;

  //double evnb_isr3pi = evnb_isr3pi;
  evnb_eeg = evnb_eeg * 2.;
  const double evnb_mcrest = evnb_kpm + evnb_rhopi + evnb_bkgrest;
  const double evnb_mcsum = evnb_eeg + evnb_omegapi + evnb_ksl + evnb_etagam + evnb_mcrest + evnb_isr3pi;
  const double evnb_bkgsum = evnb_eeg + evnb_omegapi + evnb_ksl + evnb_etagam + evnb_mcrest;
  const double sb_ratio = evnb_isr3pi / TMath::Sqrt(evnb_mcsum);
  const double s_frac = evnb_isr3pi / evnb_mcsum * 100.;
  const double b_frac = evnb_bkgsum / evnb_mcsum * 100.;

  
  ofstream myfile;
  myfile.open ("./output/compr.txt");
  
  myfile << "Number of events\n"
	 << "data = " << evnb_data << "\n"
	 << "1. eeg = " << evnb_eeg << "\n"
	 << "2. omegapi = " << evnb_omegapi << "\n"
	 << "3. ksl = " << evnb_ksl << "\n"
	 << "4. etagam = " << evnb_etagam << "\n"
	 << "5. mcrest = " << evnb_mcrest << "\n"
	 << "\tbkgrest = " << evnb_bkgrest << "\n"
	 << "\tkpm = " << evnb_kpm << "\n"
	 << "\trhopi = " << evnb_rhopi << "\n"
	 << "6. isr3pi (gen) = " << evnb_isr3pi << " (" << evnb_sig_gen << ")\n"
	 << "mcsum = " << evnb_mcsum << "\n"
	 << "s_frac [%] = " << s_frac << "\n"
	 << "b_frac [%] = " << b_frac << "\n"  
	 << "sb_ratio = " << sb_ratio << "\n";

  TH1D *hist_eeg = (TH1D *) Hlist -> At(1);
  TH1D *hist_omegapi = (TH1D *) Hlist -> At(2);
  TH1D *hist_ksl = (TH1D *) Hlist -> At(3);
  TH1D *hist_kpm = (TH1D *) Hlist -> At(4);
  TH1D *hist_rhopi = (TH1D *) Hlist -> At(5);
  TH1D *hist_etagam = (TH1D *) Hlist -> At(6);
  TH1D *hist_bkgrest = (TH1D *) Hlist -> At(7);
  TH1D *hist_isr3pi = (TH1D *) Hlist -> At(8);
  
  //checkArray(Hlist);

  
  // MC rest, merge bkgrest, kpm and rhopi
  TH1D *hist_mcrest = (TH1D*) hist_bkgrest -> Clone();
  hist_mcrest -> Add(hist_kpm, 1.);
  hist_mcrest -> Add(hist_rhopi, 1.);
  hist_mcrest -> SetName("hist_mcrest");

  Hlist -> Add(hist_mcrest);

  
  // add 2d
  //checkArray(H2dlist);

  TH2D *h2d_discrp_eeg = (TH2D *) H2dlist -> At(2);
  TH2D *h2d_discrp_omegapi = (TH2D *) H2dlist -> At(4);
  TH2D *h2d_discrp_ksl = (TH2D *) H2dlist -> At(6);
  TH2D *h2d_discrp_kpm = (TH2D *) H2dlist -> At(8);  
  TH2D *h2d_discrp_rhopi = (TH2D *) H2dlist -> At(10);  
  TH2D *h2d_discrp_etagam = (TH2D *) H2dlist -> At(12);
  TH2D *h2d_discrp_bkgrest = (TH2D *) H2dlist -> At(14); 
  TH2D *h2d_discrp_isr3pi = (TH2D *) H2dlist -> At(16);
  
  //TH2D* h2d_discrp_mcrest = (TH2D*) h2d_discrp_bkgrest -> Clone();
  //h2d_discrp_mcrest -> Add(h2d_discrp_kpm, 1.);
  //h2d_discrp_mcrest -> Add(h2d_discrp_rhopi, 1.);
  //h2d_discrp_mcrest -> SetName("h2d_discrp_mcrest");

  TH2D* h2d_discrp_mcrest = (TH2D*) h2d_discrp_rhopi -> Clone();
  
  H2dlist -> Add(h2d_discrp_mcrest);
  
  // Scale histos
  cout << "Scaling Factors (obtained from 2D fit, E_{gamma} v.s. M_{2pi})" << "\n"
       << "1: eeg     = " << sfw2d_eeg     << "\n"
       << "2: isr3pi  = " << sfw2d_isr3pi  << "\n"
       << "3: omegapi = " << sfw2d_omegapi << "\n"
       << "4: etagam  = " << sfw2d_etagam  << "\n"
       << "5: ksl     = " << sfw2d_ksl     << "\n"
       << "6: mcrest  = " << sfw2d_mcrest  << "\n\n";

  cout << "Signal scaling factor from sfw1d correction\n"
       << "isr3pi = " << sfw1d_isr3pi << "\n";

  // eeg
  TH1D *hist_eeg_sc = (TH1D*) hist_eeg -> Clone();
  hist_eeg_sc -> Scale(2.);
  hist_eeg_sc -> SetName("hist_eeg_sc");

  TH2D * h2d_discrp_eeg_sc = (TH2D*) h2d_discrp_eeg -> Clone();
  h2d_discrp_eeg_sc -> Scale(2.);
  h2d_discrp_eeg_sc -> SetName("h2d_discrp_eeg_sc");

  // isr3pi
  TH1D * hist_isr3pi_sc = (TH1D*) hist_isr3pi -> Clone();
  hist_isr3pi_sc -> Scale(1); //sfw1d_isr3pi
  hist_isr3pi_sc -> SetName("hist_isr3pi_sc");

  TH2D * h2d_discrp_isr3pi_sc = (TH2D*) h2d_discrp_isr3pi -> Clone();
  h2d_discrp_isr3pi_sc -> Scale(sfw1d_isr3pi);
  h2d_discrp_isr3pi_sc -> SetName("h2d_discrp_isr3pi_sc");

  cout << "sfw1d_isr3pi = " << sfw1d_isr3pi << endl;
  
  // omegapi
  TH1D * hist_omegapi_sc = (TH1D*) hist_omegapi -> Clone();
  hist_omegapi_sc -> Scale(1.);
  hist_omegapi_sc -> SetName("hist_omegapi_sc");

  TH2D * h2d_discrp_omegapi_sc = (TH2D*) h2d_discrp_omegapi -> Clone();
  h2d_discrp_omegapi_sc -> Scale(1.);
  h2d_discrp_omegapi_sc -> SetName("h2d_discrp_omegapi_sc");

  // etagam
  TH1D * hist_etagam_sc = (TH1D*) hist_etagam -> Clone();
  hist_etagam_sc -> Scale(1.);
  hist_etagam_sc -> SetName("hist_etagam_sc");

  TH2D * h2d_discrp_etagam_sc = (TH2D*) h2d_discrp_etagam -> Clone();
  h2d_discrp_etagam_sc -> Scale(1.);
  h2d_discrp_etagam_sc -> SetName("h2d_discrp_etagam_sc");

  // ksl
  TH1D * hist_ksl_sc = (TH1D*) hist_ksl -> Clone();
  hist_ksl_sc -> Scale(1.);
  hist_ksl_sc -> SetName("hist_ksl_sc");

  TH2D * h2d_discrp_ksl_sc = (TH2D*) h2d_discrp_ksl -> Clone();
  h2d_discrp_ksl_sc -> Scale(1.);
  h2d_discrp_ksl_sc -> SetName("h2d_discrp_ksl_sc");

  // mcrest
  TH1D * hist_mcrest_sc = (TH1D*) hist_mcrest -> Clone();
  hist_mcrest_sc -> Scale(1.);
  hist_mcrest_sc -> SetName("hist_mcrest_sc");

  TH2D * h2d_discrp_mcrest_sc = (TH2D*) h2d_discrp_mcrest -> Clone();
  h2d_discrp_mcrest_sc -> Scale(1.);
  h2d_discrp_mcrest_sc -> SetName("h2d_discrp_mcrest_sc");

  // bkgsum
  TH1D* hist_bkgsum_sc = (TH1D*) hist_eeg_sc -> Clone();
  hist_bkgsum_sc -> Add(hist_omegapi_sc, 1.);
  hist_bkgsum_sc -> Add(hist_ksl_sc, 1.);
  hist_bkgsum_sc -> Add(hist_etagam_sc, 1.);
  hist_bkgsum_sc -> Add(hist_mcrest_sc, 1.);
  hist_bkgsum_sc -> SetName("hist_bkgsum_sc");

  format_h(hist_bkgsum_sc, 6, 2);

  TH2D* h2d_discrp_bkgsum_sc = (TH2D*) h2d_discrp_eeg_sc -> Clone();
  h2d_discrp_bkgsum_sc -> Add(h2d_discrp_omegapi_sc, 1.);
  h2d_discrp_bkgsum_sc -> Add(h2d_discrp_ksl_sc, 1.);
  h2d_discrp_bkgsum_sc -> Add(h2d_discrp_etagam_sc, 1.);
  h2d_discrp_bkgsum_sc -> Add(h2d_discrp_mcrest_sc, 1.);
  h2d_discrp_bkgsum_sc -> SetName("h2d_discrp_bkgsum_sc");

  // nomalization and add to the histo list
  
  Hlist -> Add(hist_eeg_sc);
  Hlist -> Add(hist_isr3pi_sc);
  Hlist -> Add(hist_omegapi_sc);
  Hlist -> Add(hist_etagam_sc);
  Hlist -> Add(hist_ksl_sc);
  Hlist -> Add(hist_mcrest_sc);
  Hlist -> Add(hist_bkgsum_sc);

  H2dlist -> Add(h2d_discrp_eeg_sc);
  H2dlist -> Add(h2d_discrp_isr3pi_sc);
  H2dlist -> Add(h2d_discrp_omegapi_sc);
  H2dlist -> Add(h2d_discrp_etagam_sc);
  H2dlist -> Add(h2d_discrp_ksl_sc);
  H2dlist -> Add(h2d_discrp_mcrest_sc);
  H2dlist -> Add(h2d_discrp_bkgsum_sc);
  
  // save
  TFile *f_out = new TFile("./output_" + var_nm + "/hist_" + var_nm + ".root", "recreate");

  Hlist -> Write("Hlist",1);
  H2dlist -> Write("H2dlist",1);
  //hEisr_gen -> Write();
  //hangle_isr_gen -> Write();
    
  f_out -> Close();
  
  
  return 0;
  
}


