#include "compr.h"
////#include "../crx3pi/crx3pi.h"
#include "../hist.h"
#include "../method.h"
#include "plot_discrp.h"


int plot_discrp(){

  //gROOT->SetBatch(kTRUE);  
  gErrorIgnoreLevel = kError;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TGaxis::SetMaxDigits(3);
  //gROOT->SetStyle("Plain");
  //gROOT->ForceStyle();

  TFile* intree = new TFile("./output_" + var_nm + "/hist_" + var_nm + ".root");
  
  TIter next_tree(intree -> GetListOfKeys());

  TString objnm_tree, classnm_tree;

  int i = 0;
  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  TObjArray *H2dlist = (TObjArray *)intree -> Get("H2dlist");
  checkArray(H2dlist);

  TH2D *h2d_discrp_data = (TH2D *) H2dlist -> At(0);
  TH2D *h2d_discrp_eeg = (TH2D *) H2dlist -> At(2);
  TH2D *h2d_discrp_rhopi = (TH2D *) H2dlist -> At(10);
  
  TH2D *h2d_discrp_eeg_sc = (TH2D *) H2dlist -> At(19);
  TH2D *h2d_discrp_isr3pi_sc = (TH2D *) H2dlist -> At(20);
  TH2D *h2d_discrp_omegapi_sc = (TH2D *) H2dlist -> At(21);
  TH2D *h2d_discrp_etagam_sc = (TH2D *) H2dlist -> At(22);
  TH2D *h2d_discrp_ksl_sc = (TH2D *) H2dlist -> At(23);
  TH2D *h2d_discrp_mcrest_sc = (TH2D *) H2dlist -> At(24);
  TH2D *h2d_discrp_bkgsum_sc = (TH2D *) H2dlist -> At(25);

  /*
    [0]: h2d_discrp_data
    [1]: h2d_pchi2_data
    [2]: h2d_discrp_eeg
    [3]: h2d_pchi2_eeg
    [4]: h2d_discrp_omegapi
    [5]: h2d_pchi2_omegapi
    [6]: h2d_discrp_ksl
    [7]: h2d_pchi2_ksl
    [8]: h2d_discrp_kpm
    [9]: h2d_pchi2_kpm
    [10]: h2d_discrp_rhopi
    [11]: h2d_pchi2_rhopi
    [12]: h2d_discrp_etagam
    [13]: h2d_pchi2_etagam
    [14]: h2d_discrp_bkgrest
    [15]: h2d_pchi2_bkgrest
    [16]: h2d_discrp_isr3pi
    [17]: h2d_pchi2_isr3pi
    [18]: h2d_discrp_rhopi
    [19]: h2d_discrp_eeg_sc
    [20]: h2d_discrp_isr3pi_sc
    [21]: h2d_discrp_omegapi_sc
    [22]: h2d_discrp_etagam_sc
    [23]: h2d_discrp_ksl_sc
    [24]: h2d_discrp_mcrest_sc
    [25]: h2d_discrp_bkgsum_sc
    */
    
  
  TH1D * projh2X_sc = h2d_discrp_eeg_sc -> ProjectionX();
  TH1D * projh2X = h2d_discrp_eeg -> ProjectionX();

  TH1D * projh2Y_sc = h2d_discrp_eeg_sc -> ProjectionY();
  
  int xbins = projh2X_sc -> GetNbinsX(); 
  int ybins = projh2Y_sc -> GetNbinsX();

  cout << "xbins = " << xbins << ", ybins = " << ybins << endl;

  // mcsum 
  TH2D* h2d_discrp_mcsum_sc = (TH2D*) h2d_discrp_bkgsum_sc -> Clone();
  h2d_discrp_mcsum_sc -> Add(h2d_discrp_isr3pi_sc, 1.);
  h2d_discrp_mcsum_sc -> SetName("hist_mcsum_sc");

  // MCsum no eta
  TH2D* h2d_discrp_mcsum_noeta_sc = (TH2D*) h2d_discrp_eeg_sc -> Clone();
  h2d_discrp_mcsum_noeta_sc -> Add(h2d_discrp_omegapi_sc, 1.);
  h2d_discrp_mcsum_noeta_sc -> Add(h2d_discrp_ksl_sc, 1.);
  h2d_discrp_mcsum_noeta_sc -> Add(h2d_discrp_mcrest_sc, 1.);
  h2d_discrp_mcsum_noeta_sc -> SetName("h2d_discrp_mcsum_noeta_sc");
  //h2d_discrp_mcsum_noeta_sc -> Add(h2d_discrp_isr3pi_sc, 1.);

  // make a projection on y-axis
  TH1D * projh2Y = h2d_discrp_mcsum_noeta_sc -> ProjectionX();
  //projh2Y -> Draw();
    
  double evnt_mcsum_noeta = 0.;
  double evnt_data = 0.;
  double evnt_eeg = 0.;

  for (int i = 1; i <= xbins; i ++) {

    for (int j = 1; j <= ybins; j ++) {

      //cout << i << ", " << j << endl;  
      evnt_data += h2d_discrp_data -> GetBinContent(i,j);
      evnt_eeg += h2d_discrp_eeg_sc -> GetBinContent(i,j);
      evnt_mcsum_noeta += h2d_discrp_mcsum_sc -> GetBinContent(i,j);

    }
    
  }


  cout << "evnt_data = " << evnt_data << "\n"
       << "evnt_eeg = " << evnt_eeg << "\n"
       << "evnt_mcsum_noeta = " << evnt_mcsum_noeta << "\n";

  //// plot

 
  // data
  
  TPaveText *pt1;
  TCanvas *cv_data_discrp = plot_cv("cv_isr3pi_data", h2d_discrp_data, pt1, "(a)"); // Data
  //projh2X_sc -> Draw();
  //projh2X -> Draw("Same");
  // mcsum noeta

  // bkg no eta
   
  //TCanvas *cv_mcsum_noeta_discrp = new TCanvas("cv_mcsum_noeta_discrp", " ", 700, 700);
  //cv_mcsum_noeta_discrp -> SetBottomMargin(0.15);//0.007
  //cv_mcsum_noeta_discrp -> SetLeftMargin(0.15);
  //cv_mcsum_noeta_discrp -> SetRightMargin(0.15);

  TPaveText *pt2;
  TCanvas *cv_mcsum_noeta_discrp = plot_cv("cv_mcsum_noeta", h2d_discrp_mcsum_noeta_sc, pt2, "(b)"); // MC sum no etagam
  
  // etagam
  TPaveText *pt3;
  TCanvas *cv_etagam_discrp = plot_cv("cv_etagam", h2d_discrp_etagam_sc, pt3, "(c)"); // etagam
  
  // signal
  TPaveText *pt4;
  TCanvas *cv_isr3pi_discrp = plot_cv("cv_isr3pi_discrp", h2d_discrp_isr3pi_sc, pt4, "(d)"); //#pi^{+}#pi^{-}#pi^{0}#gamma
  //h2d_discrp_isr3pi_sc -> GetYaxis() -> SetRangeUser(0.6, 1.);
  //h2d_discrp_isr3pi_sc -> GetXaxis() -> SetRangeUser(0.25, 0.5);
  
  // rhopi
  //TPaveText *pt5;
  //TCanvas *cv_rhopi_discrp = new TCanvas("cv_rhopi_discrp", " ", 700, 700);
  //h2d_attri(h2d_discrp_rhopi, pt5, "#rho#pi");
  
  
  // save
  cv_isr3pi_discrp -> SaveAs("./output_discrp/cv_isr3pi_discrp.pdf");
  cv_etagam_discrp -> SaveAs("./output_discrp/cv_etagam_discrp.pdf");
  cv_mcsum_noeta_discrp -> SaveAs("./output_discrp/cv_mcsum_noeta_discrp.pdf");
  cv_data_discrp -> SaveAs("./output_discrp/cv_data_discrp.pdf");
  
  return 0;
  
}

