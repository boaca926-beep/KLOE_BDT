#include "../header_plot/plot.h"
#include "../header_method/method.h"
#include "../header_bdt/plot2d.h"

int get2Dhist() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //gROOT->SetBatch(kFALSE); //kFALSE
    
  cout << "Creating 2D histos... " << infile_nm << "\n" << endl;

  // get histos
  TFile* infile = new TFile(infile_nm);
  getObj(infile);

  TList *HSFW2D = (TList*)gROOT->FindObject("HSFW2D");
  if (!HSFW2D) { cerr << "HSFW2D not found" << endl; return 1; }
  checkList(HSFW2D);
 
  TH2D *h2d_sfw_TISR3PI_SIG = (TH2D *)HSFW2D->FindObject("h2d_sfw_TISR3PI_SIG");
  TH2D *h2d_sfw_TOMEGAPI = (TH2D *)HSFW2D->FindObject("h2d_sfw_TOMEGAPI");
  TH2D *h2d_sfw_TKPM = (TH2D *)HSFW2D->FindObject("h2d_sfw_TKPM");
  TH2D *h2d_sfw_TKSL = (TH2D *)HSFW2D->FindObject("h2d_sfw_TKSL");
  TH2D *h2d_sfw_TRHOPI = (TH2D *)HSFW2D->FindObject("h2d_sfw_TRHOPI");
  TH2D *h2d_sfw_TETAGAM = (TH2D *)HSFW2D->FindObject("h2d_sfw_TETAGAM");
  TH2D *h2d_sfw_TBKGREST = (TH2D *)HSFW2D->FindObject("h2d_sfw_TBKGREST");
  if (!h2d_sfw_TBKGREST) {
    std::cerr << "ERROR: histogram h2d_sfw_TBKGREST not found!" << std::endl;
    return 1;
}
  TH2D *h2d_sfw_TDATA = (TH2D *)HSFW2D->FindObject("h2d_sfw_TDATA");
  TH2D *h2d_sfw_TEEG = (TH2D *)HSFW2D->FindObject("h2d_sfw_TEEG");
  
  // mc rest
  TH2D *hmcrest = (TH2D*) h2d_sfw_TBKGREST -> Clone();
  hmcrest -> Add(h2d_sfw_TKPM, 1.);
  hmcrest -> Add(h2d_sfw_TRHOPI, 1.);
  hmcrest -> SetName("hmcrest");
  
  // bkg sum no etagam
  TH2D *hbkgsum_noeta = (TH2D *) h2d_sfw_TEEG -> Clone();
  //hbkgsum_noeta -> Add(h2d_sfw_TISR3PI_SIG, 1.);
  hbkgsum_noeta -> Add(h2d_sfw_TKSL, 1.);
  hbkgsum_noeta -> Add(h2d_sfw_TOMEGAPI, 1.);
  hbkgsum_noeta -> Add(hmcrest, 1.);
  hbkgsum_noeta -> SetName("hbkgsum_noeta");

  // bkg sum
  TH2D *hbkgsum = (TH2D*) hbkgsum_noeta -> Clone();
  hbkgsum -> Add(h2d_sfw_TETAGAM, 1.);
  hbkgsum -> SetName("hbkgsum");

  // mc sum
  TH2D *hmcsum = (TH2D *) hbkgsum -> Clone();
  hmcsum -> Add(h2d_sfw_TETAGAM, 1.);
  hmcsum -> Add(h2d_sfw_TISR3PI_SIG, 1.);
  hmcsum -> SetName("hmcsum");
  
  // save
  TFile *f_out = new TFile(output_path + "sfw2d_output.root", "recreate");

  h2d_sfw_TDATA -> Write();
  h2d_sfw_TISR3PI_SIG -> Write();
  h2d_sfw_TETAGAM -> Write();
    
  hmcrest -> Write();
  hbkgsum_noeta -> Write();
  hmcsum -> Write();
  hbkgsum -> Write();
  
  f_out->Close();
  infile->Close();

  cout << "2D histos are created!" << endl;
  
  return 0;
  
}
