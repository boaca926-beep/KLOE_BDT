#include "../header_plot/plot.h"
#include "../header_method/method.h"
#include "../header_bdt/plot2d.h"

gROOT->ForceStyle();

TCanvas *plot2d_cv(const TString cv_title, const TString cv_nm, TH2D *h2d, const TString pt_str)
{

  TH1D *h2d_projx = h2d -> ProjectionX();
  TH1D *h2d_projy = h2d -> ProjectionY();
  
  double binwidth_x = getbinwidth(h2d_projx);
  double binwidth_y = getbinwidth(h2d_projy);
  
  TCanvas * cv2d =  new TCanvas("cv2d_" + hist_type, cv_nm, 0, 0, 700, 700);

  cv2d -> SetBottomMargin(0.15);//0.007
  cv2d -> SetLeftMargin(0.15);
  cv2d -> SetRightMargin(0.15);

  //h2d -> SetMinimum(10);

  h2d -> GetXaxis() -> SetNdivisions(5);
  h2d -> GetXaxis() -> SetTitle("M_{2#pi} " + TString::Format("Events/[%0.2f", binwidth_x) + " MeV/c^{2}]"); //SetTitle(x_label + " " + x_unit);
  h2d -> GetXaxis() -> SetTitleOffset(1.2);
  h2d -> GetXaxis() -> SetTitleSize(0.06);
  h2d -> GetXaxis() -> CenterTitle();
  h2d -> GetXaxis() -> SetLabelSize(0.06);
  h2d -> GetXaxis() -> SetLabelOffset(0.01);
  //h2d -> GetXaxis() -> SetRangeUser(0.2, 0.6);
  
  h2d -> GetYaxis() -> SetTitle("E_{#gamma_{3}} " + TString::Format("Events/[%0.2f", binwidth_y) + " MeV]"); //SetTitle(x_label + " " + x_unit);
  h2d -> GetYaxis() -> SetLabelOffset(0.01);
  h2d -> GetYaxis() -> SetTitleOffset(1.2);
  h2d -> GetYaxis() -> SetLabelSize(0.05);
  h2d -> GetYaxis() -> SetTitleSize(0.06);
  h2d -> GetYaxis() -> CenterTitle();

  h2d -> Draw("COLZ");
  //h2d -> Draw("TEXT0COLZ");

  gPad->SetLogz(1);
  
  char display[50];

  //sprintf(display, "Data");
  
  TPaveText *pt = new TPaveText(pt1_x0, 0.7, pt1_x1, 0.85, "NDC");

  pt -> SetTextSize(0.1);
  pt -> SetFillColor(0);
  pt -> SetTextAlign(12);
  
  //pt -> AddText("Relative Error [%]");
  pt -> AddText(pt_str);
  
  pt -> Draw("Same");
  h2d -> GetZaxis() -> SetLabelSize(0.045);
  
  gPad -> SetLogz();
  //gPad -> Modified();
  //gPad -> Update();

  return cv2d;
  
}

int plot2d_ppbeta(const TString plots_dir = "") {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  //gROOT->SetBatch(kTRUE);
  gStyle->SetPaintTextFormat("4.0f");
  //gROOT->SetBatch(kFALSE); //kFALSE
  gStyle->SetTitleSize(0.2,"x");
  
  cout << "Plotting sfw2d ... " << endl;

  // get histos
  
  TFile *infile = new TFile(infile_nm);

  // Check input file
  getObj(infile);
  
  TH2D* h2d = (TH2D*)infile -> Get(hist_type);
  //cout << hist_type << endl;

  // plot

  TCanvas *cv2d = plot2d_cv("cv2d_" + hist_type, cv_nm, h2d, cv_text);

  // save
  cout << "plots_dir = " << plots_dir << endl;
  cv2d -> SaveAs(plots_dir + "sfw2d_" + cv_nm + ".pdf");
  
  return 0;

}
