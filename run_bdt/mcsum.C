//#include "../hist.h"
//#include "sfw2d.h"

int mcsum() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //gROOT->SetBatch(kFALSE); //kFALSE
    
  cout << "Plotting mcsum ... " << endl;

  //
  TFile * infile_mcsum = new TFile("./plots/hist.root");

  TIter next_tree(infile_mcsum -> GetListOfKeys());

  TString objnm_tree, classnm_tree;
  
  int i = 0;
  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    cout << " hist" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  TH2D* h2d = (TH2D*)infile_mcsum -> Get(hist_type);
  TH1D *h2d_projx = h2d -> ProjectionX();
  TH1D *h2d_projy = h2d -> ProjectionY();
    
  cout << hist_type << endl;

  // plot

  double binwidth_x = getbinwidth(h2d_projx);
  double binwidth_y = getbinwidth(h2d_projy);

  h2d -> GetXaxis() -> SetTitle("M_{2#pi} " + TString::Format("Events/%0.2f", binwidth_x) + " [MeV/c^{2}]"); //SetTitle(x_label + " " + x_unit);
  //h2d -> GetXaxis() -> SetRangeUser();
  h2d -> GetXaxis() -> SetTitleOffset(1.2);
  h2d -> GetXaxis() -> CenterTitle();
  h2d -> GetYaxis() -> SetTitle("E_{#gamma_{3}} " + TString::Format("Events/%0.2f", binwidth_y) + " [MeV]"); //SetTitle(x_label + " " + x_unit);
  h2d -> GetYaxis() -> SetTitleOffset(1.4);
  h2d -> GetYaxis() -> SetLabelSize(0.03);
  h2d -> GetYaxis() -> CenterTitle();
  
  TPaveText *pt1 = new TPaveText(0.65, 0.7, 0.8, 0.85, "NDC");

  pt1 -> SetTextSize(0.1);
  pt1 -> SetFillColor(0);
  pt1 -> SetTextAlign(12);
  
  //pt1 -> AddText("(b) MC without #eta#gamma");
  pt1 -> AddText("(b)");
    
  TCanvas * cv2d =  new TCanvas("cv2d_" + hist_type, cv_nm, 0, 0, 700, 700);

  cv2d -> SetBottomMargin(0.15);//0.007
  cv2d -> SetLeftMargin(0.15);
  cv2d -> SetRightMargin(0.15);

  h2d -> GetXaxis() -> SetTitleOffset(1.0);
  h2d -> GetXaxis() -> SetLabelOffset(0.01);
  h2d -> GetXaxis() -> CenterTitle();
  h2d -> GetXaxis() -> SetTitleSize(0.06);
  h2d -> GetXaxis() -> SetLabelSize(0.04);
  
  h2d -> GetYaxis() -> SetTitleOffset(1.1);
  h2d -> GetYaxis() -> SetLabelOffset(0.01);
  h2d -> GetYaxis() -> SetTitleSize(0.06);
  h2d -> GetYaxis() -> SetLabelSize(0.04);
  h2d -> GetYaxis() -> CenterTitle();

  h2d -> Draw("COLZ");
  //h2d -> Draw("TEXT0COLZ");

  pt1 -> Draw("Same");
  h2d -> GetZaxis() -> SetLabelSize(0.045);
  
  gPad -> SetLogz();
  
  // save
  cv2d -> SaveAs("./plots/sfw_2d_" + cv_nm + ".pdf");
  
  
  return 0;

}
