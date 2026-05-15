TF1 *ReFun_norm;

const double beta_cut = 1.98;
const double c0 = 0.11;
const double c1 = 0.8;

TCanvas *plot_cv(const TString cv_title, TH2D *h2d, TPaveText *pt, const TString pt_str) {

  TCanvas *cv_tmp = new TCanvas(cv_title, " ", 700, 700);

  cv_tmp -> SetBottomMargin(0.15);//0.007
  cv_tmp -> SetLeftMargin(0.15);
  cv_tmp -> SetRightMargin(0.15);

  h2d -> SetMinimum(10);

  h2d -> GetXaxis() -> SetNdivisions(5);
  h2d -> GetXaxis() -> SetTitle("M_{2#pi} [GeV/c^{2}]");
  h2d -> GetXaxis() -> SetTitleOffset(1.2);
  h2d -> GetXaxis() -> SetTitleSize(0.06);
  h2d -> GetXaxis() -> CenterTitle();
  h2d -> GetXaxis() -> SetLabelSize(0.06);
  h2d -> GetXaxis() -> SetLabelOffset(0.01);
  h2d -> GetXaxis() -> SetRangeUser(0.2, 0.6);
  
  h2d -> GetYaxis() -> SetTitle("#beta_{#pi}");
  h2d -> GetYaxis() -> SetLabelOffset(0.01);
  h2d -> GetYaxis() -> SetTitleOffset(1.2);
  h2d -> GetYaxis() -> SetLabelSize(0.06);
  h2d -> GetYaxis() -> SetTitleSize(0.06);
  h2d -> GetYaxis() -> CenterTitle();

  h2d -> GetZaxis() -> SetLabelSize(0.06);

  gPad->SetLogz(1);
  //cv_tmp -> SetLogx();
  //h2d -> GetXaxis() -> SetMoreLogLabels();
  
  char display[50];
  
  pt = new TPaveText(0.65, 0.2, 0.8, 0.3, "NDC");

  pt -> SetTextSize(0.1);
  pt -> SetFillColor(0);
  pt -> SetTextAlign(12);

  //sprintf(display, str);
  pt -> AddText(pt_str);
  cout << pt_str << endl;
  
  h2d -> Draw("COLZ");
  pt -> Draw("Same");

  //gPad->SetLogy(1);
  
  // beta functin
  //cout << "beta_cut = " << beta_cut << ", c0 = " << c0 << ", c1 = " << c1 << endl;

  ReFun_norm = new TF1("ReFun_norm", "[0]+1/(exp((x-[2])/[1])-1)", 0.1, 1.);
  ReFun_norm -> SetParameters(beta_cut, c0, c1);
  ReFun_norm -> SetLineColor(kRed);
  ReFun_norm -> SetLineWidth(3);

  ReFun_norm -> Draw("Same");
  
  return cv_tmp;
  
}
