const TString input_folder = "../../efficy_filfo";
  
TCanvas *plotting_efficy(const TString cv_title, const TString cv_nm, TGraphErrors *gf_sig, TGraphErrors *gf_ufo, TGraphErrors *gf_ratio){

  double x1 = 0., y1 = 0.;
  double x2 = 0., y2 = 0.;
  
  gf_sig -> GetPoint(0, x1, y1);
  gf_sig -> GetPoint(1, x2, y2);

  double Delta_m3pi = x2 - x1;
  cout << Delta_m3pi << endl;

  const double mass_min = 760., mass_max = 800.;
  
  
  TCanvas *cv = new TCanvas(cv_title, cv_nm, 1000, 800);
  
  TPad *p2 = new TPad("p2", "p2", 0., 0., 1., 0.35);
  p2 -> Draw();
  p2 -> SetBottomMargin(0.25);
  p2 -> SetTopMargin(0.04);
  p2 -> SetLeftMargin(0.1);
  //p2 -> SetGrid();

  TPad *p1 = new TPad("p1", "p1", 0., 0.35, 1., 1.);
  p1 -> Draw();
  p1 -> SetBottomMargin(0.02);//0.007
  p1 -> SetLeftMargin(0.1);
  p1 -> cd();

  gf_ufo -> GetYaxis() -> SetNdivisions(512);
  gf_ufo -> GetYaxis() -> SetTitleFont(43);
  gf_ufo -> GetYaxis() -> SetRangeUser(0., gf_ufo -> GetMaximum() * 1.2);
  //gf_ufo -> GetYaxis() -> SetRangeUser(0., 700.);
  //gf_ufo -> GetYaxis() -> SetTitle(TString::Format("Efficiency (#tilde{#varepsilon})/[%0.2f MeV/c^{2}]", Delta_m3pi));
  gf_ufo -> GetYaxis() -> SetTitle("Efficiency (#tilde{#varepsilon})");
  gf_ufo -> GetYaxis() -> SetTitleSize(33);
  gf_ufo -> GetYaxis() -> SetTitleOffset(1.5);
  gf_ufo -> GetYaxis() -> SetLabelSize(0.05);
  gf_ufo -> GetYaxis() -> CenterTitle();

  gf_ufo -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  gf_ufo -> GetXaxis() -> SetTitleOffset(1.1);
  gf_ufo -> GetXaxis() -> SetTitleSize(0.04);
  gf_ufo -> GetXaxis() -> SetLabelSize(0.04);
  gf_ufo -> GetXaxis() -> SetLabelOffset(4);
  gf_ufo -> GetXaxis() -> SetRangeUser(mass_min, mass_max);
  gf_ufo -> GetXaxis() -> SetRangeUser(720., 840.);
  gf_ufo -> GetXaxis() -> CenterTitle();
  
  //gf_ufo -> SetLineColor(1);
  //gf_ufo -> SetLineWidth(2);
  
  //gf_sig -> SetLineColor(46);
  //gf_sig -> SetLineWidth(2);
  
  gf_ufo -> Draw("APZ");
  gf_sig -> Draw("PZ");

  //
  TLegend * legd_cv_p1 = new TLegend(0.4, 0.7, 0.85, 0.85);
  
  SetLegend(legd_cv_p1);
  legd_cv_p1 -> SetTextSize(0.07);
  legd_cv_p1 -> SetNColumns(2);
  
  legd_cv_p1 -> AddEntry(gf_ufo, "#tilde{#varepsilon}_{UFO}", "lep");
  legd_cv_p1 -> AddEntry(gf_sig, "#tilde{#varepsilon}_{sig}", "lep");
  
  legd_cv_p1 -> Draw("Same");

  p2 -> cd();

  //gf_ratio -> SetLineColor(0);
  gf_ratio -> GetYaxis() -> SetNdivisions(505);
  gf_ratio -> GetYaxis() -> SetTitleSize(33);
  gf_ratio -> GetYaxis() -> SetTitleFont(43);
  gf_ratio -> GetYaxis() -> SetTitleOffset(1.4);
  ////gf_ratio -> GetYaxis() -> SetLabelFont(43); // Absolute front size in pixel (precision 3)
  gf_ratio -> GetYaxis() -> SetLabelSize(0.1);
  gf_ratio -> GetYaxis() -> SetRangeUser(0.5, 1.2);
  
  //gf_ratio -> GetYaxis() -> SetRangeUser(0., gf_ratio -> GetMaximum() * 1.2);
  gf_ratio -> GetYaxis() -> CenterTitle();
  gf_ratio -> GetYaxis() -> SetTitle("Ratio #tilde{#varepsilon}_{UFO}/#tilde{#varepsilon}_{sig}");

  gf_ratio -> GetXaxis() -> SetRangeUser(mass_min, mass_max);
  gf_ratio -> GetXaxis() -> SetTitleOffset(.8);
  gf_ratio -> GetXaxis() -> SetLabelSize(0.1);
  gf_ratio -> GetXaxis() -> SetTitleSize(0.13);
  gf_ratio -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  gf_ratio -> GetXaxis() -> CenterTitle();

  gf_ratio -> SetLineColor(kBlack);
  gf_ratio -> SetLineWidth(2);
  
  gf_ratio -> Draw("APZ");
  p2 -> SetGrid();


  TLegend *legd_cv_p2 = new TLegend(0.2, 0.65, 0.4, 0.9);
  SetLegend(legd_cv_p2);
  legd_cv_p2 -> SetTextSize(0.1);
  legd_cv_p2 -> SetNColumns(1);
  
  legd_cv_p2 -> AddEntry(gf_ratio, "#tilde{#varepsilon}_{sig}/#tilde{#varepsilon}_{UFO}", "lep");
  
  //legd_cv_p2 -> Draw("Same");
  
  return cv;



}
  
