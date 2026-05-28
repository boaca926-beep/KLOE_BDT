TGraphErrors *get_graph_syst(double xlist[], double ylist[], double xlist_err[], double ylist_err[], int length) {

  TGraphErrors *gf = new TGraphErrors(length, xlist, ylist, xlist_err, ylist_err);

  gf -> SetMarkerStyle(20);
  gf -> SetMarkerSize(0.8);
  gf -> SetLineColor(kBlack);
  gf -> SetLineWidth(2);
  gf -> SetMarkerColor(kBlack);
  gf -> GetXaxis() -> CenterTitle();
  //gf -> GetXaxis() -> SetRangeUser(-0.5., 0.5.);
  gf -> GetXaxis() -> SetTitleOffset(1.);
  gf -> GetXaxis() -> SetTitleSize(0.06);
  gf -> GetXaxis() -> SetLabelSize(0.05);  

  gf -> GetYaxis() -> SetTitleOffset(.8);
  //gf -> GetYaxis() -> SetRangeUser();
  gf -> GetYaxis() -> SetTitleSize(0.08);
  gf -> GetYaxis() -> SetLabelSize(0.05);  
  gf -> GetYaxis() -> CenterTitle();
  
  //gf -> Draw("AP");
  
  return gf;

}

//
void SetGFAttr(TGraph * gf, const TString x_title = "", const TString y_title = "Lumi_{ISR} [nb^{-1}]", const TString gf_nm = "") {

  //gf -> SetTitle("systmatics");
  gf -> SetName(gf_nm);
  gf -> SetMarkerStyle(20);
  gf -> SetMarkerSize(0.6);
  //gf -> SetLineColor(0);
  gf -> SetMarkerColor(kBlack);
  //gf -> GetXaxis() -> SetRangeUser(760., 805.);
  gf -> GetXaxis() -> SetTitle(x_title);
  gf -> GetXaxis() -> SetLabelSize(0.04);  
  gf -> GetXaxis() -> SetTitleOffset(.8);
  gf -> GetXaxis() -> SetTitleSize(.045);
  gf -> GetXaxis() -> CenterTitle();

  gf -> GetYaxis() -> SetTitle(y_title);
  gf -> GetYaxis() -> CenterTitle();
  gf -> GetYaxis() -> SetLabelSize(0.04);  
  gf -> GetYaxis() -> SetTitleOffset(.7);
  gf -> GetYaxis() -> SetTitleSize(.045);
  //gf -> GetYaxis() -> SetRangeUser(y_min, y_max);
  
}

TGraphErrors *get_graph_norm(double xlist[], double ylist[], double xlist_err[], double ylist_err[], int length) {

  TGraphErrors *gf = new TGraphErrors(length, xlist, ylist, xlist_err, ylist_err);
  
  gf -> SetName("gf");
  gf -> SetMarkerStyle(20);
  gf -> SetMarkerSize(1.1);
  gf -> SetLineColor(kBlue);
  gf -> SetMarkerColor(kBlue);
  gf -> SetLineWidth(2);

  return gf;
  
}

TGraphErrors *get_graph_band(double xlist[], double ylist[], double xlist_err[], double ylist_err[], int length) {

  TGraphErrors *gf = new TGraphErrors(length, xlist, ylist, xlist_err, ylist_err);
  
  gf -> SetName("gf");
  gf -> SetMarkerStyle(20);
  gf -> SetMarkerSize(1.1);
  gf -> SetLineColor(kBlue);
  gf -> SetLineWidth(2);
  gf -> SetFillStyle(3002);
  gf -> SetFillColor(1);
  
  return gf;

}
