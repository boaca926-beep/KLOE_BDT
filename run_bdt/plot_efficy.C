#include "../header_plot/plot.h"
#include "../header_bdt/efficy_plot.h"
#include "../header_bdt/graph.h"

TCanvas *plotting_efficy(const TString cv_title, const TString cv_nm, TGraphErrors *gf_bdt, TGraphErrors *gf_kloe, TGraphErrors *gf_ratio){

  double x1 = 0., y1 = 0.;
  double x2 = 0., y2 = 0.;
  
  gf_bdt -> GetPoint(0, x1, y1);
  gf_bdt -> GetPoint(1, x2, y2);

  double Delta_m3pi = x2 - x1;
  cout << Delta_m3pi << endl;

  const double mass_min = 760., mass_max = 800.;
  
  
  TCanvas *cv = new TCanvas(cv_title, cv_nm, 1200, 800);
  
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

  gf_kloe -> GetYaxis() -> SetNdivisions(505);
  gf_kloe -> GetYaxis() -> SetTitleFont(43);
  //gf_kloe -> GetYaxis() -> SetRangeUser(0., gf_kloe -> GetMaximum() * 1.2);
  gf_kloe -> GetYaxis() -> SetRangeUser(0., .05);
  //gf_kloe -> GetYaxis() -> SetTitle(TString::Format("Efficiency (#tilde{#varepsilon})/[%0.2f MeV/c^{2}]", Delta_m3pi));
  gf_kloe -> GetYaxis() -> SetTitle("Efficiency (#varepsilon)");
  gf_kloe -> GetYaxis() -> SetTitleSize(33);
  gf_kloe -> GetYaxis() -> SetTitleOffset(1.5);
  gf_kloe -> GetYaxis() -> SetLabelSize(0.05);
  gf_kloe -> GetYaxis() -> CenterTitle();

  gf_kloe -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  gf_kloe -> GetXaxis() -> SetTitleOffset(1.1);
  gf_kloe -> GetXaxis() -> SetTitleSize(0.04);
  gf_kloe -> GetXaxis() -> SetLabelSize(0.04);
  gf_kloe -> GetXaxis() -> SetLabelOffset(4);
  gf_kloe -> GetXaxis() -> SetRangeUser(mass_min, mass_max);
  gf_kloe -> GetXaxis() -> CenterTitle();
  
  //gf_kloe -> SetLineColor(1);
  //gf_kloe -> SetLineWidth(2);
  
  //gf_bdt -> SetLineColor(46);
  //gf_bdt -> SetLineWidth(2);
  
  gf_kloe -> Draw("APZ");
  gf_bdt -> Draw("PZ");

  //
  TLegend * legd_cv_p1 = new TLegend(0.4, 0.7, 0.85, 0.85);
  
  SetLegend(legd_cv_p1);
  legd_cv_p1 -> SetTextSize(0.07);
  legd_cv_p1 -> SetNColumns(2);
  
  legd_cv_p1 -> AddEntry(gf_kloe, "#varepsilon_{KLOE}", "lep");
  legd_cv_p1 -> AddEntry(gf_bdt, "#varepsilon_{bdt}", "lep");
  
  legd_cv_p1 -> Draw("Same");

  p2 -> cd();

  //gf_ratio -> SetLineColor(0);
  gf_ratio -> GetYaxis() -> SetNdivisions(505);
  gf_ratio -> GetYaxis() -> SetTitleSize(33);
  gf_ratio -> GetYaxis() -> SetTitleFont(43);
  gf_ratio -> GetYaxis() -> SetTitleOffset(1.4);
  ////gf_ratio -> GetYaxis() -> SetLabelFont(43); // Absolute front size in pixel (precision 3)
  gf_ratio -> GetYaxis() -> SetLabelSize(0.1);
  gf_ratio -> GetYaxis() -> SetRangeUser(0.5, 1.5);

  //gf_ratio -> GetYaxis() -> SetRangeUser(0., gf_ratio -> GetMaximum() * 1.2);
  gf_ratio -> GetYaxis() -> CenterTitle();
  gf_ratio -> GetYaxis() -> SetTitle("#varepsilon_{KLOE}/#varepsilon_{bdt}");

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
  
  legd_cv_p2 -> AddEntry(gf_ratio, "#varepsilon_{bdt}/#varepsilon_{KLOE}", "lep");
  
  //legd_cv_p2 -> Draw("Same");
  
  return cv;



}

double ratioErr(double a, double bdtma_a, double b, double bdtma_b) {// ratio = a / b

  if (b == 0) {

    cout << "Division by zero!" << endl;
    return 0;

  }

  double c = a / b;
  double rela_err = TMath::Sqrt(TMath::Power(bdtma_a / a, 2) + TMath::Power(bdtma_b / b, 2));
  double bdtma_c = c * rela_err;

  //cout << TMath::Power(bdtma_a / a, 2) + TMath::Power(bdtma_b / b, 2) << ", a = " << a << ", b = " << b << endl;
   
  return bdtma_c;
  
}

double get_ratio(double a, double b) {// ratio = a / b

  if (b == 0) {

    cout << "Division by zero!" << endl;
    return 0;

  }

  double c = a / b;

  return c;
  
}

TGraphErrors* get_graph(TH1D* h) {
    if (!h) return nullptr;
    int n = h->GetNbinsX();
    if (n <= 0) return nullptr;

    // Use dynamic arrays (std::vector) – no VLAs, no stack overflow
    std::vector<double> x(n), xerr(n), y(n), yerr(n);
    for (int i = 1; i <= n; ++i) {
        x[i-1]    = h->GetBinCenter(i);
        xerr[i-1] = h->GetBinWidth(i) * 0.5;   // half‑bin width as x error
        y[i-1]    = h->GetBinContent(i);
        yerr[i-1] = h->GetBinError(i);
    }

    TGraphErrors* g = new TGraphErrors(n, x.data(), y.data(), xerr.data(), yerr.data());
    // Optional styling (same as your get_graph_syst)
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.8);
    g->SetLineColor(kBlack);
    g->SetMarkerColor(kBlack);
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    return g;
}

int plot_efficy() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFitFormat("6.4g");

  // Plot efficiency comparsion
  
  TFile *efficy_bdt_file = TFile::Open(output_folder + "efficy_" + "bdt.root");
  TFile *efficy_kloe_file = TFile::Open(output_folder + "efficy_" + "kloe.root");
  
  if (!efficy_bdt_file || efficy_bdt_file->IsZombie()) {
    std::cerr << "ERROR: cannot open " << efficy_bdt_file << std::endl;
    return 0;
  }
  
  if (!efficy_kloe_file || efficy_kloe_file->IsZombie()) {
    std::cerr << "ERROR: cannot open " << efficy_kloe_file << std::endl;
    return 0;
  }

  TH1D *hefficy_bdt = (TH1D*)efficy_bdt_file -> Get("hefficy");
  TH1D *hefficy_kloe = (TH1D*)efficy_kloe_file -> Get("hefficy");

  //hefficy_bdt->Draw();
  //hefficy_kloe->Draw("same");

  TGraphErrors* efficy_bdt_gf = get_graph(hefficy_bdt);
  efficy_bdt_gf -> SetLineColor(kBlue);
  efficy_bdt_gf -> SetMarkerSize(.8);
  efficy_bdt_gf -> SetMarkerStyle(20);
  
  TGraphErrors* efficy_kloe_gf = get_graph(hefficy_kloe);
  efficy_kloe_gf -> SetLineColor(kBlack);
  efficy_kloe_gf -> SetMarkerSize(.8);
  efficy_kloe_gf -> SetMarkerStyle(22);
  
  //efficy_bdt_gf->Draw();
  //efficy_kloe_gf->Draw("same");

  // calcualte ratio
  int nPoints = efficy_bdt_gf -> GetN();
  double *x_efficy_bdt = efficy_bdt_gf -> GetX();
  double *x_efficy_bdt_err = efficy_bdt_gf -> GetEX();

  double *y_efficy_bdt = efficy_bdt_gf -> GetY();
  double *y_efficy_bdt_err = efficy_bdt_gf -> GetEY();
  
  double *y_efficy_kloe = efficy_kloe_gf -> GetY();
  double *y_efficy_kloe_err = efficy_kloe_gf -> GetEY();
  
  double RATIO[nPoints], RATIO_ERR[nPoints];

  cout << "nPoints = " << nPoints << endl;
  
  for (int i = 0; i < nPoints; i ++) {

    if (y_efficy_bdt[i] == 0. || y_efficy_kloe[i] == 0.) {
        RATIO[i] = 0.;
	RATIO_ERR[i] = 0.;
    }
    else {
      RATIO[i] = get_ratio(y_efficy_kloe[i], y_efficy_bdt[i]); //y_efficy_kloe[i] / y_efficy_bdt[i];
      RATIO_ERR[i] = ratioErr(y_efficy_kloe[i], y_efficy_kloe_err[i], y_efficy_bdt[i], y_efficy_bdt_err[i]);
  
    }
    
    //cout << "point " << i << ", mass =" << x_efficy_bdt[i] << ", efficy_bdt = " << y_efficy_bdt[i] << "+/-" << y_efficy_bdt_err[i] << ", efficy_kloe = " << y_efficy_kloe[i] << "+/-" << y_efficy_kloe_err[i] << ", efficy ratio (kloe/bdt) = " << RATIO[i] << "+/-" << RATIO_ERR[i] << endl;
    
  }

  TGraphErrors *gf_ratio = get_graph_syst(x_efficy_bdt, RATIO, x_efficy_bdt_err, RATIO_ERR, nPoints);
  gf_ratio -> SetName("gf_ratio");

  // plot
  TCanvas *cv_efficy = plotting_efficy("cv_efficy", "Efficiency Comparsion", efficy_bdt_gf, efficy_kloe_gf, gf_ratio);

  // save
  cv_efficy->cd();
  cv_efficy->SaveAs(output_folder + "cv_efficy_ratio.pdf");

  return 0;
  
}
