#include "resol.h"
#include "../hist.h"
#include "../fitfun.h"

///home/bo/Desktop/analysis/resol/plot_resol.sh calling plot_resol.C
//input tree: ../crx3pi/output_chains_norm/tree_cut0.root, chi2 cut = 40 gives a better fit to extract resolutions
    
int plot_resol() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptTitle(0);
  gStyle->SetStatBorderSize(0);
  //gStyle->SetOptFit(11);
  //gStyle->SetStatX(0.89);
  //gStyle->SetStatY(0.89);
  gStyle->SetFitFormat("6.4g");
  
  gStyle->SetOptStat(0);

  const TString resol_nm = "h1d_resol_" + kin_var;
  
  //const TString output = "./resol_output"; 
  cout << "Plot resolutions ... " << resol_nm << endl;
    
  TFile *intree = new TFile("histos/hist.root");
  
  TIter next_tree(intree -> GetListOfKeys());
  
  TString objnm_tree, classnm_tree;
  
  int i = 0;

  TKey *key;
  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    cout << " tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  // get histos
  TH1D* hresol = (TH1D*)intree -> Get(resol_nm);
  const double bin_size = hresol -> GetNbinsX(); 
  double norm_resol = hresol -> Integral(1, bin_size);

  cout << norm_resol << endl;
  hresol -> Scale(1. / norm_resol);
  
  //hresol -> Draw();

  // fit resolution
  TString hist_nm = hresol -> GetName();
  double xmin = hresol -> GetXaxis() -> GetXmin();
  double xmax = hresol -> GetXaxis() -> GetXmax();
  
  //cout << hist_nm << ", (" << xmin << ", " << xmax << ")" << endl;

  const int npar = 6, npar_sub = 3;
  double rms = hresol -> GetRMS();
  double mean = hresol -> GetMean();
  double peak = hresol -> GetMaximum();
  double bin_width = getbinwidth(hresol);
  double fitpara[npar] = {peak, mean, rms, 0.01 * peak, mean, 3. * rms};
  double fitgfpar[npar], fitgfpar1[npar_sub], fitgfpar2[npar_sub];
  TF1 *fitfun = new TF1("fitfun", fun_double, xmin, xmax, npar);
  fitfun -> SetParNames("p1","p2","p3","p4","p5","p6");
  fitfun -> SetParameters(fitpara);
  //fitfun -> SetParLimits(4,-.5,.5);
  fitfun -> SetLineWidth(2);
  fitfun -> SetNpx(5000);

  // fitting
  TFitResultPtr r = hresol -> Fit(fitfun, "LQ0", " ", fit_range * xmin, fit_range * xmax);

  // gauss1
  TF1 *gauss1_fit = new TF1("gauss1_fit", gauss1d, fit_range * xmin, fit_range * xmax, npar_sub);
  fitgfpar1[0] = fitfun -> GetParameter(0); // peak1
  fitgfpar1[1] = fitfun -> GetParameter(1); // mean1
  fitgfpar1[2] = fitfun -> GetParameter(2); // width1
  
  gauss1_fit -> SetParameters(fitgfpar1);
  gauss1_fit -> SetLineColor(kRed);
  gauss1_fit -> SetLineWidth(2);
  gauss1_fit -> SetLineStyle(1);
  gauss1_fit -> SetNpx(5000);

  // gauss2
  TF1 *gauss2_fit = new TF1("gauss2_fit", gauss1d, fit_range * xmin, fit_range * xmax, npar_sub);
  fitgfpar2[0] = fitfun -> GetParameter(3); // peak2
  fitgfpar2[1] = fitfun -> GetParameter(4); // mean2
  fitgfpar2[2] = fitfun -> GetParameter(5); // width2
  
  gauss2_fit -> SetParameters(fitgfpar2);
  gauss2_fit -> SetLineColor(kBlue);
  gauss2_fit -> SetLineStyle(3);
  gauss2_fit -> SetLineWidth(2);
  gauss2_fit -> SetNpx(5000); 
  
  // gauss_sum
  TF1 *gauss_sum = new TF1("gauss_sum", fun_double, fit_range * xmin, fit_range * xmax, npar);
  fitgfpar[0] = fitfun -> GetParameter(0); // peak1
  fitgfpar[1] = fitfun -> GetParameter(1); // mean1
  fitgfpar[2] = fitfun -> GetParameter(2); // width1
  fitgfpar[3] = fitfun -> GetParameter(3); // peak2
  fitgfpar[4] = fitfun -> GetParameter(4); // mean2
  fitgfpar[5] = fitfun -> GetParameter(5); // width2
  
  gauss_sum -> SetParameters(fitgfpar);
  gauss_sum -> SetLineColor(kBlack);
  //gauss_sum -> SetLineColor(kRed);
  gauss_sum -> SetLineStyle(1);
  gauss_sum -> SetLineWidth(2);
  gauss_sum -> SetNpx(5000);


  int binmax_y = hresol->GetMaximumBin();
  double binmax_x = hresol->GetXaxis()->GetBinCenter(binmax_y);
  
  cout << "fit_range = " << fit_range << "\n"
       << "residual binmax_y = " << binmax_y << ", residual = " << binmax_x << "\n";
  
  // plot
  //formatfill_h(hresol, 4, 3001);

  const double ymax = hresol -> GetMaximum();
  
  TCanvas *cv_resol = new TCanvas("cv_" + kin_var + "_resol", resol_nm, 0, 0, 700, 700);
  cv_resol -> SetBottomMargin(0.15);//0.007
  cv_resol -> SetLeftMargin(0.15);

  hresol -> SetMarkerStyle(21);
  hresol -> SetMarkerSize(0.7);
  //hresol -> GetYaxis() -> SetTitle(TString::Format("Events/[%0.2f", bin_width) + " " + unit + "]");
  hresol -> GetYaxis() -> SetTitle("Events");
  hresol -> GetYaxis() -> SetTitleOffset(1.3);
  hresol -> GetYaxis() -> CenterTitle();
  hresol -> GetYaxis() -> SetLabelSize(0.05);
  hresol -> GetYaxis() -> SetTitleSize(0.06);
  hresol -> GetYaxis() -> SetRangeUser(0., 0.1);
  cout << "ymax = " << ymax << endl;
  
  hresol -> GetXaxis() -> SetTitle(xtit + " " + unit);
  hresol -> GetXaxis() -> SetRangeUser(xmin, xmax);
  hresol -> GetXaxis() -> SetTitleOffset(1.);
  hresol -> GetXaxis() -> CenterTitle();
  hresol -> GetXaxis() -> SetLabelSize(0.045);
  hresol -> GetXaxis() -> SetTitleSize(0.06);

  TPaveText *pt1 = new TPaveText(0.3, 0.82, 0.8, 0.83, "NDC");
  
  pt1 -> SetTextSize(0.07);
  pt1 -> SetFillColor(0);
  pt1 -> SetTextAlign(12);

  TPaveText *pt2 = new TPaveText(0.2, 0.72, 0.5, 0.73, "NDC");
  
  pt2 -> SetTextSize(0.10);
  pt2 -> SetFillColor(0);
  pt2 -> SetTextAlign(12);

 
  if (kin_var == "Angle12") {
    pt1 -> AddText(Form("#delta#angle_{#gamma#gamma}=%0.2e^{#circ}", fitfun -> GetParameter(2)));
    pt2 -> AddText("(a)");
  }
  else if (kin_var == "betapi0") {
    //pt1 -> AddText(Form("#sigma=%0.4g#pm%0.1e", fitfun -> GetParameter(2), fitfun -> GetParError(2)));
    pt1 -> AddText(Form("#delta#beta_{#pi}=%0.1e", fitfun -> GetParameter(2)));
    pt2 -> AddText("(c)");
  }
  else if (kin_var == "deltaE") {
    pt1 -> AddText(Form("#deltaE_{diff}=%0.3g MeV", fitfun -> GetParameter(2)));
    pt2 -> AddText("(b)");
  }
  else if (kin_var == "Eisr") {
    pt1 -> AddText(Form("#deltaE_{#gamma_{3}}=%0.3g MeV", fitfun -> GetParameter(2)));
  }
  else if (kin_var == "ppIM") {
    pt1 -> AddText(Form("#deltaM_{2#pi}=%0.3g MeV/c^{2}", fitfun -> GetParameter(2)));
  }
  else if (kin_var == "IM3pi") {
    pt1 -> AddText(Form("#deltaM_{3#pi}=%0.3g MeV/c^{2}", fitfun -> GetParameter(2)));
  }
  
  
  hresol -> Draw();
  gauss1_fit -> Draw("Same");
  gauss2_fit -> Draw("Same");
  gauss_sum -> Draw("Same");
  pt1 -> Draw("Same");
  pt2 -> Draw("Same");

  //cv_resol -> Update();

  //gStyle->SetOptStat(0);

  /*
  
  TPaveStats *p1 = (TPaveStats*)hresol -> GetListOfFunctions() -> FindObject("stats");
  p1 -> SetName("mystats");
  hresol -> SetStats(0);
  p1 -> SetX1NDC(0.3);
  p1 -> SetY1NDC(0.65);
  p1 -> SetTextSize(0.05);
  //p1 -> SetOptFit(11);
  p1 -> GetLineWith("p3") -> SetTextColor(kBlue);
  p1 -> Draw("Same");
  */

  // save
  cv_resol -> SaveAs("plots/" + kin_var + "_resol.pdf");
  
  
  return 0;

}
