// Determine Z-value convergency.
// Z-value (significance) of pi0 photon energy scale correction
// Z-value results in Thesis_Syst
#include "../header/graph.h"

void Z_plot() {
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  TString tuning_type = "scaled";   // FIXED: TSTring -> TString

  int nb_points = 0;
  double iter_max = 5.0;            // FIXED: double *iter_max -> double iter_max
  double *bias = nullptr;
  double *bias_err = nullptr;
  double *Z = nullptr;
  double *iter = nullptr;
  
  // Define arrays outside if blocks (static to persist)
  static double bias_tuning[]     = {6.136, 3.88664};
  static double bias_err_tuning[] = {0.045, 0.0448303};
  static double Z_tuning[]        = {0., 0.};
  static double iter_tuning[]     = {0., 1.};
  
  static double bias_scaled[]     = {-0.580, -0.321, -0.169, -0.0898698, -0.0488806};
  static double bias_err_scaled[] = {0.045, 0.045, 0.045, 0.044744, 0.0447478};
  static double Z_scaled[]        = {0., 0., 0., 0., 0.};
  static double iter_scaled[]     = {0., 1., 2., 3., 4.};
  
  if (tuning_type == "tuning") {  // energy pull tuning + energy scaling
    nb_points = 2;
    bias = bias_tuning;
    bias_err = bias_err_tuning;
    Z = Z_tuning;
    iter = iter_tuning;
  }
  else if (tuning_type == "scaled") {  // test of energy scaling
    nb_points = 5;
    bias = bias_scaled;
    bias_err = bias_err_scaled;
    Z = Z_scaled;
    iter = iter_scaled;
  }
  
  for (int i = 0; i < nb_points; i++) {
    Z[i] = TMath::Abs(bias[i]) / bias_err[i];
    std::cout << "Iteration " << iter[i] 
              << ", bias = " << bias[i] << " +/- " << bias_err[i] 
              << ", Z = " << Z[i] << std::endl;
  }

  TCanvas *c1 = new TCanvas("c1", "Z-value convergence", 900, 600);
  gPad->SetRightMargin(0.12);
  
  TGraph *gf_Z = new TGraph(nb_points, iter, Z);
  TGraph *gf_plot = (TGraph*)gf_Z->Clone("gf_plot");
  
  // ---------- Fit ----------
  TF1 *fit_exp = new TF1("fit_exp", "[0] * exp(-[1] * x)", -0.5, iter_max);
  fit_exp->SetParameter(0, Z[0]);
  fit_exp->SetParameter(1, 0.5);
  gf_Z->Fit(fit_exp, "R");
  fit_exp->SetLineColor(kRed);
  fit_exp->SetLineStyle(2);
  fit_exp->SetLineWidth(2);
  
  // ---------- Style ----------
  gf_plot->SetTitle("Z-value convergence;Iteration;Z (significance)");
  gf_plot->SetMarkerStyle(20);
  gf_plot->SetMarkerSize(1.8);
  gf_plot->SetMarkerColor(kBlue);
  gf_plot->SetLineColor(kBlue);
  gf_plot->SetLineWidth(3);
  
  // ---------- X-axis ----------
  gf_plot->GetXaxis()->SetLimits(0., iter_max);
  gf_plot->GetXaxis()->SetNdivisions(3);
  gf_plot->GetXaxis()->CenterTitle(true);
  gf_plot->GetXaxis()->SetLabelSize(0.045);
  gf_plot->GetXaxis()->SetTitleSize(0.05);
  
  // ---------- Y-axis ----------
  gf_plot->GetYaxis()->CenterTitle(true);
  gf_plot->GetYaxis()->SetLabelSize(0.045);
  gf_plot->GetYaxis()->SetTitleSize(0.05);
  gf_plot->GetYaxis()->SetNdivisions(505);
  gf_plot->GetYaxis()->SetRangeUser(0., 20.);
  
  // ---------- Threshold Line ----------
  TLine *lineZ2 = new TLine(0., 2, iter_max, 2.);
  lineZ2->SetLineColor(kGray + 2);
  lineZ2->SetLineStyle(2);
  lineZ2->SetLineWidth(2);
  
  // ---------- Draw ----------
  gf_plot->Draw("AP");
  fit_exp->Draw("SAME");
  lineZ2->Draw("SAME");

  // ---------- Legend ----------
  TLegend *leg = new TLegend(0.15, 0.72, 0.45, 0.88);
  leg->SetHeader("Convergence", "C");
  leg->SetTextSize(0.04);
  leg->AddEntry(gf_plot, "Z-value (data)", "PL");
  leg->AddEntry(fit_exp, "Exp. fit", "L");
  leg->AddEntry(lineZ2, "Stop criterion", "L");
  leg->Draw();

  // Create directory if needed
  gSystem->mkdir("../massBias_" + tuning_type, kTRUE);
  c1->SaveAs("../massBias_" + tuning_type + "/Z_plot.pdf");

  std::cout << "\n=== Fit Results ===" << std::endl;
  std::cout << "Amplitude (Z0) = " << fit_exp->GetParameter(0) 
            << " +/- " << fit_exp->GetParError(0) << std::endl;
  std::cout << "Decay rate (k)  = " << fit_exp->GetParameter(1) 
            << " +/- " << fit_exp->GetParError(1) << std::endl;

  // Calculate width Bias
  const double width_data = 7.097;
  const double width_data_err = 0.310;

  const double width_mc_raw = 6.372;
  const double width_mc_raw_err = 0.064;
  
  const double width_mc_tuning = 6.466;
  const double width_mc_tuning_err = 0.061;

  double width_bias_raw = width_mc_raw - width_data;
  double width_bias_raw_err = TMath::Sqrt(TMath::Power(width_mc_raw_err, 2) + TMath::Power(width_data_err, 2));

  double width_bias_tuning = width_mc_tuning - width_data;
  double width_bias_tuning_err = TMath::Sqrt(TMath::Power(width_mc_tuning_err, 2) + TMath::Power(width_data_err, 2));

  cout << "width_bias_raw = " << width_bias_raw << " +/- " << width_bias_raw_err << endl;
  cout << "width_bias_tuning = " << width_bias_tuning << " +/- " << width_bias_tuning_err << endl;
}
