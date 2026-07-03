// Determine Z-value convergency.
// Z-value (significance) of pi0 photon energy scale correction
// Z-value results in Thesis_Syst
#include "../header/graph.h"
void Z_plot() {

  gStyle->SetOptFit(0);
  
  const int nb_points = 3;
  double bias[nb_points]     = {-0.580, -0.321, -0.169};
  double bias_err[nb_points] = {0.045, 0.045, 0.045};
  double Z[nb_points]        = {0., 0., 0.};
  double iter[nb_points]     = {0., 1., 2.};
 
  for (int i = 0; i < nb_points; i++) {
    Z[i] = TMath::Abs(bias[i]) / bias_err[i];
    std::cout << "Iteration " << iter[i] 
              << ", bias = " << bias[i] << " +/- " << bias_err[i] 
              << ", Z = " << Z[i] << std::endl;
  }

  TCanvas *c1 = new TCanvas("c1", "Z-value Convergence", 800, 600);
  gPad->SetRightMargin(0.12);
  
  TGraph *gf_Z = new TGraph(nb_points, iter, Z);
  TGraph *gf_plot = (TGraph*)gf_Z->Clone("gf_plot");
  
  // ---------- Fit ----------
  TF1 *fit_exp = new TF1("fit_exp", "[0] * exp(-[1] * x)", -0.5, 3.0);
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
  gf_plot->GetXaxis()->SetLimits(0., 2.5);   // Keep it tight to data
  gf_plot->GetXaxis()->SetNdivisions(3);          // 0,1,2
  gf_plot->GetXaxis()->CenterTitle(true);
  gf_plot->GetXaxis()->SetLabelSize(0.045);
  gf_plot->GetXaxis()->SetTitleSize(0.05);
  
  // ---------- Y-axis ----------
  gf_plot->GetYaxis()->CenterTitle(true);
  gf_plot->GetYaxis()->SetLabelSize(0.045);
  gf_plot->GetYaxis()->SetTitleSize(0.05);
  gf_plot->GetYaxis()->SetRangeUser(0., 15.);
  
  // ---------- Threshold Line ----------
  // Now spans full x-axis range (-0.5 to 2.5)
  TLine *lineZ2 = new TLine(0., 2, 2.5, 2);
  lineZ2->SetLineColor(kGray + 2);
  lineZ2->SetLineStyle(2);
  lineZ2->SetLineWidth(2);
  
  // ---------- Draw ----------
  gf_plot->Draw("AP");       // Line connects points now
  fit_exp->Draw("SAME");
  lineZ2->Draw("SAME");

  // ---------- Legend ----------
  TLegend *leg = new TLegend(0.15, 0.72, 0.45, 0.88);
  leg->SetHeader("Convergence", "C");
  leg->SetTextSize(0.04);
  leg->AddEntry(gf_plot, "Z-value (data)", "PL");
  leg->AddEntry(fit_exp, "Exp. fit", "L");
  leg->Draw();

  c1->SaveAs("Z_plot.pdf");
  c1->SaveAs("Z_plot.png");

  std::cout << "\n=== Fit Results ===" << std::endl;
  std::cout << "Amplitude (Z0) = " << fit_exp->GetParameter(0) 
            << " +/- " << fit_exp->GetParError(0) << std::endl;
  std::cout << "Decay rate (k)  = " << fit_exp->GetParameter(1) 
            << " +/- " << fit_exp->GetParError(1) << std::endl;
}
