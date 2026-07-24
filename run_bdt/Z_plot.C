// ============================================================================
// Z_plot.C
//
// Reads mass bias results from a file (one line per iteration):
//   iteration  mpi0_data  mpi0_data_err  mpi0_mc  mpi0_mc_err
// Computes Z = |m_data - m_mc| / sqrt( sigma_data^2 + sigma_mc^2 )
// Plots Z vs iteration with exponential fit.
// Stop criterion: Z < 2.
// Usage:
//   root -l -q 'Z_plot.C("massbias_iterations.txt")'
// ============================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLine.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TString.h"

void Z_plot(const TString filename = "../pull_scan/massbias_iterations.txt") {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  // ------------------------------------------------------------
  // Read the data file
  // ------------------------------------------------------------
  std::ifstream infile(filename);
  if (!infile.is_open()) {
    std::cerr << "ERROR: Cannot open file " << filename << std::endl;
    return;
  }

  std::vector<double> iter_vec, Z_vec, bias_vec, bias_err_vec;
  double iter, m_data, m_data_err, m_mc, m_mc_err;

  while (infile >> iter >> m_data >> m_data_err >> m_mc >> m_mc_err) {
    double bias = m_data - m_mc;
    double bias_err = std::sqrt(m_data_err * m_data_err + m_mc_err * m_mc_err);
    double Z = std::abs(bias) / bias_err;

    iter_vec.push_back(iter);
    bias_vec.push_back(bias);
    bias_err_vec.push_back(bias_err);
    Z_vec.push_back(Z);

    std::cout << "Iter " << iter 
              << ": bias = " << bias << " +/- " << bias_err 
              << ", Z = " << Z << std::endl;
  }
  infile.close();

  int nb_points = iter_vec.size();
  if (nb_points == 0) {
    std::cerr << "ERROR: No data read from file." << std::endl;
    return;
  }

  // ------------------------------------------------------------
  // Prepare TGraphs
  // ------------------------------------------------------------
  TGraph *gZ = new TGraph(nb_points, iter_vec.data(), Z_vec.data());
  gZ->SetName("gZ");

  // Also create a graph with error bars (optional)
  TGraphErrors *gZ_err = new TGraphErrors(nb_points, iter_vec.data(), Z_vec.data(),
                                          0, bias_err_vec.data());
  gZ_err->SetName("gZ_err");

  // ------------------------------------------------------------
  // Fit an exponential decay: Z(iter) = Z0 * exp(-k * iter)
  // ------------------------------------------------------------
  double max_iter = *std::max_element(iter_vec.begin(), iter_vec.end());
  double min_iter = *std::min_element(iter_vec.begin(), iter_vec.end());

  TF1 *fit_exp = new TF1("fit_exp", "[0] * exp(-[1] * x)", min_iter - 0.5, max_iter + 0.5);
  fit_exp->SetParameter(0, Z_vec[0]);
  fit_exp->SetParameter(1, 0.5);
  gZ->Fit(fit_exp, "R");
  fit_exp->SetLineColor(kRed);
  fit_exp->SetLineStyle(2);
  fit_exp->SetLineWidth(2);

  // ------------------------------------------------------------
  // Plotting
  // ------------------------------------------------------------
  TCanvas *c1 = new TCanvas("c1", "Z-value convergence", 1000, 700);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);

  // Style the graph
  gZ->SetTitle(";Iteration;Z (significance)");
  gZ->SetMarkerStyle(20);
  gZ->SetMarkerSize(1.8);
  gZ->SetMarkerColor(kBlue);
  gZ->SetLineColor(kBlue);
  gZ->SetLineWidth(3);

  // X axis
  gZ->GetXaxis()->SetLimits(min_iter, max_iter + 0.5);
  gZ->GetXaxis()->SetNdivisions(15);
  gZ->GetXaxis()->CenterTitle(true);
  gZ->GetXaxis()->SetLabelSize(0.045);
  gZ->GetXaxis()->SetTitleSize(0.05);
  
  // Y axis
  double zmax = *std::max_element(Z_vec.begin(), Z_vec.end());
  gZ->GetYaxis()->CenterTitle(true);
  gZ->GetYaxis()->SetLabelSize(0.045);
  gZ->GetYaxis()->SetTitleSize(0.05);
  gZ->GetYaxis()->SetNdivisions(505);
  // Use linear scale; if you want log, uncomment and adjust range
  // gZ->GetYaxis()->SetRangeUser(0.1, zmax * 1.2);
  gZ->GetYaxis()->SetRangeUser(0, zmax * 1.2);

  // Draw
  gZ->Draw("AP");
  fit_exp->Draw("SAME");

  // Stop criterion line (Z = 2)
  TLine *lineZ2 = new TLine(min_iter, 2.0, max_iter + 0.5, 2.0);
  lineZ2->SetLineColor(kGray + 2);
  lineZ2->SetLineStyle(2);
  lineZ2->SetLineWidth(2);
  lineZ2->Draw("SAME");

  // Legend
  TLegend *leg = new TLegend(0.55, 0.6, 0.85, 0.88);
  leg->SetTextFont(132);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(gZ, "Z-value (data)", "PL");
  leg->AddEntry(fit_exp, "Exp. fit", "L");
  leg->AddEntry(lineZ2, "Stop criterion (Z = 2)", "L");
  leg->Draw();

  // ------------------------------------------------------------
  // Save the plot
  // ------------------------------------------------------------
  gSystem->mkdir("../massBias_tuning", kTRUE);
  c1->SaveAs("../massBias_tuning/Z_plot.pdf");

  // ------------------------------------------------------------
  // Print fit results
  // ------------------------------------------------------------
  std::cout << "\n=== Fit Results ===" << std::endl;
  std::cout << "Z0 (amplitude) = " << fit_exp->GetParameter(0)
            << " +/- " << fit_exp->GetParError(0) << std::endl;
  std::cout << "Decay rate k   = " << fit_exp->GetParameter(1)
            << " +/- " << fit_exp->GetParError(1) << std::endl;
  std::cout << "Convergence reached when Z < 2." << std::endl;
}
