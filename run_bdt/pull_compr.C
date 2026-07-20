// ============================================================================
// pull_compr.C
//
// Compare bias (mean of pull) vs energy between signal and data.
// Reads:
//   ../pull_scan/pull_scan_TISR3PI_SIG_PEAK.root
//   ../pull_scan/pull_scan_TDATA.root
// Plots:
//   - Top pad: Bias vs energy for signal and data (real points only).
//   - Bottom pad: Difference (data - signal) for all 28 signal bins,
//                 with zeros for missing data bins, shown in gray.
// Saves:
//   - ../pull_scan/bias_comparison.pdf
//   - ../pull_scan/bias_comparison.root
// ============================================================================

#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TGaxis.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TLine.h>
#include <TPad.h>
#include <iostream>
#include <vector>

using namespace std;

void pull_compr() {

  cout << "\n========================================" << endl;
  cout << "  COMPARE BIAS: DATA vs SIGNAL" << endl;
  cout << "========================================\n" << endl;

  gROOT->GetListOfCanvases()->Delete();
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);

  // ----------------------------------------------------------------------
  // Input files
  // ----------------------------------------------------------------------
  TFile *f_sig = TFile::Open("../pull_scan/pull_scan_TISR3PI_SIG_PEAK.root");
  if (!f_sig || f_sig->IsZombie()) {
    cerr << "ERROR: Cannot open signal file." << endl;
    return;
  }
  TFile *f_data = TFile::Open("../pull_scan/pull_scan_TDATA.root");
  if (!f_data || f_data->IsZombie()) {
    cerr << "ERROR: Cannot open data file." << endl;
    f_sig->Close();
    return;
  }

  TGraphErrors *g_bias_sig = (TGraphErrors*)f_sig->Get("g_bias_vs_E");
  TGraphErrors *g_bias_data = (TGraphErrors*)f_data->Get("g_bias_vs_E");
  if (!g_bias_sig || !g_bias_data) {
    cerr << "ERROR: Cannot find g_bias_vs_E in one of the files." << endl;
    f_sig->Close(); f_data->Close();
    return;
  }

  int n_sig = g_bias_sig->GetN();
  int n_data = g_bias_data->GetN();
  cout << "✓ Loaded signal bias with " << n_sig << " points" << endl;
  cout << "✓ Loaded data bias with " << n_data << " points" << endl;

  // ----------------------------------------------------------------------
  // Build the padded difference graphs
  // ----------------------------------------------------------------------
  vector<double> E_common, diff_common, err_common;
  vector<double> E_zero, diff_zero, err_zero;
  vector<double> E_pad, diff_pad, err_pad;

  for (int j = 0; j < n_sig; ++j) {
    double x_sig, y_sig;
    g_bias_sig->GetPoint(j, x_sig, y_sig);
    double ey_sig = g_bias_sig->GetErrorY(j);

    bool found = false;
    for (int i = 0; i < n_data; ++i) {
      double x_data, y_data;
      g_bias_data->GetPoint(i, x_data, y_data);
      double ey_data = g_bias_data->GetErrorY(i);
      if (x_sig == x_data) {
        double diff_val = y_data - y_sig;
        double err = TMath::Sqrt(ey_data*ey_data + ey_sig*ey_sig);
        E_common.push_back(x_sig);
        diff_common.push_back(diff_val);
        err_common.push_back(err);
        E_pad.push_back(x_sig);
        diff_pad.push_back(diff_val);
        err_pad.push_back(err);
        found = true;
        break;
      }
    }
    if (!found) {
      E_zero.push_back(x_sig);
      diff_zero.push_back(0.0);
      err_zero.push_back(0.0);
      E_pad.push_back(x_sig);
      diff_pad.push_back(0.0);
      err_pad.push_back(0.0);
    }
  }

  // Remove energy points outside [25, 270] MeV
  vector<double> E_common_clean, diff_common_clean, err_common_clean;
  for (size_t i = 0; i < E_common.size(); ++i) {
    if (E_common[i] >= 25.0 && E_common[i] <= 270.0) {
      E_common_clean.push_back(E_common[i]);
      diff_common_clean.push_back(diff_common[i]);
      err_common_clean.push_back(err_common[i]);
    }
  }

  // Replace the original vectors
  //E_common = E_common_clean;
  //diff_common = diff_common_clean;
  //err_common = err_common_clean;

  cout << "✓ Kept " << E_common.size() << " points with energy between 25 and 270 MeV." << endl;
  
  TGraphErrors *g_bias_diff_common = new TGraphErrors(E_common.size(), &E_common[0],
                                                      &diff_common[0], 0, &err_common[0]);
  g_bias_diff_common->SetName("g_bias_diff_common");

  TGraphErrors *g_bias_diff_common_clean = new TGraphErrors(E_common_clean.size(), &E_common_clean[0],
                                                      &diff_common_clean[0], 0, &err_common_clean[0]);
  g_bias_diff_common_clean->SetName("g_bias_diff_common_clean");
  
  TGraphErrors *g_bias_diff_zero = new TGraphErrors(E_zero.size(), &E_zero[0],
                                                    &diff_zero[0], 0, &err_zero[0]);
  g_bias_diff_zero->SetName("g_bias_diff_zero");

  TGraphErrors *g_bias_diff_padded = new TGraphErrors(E_pad.size(), &E_pad[0],
                                                      &diff_pad[0], 0, &err_pad[0]);
  g_bias_diff_padded->SetName("g_bias_diff_padded");

  // ----------------------------------------------------------------------
  // Define the common X-axis range
  // ----------------------------------------------------------------------
  const double xMin = 0.0;
  const double xMax = 350.0;

  // ----------------------------------------------------------------------
  // Draw plots: one canvas with two pads
  // ----------------------------------------------------------------------
  TCanvas *c1 = new TCanvas("c1", "Bias Comparison", 1400, 900);
  c1->SetBottomMargin(0.12);
  c1->SetLeftMargin(0.12);

  // ---- Top pad ----
  TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
  pad1->SetBottomMargin(0.01);
  pad1->SetLeftMargin(0.12);   // extra space to accommodate varying label widths
  //pad1->SetRightMargin(0.05);
  pad1->Draw();
  pad1->cd();

  // Set uniform axis styles
  g_bias_sig->GetYaxis()->SetTitleSize(0.08);
  g_bias_sig->GetYaxis()->SetTitleOffset(.5);
  g_bias_sig->GetYaxis()->SetLabelSize(0.06);
  g_bias_sig->GetYaxis()->SetNdivisions(505);
 
  g_bias_sig->SetMarkerStyle(20);
  g_bias_sig->SetMarkerSize(1.2);
  g_bias_sig->SetMarkerColor(kBlue);
  g_bias_sig->SetLineColor(kBlue);
  g_bias_sig->GetXaxis()->SetTitle("");
  g_bias_sig->GetYaxis()->SetTitle("<pull> (bias)");
  g_bias_sig->GetYaxis()->CenterTitle();
  g_bias_sig->Draw("AP");
  // Set X range on the histogram
  g_bias_sig->GetHistogram()->GetXaxis()->SetRangeUser(xMin, xMax);

  g_bias_data->SetMarkerStyle(21);
  g_bias_data->SetMarkerSize(1.2);
  g_bias_data->SetMarkerColor(kRed);
  g_bias_data->SetLineColor(kRed);
  g_bias_data->Draw("P same");

  // Hide x-axis on top pad by modifying the histogram's x-axis
  g_bias_sig->GetHistogram()->GetXaxis()->SetLabelOffset(999);
  g_bias_sig->GetHistogram()->GetXaxis()->SetTickLength(0);
  g_bias_sig->GetHistogram()->GetXaxis()->SetTitleOffset(999);
  gPad->Modified(); gPad->Update();

  TLegend *leg_top = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg_top->SetFillStyle(0);
  leg_top->SetBorderSize(0);
  leg_top->AddEntry(g_bias_sig, "Signal", "lp");
  leg_top->AddEntry(g_bias_data, "Data", "lp");
  leg_top->Draw();

  TLine *line0_top = new TLine(xMin, 0, xMax, 0);
  line0_top->SetLineStyle(2);
  line0_top->SetLineColor(kGray+2);
  line0_top->Draw();

  // ---- Bottom pad ----
  c1->cd();
  TPad *pad2 = new TPad("pad2", "", 0, 0, 1, 0.3);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.3);
  pad2->SetLeftMargin(0.12);   // same as top
  //pad2->SetRightMargin(0.05);
  pad2->Draw();
  pad2->cd();

  g_bias_diff_padded->GetYaxis()->SetTitleSize(0.1);
  g_bias_diff_padded->GetYaxis()->SetTitleOffset(.4);
  g_bias_diff_padded->GetYaxis()->SetLabelSize(0.1);

  g_bias_diff_padded->GetXaxis()->SetTitleSize(0.15);
  g_bias_diff_padded->GetXaxis()->SetTitleOffset(.8);
  g_bias_diff_padded->GetXaxis()->SetLabelSize(0.1);
  
  g_bias_diff_padded->SetMarkerStyle(20);
  g_bias_diff_padded->SetMarkerSize(1.2);
  g_bias_diff_padded->SetMarkerColor(kBlack);
  g_bias_diff_padded->SetLineColor(kBlack);
  g_bias_diff_padded->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
  g_bias_diff_padded->GetYaxis()->SetTitle("Data - Signal bias");
  g_bias_diff_padded->GetXaxis()->CenterTitle();
  g_bias_diff_padded->GetYaxis()->CenterTitle();
  g_bias_diff_padded->Draw("AP");
  g_bias_diff_padded->GetHistogram()->GetXaxis()->SetRangeUser(xMin, xMax);
  g_bias_diff_padded->SetMarkerColor(kGray);
  g_bias_diff_padded->SetLineColor(kGray);
  g_bias_diff_padded->GetYaxis()->SetNdivisions(505);
  
  
  g_bias_diff_zero->SetMarkerStyle(24);
  g_bias_diff_zero->SetMarkerSize(1.2);
  g_bias_diff_zero->SetMarkerColor(kBlue);
  g_bias_diff_zero->SetLineColor(kBlue);
  g_bias_diff_zero->Draw("P same");

  g_bias_diff_common->SetMarkerStyle(24);
  g_bias_diff_common->SetMarkerSize(1.2);
  g_bias_diff_common->SetMarkerColor(kBlack);
  g_bias_diff_common->SetLineColor(kBlack);
  g_bias_diff_common->Draw("P same");

  TLegend *leg_bot = new TLegend(0.7, 0.7, 0.9, 0.95);
  leg_bot->SetFillStyle(0);
  leg_bot->SetBorderSize(0);
  leg_bot->AddEntry(g_bias_diff_common, "Data - Signal (real)", "lp");
  leg_bot->AddEntry(g_bias_diff_zero, "No data (set to 0)", "lp");
  leg_bot->Draw();

  TLine *line0_bot = new TLine(xMin, 0, xMax, 0);
  line0_bot->SetLineStyle(2);
  line0_bot->SetLineColor(kGray+2);
  line0_bot->Draw();

  //
  TCanvas *c2 = new TCanvas("c2", "(Data - Signal) Bias", 1400, 900);
  c2->SetBottomMargin(0.12);
  c2->SetLeftMargin(0.15);   // extra space to accommodate varying label widths

  g_bias_diff_common_clean->GetYaxis()->SetTitleSize(0.05);
  g_bias_diff_common_clean->GetYaxis()->SetTitleOffset(1.4);
  g_bias_diff_common_clean->GetYaxis()->SetLabelSize(0.05);
  g_bias_diff_common_clean->GetHistogram()->GetYaxis()->SetRangeUser(-0.3, 0.3);
  
  g_bias_diff_common_clean->GetXaxis()->SetTitleSize(0.05);
  g_bias_diff_common_clean->GetXaxis()->SetTitleOffset(1.2);
  g_bias_diff_common_clean->GetXaxis()->SetLabelSize(0.05);
  
  g_bias_diff_common_clean->SetMarkerStyle(20);
  g_bias_diff_common_clean->SetMarkerSize(1.2);
  g_bias_diff_common_clean->SetMarkerColor(kBlack);
  g_bias_diff_common_clean->SetLineColor(kBlack);
  g_bias_diff_common_clean->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
  g_bias_diff_common_clean->GetYaxis()->SetTitle("(Data - Signal) [bias]");
  g_bias_diff_common_clean->GetXaxis()->CenterTitle();
  g_bias_diff_common_clean->GetYaxis()->CenterTitle();
  g_bias_diff_common_clean->Draw("AP");
  //g_bias_diff_common_clean->GetHistogram()->GetXaxis()->SetRangeUser(xMin, xMax);
  g_bias_diff_common_clean->SetMarkerColor(kGray+5);
  g_bias_diff_common_clean->SetLineColor(kGray+5);
  g_bias_diff_common_clean->GetYaxis()->SetNdivisions(505);
  
  g_bias_diff_common_clean->Draw("AP");

  double xMinClean = 1e9, xMaxClean = -1e9;
  for (int i = 0; i < g_bias_diff_common_clean->GetN(); ++i) {
    double x, y;
    g_bias_diff_common_clean->GetPoint(i, x, y);
    if (x < xMinClean) xMinClean = x;
  if (x > xMaxClean) xMaxClean = x;
  }
  cout << "✓ Clean difference graph x-axis limits: " << xMinClean << " - " << xMaxClean << " MeV" << endl;
  
  TLine *line0_fit = new TLine(xMinClean, 0, xMaxClean, 0);
  line0_fit->SetLineStyle(2);
  line0_fit->SetLineColor(kGray+2);
  line0_fit->Draw();

  // Linear fit
  TF1 *linFit = new TF1("linFit", "pol0", 25, 270);
  g_bias_diff_common_clean->Fit(linFit, "RQ");
  linFit->SetLineColor(kRed);
  linFit->SetLineWidth(2);
  linFit->Draw("same");

  // Display fit results
  TPaveText *pt = new TPaveText(0.2, 0.8, 0.7, 0.85, "NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  //pt->AddText(Form("Slope = %.3f #pm %.3f", linFit->GetParameter(1), linFit->GetParError(1)));
  pt->AddText(Form("Constant = %.3f #pm %.3f    #chi^{2}/NDF = %.2f", linFit->GetParameter(0), linFit->GetParError(0), linFit->GetChisquare()/linFit->GetNDF()));
  pt->Draw();

  TLegend *leg_fit = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg_fit->SetFillStyle(0);
  leg_fit->SetBorderSize(0);
  leg_fit->AddEntry(g_bias_diff_common_clean, "Pull Bias", "lep");
  leg_fit->AddEntry(linFit, "Fit", "l");
  leg_fit->Draw();

  // ----------------------------------------------------------------------
  // Save outputs
  // ----------------------------------------------------------------------
  TFile *fout = new TFile("../pull_scan/bias_comparison.root", "RECREATE");
  fout->cd();
  g_bias_sig->Write("g_bias_sig");
  g_bias_data->Write("g_bias_data");
  g_bias_diff_common->Write("g_bias_diff_common");
  g_bias_diff_zero->Write("g_bias_diff_zero");
  g_bias_diff_padded->Write("g_bias_diff_padded");
  fout->Close();

  c1->Print("../pull_scan/bias_comparison.pdf");
  c2->Print("../pull_scan/bias_average.pdf");

  f_sig->Close();
  f_data->Close();

  cout << "\n✅ Output written to ../pull_scan/bias_comparison.pdf and .root" << endl;
}
