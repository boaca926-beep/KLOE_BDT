// ============================================================================
// resol_ratio_compr.C
//
// Compare resolution ratio (Data/Signal) vs energy.
// Reads:
//   ../pull_scan/pull_scan_TISR3PI_SIG_PEAK.root
//   ../pull_scan/pull_scan_TDATA.root
// Plots:
//   - Signal and data resolution points.
//   - Ratio (Data/Signal) with zero for missing bins.
//   - Constant fit to the ratio over [15, 350] MeV.
// Saves:
//   - ../pull_scan/resolution_ratio.pdf
//   - ../pull_scan/resolution_ratio.root
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
#include <TF1.h>
#include <iostream>
#include <vector>

using namespace std;

void resol_compr() {

  cout << "\n========================================" << endl;
  cout << "  COMPARE RESOLUTION RATIO: DATA / SIGNAL" << endl;
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

  TGraphErrors *g_res_sig = (TGraphErrors*)f_sig->Get("g_resolution_vs_E");
  TGraphErrors *g_res_data = (TGraphErrors*)f_data->Get("g_resolution_vs_E");
  if (!g_res_sig || !g_res_data) {
    cerr << "ERROR: Cannot find g_resolution_vs_E in one of the files." << endl;
    f_sig->Close(); f_data->Close();
    return;
  }

  int n_sig = g_res_sig->GetN();
  int n_data = g_res_data->GetN();
  cout << "✓ Loaded signal resolution with " << n_sig << " points" << endl;
  cout << "✓ Loaded data resolution with " << n_data << " points" << endl;

  // ----------------------------------------------------------------------
  // Build the padded ratio graphs
  // ----------------------------------------------------------------------
  vector<double> E_common, ratio_common, err_common;
  vector<double> E_zero, ratio_zero, err_zero;
  vector<double> E_pad, ratio_pad, err_pad;

  const double eps = 1e-6;

  for (int j = 0; j < n_sig; ++j) {
    double x_sig, y_sig;
    g_res_sig->GetPoint(j, x_sig, y_sig);
    double ey_sig = g_res_sig->GetErrorY(j);

    bool found = false;
    for (int i = 0; i < n_data; ++i) {
      double x_data, y_data;
      g_res_data->GetPoint(i, x_data, y_data);
      double ey_data = g_res_data->GetErrorY(i);
      if (TMath::Abs(x_sig - x_data) < eps) {
        // Compute ratio Data / Signal
        if (y_sig == 0) {
          cerr << "WARNING: Signal resolution is zero at E = " << x_sig << " – skipping ratio." << endl;
          continue;
        }
        double ratio_val = y_data / y_sig;
        // Error propagation for ratio: sigma_ratio = ratio * sqrt( (ey_data/y_data)^2 + (ey_sig/y_sig)^2 )
        // Assume y_data > 0 (resolution is positive)
        double err_ratio = ratio_val * TMath::Sqrt( TMath::Power(ey_data / y_data, 2) +
                                                    TMath::Power(ey_sig / y_sig, 2) );
        E_common.push_back(x_sig);
        ratio_common.push_back(ratio_val);
        err_common.push_back(err_ratio);
        E_pad.push_back(x_sig);
        ratio_pad.push_back(ratio_val);
        err_pad.push_back(err_ratio);
        cout << i+1 << ": y_data = " << y_data << ", y_sig = " << y_sig << ": ratio = " << ratio_val << endl;
        found = true;
        break;
      }
    }
    if (!found) {
      // No matching data point: set ratio to 0 (indicates missing)
      E_zero.push_back(x_sig);
      ratio_zero.push_back(0.0);
      err_zero.push_back(0.0);
      E_pad.push_back(x_sig);
      ratio_pad.push_back(0.0);
      err_pad.push_back(0.0);
    }
  }

  if (E_common.empty()) {
    cerr << "WARNING: No common energy points found between signal and data!" << endl;
  }

  // Remove energy points outside [15, 350] MeV for the fit
  vector<double> E_common_clean, ratio_common_clean, err_common_clean;
  for (size_t i = 0; i < E_common.size(); ++i) {
    if (E_common[i] >= 15.0 && E_common[i] <= 350.0) {
      E_common_clean.push_back(E_common[i]);
      ratio_common_clean.push_back(ratio_common[i]);
      err_common_clean.push_back(err_common[i]);
    }
  }

  TGraphErrors *g_res_ratio_common = new TGraphErrors(E_common.size(), &E_common[0],
                                                      &ratio_common[0], 0, &err_common[0]);
  g_res_ratio_common->SetName("g_res_ratio_common");

  TGraphErrors *g_res_ratio_common_clean = new TGraphErrors(E_common_clean.size(), &E_common_clean[0],
                                                            &ratio_common_clean[0], 0, &err_common_clean[0]);
  g_res_ratio_common_clean->SetName("g_res_ratio_common_clean");

  TGraphErrors *g_res_ratio_zero = new TGraphErrors(E_zero.size(), &E_zero[0],
                                                    &ratio_zero[0], 0, &err_zero[0]);
  g_res_ratio_zero->SetName("g_res_ratio_zero");

  TGraphErrors *g_res_ratio_padded = new TGraphErrors(E_pad.size(), &E_pad[0],
                                                      &ratio_pad[0], 0, &err_pad[0]);
  g_res_ratio_padded->SetName("g_res_ratio_padded");

  double xMinClean = 1e9, xMaxClean = -1e9;
  for (int i = 0; i < g_res_ratio_common_clean->GetN(); ++i) {
    double x, y;
    g_res_ratio_common_clean->GetPoint(i, x, y);
    if (x < xMinClean) xMinClean = x;
    if (x > xMaxClean) xMaxClean = x;
  }

  // Fit a constant to the cleaned ratio
  TF1 *linFit = new TF1("linFit", "pol0", xMinClean, xMaxClean);
  g_res_ratio_common_clean->Fit(linFit, "RQ");

  cout << "✓ Kept " << E_common.size() << " points with energy between "
       << xMinClean << " and " << xMaxClean << " MeV" << endl << endl;

  // ----------------------------------------------------------------------
  // Draw plot: one canvas showing signal, data, and ratio
  // ----------------------------------------------------------------------
  const double xMin = 0.0;
  const double xMax = 350.0;

  TCanvas *c1 = new TCanvas("c1", "Resolution Ratio Comparison", 1400, 900);
  c1->SetBottomMargin(0.12);
  c1->SetLeftMargin(0.12);

  // Set axis styles for the signal graph (use it to define the pad)
  g_res_sig->GetXaxis()->SetTitleSize(0.06);
  g_res_sig->GetXaxis()->SetTitleOffset(.9);
  g_res_sig->GetXaxis()->SetLabelSize(0.05);
  g_res_sig->GetXaxis()->SetNdivisions(505);

  g_res_sig->GetYaxis()->SetTitleSize(0.06);
  g_res_sig->GetYaxis()->SetTitleOffset(0.8);
  g_res_sig->GetYaxis()->SetLabelSize(0.05);
  g_res_sig->GetYaxis()->SetNdivisions(505);
  // Adjust Y range to show both resolution and ratio (ratio might be around 1)
  g_res_sig->GetHistogram()->GetYaxis()->SetRangeUser(0.8, 1.4); // adjust if needed

  g_res_sig->SetMarkerStyle(20);
  g_res_sig->SetMarkerSize(1.2);
  g_res_sig->SetMarkerColor(kBlue);
  g_res_sig->SetLineColor(kBlue);
  g_res_sig->GetXaxis()->SetTitle("E_{#gamma} (MeV)");
  g_res_sig->GetYaxis()->SetTitle("#sigma (pull)");
  g_res_sig->GetYaxis()->CenterTitle();
  g_res_sig->Draw("AP");

  g_res_data->SetMarkerStyle(21);
  g_res_data->SetMarkerSize(1.2);
  g_res_data->SetMarkerColor(kRed);
  g_res_data->SetLineColor(kRed);
  g_res_data->Draw("P same");

  // Draw the ratio (data/signal) for common points (black)
  g_res_ratio_common->SetMarkerStyle(24);
  g_res_ratio_common->SetMarkerSize(1.2);
  g_res_ratio_common->SetMarkerColor(kBlack);
  g_res_ratio_common->SetLineColor(kBlack);
  g_res_ratio_common->Draw("P same");

  // Zero-filled missing bins (blue) – optional
  g_res_ratio_zero->SetMarkerStyle(24);
  g_res_ratio_zero->SetMarkerSize(1.2);
  g_res_ratio_zero->SetMarkerColor(kBlue);
  g_res_ratio_zero->SetLineColor(kBlue);
  g_res_ratio_zero->Draw("P same");   // uncomment if you want to show missing bins

  // Line at ratio = 1 (if data=signal)
  TLine *line1 = new TLine(xMin, 1.0, xMax, 1.0);
  line1->SetLineStyle(2);
  line1->SetLineColor(kGray+2);
  line1->Draw();

  // Display fit results
  TPaveText *pt = new TPaveText(0.15, 0.8, 0.7, 0.85, "NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->AddText(Form("<Ratio> = %.3f #pm %.3f    #chi^{2}/NDF = %.2f",
                   linFit->GetParameter(0), linFit->GetParError(0),
                   linFit->GetChisquare()/linFit->GetNDF()));
  pt->Draw();

  double ratio_diff_avg = linFit->GetParameter(0);
  double ratio_diff_err = linFit->GetParError(0);
  double Z_value = (1 - ratio_diff_avg) / ratio_diff_err;
  cout << "Z_value = " << Z_value << endl;

  TPaveText *pt1 = new TPaveText(0.15, 0.7, 0.7, 0.75, "NDC");
  pt1->SetFillColor(0);
  pt1->SetBorderSize(0);
  pt1->AddText(Form("Significance = %.2f#sigma", Z_value));
  pt1->Draw();

  linFit->SetLineColor(kRed);
  linFit->SetLineWidth(2);
  linFit->Draw("same");

  // Legend
  TLegend *leg = new TLegend(0.7, 0.65, 0.9, 0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(g_res_sig, "Signal", "lp");
  leg->AddEntry(g_res_data, "Data", "lp");
  leg->AddEntry(g_res_ratio_common, "Ratio (Data/Signal)", "lp");
  leg->AddEntry(linFit, "Fit", "l");
  leg->Draw();

  cout << "SUMMARY" << "\n"
       << "Average resolution ratio (Data/Signal) = " << ratio_diff_avg << " +/- " << ratio_diff_err << endl;

  // ----------------------------------------------------------------------
  // Save outputs
  // ----------------------------------------------------------------------
  TFile *fout = new TFile("../pull_scan/resolution_ratio.root", "RECREATE");
  fout->cd();
  g_res_sig->Write("g_res_sig");
  g_res_data->Write("g_res_data");
  g_res_ratio_common->Write("g_res_ratio_common");
  g_res_ratio_zero->Write("g_res_ratio_zero");
  g_res_ratio_padded->Write("g_res_ratio_padded");
  fout->Close();

  c1->Print("../pull_scan/resolution_ratio.pdf");

  f_sig->Close();
  f_data->Close();

  cout << "\n✅ Output written to ../pull_scan/resolution_ratio.pdf and .root" << endl;
}
