// omega_fit_bdt.C – adapted for merged non-resonant background
// All masses in MeV
// Now uses Breit-Wigner fit on pure ω signal.

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TFitResult.h>
#include <iostream>
#include <cmath>
#include <vector>

#include "../header_bdt/sfw2d.txt"
#include "../header_bdt/omega_fit_bdt.h"

Double_t breitwigner(Double_t *x, Double_t *par) {
    return par[0] / ((x[0] - par[1]) * (x[0] - par[1]) + par[2] * par[2]);
}

void omega_fit_bdt_merged() {
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  // ------------------------------------------------------------------
  // 1. Open tree file and sfw2d file
  // ------------------------------------------------------------------
  TFile *ftree = TFile::Open(treeFile);
  if (!ftree || ftree->IsZombie()) {
    std::cerr << "ERROR: cannot open " << treeFile << std::endl;
    return;
  }

  TFile *fsfw2d = TFile::Open(sfw2dFile);
  if (!fsfw2d || fsfw2d->IsZombie()) {
    std::cerr << "ERROR: cannot open " << sfw2dFile << std::endl;
    return;
  }

  TTree *tdata = (TTree*) ftree->Get("TDATA");
  if (!tdata) { std::cerr << "ERROR: TDATA not found." << std::endl; return; }
  
  double low, high;
  int nbins = NBINS;
  
  low  = MASS_MIN;          // e.g., 600.0 MeV
  high = MASS_MAX;          // e.g., 1000.0 MeV

  std::cout << "Mass range [" << low << ", " << high << "] MeV/c²\n";

  // ------------------------------------------------------------------
  // 2. Data histogram
  // ------------------------------------------------------------------
  TH1D *h_data = new TH1D("h_data", "", nbins, low, high);
  h_data->Sumw2();
  h_data->SetDirectory(0);
  double mtest;
  tdata->SetBranchAddress("Br_m3pi_bdt", &mtest);
  for (Long64_t i = 0; i < tdata->GetEntries(); ++i) {
    tdata->GetEntry(i);
    h_data->Fill(mtest);
  }

  // ------------------------------------------------------------------
  // 3. Load scaling factors (no separate non-resonant)
  // ------------------------------------------------------------------
  double eeg_sfw, isr3pi_sfw, omegapi_sfw, ksl_sfw, mcrest_sfw;
      
  if (fsfw2d && !fsfw2d->IsZombie()) {
    TTree *fitTree = (TTree*)fsfw2d->Get("TRESULT");
    if (fitTree) {
      fitTree->SetBranchAddress("Br_eeg_sfw", &eeg_sfw);
      fitTree->SetBranchAddress("Br_isr3pi_sfw", &isr3pi_sfw);
      fitTree->SetBranchAddress("Br_omegapi_sfw", &omegapi_sfw);
      fitTree->SetBranchAddress("Br_ksl_sfw", &ksl_sfw);
      fitTree->SetBranchAddress("Br_mcrest_sfw", &mcrest_sfw);  // includes non-resonant
      fitTree->GetEntry(0);
      
      cout << "\n=== Scaling factors from SFW2D (merged) ===" << endl;
      cout << "eeg_sfw = " << eeg_sfw << endl;
      cout << "isr3pi_sfw = " << isr3pi_sfw << endl;
      cout << "omegapi_sfw = " << omegapi_sfw << endl;
      cout << "ksl_sfw = " << ksl_sfw << endl;
      cout << "mcrest_sfw = " << mcrest_sfw << " (includes non-resonant)" << endl;
    } else {
      cout << "WARNING: TRESULT tree not found in sfw2d.root" << endl;
    }
    fsfw2d->Close();
  } else {
    cout << "WARNING: sfw2d.root not found, using unscaled MC" << endl;
  }
  
  // ------------------------------------------------------------------
  // 4. Load scaled MC components
  // ------------------------------------------------------------------
  auto makeScaledHist = [&](const char* tname, double scale, int color = 1, int style = 1) -> TH1D* {
    TTree *t = (TTree*) ftree->Get(tname);
    if (!t) return nullptr;
    if (!t->GetBranch("Br_m3pi_bdt")) return nullptr;
    TH1D *h = new TH1D(Form("h_%s", tname), "", nbins, low, high);
    h->Sumw2();
    h->SetDirectory(0);
    double val;
    t->SetBranchAddress("Br_m3pi_bdt", &val);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      h->Fill(val);
    }
    h->Scale(scale);
    h->SetLineColor(color);
    h->SetLineStyle(style);
    h->SetLineWidth(2);
    return h;
  };

  TH1D *h_eeg       = makeScaledHist("TEEG", eeg_sfw, 6, 7);
  TH1D *h_isr3pi    = makeScaledHist("TISR3PI_SIG_PEAK", isr3pi_sfw, 4, 2);
  TH1D *h_omegapi   = makeScaledHist("TOMEGAPI", omegapi_sfw, 7, 5);
  TH1D *h_ksl       = makeScaledHist("TKSL", ksl_sfw, 28, 4);

  // MC Rest (includes non-resonant via MCREST)
  TH1D *h_mcrest = new TH1D("h_mcrest", "", nbins, low, high);
  h_mcrest->Sumw2();
  h_mcrest->SetDirectory(0);
  // Build from TBKGREST + TKPM + TRHOPI + TETAGAM + TISR3PI_SIG_NON_RESON
  // But since non-resonant is already merged into MCREST scaling factor,
  // we just load the combined MCREST histogram from the merged version.
  // However, we don't have a single tree for MCREST; we build it.
  // Use the same components as in compr_bdt_merged.C
  for (const char* name : {"TKPM", "TRHOPI", "TBKGREST", "TETAGAM"}) {
    TTree *t = (TTree*) ftree->Get(name);
    if (!t) continue;
    if (!t->GetBranch("Br_m3pi_bdt")) continue;
    double val;
    t->SetBranchAddress("Br_m3pi_bdt", &val);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      h_mcrest->Fill(val);
    }
  }
  // Add non-resonant explicitly (since it's not in the above list)
  TTree *t_nonres = (TTree*) ftree->Get("TISR3PI_SIG_NON_RESON");
  if (t_nonres && t_nonres->GetBranch("Br_m3pi_bdt")) {
    double val;
    t_nonres->SetBranchAddress("Br_m3pi_bdt", &val);
    for (Long64_t i = 0; i < t_nonres->GetEntries(); ++i) {
      t_nonres->GetEntry(i);
      h_mcrest->Fill(val);
    }
  }
  h_mcrest->Scale(mcrest_sfw);
  h_mcrest->SetLineColor(37);
  h_mcrest->SetLineStyle(6);
  h_mcrest->SetLineWidth(2);

  if (!h_isr3pi) { std::cerr << "ERROR: No ISR3pi histogram." << std::endl; return; }

  // ------------------------------------------------------------------
  // 5. Build total MC for display (all components)
  // ------------------------------------------------------------------
  std::vector<TH1D*> comps = {h_eeg, h_omegapi, h_ksl, h_mcrest, h_isr3pi};
  TH1D *h_mc_total = (TH1D*) h_mcrest->Clone("h_mc_total");
  h_mc_total->Reset();
  h_mc_total->Sumw2();
  for (auto h : comps) if (h) h_mc_total->Add(h);
  h_mc_total->SetLineColor(kRed);
  h_mc_total->SetLineStyle(1);
  h_mc_total->SetLineWidth(2);

  // ------------------------------------------------------------------
  // 6. Subtract all backgrounds to get pure ω signal (data)
  // ------------------------------------------------------------------
  TH1D *h_data_isr = (TH1D*) h_data->Clone("h_data_isr");
  h_data_isr->SetDirectory(0);
  
  for (int bin = 1; bin <= h_data_isr->GetNbinsX(); ++bin) {
    double data_val = h_data->GetBinContent(bin);
    double data_err = h_data->GetBinError(bin);
    double bkg_total = 0;
    double bkg_err_sq = 0;
    
    if (h_eeg) {
      bkg_total += h_eeg->GetBinContent(bin);
      bkg_err_sq += pow(h_eeg->GetBinError(bin), 2);
    }
    if (h_omegapi) {
      bkg_total += h_omegapi->GetBinContent(bin);
      bkg_err_sq += pow(h_omegapi->GetBinError(bin), 2);
    }
    if (h_ksl) {
      bkg_total += h_ksl->GetBinContent(bin);
      bkg_err_sq += pow(h_ksl->GetBinError(bin), 2);
    }
    if (h_mcrest) {
      bkg_total += h_mcrest->GetBinContent(bin);
      bkg_err_sq += pow(h_mcrest->GetBinError(bin), 2);
    }
    
    double new_val = data_val - bkg_total;
    double new_err = sqrt(data_err*data_err + bkg_err_sq);
    
    h_data_isr->SetBinContent(bin, new_val);
    h_data_isr->SetBinError(bin, new_err);
    
    if (h_data_isr->GetBinContent(bin) < 0) {
      h_data_isr->SetBinContent(bin, 0);
    }
  }

  // ------------------------------------------------------------------
  // 7. Breit-Wigner fit to pure ω peak (data and MC)
  // ------------------------------------------------------------------
  double fit_low = 760;   // MeV
  double fit_high = 810;  // MeV

  // MC peak
  TH1D *h_mc_peak = makeScaledHist("TISR3PI_SIG_PEAK", isr3pi_sfw);
  h_mc_peak->SetDirectory(0);

  // Fit data
  TF1 *f_data = new TF1("f_data", breitwigner, fit_low, fit_high, 3);
  f_data->SetParameters(h_data_isr->GetMaximum()*4, 782, 5);
  f_data->SetParLimits(1, 760, 800);
  f_data->SetParLimits(2, 1, 20);
  TFitResultPtr r_data = h_data_isr->Fit(f_data, "RQS", "", fit_low, fit_high);
  double data_mean = f_data->GetParameter(1);
  double data_mean_err = f_data->GetParError(1);
  double data_width = f_data->GetParameter(2);
  double data_width_err = f_data->GetParError(2);
  double data_chi2 = f_data->GetChisquare() / f_data->GetNDF();

  // Fit MC peak
  TF1 *f_mc = new TF1("f_mc", breitwigner, fit_low, fit_high, 3);
  f_mc->SetParameters(h_mc_peak->GetMaximum()*4, 782, 5);
  f_mc->SetParLimits(1, 760, 800);
  f_mc->SetParLimits(2, 1, 20);
  TFitResultPtr r_mc = h_mc_peak->Fit(f_mc, "RQS", "", fit_low, fit_high);
  double mc_mean = f_mc->GetParameter(1);
  double mc_mean_err = f_mc->GetParError(1);
  double mc_width = f_mc->GetParameter(2);
  double mc_width_err = f_mc->GetParError(2);
  double mc_chi2 = f_mc->GetChisquare() / f_mc->GetNDF();

  // Compute bias and width ratio
  double mass_bias = mc_mean - data_mean;
  double mass_bias_err = sqrt(mc_mean_err*mc_mean_err + data_mean_err*data_mean_err);
  double width_ratio = mc_width / data_width;
  double width_ratio_err = width_ratio * sqrt( pow(mc_width_err/mc_width,2) + pow(data_width_err/data_width,2) );

  cout << "\n========================================" << endl;
  cout << "Omega fit results (Breit-Wigner)" << endl;
  cout << "========================================" << endl;
  cout << "MC mean    = " << mc_mean << " +/- " << mc_mean_err << " MeV" << endl;
  cout << "Data mean  = " << data_mean << " +/- " << data_mean_err << " MeV" << endl;
  cout << "Mass bias  = " << mass_bias << " +/- " << mass_bias_err << " MeV" << endl;
  cout << "MC width   = " << mc_width << " +/- " << mc_width_err << " MeV" << endl;
  cout << "Data width = " << data_width << " +/- " << data_width_err << " MeV" << endl;
  cout << "Width ratio (MC/Data) = " << width_ratio << " +/- " << width_ratio_err << endl;
  cout << "χ²/ndf (data) = " << data_chi2 << ", χ²/ndf (MC) = " << mc_chi2 << endl;

  // ---- Write results to header file ----
  std::ofstream myfile;
  TString myfile_nm = "../header_bdt/omega_fit_params.h";
  myfile.open(myfile_nm.Data());
  myfile << "const double OMEGA_MEAN_DATA = " << data_mean << ";\n";
  myfile << "const double OMEGA_MEAN_DATA_ERR = " << data_mean_err << ";\n";
  myfile << "const double OMEGA_MEAN_MC = " << mc_mean << ";\n";
  myfile << "const double OMEGA_MEAN_MC_ERR = " << mc_mean_err << ";\n";
  myfile << "const double OMEGA_MASS_BIAS = " << mass_bias << ";\n";
  myfile << "const double OMEGA_MASS_BIAS_ERR = " << mass_bias_err << ";\n";
  myfile << "const double OMEGA_WIDTH_DATA = " << data_width << ";\n";
  myfile << "const double OMEGA_WIDTH_DATA_ERR = " << data_width_err << ";\n";
  myfile << "const double OMEGA_WIDTH_MC = " << mc_width << ";\n";
  myfile << "const double OMEGA_WIDTH_MC_ERR = " << mc_width_err << ";\n";
  myfile << "const double OMEGA_WIDTH_RATIO = " << width_ratio << ";\n";
  myfile << "const double OMEGA_WIDTH_RATIO_ERR = " << width_ratio_err << ";\n";
  myfile.close();

  cout << "Saved omega fit params to " << myfile_nm << endl;

  // ------------------------------------------------------------------
  // 8. Plot results
  // ------------------------------------------------------------------
  TCanvas *c = new TCanvas("c", "Omega peak (merged background)", 1200, 700);
  c->SetBottomMargin(0.12);
  c->SetLeftMargin(0.12);

  // Normalize MC peak to data integral in fit range for display
  TH1D *h_mc_peak_norm = (TH1D*)h_mc_peak->Clone("h_mc_peak_norm");
  double mc_int = h_mc_peak_norm->Integral(h_mc_peak_norm->FindBin(fit_low), h_mc_peak_norm->FindBin(fit_high));
  double data_int = h_data_isr->Integral(h_data_isr->FindBin(fit_low), h_data_isr->FindBin(fit_high));
  if (mc_int > 0) h_mc_peak_norm->Scale(data_int / mc_int);

  h_data_isr->SetMarkerStyle(20);
  h_data_isr->SetMarkerSize(0.6);
  h_data_isr->SetLineColor(1);
  h_data_isr->GetYaxis()->SetRangeUser(0, h_data_isr->GetMaximum() * 1.6);
  h_data_isr->Draw("E1");
  h_mc_peak_norm->SetLineColor(kBlue);
  h_mc_peak_norm->SetLineWidth(2);
  h_mc_peak_norm->Draw("hist same");

  f_data->SetLineColor(kRed);
  f_data->SetLineStyle(2);
  f_data->Draw("same");

  TLegend *leg = new TLegend(0.6, 0.7, 0.85, 0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_data_isr, "Data (bkg sub)", "lep");
  leg->AddEntry(h_mc_peak_norm, "MC peak (norm)", "l");
  leg->AddEntry(f_data, "BW fit to data", "l");
  leg->Draw();

  TPaveText *pt = new TPaveText(0.15, 0.8, 0.5, 0.9, "NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.04);
  pt->AddText(Form("MC mean = %.2f MeV", mc_mean));
  pt->AddText(Form("Data mean = %.2f MeV", data_mean));
  pt->AddText(Form("Mass bias = %.2f MeV", mass_bias));
  pt->Draw();

  c->SaveAs(output_path + "omega_peak_fit_merged.pdf");
  delete c;
  delete h_mc_peak_norm;

  // ------------------------------------------------------------------
  // 9. Clean up
  // ------------------------------------------------------------------
  if (ftree) { ftree->Close(); delete ftree; }
  if (fsfw2d) { fsfw2d->Close(); delete fsfw2d; }
  for (auto h : comps) if (h) delete h;
  delete h_data;
  delete h_data_isr;
  delete h_mc_peak;
  delete h_mc_total;
  delete f_data;
  delete f_mc;

  std::cout << "\n=== Summary ===" << std::endl;
  std::cout << "Saved " << output_path << "omega_peak_fit_merged.pdf" << std::endl;
  std::cout << "Fit results written to " << myfile_nm << std::endl;
}
