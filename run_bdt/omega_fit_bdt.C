// omega_fit_bdt.C – template fit with proper memory management
// All masses in MeV
// Adapted plotting style from correct_and_plot.C

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

// Global pointers for template histograms (detached from file)
TH1D *gSigTemplate = nullptr;
TH1D *gBkgTemplate = nullptr;

// Fit function
Double_t template_sum(Double_t *x, Double_t *par) {

  if (!gSigTemplate || !gBkgTemplate) {
    return 0.0;
  }

  int bin = gSigTemplate->FindBin(x[0]);
  Double_t sig = gSigTemplate->GetBinContent(bin);
  Double_t bkg = gBkgTemplate->GetBinContent(bin);
  return par[0] * sig + par[1] * bkg;
}

void omega_fit_bdt() {
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  // ------------------------------------------------------------------
  // 1. Open tree file
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
  // 2. Data histogram (detached from file)
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

  // Scaling factors
  double eeg_sfw, isr3pi_sfw, omegapi_sfw, nonReson_sfw, ksl_sfw, mcrest_sfw;
      
  if (fsfw2d && !fsfw2d->IsZombie()) {
    TTree *fitTree = (TTree*)fsfw2d->Get("TRESULT");
    if (fitTree) {
      fitTree->SetBranchAddress("Br_eeg_sfw", &eeg_sfw);
      fitTree->SetBranchAddress("Br_isr3pi_sfw", &isr3pi_sfw);
      fitTree->SetBranchAddress("Br_omegapi_sfw", &omegapi_sfw);
      fitTree->SetBranchAddress("Br_nonReson_sfw", &nonReson_sfw);
      fitTree->SetBranchAddress("Br_ksl_sfw", &ksl_sfw);
      fitTree->SetBranchAddress("Br_mcrest_sfw", &mcrest_sfw);
      fitTree->GetEntry(0);
      
      cout << "\n=== Scaling factors from SFW2D ===" << endl;
      cout << "eeg_sfw = " << eeg_sfw * 2. << endl;
      cout << "isr3pi_sfw = " << isr3pi_sfw << endl;
      cout << "omegapi_sfw = " << omegapi_sfw << endl;
      cout << "nonReson_sfw = " << nonReson_sfw << endl;
      cout << "ksl_sfw = " << ksl_sfw << endl;
      cout << "mcrest_sfw = " << mcrest_sfw << endl;
    } else {
      cout << "WARNING: TRESULT tree not found in sfw2d.root" << endl;
    }
    fsfw2d->Close();
  } else {
    cout << "WARNING: sfw2d.root not found, using unscaled MC" << endl;
  }
  
  // ------------------------------------------------------------------
  // 3. Load scaled MC components (detach each histogram)
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

  // Same color codes and line styles as correct_and_plot.C
  TH1D *h_eeg       = makeScaledHist("TEEG", eeg_sfw * 2., 6, 7);
  TH1D *h_isr3pi    = makeScaledHist("TISR3PI_SIG_PEAK", isr3pi_sfw, 4, 2);
  TH1D *h_nonReson  = makeScaledHist("TISR3PI_SIG_NON_RESON", nonReson_sfw, 2, 3);
  TH1D *h_omegapi   = makeScaledHist("TOMEGAPI", omegapi_sfw, 7, 5);
  TH1D *h_ksl       = makeScaledHist("TKSL", ksl_sfw, 28, 4);

  // Save unnormalized signal for later use
  TH1D *h_isr3pi_unnorm = (TH1D*)h_isr3pi->Clone("h_isr3pi_unnorm");
  h_isr3pi_unnorm->SetLineWidth(2);
  h_isr3pi_unnorm->SetLineStyle(1);
  h_isr3pi_unnorm->SetLineColor(kBlue);
  h_isr3pi_unnorm->SetDirectory(0);
  
  // MC Rest
  TH1D *h_mcrest = new TH1D("h_mcrest", "", nbins, low, high);
  h_mcrest->Sumw2();
  h_mcrest->SetDirectory(0);
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
  h_mcrest->Scale(mcrest_sfw);
  h_mcrest->SetLineColor(37);
  h_mcrest->SetLineStyle(6);
  h_mcrest->SetLineWidth(2);

  if (!h_isr3pi) { std::cerr << "ERROR: No ISR3pi histogram." << std::endl; return; }

  // Summary
  cout << "\n=== Summary ===" << endl;
  double peak_nb = h_isr3pi->Integral();
  double distorted_nb = h_nonReson->Integral();
  double signal_sum = peak_nb + distorted_nb;
  double purity = peak_nb / signal_sum;
  double h_low = h_data->GetXaxis()->GetXmin();
  double h_max = h_data->GetXaxis()->GetXmax();
  
  std::cout << "Data: " << h_data->Integral() << "\n"
	    << "SIGNAL: " << signal_sum << ", sfw = " << isr3pi_sfw << ", purity = " << purity * 100. << "\n"
	    << "\tpeak: " << peak_nb  << ", sfw = " << isr3pi_sfw << "\n"
    	    << "\tdistorted: " << distorted_nb << ", sfw = " << nonReson_sfw << "\n"
    	    << "EEG: " << h_eeg->Integral() << ", sfw = " << eeg_sfw * 2. << "\n"
	    << "OMEGAPI: " << h_omegapi->Integral() << ", sfw = " << omegapi_sfw << "\n"
	    << "KSL: " << h_ksl->Integral() << ", sfw = " << ksl_sfw << "\n"
	    << "MCREST: " << h_mcrest->Integral() << ", sfw = " << mcrest_sfw << "\n"
    	    << std::endl;

  // ------------------------------------------------------------------
  // 4. Subtract backgrounds -> h_data_isr (detached)
  // ------------------------------------------------------------------
  TH1D *h_data_isr = (TH1D*) h_data->Clone("h_data_isr");
  h_data_isr->SetDirectory(0);
  
  // Propagate errors properly
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
    // nonReson is NOT subtracted (it's used in the fit)
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
  // 5. Reload original MC for correction (before normalization)
  // ------------------------------------------------------------------
  TH1D *h_isr3pi_orig = makeScaledHist("TISR3PI_SIG_PEAK", isr3pi_sfw);
  if (!h_isr3pi_orig) {
    std::cerr << "ERROR: Cannot reload ISR3pi for correction" << std::endl;
    return;
  }
  h_isr3pi_orig->SetDirectory(0);

  // ------------------------------------------------------------------
  // 6. Compute correction weights using Data/MC ratio (OPTION 1)
  // ------------------------------------------------------------------
  double peak_low = 700;   // MeV
  double peak_high = 850;  // MeV
  
  TH1D *h_weight = (TH1D*) h_data_isr->Clone("h_weight");
  h_weight->SetDirectory(0);
  h_weight->Reset();

  double orig_integral = h_isr3pi_orig->Integral();
  std::cout << "\n=== Computing Correction Weights ===" << std::endl;
  std::cout << "Original MC integral: " << orig_integral << std::endl;

  int n_bins = h_weight->GetNbinsX();
  int n_used = 0;

  for (int bin = 1; bin <= n_bins; ++bin) {
    double x = h_weight->GetBinCenter(bin);
    double data_val = h_data_isr->GetBinContent(bin);
    double mc_val = h_isr3pi_orig->GetBinContent(bin);
    
    if (x >= peak_low && x <= peak_high && mc_val > 0 && data_val > 0) {
      double ratio = data_val / mc_val;
      if (ratio > 2.5) ratio = 2.5;
      if (ratio < 0.4) ratio = 0.4;
      h_weight->SetBinContent(bin, ratio);
      n_used++;
    } else {
      h_weight->SetBinContent(bin, 1.0);
    }
  }

  std::cout << "Used " << n_used << " bins for weight calculation" << std::endl;

  // Smooth weights
  TH1D *h_weight_smooth = (TH1D*) h_weight->Clone("h_weight_smooth");
  h_weight_smooth->SetDirectory(0);
  for (int bin = 2; bin <= n_bins - 1; ++bin) {
    double w_avg = (h_weight->GetBinContent(bin-1) +
                    h_weight->GetBinContent(bin) +
                    h_weight->GetBinContent(bin+1)) / 3.0;
    h_weight_smooth->SetBinContent(bin, w_avg);
  }
  h_weight_smooth->SetBinContent(1, h_weight->GetBinContent(1));
  h_weight_smooth->SetBinContent(n_bins, h_weight->GetBinContent(n_bins));

  h_weight_smooth->Draw();
  
  // ------------------------------------------------------------------
  // 7. Apply correction to ISR3pi MC (preserve yield)
  // ------------------------------------------------------------------
  TH1D *h_isr3pi_corrected = (TH1D*) h_isr3pi_orig->Clone("h_isr3pi_corrected");
  h_isr3pi_corrected->SetDirectory(0);

  for (int bin = 1; bin <= n_bins; ++bin) {
    double w = h_weight_smooth->GetBinContent(bin);
    double old = h_isr3pi_corrected->GetBinContent(bin);
    h_isr3pi_corrected->SetBinContent(bin, old * w);
    double err = h_isr3pi_corrected->GetBinError(bin);
    h_isr3pi_corrected->SetBinError(bin, err * w);
  }

  double new_integral = h_isr3pi_corrected->Integral();
  if (new_integral > 0 && orig_integral > 0) {
    double renorm = orig_integral / new_integral;
    h_isr3pi_corrected->Scale(renorm);
    std::cout << "Renormalisation factor: " << renorm << std::endl;
    std::cout << "Original integral: " << orig_integral << std::endl;
    std::cout << "Corrected integral: " << h_isr3pi_corrected->Integral() << std::endl;
  }

  h_isr3pi_corrected->SetLineColor(kGreen);
  h_isr3pi_corrected->SetLineWidth(2);
  h_isr3pi_corrected->SetLineStyle(1);

  // ------------------------------------------------------------------
  // 8. Normalize templates to unit area for the fit
  // ------------------------------------------------------------------
  // Use corrected signal and original non-resonant
  TH1D *h_signal_template = (TH1D*) h_isr3pi_corrected->Clone("h_signal_template");
  h_signal_template->SetDirectory(0);
  
  TH1D *h_bkg_template = (TH1D*) h_nonReson->Clone("h_bkg_template");
  h_bkg_template->SetDirectory(0);
  
  double sig_int = h_signal_template->Integral();
  double bkg_int = h_bkg_template->Integral();
  if (sig_int > 0) h_signal_template->Scale(1.0 / sig_int);
  if (bkg_int > 0) h_bkg_template->Scale(1.0 / bkg_int);
  
  // Set global pointers for the fit function
  gSigTemplate = h_signal_template;
  gBkgTemplate = h_bkg_template;

  // ------------------------------------------------------------------
  // 9. Template fit
  // ------------------------------------------------------------------
  TF1 *total_func = new TF1("total_func", template_sum, low, high, 2);
  total_func->SetParameters(1000, 1000);
  total_func->SetParNames("alpha", "beta");
  
  TFitResultPtr r = h_data_isr->Fit(total_func, "RQS", "", peak_low, peak_high);
  
  if (!r->IsValid()) {
    std::cerr << "WARNING: Fit did not converge!" << std::endl;
  }
  
  double alpha = total_func->GetParameter(0);
  double beta  = total_func->GetParameter(1);
  double chi2 = r->Chi2();
  int ndf = r->Ndf();
  double chi2_ndf = chi2 / ndf;
  
  std::cout << "Fit results: α = " << alpha << ", β = " << beta << std::endl;
  std::cout << "Fit quality: χ² = " << chi2 << ", ndf = " << ndf << ", #chi^{2}/ndf = " << chi2_ndf << std::endl;

  // ------------------------------------------------------------------
  // 9b. Calculate purity from fit results
  // ------------------------------------------------------------------
  int lowBin  = gSigTemplate->FindBin(peak_low);
  int highBin = gSigTemplate->FindBin(peak_high);
  
  double sig_frac_peak = gSigTemplate->Integral(lowBin, highBin);
  double bkg_frac_peak = gBkgTemplate->Integral(lowBin, highBin);
  
  double fitted_signal_yield = alpha * sig_frac_peak;
  double fitted_background_yield = beta * bkg_frac_peak;
  double fitted_signal_sum = fitted_signal_yield + fitted_background_yield;
  
  double updated_purity = fitted_signal_yield / (fitted_signal_yield + fitted_background_yield);
  
  std::cout << "Fitted signal = " << fitted_signal_sum << "\n"
            << "\tpeak = " << fitted_signal_yield << "\n"
            << "\tdistorted = " << fitted_background_yield << "\n"
            << "Updated purity (from fit, in peak region): " << updated_purity * 100. << "%" << std::endl;

  // ------------------------------------------------------------------
  // 10. Create signal & background histograms (scaled)
  // ------------------------------------------------------------------
  TH1D *h_signal = (TH1D*) h_signal_template->Clone("h_signal");
  TH1D *h_background = (TH1D*) h_bkg_template->Clone("h_background");
  h_signal->SetDirectory(0); 
  h_background->SetDirectory(0);
  h_signal->Scale(alpha);
  h_background->Scale(beta);
  h_signal->SetLineColor(kBlue);
  h_signal->SetLineWidth(2);
  h_background->SetLineColor(kBlue);
  h_background->SetLineStyle(3);
  h_background->SetLineWidth(2);
  
  // ------------------------------------------------------------------
  // 11. Build total MC sum for plotting
  // ------------------------------------------------------------------
  std::vector<TH1D*> comps;
  comps.push_back(h_eeg);
  comps.push_back(h_omegapi);
  comps.push_back(h_ksl);
  comps.push_back(h_mcrest);
  comps.push_back(h_signal);
  comps.push_back(h_background);

  TH1D *h_mc_total = (TH1D*) h_mcrest->Clone("h_mc_total");
  h_mc_total->Reset();
  h_mc_total->Sumw2();
  for (auto h : comps) if (h) h_mc_total->Add(h);
  h_mc_total->SetLineColor(kRed);
  h_mc_total->SetLineStyle(1);
  h_mc_total->SetLineWidth(2);

  // ================================================================
  // 12. Calculate mass bias from mean difference (FIXED)
  // ================================================================
  // Use PURE PEAK MC (no distorted component) and data after
  // subtracting the distorted component – this gives the true
  // detector energy shift for the ω resonance.
  // ------------------------------------------------------------------
  TH1D *h_mc_peak = makeScaledHist("TISR3PI_SIG_PEAK", isr3pi_sfw);
  if (!h_mc_peak) {
    std::cerr << "ERROR: Cannot load peak MC for mass bias" << std::endl;
    return;
  }
  h_mc_peak->SetDirectory(0);

  // Data with distorted signal removed (already created as h_data_bw)
  TH1D *h_data_bw = (TH1D*)h_data_isr->Clone("h_data_bw");
  h_data_bw->Add(h_background, -1.0);   // subtract distorted from data
  h_data_bw->SetDirectory(0);

  // Define peak range for mean calculation
  double mass_low = 760;
  double mass_high = 810;

  // Set axis range for mean/RMS calculation
  h_mc_peak->GetXaxis()->SetRangeUser(mass_low, mass_high);
  h_data_bw->GetXaxis()->SetRangeUser(mass_low, mass_high);
  
  double mc_mean = h_mc_peak->GetMean();
  double data_mean = h_data_bw->GetMean();
  double mc_rms = h_mc_peak->GetRMS();
  double data_rms = h_data_bw->GetRMS();
  
  // Reset axis range after calculation
  h_mc_peak->GetXaxis()->SetRange();
  h_data_bw->GetXaxis()->SetRange();

  // Mass bias (energy shift): data - MC? Usually shift = MC - Data,
  // but we want the correction to apply to MC to match data.
  // The sign: if MC mean > data mean, MC needs to be shifted down.
  double mass_bias = -(mc_mean - data_mean);  // same as data_mean - mc_mean

  cout << "\n========================================" << endl;
  cout << "Mass Bias from Mean Difference (FIXED)" << endl;
  cout << "========================================" << endl;
  cout << "MC mean    = " << mc_mean << " +/- " << mc_rms << " MeV/c^2" << endl;
  cout << "Data mean  = " << data_mean << " +/- " << data_rms << " MeV/c^2" << endl;
  cout << "MC - Data  = " << mc_mean - data_mean << " MeV/c^2" << endl;
  cout << "mass bias  = " << mass_bias << " MeV/c^2" << endl;
  cout << "========================================" << endl;

  // ---- Plot MC vs Data comparison ----
  TCanvas *c_bias = new TCanvas("c_bias", "Mass Bias from Mean Difference", 800, 600);
  c_bias->cd();

  TH1D *h_mc_norm = (TH1D*)h_mc_peak->Clone("h_mc_norm");
  TH1D *h_data_norm = (TH1D*)h_data_bw->Clone("h_data_norm");
  h_mc_norm->Scale(1.0 / h_mc_norm->Integral());
  h_data_norm->Scale(1.0 / h_data_norm->Integral());

  h_mc_norm->SetLineColor(kRed);
  h_mc_norm->SetLineWidth(2);
  h_data_norm->SetMarkerStyle(20);
  h_data_norm->SetMarkerSize(0.6);
  h_data_norm->SetLineColor(kBlack);

  h_data_norm->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
  h_data_norm->GetYaxis()->SetTitle("Normalized Events");
  h_data_norm->GetYaxis()->SetRangeUser(0, h_data_norm->GetMaximum() * 1.4);
  h_data_norm->Draw("E1");
  h_mc_norm->Draw("hist same");

  TLegend *leg_bias = new TLegend(0.6, 0.7, 0.85, 0.85);
  leg_bias->SetFillStyle(0);
  leg_bias->SetBorderSize(0);
  leg_bias->AddEntry(h_data_norm, "Data (bkg sub)", "lep");
  leg_bias->AddEntry(h_mc_norm, "MC (peak only)", "l");
  leg_bias->Draw();

  TPaveText *pt_bias = new TPaveText(0.15, 0.8, 0.5, 0.9, "NDC");
  pt_bias->SetFillColor(0);
  pt_bias->SetBorderSize(0);
  pt_bias->SetTextAlign(12);
  pt_bias->SetTextSize(0.04);
  pt_bias->AddText(Form("MC mean = %.2f MeV", mc_mean));
  pt_bias->AddText(Form("Data mean = %.2f MeV", data_mean));
  pt_bias->AddText(Form("Mass bias = %.2f MeV", mass_bias));
  pt_bias->Draw();

  c_bias->SaveAs(output_path + "mass_bias_mean_diff.pdf");
  delete c_bias;
  delete h_mc_norm;
  delete h_data_norm;
  delete h_data_bw;

  // ---- Write bias to header ----
  std::ofstream myfile;
  TString myfile_nm = "../header_bdt/bias.h";
  myfile.open(myfile_nm.Data());
  myfile << "const double energy_shift = " << mass_bias << ";\n";
  myfile.close();

  cout << "Saved mass bias to: " << myfile_nm << endl;
  // ================================================================

  // ------------------------------------------------------------------
  // 13. Pull distribution (full range)
  // ------------------------------------------------------------------
  TH1D *h_pull = new TH1D("h_pull", "", nbins, low, high);
  for (int bin = 1; bin <= nbins; ++bin) {
    double d = h_data->GetBinContent(bin);
    double m = h_mc_total->GetBinContent(bin);
    double err = std::sqrt(d + m);
    if (err > 0) h_pull->SetBinContent(bin, (d - m)/err);
    else h_pull->SetBinContent(bin, 0);
  }
  h_pull->SetMarkerStyle(20);
  h_pull->SetMarkerSize(0.6);
  h_pull->SetLineWidth(0);

  // ------------------------------------------------------------------
  // 14. Main plotting
  // ------------------------------------------------------------------
  TCanvas *c = new TCanvas("c", "3π mass projection (combined)", 1200, 700);
  c->SetBottomMargin(0.12);
  c->SetLeftMargin(0.12);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
  pad1->SetBottomMargin(0.02);
  pad1->SetLeftMargin(0.12);
  pad1->Draw();
  pad1->cd();

  double max_val = h_data->GetMaximum();
  max_val = std::max(max_val, h_mc_total->GetMaximum());
  for (auto h : comps) if (h) max_val = std::max(max_val, h->GetMaximum());
  h_data->GetYaxis()->SetRangeUser(0, max_val * 1.2);
  double bin_width = h_data->GetBinWidth(1);

  TH1D *h_signal_data = (TH1D*) h_data->Clone("h_signal_data");
  
  h_data->SetMarkerStyle(20);
  h_data->SetMarkerSize(0.6);
  h_data->Draw("E1");
  h_mc_total->Draw("hist same");
  for (auto h : comps) {
    if (!h) continue;
    TString h_nm = h->GetName(); 
    if (h_nm != "h_TEEG") {
      h->Draw("hist same");
    }
  }
  
  h_data->GetYaxis()->SetTitle(Form("Events / [%.1f MeV/c^{2}]", bin_width));
  h_data->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
  h_data->GetYaxis()->CenterTitle();
  h_data->GetXaxis()->SetTitleSize(0.06);
  h_data->GetXaxis()->SetLabelOffset(0.1);
  h_data->GetYaxis()->SetTitleSize(0.07);
  h_data->GetYaxis()->SetTitleOffset(0.7);
  h_data->GetYaxis()->SetLabelSize(0.04);
  h_data->GetYaxis()->SetNdivisions(505);

  TPaveText *pt0 = new TPaveText(0.7, 0.62, 0.85, 0.85, "NDC");
  pt0->SetFillColor(0);
  pt0->SetBorderSize(0);
  pt0->SetTextAlign(12);
  pt0->SetTextSize(0.04);
  pt0->SetTextFont(42);
  pt0->AddText(Form("#chi^{2}/ndf = %.2f", chi2_ndf));
  //pt0->AddText(Form("Mass bias = %.2f [MeV/c^{2}]", TMath::Abs(mass_bias)));
  pt0->AddText(Form("Purity = %.1f%%", updated_purity * 100.));
  pt0->Draw();

  TLegend *leg = new TLegend(0.15, 0.35, 0.6, 0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_data, "Data", "lep");
  leg->AddEntry(h_mc_total, "Total MC", "l");
  leg->AddEntry(h_signal, "#omega peak (signal)", "l");
  leg->AddEntry(h_background, "Distorted signal", "l");
  leg->AddEntry(h_omegapi, "#omega#pi^{0}", "l");
  leg->AddEntry(h_ksl, "K_{S}K_{L}", "l");
  leg->AddEntry(h_mcrest, "Others", "l");
  leg->Draw();

  // Pull pad
  c->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.3);
  pad2->SetLeftMargin(0.12);
  pad2->Draw();
  pad2->cd();
  gPad->SetGrid();

  h_pull->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
  h_pull->GetXaxis()->SetTitleSize(0.12);
  h_pull->GetXaxis()->SetTitleOffset(1.0);
  h_pull->GetXaxis()->SetLabelSize(0.1);
  h_pull->GetYaxis()->SetTitle("Pull");
  h_pull->GetYaxis()->SetTitleSize(0.2);
  h_pull->GetYaxis()->SetTitleOffset(0.2);
  h_pull->GetYaxis()->SetLabelSize(0.1);
  h_pull->GetYaxis()->SetRangeUser(-2, 2);
  h_pull->GetYaxis()->SetNdivisions(505);
  h_pull->GetXaxis()->CenterTitle();
  h_pull->GetYaxis()->CenterTitle();
  h_pull->Draw("P");

  TLine *line = new TLine(low, 0, high, 0);
  line->SetLineStyle(2);
  line->Draw();

  c->Update();
  c->Modified();
  c->SaveAs(output_path + "omega_combined_fit.pdf");

  // ------------------------------------------------------------------
  // 15. Background-subtracted ω signal (with h_signal_final drawn)
  // ------------------------------------------------------------------
  // Subtract fixed physics backgrounds and distorted signal
  h_signal_data->Add(h_eeg, -1.0);
  h_signal_data->Add(h_omegapi, -1.0);
  h_signal_data->Add(h_ksl, -1.0);
  h_signal_data->Add(h_mcrest, -1.0);
  h_signal_data->Add(h_background, -1.0);
  // Note: h_background is the fitted tail, so data points do NOT contain the tail
  
  for (int bin = 1; bin <= h_signal_data->GetNbinsX(); ++bin) {
    if (h_signal_data->GetBinContent(bin) < 0) {
      h_signal_data->SetBinContent(bin, 0);
    }
  }
  
  // ----- Create h_signal_final (peak + distorted) -----
  TH1D *h_signal_final = (TH1D*) h_signal->Clone("h_signal_final");
  h_signal_final->Add(h_background, 1.0);
  h_signal_final->SetLineColor(kRed);
  h_signal_final->SetLineWidth(2);
  h_signal_final->SetLineStyle(1);
  // -----------------------------------------------------

  TCanvas *c2 = new TCanvas("c2", "Pure ω peak (all backgrounds subtracted)", 1200, 700);
  c2->cd();
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.12);
 
  h_signal_data->SetMarkerStyle(20);
  h_signal_data->SetMarkerSize(0.6);
  h_signal_data->SetLineColor(1);
  
  h_signal_data->GetXaxis()->SetNdivisions(505);
  h_signal_data->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
  h_signal_data->GetXaxis()->CenterTitle();
  h_signal_data->GetXaxis()->SetTitleSize(0.05);
  h_signal_data->GetXaxis()->SetTitleOffset(0.8);
  h_signal_data->GetXaxis()->SetLabelSize(0.04);

  const double ymax = h_signal_data->GetMaximum();
  h_signal_data->GetYaxis()->SetTitle(Form("Events / [%.1f MeV/c^{2}]", bin_width));
  h_signal_data->GetYaxis()->SetRangeUser(0.01, ymax * 1.8);
  h_signal_data->GetYaxis()->CenterTitle();
  h_signal_data->GetYaxis()->SetTitleSize(0.05);
  h_signal_data->GetYaxis()->SetTitleOffset(1.2);
  h_signal_data->GetYaxis()->SetLabelSize(0.04);
  h_signal_data->GetYaxis()->SetNdivisions(505);
  
  h_signal_data->Draw("E1");
  
  // Draw the total fitted signal (peak + distorted)
  //h_signal_final->Draw("hist same");
  
  // Optionally draw the individual components (peak and distorted)
  h_signal->Draw("hist same");
  // h_background->Draw("hist same"); // optional
  
  TPaveText *pt = new TPaveText(0.65, 0.75, 0.85, 0.85, "NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.04);
  pt->SetTextFont(42);
  pt->AddText(Form("Purity = %.1f%%", updated_purity * 100.));
  pt->AddText(Form("#chi^{2}/ndf = %.2f", chi2_ndf));
  pt->Draw();

  TLegend *leg2 = new TLegend(0.15, 0.65, 0.5, 0.9);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.04);
  leg2->AddEntry(h_signal_data, "Data - all background (pure peak)", "lep");
  //leg2->AddEntry(h_signal_final, "Fitted total signal (peak + distort)", "l");
  leg2->AddEntry(h_signal, "Fitted peak only", "l");
  // leg2->AddEntry(h_background, "Fitted distorted", "l");
  leg2->Draw();
  
  c2->Update();
  c2->Modified();
  c2->SaveAs(output_path + "omega_background_subtracted.pdf");
  
  // ------------------------------------------------------------------
  // 16. Correction weights canvas
  // ------------------------------------------------------------------
  TAxis* xAxis = h_weight_smooth->GetXaxis();
  double xMin = xAxis->GetXmin();
  double xMax = xAxis->GetXmax();
 
  TLine *line_weight = new TLine(xMin, 1, xMax, 1.);
  line_weight->SetLineColor(kGray + 2);
  line_weight->SetLineStyle(2);
  line_weight->SetLineWidth(2);
  
  TCanvas *c_weight = new TCanvas("c_weight", "Correction weights", 800, 600);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.12);

  h_weight_smooth->GetXaxis()->SetNdivisions(505);
  h_weight_smooth->SetTitle("Correction weight for ISR3pi");
  h_weight_smooth->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
  h_weight_smooth->GetYaxis()->SetTitle("Weight");
  h_weight_smooth->GetXaxis()->CenterTitle();
  h_weight_smooth->GetXaxis()->SetTitleSize(0.05);
  h_weight_smooth->GetXaxis()->SetTitleOffset(0.8);
  h_weight_smooth->GetXaxis()->SetLabelSize(0.04);
  h_weight_smooth->SetLineColor(kBlue);
  h_weight_smooth->SetLineWidth(2);

  h_weight_smooth->Draw();
  line_weight->Draw("same");
  c_weight->SaveAs(output_path + "omega_correction_weights.pdf");

  // ------------------------------------------------------------------
  // 17. Save results
  // ------------------------------------------------------------------
  TFile *fout = new TFile(output_path + "omega_fit_results.root", "RECREATE");
  h_data_isr->Write();
  h_isr3pi_corrected->Write("h_isr3pi_corrected");
  h_isr3pi_orig->Write("h_isr3pi_original");
  h_nonReson->Write();
  h_signal->Write();
  h_background->Write();
  h_signal_final->Write();
  h_weight_smooth->Write();
  h_mc_total->Write();
  h_pull->Write();
  total_func->Write();
  fout->Close();

  // ------------------------------------------------------------------
  // 18. Clean up
  // ------------------------------------------------------------------
  gSigTemplate = nullptr;
  gBkgTemplate = nullptr;
  
  if (total_func) delete total_func;
  if (c) delete c;
  if (c2) delete c2;
  if (c_weight) delete c_weight;
  if (h_pull) delete h_pull;
  if (h_signal_data) delete h_signal_data;
  if (h_mc_total) delete h_mc_total;
  if (h_signal_final) delete h_signal_final;

  if (fsfw2d) {
    fsfw2d->Close();
    delete fsfw2d;
  }
  
  if (ftree) {
    ftree->Close();
    delete ftree;
  }

  if (fout) {
    fout->Close();
    delete fout;
  }
  
  std::cout << "\n=== Summary ===" << std::endl;
  std::cout << "Saved " << output_path << "omega_fit_results.root" << std::endl;
  std::cout << "Saved " << output_path << "omega_combined_fit.pdf" << std::endl;
  std::cout << "Saved " << output_path << "omega_background_subtracted.pdf" << std::endl;
  std::cout << "Saved " << output_path << "omega_correction_weights.pdf" << std::endl;
  std::cout << "Fit: α = " << alpha << ", β = " << beta << std::endl;
  std::cout << "χ²/ndf = " << chi2_ndf << std::endl;
}
