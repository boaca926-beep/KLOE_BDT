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
  std::cout << "Data integral: " << h_data->Integral() << std::endl;

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
  std::cout << "EEG integral: " << h_eeg->Integral() << ", sfw = " << eeg_sfw * 2. << "\n"
	    << "ISR3PI integral: " << h_isr3pi->Integral() << ", sfw = " << isr3pi_sfw << "\n"
	    << "OMEGAPI integral: " << h_omegapi->Integral() << ", sfw = " << omegapi_sfw << "\n"
	    << "NonReson integral: " << h_nonReson->Integral() << ", sfw = " << nonReson_sfw << "\n"
	    << "KSL integral: " << h_ksl->Integral() << ", sfw = " << ksl_sfw << "\n"
	    << "MCREST integral: " << h_mcrest->Integral() << ", sfw = " << mcrest_sfw << "\n"
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
  // 5. Normalize templates to unit area
  // ------------------------------------------------------------------
  double sig_int = h_isr3pi->Integral();
  double bkg_int = h_nonReson->Integral();
  if (sig_int > 0) h_isr3pi->Scale(1.0 / sig_int);
  if (bkg_int > 0) h_nonReson->Scale(1.0 / bkg_int);
  
  // Set global pointers for the fit function
  gSigTemplate = h_isr3pi;
  gBkgTemplate = h_nonReson;

  // ------------------------------------------------------------------
  // 6. Template fit - MINIMAL CHANGE: expanded fit range
  // ------------------------------------------------------------------
  double peak_low = 700;   // MeV - was 740
  double peak_high = 850;  // MeV - was 820
  
  TF1 *total_func = new TF1("total_func", template_sum, low, high, 2);
  total_func->SetParameters(1000, 1000);
  total_func->SetParNames("alpha", "beta");
  
  // Single fit with quality assessment
  TFitResultPtr r = h_data_isr->Fit(total_func, "RQS", "", peak_low, peak_high);
  
  // Check if fit converged
  if (!r->IsValid()) {
    std::cerr << "WARNING: Fit did not converge!" << std::endl;
  }
  
  double alpha = total_func->GetParameter(0);
  double beta  = total_func->GetParameter(1);
  double chi2 = r->Chi2();
  int ndf = r->Ndf();
  double chi2_ndf = chi2 / ndf;
  
  std::cout << "Fit results: α = " << alpha << ", β = " << beta << std::endl;
  std::cout << "Fit quality: χ² = " << chi2 << ", ndf = " << ndf << ", χ²/ndf = " << chi2_ndf << std::endl;

  // ------------------------------------------------------------------
  // 7. Create signal & background histograms (scaled)
  // ------------------------------------------------------------------
  TH1D *h_signal = (TH1D*) h_isr3pi->Clone("h_signal");
  TH1D *h_background = (TH1D*) h_nonReson->Clone("h_background");
  h_signal->SetDirectory(0); 
  h_background->SetDirectory(0);
  h_signal->Scale(alpha);
  h_background->Scale(beta);
  h_signal->SetLineColor(kBlue);
  h_signal->SetLineWidth(2);
  h_background->SetLineColor(kRed);
  h_background->SetLineStyle(3);
  h_background->SetLineWidth(2);
  
  // ------------------------------------------------------------------
  // 8. Compute correction weights (using h_signal as desired shape)
  // ------------------------------------------------------------------
  TH1D *h_weight = (TH1D*) h_data_isr->Clone("h_weight");
  h_weight->SetDirectory(0);
  h_weight->Reset();
  for (int bin = 1; bin <= h_weight->GetNbinsX(); ++bin) {
    double x = h_weight->GetBinCenter(bin);
    double mc_val = h_isr3pi->GetBinContent(bin);
    double desired = h_signal->GetBinContent(bin);
    if (x >= peak_low && x <= peak_high && mc_val > 0 && desired > 0) {
      double ratio = desired / mc_val;
      if (ratio > 2.5) ratio = 2.5;
      if (ratio < 0.4) ratio = 0.4;
      h_weight->SetBinContent(bin, ratio);
    } else {
      h_weight->SetBinContent(bin, 1.0);
    }
  }
  // Smooth weights
  TH1D *h_weight_smooth = (TH1D*) h_weight->Clone("h_weight_smooth");
  h_weight_smooth->SetDirectory(0);
  for (int bin = 2; bin <= h_weight_smooth->GetNbinsX()-1; ++bin) {
    double w_avg = (h_weight->GetBinContent(bin-1) +
		    h_weight->GetBinContent(bin) +
		    h_weight->GetBinContent(bin+1)) / 3.0;
    h_weight_smooth->SetBinContent(bin, w_avg);
  }

  // ------------------------------------------------------------------
  // 9. Apply correction to ISR3pi MC (use original scaled MC)
  // ------------------------------------------------------------------
  // Reload original ISR3pi with proper normalization for correction
  TH1D *h_isr3pi_orig = makeScaledHist("TISR3PI_SIG_PEAK", isr3pi_sfw);
  if (!h_isr3pi_orig) {
    std::cerr << "ERROR: Cannot reload ISR3pi for correction" << std::endl;
    return;
  }
  
  TH1D *h_isr3pi_corrected = (TH1D*) h_isr3pi_orig->Clone("h_isr3pi_corrected");
  h_isr3pi_corrected->SetDirectory(0);
  for (int bin = 1; bin <= h_isr3pi_corrected->GetNbinsX(); ++bin) {
    double w = h_weight_smooth->GetBinContent(bin);
    double old = h_isr3pi_corrected->GetBinContent(bin);
    h_isr3pi_corrected->SetBinContent(bin, old * w);
    double err = h_isr3pi_corrected->GetBinError(bin);
    h_isr3pi_corrected->SetBinError(bin, err * w);
  }
  // Renormalise
  double orig_int = h_isr3pi_orig->Integral(peak_low, peak_high);
  double new_int = h_isr3pi_corrected->Integral(peak_low, peak_high);
  if (new_int > 0 && orig_int > 0) {
    double renorm = orig_int / new_int;
    h_isr3pi_corrected->Scale(renorm);
    std::cout << "Renormalisation factor: " << renorm << std::endl;
  }
  h_isr3pi_corrected->SetLineColor(kGreen);
  h_isr3pi_corrected->SetLineWidth(2);

  // ------------------------------------------------------------------
  // 10. Build total MC sum for plotting
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

  // ------------------------------------------------------------------
  // 11. Pull distribution (full range)
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
  // 12. Main plotting (like correct_and_plot.C)
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

  h_signal->SetLineColor(kBlue);
  h_signal->SetLineStyle(1);
  h_signal->SetLineWidth(2);
  
  TH1D *h_signal_data = (TH1D*) h_data->Clone("h_signal_data");
  
  h_data->SetMarkerStyle(20);
  h_data->SetMarkerSize(0.6);
  h_data->Draw("E0");
  h_mc_total->Draw("hist same");
  for (auto h : comps) if (h) h->Draw("hist same");
  
  h_data->GetYaxis()->SetTitle(Form("Events / [%.1f MeV/c^{2}]", bin_width));
  h_data->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
  h_data->GetYaxis()->CenterTitle();
  h_data->GetXaxis()->SetTitleSize(0.06);
  h_data->GetXaxis()->SetLabelOffset(0.1);
  h_data->GetYaxis()->SetTitleSize(0.07);
  h_data->GetYaxis()->SetTitleOffset(0.7);
  h_data->GetYaxis()->SetLabelSize(0.04);
  h_data->GetYaxis()->SetNdivisions(505);

  // Legend - same order as correct_and_plot.C
  TLegend *leg = new TLegend(0.15, 0.35, 0.6, 0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_data, "Data", "lep");
  leg->AddEntry(h_mc_total, "Total MC", "l");
  leg->AddEntry(h_signal, "#omega peak (signal)", "l");
  leg->AddEntry(h_background, "Distorted signal", "l");
  leg->AddEntry(h_eeg, "e^{+}e^{-}#gamma", "l");
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
  h_pull->GetYaxis()->SetRangeUser(-5, 5);
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
  // 13. Background-subtracted ω signal (ADAPTED: like linear version)
  // ------------------------------------------------------------------
  // Start from raw data and subtract ALL backgrounds explicitly
  //TH1D *h_signal_data = (TH1D*) h_data->Clone("h_signal_data");
  h_signal_data->Add(h_eeg, -1.0);
  h_signal_data->Add(h_omegapi, -1.0);
  h_signal_data->Add(h_ksl, -1.0);
  h_signal_data->Add(h_mcrest, -1.0);
  h_signal_data->Add(h_background, -1.0);  // non-resonant ISR
  
  // Set negative bins to zero
  for (int bin = 1; bin <= h_signal_data->GetNbinsX(); ++bin) {
    if (h_signal_data->GetBinContent(bin) < 0) {
      h_signal_data->SetBinContent(bin, 0);
    }
  }
  
  TCanvas *c2 = new TCanvas("c2", "Background-subtracted ω signal", 1200, 700);
  
  c2->cd();
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);
 
  h_signal_data->SetMarkerStyle(20);
  h_signal_data->SetMarkerSize(0.6);
  
  // Apply ALL axis settings BEFORE drawing
  // X-axis settings - FIXED
  h_signal_data->GetXaxis()->SetNdivisions(505);
  h_signal_data->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
  h_signal_data->GetXaxis()->CenterTitle();
  h_signal_data->GetXaxis()->SetTitleSize(0.05);
  h_signal_data->GetXaxis()->SetTitleOffset(0.8);   // SMALLER = closer to axis
  h_signal_data->GetXaxis()->SetLabelSize(0.04);
  //h_signal_data->GetXaxis()->SetRangeUser(700, 900);
 
  h_signal_data->GetYaxis()->SetTitle(Form("Events / [%.1f MeV/c^{2}]", bin_width));
  h_signal_data->GetYaxis()->CenterTitle();
  h_signal_data->GetYaxis()->SetTitleSize(0.05);
  h_signal_data->GetYaxis()->SetTitleOffset(1.2);
  h_signal_data->GetYaxis()->SetLabelSize(0.04);
  h_signal_data->GetYaxis()->SetNdivisions(505);
  
  // Draw with full options
  h_signal_data->Draw("E1");
  
  // Overlay fitted signal for comparison
  h_signal->Draw("hist same");
  
  TLegend *leg2 = new TLegend(0.2, 0.7, 0.5, 0.9);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.04);
  leg2->AddEntry(h_signal_data, "Data - all backgrounds", "lep");
  leg2->AddEntry(h_signal, "Fitted #omega peak (signal)", "l");
  leg2->Draw();
  
  // Force canvas to update
  c2->Update();
  c2->Modified();
  
  c2->SaveAs(output_path + "omega_background_subtracted.pdf");
  
  // ------------------------------------------------------------------
  // 14. Correction weights canvas
  // ------------------------------------------------------------------
  TCanvas *c_weight = new TCanvas("c_weight", "Correction weights", 800, 600);
  h_weight_smooth->SetTitle("Correction weight for ISR3pi");
  h_weight_smooth->GetXaxis()->SetTitle("M_{3π} [MeV/c^{2}]");
  h_weight_smooth->GetYaxis()->SetTitle("Weight");
  h_weight_smooth->SetLineColor(kBlue);
  h_weight_smooth->SetLineWidth(2);
  h_weight_smooth->Draw();
  c_weight->SaveAs(output_path + "omega_correction_weights.pdf");

  // ------------------------------------------------------------------
  // 15. Save results
  // ------------------------------------------------------------------
  TFile *fout = new TFile(output_path + "omega_fit_results.root", "RECREATE");
  h_data_isr->Write();
  h_isr3pi->Write();
  h_nonReson->Write();
  h_signal->Write();
  h_background->Write();
  h_weight_smooth->Write();
  h_isr3pi_corrected->Write();
  h_mc_total->Write();
  h_pull->Write();
  total_func->Write();
  fout->Close();

  // ------------------------------------------------------------------
  // 16. Clean up
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
