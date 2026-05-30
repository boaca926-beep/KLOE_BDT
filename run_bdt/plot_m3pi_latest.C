// plot_m3pi_latest.C
// Uses the latest scaling factors from the RooFit 2D fit (IMpp vs deltaE)
// to produce the 3π mass projection, with optional corrected ISR3π shape.

// For validation only: Keep the original MC for efficiency and yields; use the corrected shape only for a “post‑fit” comparison plot. This is what you have done – it improves the visual agreement without altering the physics results.

// In your current approach, the correction does not affect bin‑by‑bin efficiency because you are only using it for plotting. The original MC (unweighted) remains the source for efficiency. This is perfectly valid as long as you do not claim that the corrected shape is part of the analysis model – it is a diagnostic tool.

// plot_m3pi_latest.C
// Uses the latest scaling factors from the RooFit 2D fit (IMpp vs deltaE)
// to produce the 3π mass projection, with optional corrected ISR3π shape.
// Follows line styles from plot_compr.C (numeric styles: 2,3,4,5,6,7)
// Colors remain as originally defined.

// plot_m3pi_latest.C
// Uses the latest scaling factors from the RooFit 2D fit (IMpp vs deltaE)
// to produce the 3π mass projection, with corrected ISR3π shape from either
// corrected_isr3pi_sample.root or corrected_isr3pi_tmp.root.

/*
The template fit is generally better for the reasons discussed earlier:

1.    It preserves the shape correlations from MC (important for efficiency).

2.    It uses realistic, data‑driven templates from the BDT classification.

3.    It has only two free parameters, making the fit more stable.

4.    The side‑band fit with a polynomial background can introduce biases depending on the chosen sideband range.
*/    

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "../header_bdt/sfw2d_bdt.txt"
#include "../header_bdt/correct_omega.h"
#include "../header_bdt/path.h"

void pt_style(TPaveText *pt, TString text) {
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.04);
  pt->SetTextFont(42);
  pt->AddText(text);
  //pt->AddText(line4);
  
}

void plot_m3pi_latest() {
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  
  // Latest scaling factors
  std::cout << "Using latest scaling factors from RooFit 2D fit:\n"
            << "  EEG     : " << eeg_sfw << "\n"
            << "  ISR3pi  : " << isr3pi_sfw << "\n"
            << "  OmegaPi : " << omegapi_sfw << "\n"
            << "  EtaGamma: " << etagam_sfw << "\n"
            << "  KSL     : " << ksl_sfw << "\n"
            << "  MC Rest : " << mcrest_sfw << std::endl;

  // Open tree file
  TString treeFile = "/home/kloe/Desktop/input_bdt_TDATA_chain/cut/tree_pre_bdt.root";
  TFile *ftree = TFile::Open(treeFile);
  if (!ftree || ftree->IsZombie()) {
    std::cerr << "ERROR: cannot open " << treeFile << std::endl;
    return;
  }

  // 3π mass range (MeV)
  TTree *tdata = (TTree*) ftree->Get("TDATA");
  if (!tdata) { std::cerr << "No TDATA tree.\n"; return; }
  double mass_min = 600.0;
  double mass_max = 1000.0;

  double omega_min = mass_min; // 760.0
  double omega_max = mass_max; // 1000.0
  
  const int nBins = 200;

  // Helper to fill a histogram from a tree
  auto fillHist = [&](const char* treeName, double scale, int color, int style, const char* title) -> TH1D* {
    TTree *t = (TTree*) ftree->Get(treeName);
    if (!t) return nullptr;
    TH1D *h = new TH1D(title, "", nBins, mass_min, mass_max);
    h->Sumw2();
    double m;
    t->SetBranchAddress("Br_m3pi_bdt", &m);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      h->Fill(m);
    }
    h->Scale(scale);
    h->SetLineColor(color);
    h->SetLineStyle(style);
    return h;
  };

  // Data
  TH1D *h_data = new TH1D("h_data", "", nBins, mass_min, mass_max);
  h_data->SetMarkerStyle(21);
  h_data->SetMarkerSize(0.7);
  h_data->SetLineColor(1);
  h_data->Sumw2();
  double m3pi;
  tdata->SetBranchAddress("Br_m3pi_bdt", &m3pi);
  for (Long64_t i = 0; i < tdata->GetEntries(); ++i) {
    tdata->GetEntry(i);
    h_data->Fill(m3pi);
  }

  // --- Declare pointers for components that will be used later ---
  TH1D *h_eeg = nullptr;
  TH1D *h_omegapi = nullptr;
  TH1D *h_ksl = nullptr;
  TH1D *h_etagam = nullptr;

  // MC components (fixed backgrounds)
  std::vector<TH1D*> comps;
  h_eeg = fillHist("TEEG",        eeg_sfw,     6, 7, "EEG");
  comps.push_back(h_eeg);
  h_omegapi = fillHist("TOMEGAPI",    omegapi_sfw, 7, 5, "#omega#pi^{0}");
  comps.push_back(h_omegapi);
  h_ksl = fillHist("TKSL",        ksl_sfw,     28, 4, "K_{S}K_{L}");
  comps.push_back(h_ksl);
  h_etagam = fillHist("TETAGAM",     etagam_sfw,  3, 3, "#eta#gamma");
  comps.push_back(h_etagam);

  // ---------- ISR3π: try corrected histograms from either _sample.root or _tmp.root ----------
  TH1D *h_isr_peak = nullptr;
  TH1D *h_isr_bkg  = nullptr;
  
  const char* corrFiles[2];
  TString file1 = output_path + "corrected_isr3pi_hybrid.root";
  //TString file1 = output_path + "corrected_isr3pi_hybrid_lower_linear.root";
  TString file2 = output_path + "corrected_isr3pi_tmp.root";
  corrFiles[0] = file1.Data();
  corrFiles[1] = file2.Data();

  TFile *fcorr = nullptr;   // <-- DECLARE and initialise
  
  for (int i = 0; i < 1; ++i) {
    fcorr = TFile::Open(corrFiles[i]);
    if (fcorr && !fcorr->IsZombie()) {
      TH1D *h_tmp_peak = (TH1D*) fcorr->Get("h_signal");
      TH1D *h_tmp_bkg  = (TH1D*) fcorr->Get("h_background");
      if (h_tmp_peak && h_tmp_bkg) {
        std::cout << "Using separated ISR3π peak and background from " << corrFiles[i] << std::endl;
        h_isr_peak = (TH1D*) h_tmp_peak->Clone("h_isr_peak");
        h_isr_bkg  = (TH1D*) h_tmp_bkg->Clone("h_isr_bkg");
        h_isr_peak->SetTitle("3#pi (peak)");
        h_isr_bkg->SetTitle("Combinatorial");
        h_isr_peak->SetLineColor(kBlue);
        h_isr_peak->SetLineStyle(1);
        h_isr_peak->SetLineWidth(2);
        h_isr_bkg->SetLineColor(kRed);
        h_isr_bkg->SetLineStyle(3);
        h_isr_bkg->SetLineWidth(2);
        break; // success
      } else {
        std::cerr << "WARNING: " << corrFiles[i] << " exists but missing h_signal/h_background." << std::endl;
        fcorr->Close();
        fcorr = nullptr;
      }
    }
  }

  // Fallback if no correction file worked
  if (!h_isr_peak || !h_isr_bkg) {
    h_isr_peak = fillHist("TISR3PI_SIG", isr3pi_sfw, 4, 2, "ISR3#pi");
    h_isr_peak->SetTitle("ISR3#pi");
    h_isr_peak->SetLineColor(kBlue);
    h_isr_peak->SetLineStyle(1);
    h_isr_peak->SetLineWidth(2);
    comps.push_back(h_isr_peak);
    std::cout << "Using original ISR3π histogram from tree (no separation).\n";
    if (fcorr) { fcorr->Close(); fcorr = nullptr; }
  } else {
    comps.push_back(h_isr_peak);
    comps.push_back(h_isr_bkg);
  }

  // MC Rest (combination of TKPM, TRHOPI, TBKGREST)
  TH1D *h_mcrest = new TH1D("h_mcrest", "", nBins, mass_min, mass_max);
  h_mcrest->Sumw2();
  for (const char* name : {"TKPM", "TRHOPI", "TBKGREST"}) {
    TTree *t = (TTree*) ftree->Get(name);
    if (!t) continue;
    double m;
    t->SetBranchAddress("Br_m3pi_bdt", &m);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      h_mcrest->Fill(m);
    }
  }
  h_mcrest->Scale(mcrest_sfw);
  h_mcrest->SetLineColor(37);
  h_mcrest->SetLineStyle(6);
  comps.push_back(h_mcrest);

  // Build total MC sum
  TH1D *h_mc_total = (TH1D*) h_mcrest->Clone("h_mc_total");
  h_mc_total->Reset();
  h_mc_total->Sumw2();
  for (auto h : comps) if (h) h_mc_total->Add(h);
  h_mc_total->SetLineColor(kRed);
  h_mc_total->SetLineWidth(2);
  h_mc_total->SetLineStyle(1);

  // Pulls
  TH1D *h_pull = new TH1D("h_pull", "", nBins, mass_min, mass_max);
  for (int bin = 1; bin <= nBins; ++bin) {
    double d = h_data->GetBinContent(bin);
    double m = h_mc_total->GetBinContent(bin);
    double err = std::sqrt(d + m);
    if (err > 0) h_pull->SetBinContent(bin, (d - m)/err);
    else h_pull->SetBinContent(bin, 0);
  }
  h_pull->SetMarkerStyle(20);
  h_pull->SetMarkerSize(0.6);
  h_pull->SetLineWidth(0);

  // Canvas and pads
  TCanvas *c = new TCanvas("c", "3π mass projection with pulls", 1200, 700);
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
  
  h_data->Draw("E0");
  h_mc_total->Draw("hist same");
  for (auto h : comps) if (h) h->Draw("hist same");
  h_data->GetYaxis()->SetTitle(Form("Events / [%.1f MeV/c^{2}]", bin_width));
  h_data->GetXaxis()->SetLabelOffset(0.1);
  //h_data->GetYaxis()->SetLabelOffset(0.005);
  h_data->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
  h_data->GetXaxis()->SetRangeUser(omega_min, omega_max);
  h_data->GetYaxis()->CenterTitle();
  h_data->GetXaxis()->SetTitleSize(0.06);
  h_data->GetYaxis()->SetTitleSize(0.07);
  h_data->GetYaxis()->SetTitleOffset(.7);
  h_data->GetYaxis()->SetLabelSize(0.04);
  h_data->GetYaxis()->SetNdivisions(505);
  
  // Legend – order must match comps vector
  TLegend *leg = new TLegend(0.65, 0.25, 0.9, 0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);           
  leg->AddEntry(h_data, "Data", "lep");
  leg->AddEntry(h_mc_total, "Total MC", "l");

  int idx = 0;
  // EEG
  if (comps.size() > idx) leg->AddEntry(comps[idx++], "EEG", "l");
  // OmegaPi
  if (comps.size() > idx) leg->AddEntry(comps[idx++], "#omega#pi^{0}", "l");
  // KSL
  if (comps.size() > idx) leg->AddEntry(comps[idx++], "K_{S}K_{L}", "l");
  // EtaGamma
  if (comps.size() > idx) leg->AddEntry(comps[idx++], "#eta#gamma", "l");

  // ISR3π part – conditional on whether we have separated peak+background
  if (h_isr_peak && h_isr_bkg) {

    // Two entries: peak (corrected) and non‑resonant background
    if (comps.size() > idx) leg->AddEntry(comps[idx++], "3#pi (peak) corrected", "l");
    if (comps.size() > idx) leg->AddEntry(comps[idx++], "Combinatorial", "l");
  } else if (h_isr_peak && !h_isr_bkg) {
    // Fallback: single ISR3π component
    if (comps.size() > idx) leg->AddEntry(comps[idx++], "ISR3#pi", "l");
  }
  // MC Rest (Others)
  if (comps.size() > idx) leg->AddEntry(comps[idx++], "Others", "l");

  leg->Draw();

  // Lower pad (pulls)
  c->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.3);
  pad2->SetLeftMargin(0.12);
  pad2->Draw();
  pad2->cd();
  gPad->SetGrid();   // adds grid lines to both X and Y axes

  h_pull->GetXaxis()->SetTitle("M_{3#pi} [MeV]");
  h_pull->GetXaxis()->SetTitleSize(0.12);
  h_pull->GetXaxis()->SetTitleOffset(1.0);
  h_pull->GetXaxis()->SetLabelSize(0.1);
  h_pull->GetYaxis()->SetTitle("Pull");
  h_pull->GetYaxis()->SetTitleSize(0.2);       // increased for readability
  h_pull->GetYaxis()->SetTitleOffset(0.2);      // avoid overlap
  h_pull->GetYaxis()->SetLabelSize(0.1);
  //h_pull->GetXaxis()->SetRangeUser(760, 820);   // optional: zoom to ω peak
  h_pull->GetYaxis()->SetRangeUser(-50, 50);
  h_pull->GetYaxis()->SetNdivisions(505);
  h_pull->GetXaxis()->CenterTitle();
  h_pull->GetYaxis()->CenterTitle();
  h_pull->Draw("P");
 
  TLine *line = new TLine(omega_min, 0, omega_max, 0);
  line->SetLineStyle(2);
  line->Draw();

  c->SaveAs(output_path + "m3pi_projection_with_pulls_sample.pdf");
  std::cout << "\nSaved " + output_path + "m3pi_projection_with_pulls_sample.pdf\n";
  delete c;
  
  // --- Minimal addition: background‑subtracted ω signal ---
  TH1D *h_signal_data = (TH1D*) h_data->Clone("h_signal_data");
  h_signal_data->Add(h_eeg, -1.0);
  h_signal_data->Add(h_omegapi, -1.0);
  h_signal_data->Add(h_ksl, -1.0);
  h_signal_data->Add(h_etagam, -1.0);
  h_signal_data->Add(h_mcrest, -1.0);
  if (h_isr_bkg) h_signal_data->Add(h_isr_bkg, -1.0);
  // Set negative bins to zero
  for (int bin = 1; bin <= h_signal_data->GetNbinsX(); ++bin)
    if (h_signal_data->GetBinContent(bin) < 0) h_signal_data->SetBinContent(bin, 0);

  TCanvas *c2 = new TCanvas("c2", "Background‑subtracted ω signal", 1200, 700);
  c2->SetBottomMargin(0.12);
  c2->SetLeftMargin(0.12);

  h_signal_data->SetMarkerStyle(20);
  h_signal_data->SetMarkerSize(0.6);
  h_signal_data->GetXaxis()->SetLabelOffset(0.007);   // keep default small value
  h_signal_data->GetYaxis()->SetTitle(Form("Events / [%.1f MeV/c^{2}]", bin_width));
  h_signal_data->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
  h_signal_data->GetXaxis()->SetRangeUser(omega_min, omega_max - 50.0);
  h_signal_data->GetXaxis()->CenterTitle();
  h_signal_data->GetYaxis()->CenterTitle();
  h_signal_data->GetXaxis()->SetTitleSize(0.05);
  h_signal_data->GetYaxis()->SetTitleSize(0.05);
  h_signal_data->GetYaxis()->SetTitleOffset(1.2);
  h_signal_data->GetYaxis()->SetLabelSize(0.04);

  h_signal_data->Draw("E0");
  if (h_isr_peak && h_isr_bkg) {
    h_isr_peak->SetLineColor(kBlue);
    h_isr_peak->Draw("hist same");
    h_isr_bkg->SetLineColor(kRed);
    h_isr_bkg->Draw("hist same");

    // Define the mass window (MeV/c²)
    double low = 650.0, high = 900.0;
    int bin_low = h_isr_peak->FindBin(low);
    int bin_high = h_isr_peak->FindBin(high);
    
    double peak_entries = h_isr_peak->Integral(bin_low, bin_high);
    double comb_entries = h_isr_bkg->Integral(bin_low, bin_high);
    double total = peak_entries + comb_entries;
    double purity = (total > 0) ? (peak_entries / total) * 100.0 : 0.0;
    
    cout << "Data-driven peak entries in mass region m3pi = [" << low << ", " << high << "] MeV/c^2:\n"
         << "peak entries: " << peak_entries << "\n"
         << "combinatorial entries: " << comb_entries << "\n"
         << "purity: " << purity << "%\n";

    // mass range for purity estimation
    TLine *line1 = new TLine(650, 0, 650, 400);
    line1->SetLineColor(kBlack);
    line1->SetLineWidth(3);

    TLine *line2 = new TLine(900, 0, 900, 400);
    line2->SetLineColor(kBlack);
    line2->SetLineWidth(3);
    line2->Draw();
  
    line1->Draw();
    line2->Draw("same");
  
    // Create a transparent TPaveText
    TPaveText *pt = new TPaveText(0.6, 0.6, 0.7, 0.68, "NDC");
    TString line = Form("Data-driven purity = %.2f", purity* 1e-2);
    pt_style(pt, line);
    pt->SetTextColor(kRed);
    pt->Draw("same");
    
  }

  TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.04);           
  leg2->AddEntry(h_signal_data, "Data - backgrounds", "lep");
  leg2->AddEntry(h_isr_peak, "Corrected #omega peak", "l");
  leg2->AddEntry(h_isr_bkg, "Comb. background", "l");
  leg2->Draw("same");
  
  c2->SaveAs(output_path + "background_subtracted_omega_sample.pdf");
  delete c2;

  if (h_isr_peak && h_isr_bkg) {
    TString outFileName = output_path + "hist_sample.root";   // use same directory as PDFs
    TFile *f_output = new TFile(outFileName, "RECREATE");
    if (f_output && !f_output->IsZombie()) {
      h_signal_data->Write();
      h_isr_peak->Write();
      h_isr_bkg->Write();
      f_output->Close();
    }
    delete f_output;
  }

  ftree->Close();
  if (fcorr) fcorr->Close();

  gSystem->Exit(0);
  
}
