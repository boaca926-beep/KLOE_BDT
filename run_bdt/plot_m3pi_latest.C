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
// to produce the 3π mass projection, with separated ISR3π peak and background.

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

void plot_m3pi_latest() {
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
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

  // MC components (fixed backgrounds)
  std::vector<TH1D*> comps;
  comps.push_back(fillHist("TEEG",        eeg_sfw,     6, 7, "EEG"));
  comps.push_back(fillHist("TOMEGAPI",    omegapi_sfw, 7, 5, "#omega#pi^{0}"));
  comps.push_back(fillHist("TKSL",        ksl_sfw,     28, 4, "K_{S}K_{L}"));
  comps.push_back(fillHist("TETAGAM",     etagam_sfw,  3, 3, "#eta#gamma"));

  // ---------- ISR3π: try to use separated peak and background ----------
  TH1D *h_isr_peak = nullptr;
  TH1D *h_isr_bkg  = nullptr;
  TFile *fcorr = TFile::Open("corrected_isr3pi_tmp.root");
  if (fcorr && !fcorr->IsZombie()) {
    TH1D *h_tmp_peak = (TH1D*) fcorr->Get("h_signal");
    TH1D *h_tmp_bkg  = (TH1D*) fcorr->Get("h_background");
    if (h_tmp_peak && h_tmp_bkg) {
      std::cout << "Using separated ISR3π peak and background from corrected_isr3pi.root\n";
      h_isr_peak = (TH1D*) h_tmp_peak->Clone("h_isr_peak");
      h_isr_bkg  = (TH1D*) h_tmp_bkg->Clone("h_isr_bkg");
      h_isr_peak->SetTitle("3#pi (peak)");
      h_isr_bkg->SetTitle("Non-resonant 3#pi");
      // Styles: peak solid, background dotted
      h_isr_peak->SetLineColor(kBlue);
      h_isr_peak->SetLineStyle(1);
      h_isr_peak->SetLineWidth(2);
      h_isr_bkg->SetLineColor(9);
      h_isr_bkg->SetLineStyle(3);
      h_isr_bkg->SetLineWidth(2);
    } else {
      std::cerr << "WARNING: Separated histograms not found. Falling back to original.\n";
    }
  }
  if (!h_isr_peak || !h_isr_bkg) {
    // Fallback: use original tree-based histogram (scaled)
    h_isr_peak = fillHist("TISR3PI_SIG", isr3pi_sfw, 4, 2, "ISR3#pi");
    h_isr_peak->SetTitle("ISR3#pi");
    h_isr_peak->SetLineColor(kBlue);
    h_isr_peak->SetLineStyle(1);
    h_isr_peak->SetLineWidth(2);
    // No separate background; we just push one component.
    comps.push_back(h_isr_peak);
    std::cout << "Using original ISR3π histogram from tree (no separation).\n";
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
  TCanvas *c = new TCanvas("c", "3π mass projection with pulls", 800, 800);
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

  h_data->Draw("E1");
  h_mc_total->Draw("hist same");
  for (auto h : comps) if (h) h->Draw("hist same");

  h_data->GetXaxis()->SetTitle("M_{3#pi} [MeV]");
  h_data->GetYaxis()->SetTitle("Events");
  h_data->GetYaxis()->CenterTitle();
  h_data->GetYaxis()->SetTitleSize(0.05);
  h_data->GetYaxis()->SetTitleOffset(1.2);
  h_data->GetYaxis()->SetLabelSize(0.04);

  // Legend
  TLegend *leg = new TLegend(0.65, 0.25, 0.9, 0.9);
  leg->SetTextFont(132);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h_data, "Data", "lep");
  leg->AddEntry(h_mc_total, "Total MC", "l");
  // Add fixed components (order as in comps)
  int idx = 0;
  if (comps.size() > idx) leg->AddEntry(comps[idx++], "EEG", "l");
  if (comps.size() > idx) leg->AddEntry(comps[idx++], "#omega#pi^{0}", "l");
  if (comps.size() > idx) leg->AddEntry(comps[idx++], "K_{S}K_{L}", "l");
  if (comps.size() > idx) leg->AddEntry(comps[idx++], "#eta#gamma", "l");
  // For ISR3π: if separated, two entries; else one
  if (h_isr_peak && h_isr_bkg) {
    if (comps.size() > idx) leg->AddEntry(comps[idx++], "3#pi (peak)", "l");
    if (comps.size() > idx) leg->AddEntry(comps[idx++], "Non-resonant 3#pi", "l");
  } else if (h_isr_peak && !h_isr_bkg) {
    if (comps.size() > idx) leg->AddEntry(comps[idx++], "ISR3#pi", "l");
  }
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

  h_pull->GetXaxis()->SetTitle("M_{3#pi} [MeV]");
  h_pull->GetXaxis()->SetTitleSize(0.12);
  h_pull->GetXaxis()->SetTitleOffset(1.0);
  h_pull->GetXaxis()->SetLabelSize(0.1);
  h_pull->GetYaxis()->SetTitle("Pull");
  h_pull->GetYaxis()->SetTitleSize(0.12);
  h_pull->GetYaxis()->SetTitleOffset(0.5);
  h_pull->GetYaxis()->SetLabelSize(0.1);
  h_pull->GetYaxis()->SetRangeUser(-10, 10);
  h_pull->GetYaxis()->SetNdivisions(505);
  h_pull->GetXaxis()->CenterTitle();
  h_pull->GetYaxis()->CenterTitle();
  h_pull->Draw("P");

  TLine *line = new TLine(mass_min, 0, mass_max, 0);
  line->SetLineStyle(2);
  line->Draw();

  c->SaveAs("m3pi_projection_with_pulls.pdf");
  std::cout << "\nSaved m3pi_projection_with_pulls.pdf\n";

  delete c;
  ftree->Close();
  if (fcorr) fcorr->Close();
}
