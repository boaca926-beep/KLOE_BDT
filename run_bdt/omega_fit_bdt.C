// correct_omega_peak_sample.C – template fit with proper memory management
// Adapted to follow sfw2d.C pattern: close input file after using globals.

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <cmath>

#include "../header_bdt/sfw2d.txt"
#include "../header_bdt/omega_fit_bdt.h"

// Global pointers for template histograms (detached from file)
TH1D *gSigTemplate = nullptr;
TH1D *gNonResonTemplate = nullptr;
TH1D *gBkgTemplate = nullptr;

// Fit function
Double_t template_sum(Double_t *x, Double_t *par) {
  int bin = gSigTemplate->FindBin(x[0]);
  Double_t sig = gSigTemplate->GetBinContent(bin);
  Double_t non_reson = gNonResonTemplate->GetBinContent(bin);
  Double_t bkg = gBkgTemplate->GetBinContent(bin); // bkg withoug non-reson
  return par[0] * sig + par[1] * bkg + non_reson;
}

void omega_fit_bdt() {
  // ------------------------------------------------------------------
  // 1. Open tree file
  // ------------------------------------------------------------------
  //TString treeFile = "/home/kloe/Desktop/input_bdt_TDATA_chain/cut/tree_pre_bdt.root";
  TFile *ftree = TFile::Open(treeFile);
  if (!ftree || ftree->IsZombie()) {
    std::cerr << "ERROR: cannot open " << treeFile << std::endl;
    return;
  }

  TTree *tdata = (TTree*) ftree->Get("TDATA");
  if (!tdata) { std::cerr << "ERROR: TDATA not found." << std::endl; return; }
  
  // Determine unit
  double mtest;
  tdata->SetBranchAddress("Br_m3pi_bdt", &mtest);
  tdata->GetEntry(0);
  bool is_mev = (mtest > 10);
  double low, high;
  int nbins = NBINS;
  
  low  = MASS_MIN;          // e.g., 600.0 MeV
  high = MASS_MAX;          // e.g., 1000.0 MeV

  std::cout << "Mass unit: " << " range [" << low << ", " << high << "] MeV/c² \n";

  // ------------------------------------------------------------------
  // 2. Data histogram (detached from file)
  // ------------------------------------------------------------------
  TH1D *h_data = new TH1D("h_data", "", nbins, low, high);
  h_data->Sumw2();
  h_data->SetDirectory(0);                    // DETACH
  tdata->SetBranchAddress("Br_m3pi_bdt", &mtest);
  for (Long64_t i = 0; i < tdata->GetEntries(); ++i) {
    tdata->GetEntry(i);
    h_data->Fill(mtest);
  }
  std::cout << "Data integral: " << h_data->Integral() << std::endl;

  // ------------------------------------------------------------------
  // 3. Load scaled MC components (detach each histogram)
  // ------------------------------------------------------------------
  auto makeScaledHist = [&](const char* tname, double scale) -> TH1D* {
    TTree *t = (TTree*) ftree->Get(tname);
    if (!t) return nullptr;
    if (!t->GetBranch("Br_m3pi_bdt")) return nullptr;
    TH1D *h = new TH1D(Form("h_%s", tname), "", nbins, low, high);
    h->Sumw2();
    h->SetDirectory(0);                     // DETACH
    double val;
    t->SetBranchAddress("Br_m3pi_bdt", &val);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      h->Fill(val);
    }
    h->Scale(scale);
    return h;
  };

  TH1D *h_eeg       = makeScaledHist("TEEG", eeg_sfw * 2.);
  TH1D *h_isr3pi    = makeScaledHist("TISR3PI_SIG_PEAK", isr3pi_sfw);
  TH1D *h_nonReson  = makeScaledHist("TISR3PI_SIG_NON_RESON", nonReson_sfw);
  TH1D *h_omegapi   = makeScaledHist("TOMEGAPI", omegapi_sfw);
  TH1D *h_ksl       = makeScaledHist("TKSL", ksl_sfw);

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

  if (!h_isr3pi) { std::cerr << "ERROR: No ISR3pi histogram." << std::endl; return; }

  // Summary
  std::cout << "EEG integral: " << h_eeg->Integral() << std::endl;

}
