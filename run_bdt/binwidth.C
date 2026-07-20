// ============================================================================
// binwidth.C
// Determine core resolution using a double Gaussian fit in a restricted range.
// Extracts the narrower Gaussian as the core resolution.
// Always uses the uncorrected fitted energies (Br_e1_fit, Br_e2_fit).
// Improved initial parameters from histogram FWHM.
// ============================================================================

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TGaxis.h>
#include <TStyle.h>
#include <TKey.h>
#include <iostream>

using namespace std;

// ----------------------------------------------------------------------------
// Main macro
// ----------------------------------------------------------------------------
void binwidth(const TString tree_type = "TISR3PI_SIG_PEAK",
              const TString sample_type = "Signal") {

  cout << "\n========================================" << endl;
  cout << "  DOUBLE GAUSSIAN CORE RESOLUTION" << endl;
  cout << "  Sample: " << sample_type << endl;
  cout << "  Tree:   " << tree_type << endl;
  cout << "  Using uncorrected fitted energies (Br_e1_fit, etc.)" << endl;
  cout << "========================================\n" << endl;

  gROOT->GetListOfCanvases()->Delete();
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();
  TH1::AddDirectory(kFALSE);

  gSystem->Exec("mkdir -p ../pull_scan");

  // Input file
  const TString input_file_nm = "/home/bo/Desktop/bdt_tuning_TDATA_norm/cut/tree_pre.root";
  
  TFile* tree_file = new TFile(input_file_nm);
  if (!tree_file || tree_file->IsZombie()) {
    cerr << "ERROR: Cannot open " << input_file_nm << endl;
    return;
  }

  TTree *SIG_TREE = (TTree*)tree_file->Get(tree_type);
  if (!SIG_TREE) {
    cerr << "ERROR: Cannot find " << tree_type << endl;
    tree_file->Close();
    return;
  }
  cout << "✓ Loaded " << tree_type << " with " << SIG_TREE->GetEntries() << " entries" << endl;

  // Histograms
  TH1D *h_phoE = new TH1D("hist_phoE_sig", "", 300, 0., 350.);
  TH1D *h_diff = new TH1D("hist_phoE_diff_sig", "", 400, -50., 50.);

  double e1=0, e2=0, e3=0;
  double e1_true=0, e2_true=0, e3_true=0;

  // Always use uncorrected fitted energies
  SIG_TREE->SetBranchAddress("Br_e1_fit", &e1);
  SIG_TREE->SetBranchAddress("Br_e2_fit", &e2);
  SIG_TREE->SetBranchAddress("Br_e3_fit", &e3);
  SIG_TREE->SetBranchAddress("Br_e1_bdt_true", &e1_true);
  SIG_TREE->SetBranchAddress("Br_e2_bdt_true", &e2_true);
  SIG_TREE->SetBranchAddress("Br_e3_bdt_true", &e3_true);

  Long64_t nentries = SIG_TREE->GetEntries();
  for (Long64_t irow = 0; irow < nentries; ++irow) {
    SIG_TREE->GetEntry(irow);
    // Only use π⁰ photons (e1 and e2) for resolution
    h_diff->Fill(e1 - e1_true);
    h_diff->Fill(e2 - e2_true);
    // Optionally include e3 (ISR) – currently commented out
    //h_diff->Fill(e3 - e3_true);
    h_phoE->Fill(e1);
    h_phoE->Fill(e2);
    //h_phoE->Fill(e3);
  }

  cout << "Histograms filled with " << h_diff->GetEntries() << " entries" << endl;

  // Normalize
  double integral = h_diff->Integral();
  if (integral <= 0) {
    cerr << "ERROR: pull histogram has zero entries." << endl;
    tree_file->Close();
    return;
  }
  h_diff->Scale(1.0 / integral);

  // Estimate RMS and mean
  double rms = h_diff->GetRMS();
  double mean = 0.;
  double peak = h_diff->GetMaximum();

  // Fit range: restrict to core (±1.5 sigma) to avoid tails
  const double fit_factor = 1.5;   // adjust if needed (1.0–2.0)
  double fit_min = mean - fit_factor * rms;
  double fit_max = mean + fit_factor * rms;
  if (fit_min < h_diff->GetXaxis()->GetXmin()) fit_min = h_diff->GetXaxis()->GetXmin();
  if (fit_max > h_diff->GetXaxis()->GetXmax()) fit_max = h_diff->GetXaxis()->GetXmax();

  cout << "\nFit range (core): [" << fit_min << ", " << fit_max << "] MeV" << endl;

  // ----------------------------------------------------------------
  // Compute better initial guesses from histogram properties
  // ----------------------------------------------------------------
  // Find FWHM
  int bin_peak = h_diff->GetMaximumBin();
  double x_peak = h_diff->GetBinCenter(bin_peak);
  double half_max = peak / 2.0;

  int bin_left = bin_peak, bin_right = bin_peak;
  while (bin_left > 1 && h_diff->GetBinContent(bin_left) > half_max) bin_left--;
  while (bin_right < h_diff->GetNbinsX() && h_diff->GetBinContent(bin_right) > half_max) bin_right++;

  double fwhm = h_diff->GetBinCenter(bin_right) - h_diff->GetBinCenter(bin_left);
  double sigma_guess_core = fwhm / 2.3548; // Gaussian FWHM = 2.3548 * sigma
  if (sigma_guess_core <= 0) sigma_guess_core = rms * 0.5; // fallback

  double sigma_guess_tail = rms * 1.2; // tail ~ RMS

  // ----------------------------------------------------------------
  // Double Gaussian fit (always used, no fallback to single Gaussian)
  // ----------------------------------------------------------------
  TF1 *doubleGaus = new TF1("doubleGaus", "gaus(0)+gaus(3)", fit_min, fit_max);
  
  // Set initial parameters
  doubleGaus->SetParameters(peak * 0.8, mean, sigma_guess_core,
                            peak * 0.2, mean, sigma_guess_tail);
  
  // Set parameter limits
  doubleGaus->SetParLimits(0, 0, peak * 2);
  doubleGaus->SetParLimits(1, mean - 1.0, mean + 1.0);
  doubleGaus->SetParLimits(2, 0.1, 5.0);
  doubleGaus->SetParLimits(3, 0, peak * 2);
  doubleGaus->SetParLimits(4, mean - 2.0, mean + 2.0);
  doubleGaus->SetParLimits(5, 0.5, 10.0);

  Int_t status = h_diff->Fit(doubleGaus, "RQS");
  double chi2 = doubleGaus->GetChisquare();
  int ndf = doubleGaus->GetNDF();
  double chi2ndf = (ndf > 0) ? chi2 / ndf : 0;

  // Extract both sigmas
  double s1 = doubleGaus->GetParameter(2);
  double s2 = doubleGaus->GetParameter(5);
  double e1_err = doubleGaus->GetParError(2);
  double e2_err = doubleGaus->GetParError(5);

  // The core is the narrower Gaussian
  bool useFirst = (s1 < s2);
  double core_sigma = useFirst ? s1 : s2;
  double core_sigma_err = useFirst ? e1_err : e2_err;

  // Warn if fit quality is poor, but keep the result
  if (status != 0 || chi2ndf > 3.0) {
    cout << "\n⚠️  Double Gaussian fit status = " << status 
         << ", χ²/NDF = " << chi2ndf << " (fit may be poor)" << endl;
  }

  // Always use the double Gaussian result
  double final_sigma = core_sigma;
  double final_sigma_err = core_sigma_err;
  double final_chi2ndf = chi2ndf;
  TF1 *finalFit = doubleGaus;

  cout << "\nDouble Gaussian fit results:" << endl;
  cout << "  sigma1 = " << s1 << " ± " << e1_err << " MeV" << endl;
  cout << "  sigma2 = " << s2 << " ± " << e2_err << " MeV" << endl;
  cout << "  Core sigma (narrower) = " << final_sigma << " ± " << final_sigma_err << " MeV" << endl;
  cout << "  χ²/NDF = " << final_chi2ndf << endl;

  // Bin width recommendation
  double binwidth = 2.5 * final_sigma;   // fine scan; use 3*sigma for coarser

  cout << "\n========================================" << endl;
  cout << "  BIN WIDTH RECOMMENDATION (double Gaus core)" << endl;
  cout << "========================================" << endl;
  cout << "  core sigma = " << final_sigma << " ± " << final_sigma_err << " MeV" << endl;
  cout << "  binwidth   = " << binwidth << " MeV (2.5 × sigma)" << endl;
  cout << "========================================" << endl;

  // ----------------------------------------------------------------------
  // Draw plots
  // ----------------------------------------------------------------------
  TCanvas *c1 = new TCanvas("c1", "Photon energy spectrum " + sample_type, 900, 700);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.15);
  h_phoE->SetLineWidth(2);
  h_phoE->SetLineColor(4);
  h_phoE->GetXaxis()->SetTitle("Photon energy E_{rec} (MeV)");
  h_phoE->GetYaxis()->SetTitle("Entries");
  h_phoE->GetXaxis()->CenterTitle();
  h_phoE->GetYaxis()->CenterTitle();
  h_phoE->GetYaxis()->SetRangeUser(0.01, 1.6 * h_phoE->GetMaximum());
  h_phoE->Draw("hist");

  TLegend *leg1 = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->AddEntry(h_phoE, sample_type, "l");
  leg1->Draw();

  // Pull distribution with fit
  TCanvas *c2 = new TCanvas("c2", "Core pull distribution", 900, 700);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.15);
  h_diff->SetLineWidth(2);
  h_diff->SetLineColor(4);
  h_diff->GetXaxis()->SetTitle("E_{fit} - E_{true} (MeV)");
  h_diff->GetYaxis()->SetTitle("Normalized Entries");
  h_diff->GetXaxis()->CenterTitle();
  h_diff->GetYaxis()->CenterTitle();
  h_diff->GetYaxis()->SetRangeUser(0., 1.6 * h_diff->GetMaximum());
  h_diff->Draw("hist");

  // Draw the double Gaussian fit
  finalFit->SetLineColor(kRed);
  finalFit->SetLineWidth(2);
  finalFit->Draw("same");

  // Draw vertical lines for fit range
  TLine *lmin = new TLine(fit_min, 0, fit_min, h_diff->GetMaximum());
  lmin->SetLineColor(kGreen);
  lmin->SetLineStyle(2);
  lmin->Draw();
  TLine *lmax = new TLine(fit_max, 0, fit_max, h_diff->GetMaximum());
  lmax->SetLineColor(kGreen);
  lmax->SetLineStyle(2);
  lmax->Draw();

  TPaveText *pt1 = new TPaveText(0.2, 0.75, 0.5, 0.85, "NDC");
  pt1->SetFillColor(0);
  pt1->SetBorderSize(0);
  pt1->SetTextAlign(12);
  pt1->SetTextSize(0.05);
  pt1->AddText(Form("#sigma = %.2f MeV", final_sigma));
  pt1->AddText(Form("#chi^{2}/NDF = %.2f", final_chi2ndf));
  pt1->Draw();

  TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->AddEntry(h_diff, "Data", "l");
  leg2->AddEntry(finalFit, "Double Gaussian (core)", "l");
  leg2->Draw("same");

  // Save histograms
  TFile *fout = new TFile("binwidth_output.root", "RECREATE");
  h_diff->Write();
  finalFit->Write();
  h_phoE->Write();
  fout->Close();

  tree_file->Close();

  cout << "\n✅ Histograms saved to binwidth_output.root" << endl;
}
