// ============================================================================
// track_scan.C
//
// Fits residual distributions for invaraint mass of pi+p- system in mass bins.
// Options for fit_type:
//   "pull"     : ppIM_rec - ppIM_true)   
// Options for fit_model:
//   "gausPoly"    : Gaussian + linear polynomial (default)
//   "doubleGaus"  : Double Gaussian (core + tail)
//   "crystalBall" : Crystal Ball (asymmetric tail)
//   "gausCheb2"   : Gaussian + 2nd-order polynomial (Chebyshev-like)
// Usage:
//   .x track_scan.C("TISR3PI_SIG_PEAK", "Signal", "pull", true, "gausCheb2")
// PDF output: side‑by‑side (16 histograms per page) in a 4x4 grid.
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
#include <TList.h>
#include <TKey.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <iostream>
#include <vector>

using namespace std;

// ----------------------------------------------------------------------------
// Fit a single histogram with the chosen model
// ----------------------------------------------------------------------------
bool fitHist(TH1D *h, double &mean, double &sigma,
             double &mean_err, double &sigma_err,
             double &core_amp, double &fit_min, double &fit_max,
             double &chi2ndf, const TString &model = "gausPoly") {
  if (!h || h->GetEntries() < 500) return false;

  double xmin = h->GetXaxis()->GetXmin();
  double xmax = h->GetXaxis()->GetXmax();
  double rms = h->GetRMS();
  double peak = h->GetMaximum();
  double mean0 = h->GetMean();

  const double fit_factor = 1.5;
  fit_min = mean0 - fit_factor * rms;
  fit_max = mean0 + fit_factor * rms;
  if (fit_min < xmin) fit_min = xmin;
  if (fit_max > xmax) fit_max = xmax;

  int bin_peak = h->GetMaximumBin();
  double half_max = peak / 2.0;
  int bin_left = bin_peak, bin_right = bin_peak;
  while (bin_left > 1 && h->GetBinContent(bin_left) > half_max) bin_left--;
  while (bin_right < h->GetNbinsX() && h->GetBinContent(bin_right) > half_max) bin_right++;
  double fwhm = h->GetBinCenter(bin_right) - h->GetBinCenter(bin_left);
  double sigma_guess_core = fwhm / 2.3548;
  if (sigma_guess_core <= 0 || sigma_guess_core > rms) sigma_guess_core = rms * 0.7;

  TF1 *fit = nullptr;
  bool ok = false;

  // ------------------------------------------------
  // Choose model
  // ------------------------------------------------
  if (model == "gausCheb2") {
    // Gaussian + 2nd-order polynomial (equivalent to Chebyshev of order 2)
    fit = new TF1("gausCheb2", "gaus(0)+pol2(3)", fit_min, fit_max);
    // Parameters: [0]=amp, [1]=mean, [2]=sigma, [3]=p0, [4]=p1, [5]=p2
    double amp_g = peak * 0.9;
    fit->SetParameters(amp_g, mean0, sigma_guess_core, 0.0, 0.0, 0.0);
    fit->SetParLimits(0, 0, peak * 3);
    fit->SetParLimits(1, mean0 - 1.0, mean0 + 1.0);
    fit->SetParLimits(2, 0.005, 5);
    // Polynomial coefficients: allow modest values
    fit->SetParLimits(3, -peak*0.5, peak*0.5);
    fit->SetParLimits(4, -peak*0.1, peak*0.1);
    fit->SetParLimits(5, -peak*0.01, peak*0.01);

    Int_t status = h->Fit(fit, "RQS");
    double chi2 = fit->GetChisquare();
    int ndf = fit->GetNDF();
    chi2ndf = (ndf > 0) ? chi2 / ndf : 0;

    if (status == 0 && chi2ndf < 5.0) {
      mean = fit->GetParameter(1);
      sigma = fit->GetParameter(2);
      mean_err = fit->GetParError(1);
      sigma_err = fit->GetParError(2);
      core_amp = fit->GetParameter(0);
      ok = (sigma > 0.005 && sigma_err > 0);
    }
    delete fit;
    if (ok) return true;
  }
  else if (model == "crystalBall") {
    fit = new TF1("crystalBall", "crystalball", fit_min, fit_max);
    fit->SetParameters(peak, mean0, sigma_guess_core, 1.5, 2.0);
    fit->SetParLimits(0, 0, peak * 3);
    fit->SetParLimits(1, mean0 - 1.0, mean0 + 1.0);
    fit->SetParLimits(2, 0.005, 5);
    fit->SetParLimits(3, 0.1, 10);
    fit->SetParLimits(4, 1.0, 10);

    Int_t status = h->Fit(fit, "RQS");
    double chi2 = fit->GetChisquare();
    int ndf = fit->GetNDF();
    chi2ndf = (ndf > 0) ? chi2 / ndf : 0;

    if (status == 0 && chi2ndf < 5.0) {
      mean = fit->GetParameter(1);
      sigma = fit->GetParameter(2);
      mean_err = fit->GetParError(1);
      sigma_err = fit->GetParError(2);
      core_amp = fit->GetParameter(0);
      ok = (sigma > 0.005 && sigma_err > 0);
    }
    delete fit;
    if (ok) return true;
  }
  else if (model == "doubleGaus") {
    fit = new TF1("doubleGaus", "gaus(0)+gaus(3)", fit_min, fit_max);
    double amp1 = peak * 0.8;
    double amp2 = peak * 0.2;
    fit->SetParameters(amp1, mean0, sigma_guess_core,
                       amp2, mean0, sigma_guess_core * 1.5);
    fit->SetParLimits(0, 0, peak * 3);
    fit->SetParLimits(1, mean0 - 1.0, mean0 + 1.0);
    fit->SetParLimits(2, 0.005, 5);
    fit->SetParLimits(3, 0, peak);
    fit->SetParLimits(4, mean0 - 2.0, mean0 + 2.0);
    fit->SetParLimits(5, 0.005, 10);

    Int_t status = h->Fit(fit, "RQS");
    double chi2 = fit->GetChisquare();
    int ndf = fit->GetNDF();
    chi2ndf = (ndf > 0) ? chi2 / ndf : 0;

    if (status == 0 && chi2ndf < 5.0) {
      mean = fit->GetParameter(1);
      sigma = fit->GetParameter(2);
      mean_err = fit->GetParError(1);
      sigma_err = fit->GetParError(2);
      core_amp = fit->GetParameter(0);
      ok = (sigma > 0.005 && sigma_err > 0);
    }
    delete fit;
    if (ok) return true;
  }
  else { // "gausPoly" (default)
    fit = new TF1("gausPoly", "gaus(0)+pol1(3)", fit_min, fit_max);
    double amp_g = peak * 0.9;
    fit->SetParameters(amp_g, mean0, sigma_guess_core, 0.0, 0.0);
    fit->SetParLimits(0, 0, peak * 3);
    fit->SetParLimits(1, mean0 - 1.0, mean0 + 1.0);
    fit->SetParLimits(2, 0.005, 5);
    fit->SetParLimits(3, -peak, peak);
    fit->SetParLimits(4, -peak, peak);

    Int_t status = h->Fit(fit, "RQS");
    double chi2 = fit->GetChisquare();
    int ndf = fit->GetNDF();
    chi2ndf = (ndf > 0) ? chi2 / ndf : 0;

    if (status == 0 && chi2ndf < 5.0) {
      mean = fit->GetParameter(1);
      sigma = fit->GetParameter(2);
      mean_err = fit->GetParError(1);
      sigma_err = fit->GetParError(2);
      core_amp = fit->GetParameter(0);
      ok = (sigma > 0.005 && sigma_err > 0);
    }
    delete fit;
    if (ok) return true;
  }

  // Fallback: single Gaussian
  TF1 *sgaus = new TF1("singleGaus", "gaus", fit_min, fit_max);
  sgaus->SetParameters(peak, mean0, sigma_guess_core);
  sgaus->SetParLimits(0, 0, peak * 3);
  sgaus->SetParLimits(1, mean0 - 1.0, mean0 + 1.0);
  sgaus->SetParLimits(2, 0.005, 5);

  Int_t status = h->Fit(sgaus, "RQS");
  double chi2 = sgaus->GetChisquare();
  int ndf = sgaus->GetNDF();
  chi2ndf = (ndf > 0) ? chi2 / ndf : 0;

  if (status == 0 && chi2ndf < 5.0) {
    mean = sgaus->GetParameter(1);
    sigma = sgaus->GetParameter(2);
    mean_err = sgaus->GetParError(1);
    sigma_err = sgaus->GetParError(2);
    core_amp = sgaus->GetParameter(0);
    ok = (sigma > 0.005 && sigma_err > 0);
  }
  delete sgaus;
  return ok;
}

// ----------------------------------------------------------------------------
// Draw a single bin histogram into a given TPad (for 4x4 grid)
// ----------------------------------------------------------------------------
void drawBinHistInPad(TH1D *h, int bin, double mean, double sigma,
                      double core_amp, double fit_min, double fit_max,
                      double mass_min, double mass_max,
                      double chi2ndf, const TString &fit_type, const TString &sample_type,
                      TPad *pad) {
  if (!h || !pad) return;

  pad->cd();
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.02);

  h->SetLineWidth(1);
  h->SetLineColor(kBlue);

  h->GetXaxis()->SetNdivisions(505);
  h->GetXaxis()->SetTitle(Form("E_{#gamma} %s", fit_type.Data()));
  //h->GetYaxis()->SetTitle("Entries");
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelSize(0.07);
  h->GetYaxis()->SetLabelSize(0.07);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetRangeUser(0.01, 1.4 * h->GetMaximum());
  h->GetYaxis()->SetNdivisions(505);
  
  //cout << sample_type << endl;
  h->Draw("E0");
  
  TF1 *core = new TF1("core", "gaus", fit_min, fit_max);
  core->SetParameters(core_amp, mean, sigma);
  core->SetLineColor(kRed);
  core->SetLineWidth(1);
  //core->Draw("same");

  TPaveText *pt = new TPaveText(0.65, 0.6, 0.95, 0.9, "NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.05);
  pt->AddText(Form("Bias = %.3f", mean));
  pt->AddText(Form("#sigma = %.3f", sigma));
  pt->AddText(Form("#chi^{2}/NDF = %.2f", chi2ndf));
  pt->AddText(Form("Entries = %d", (int)h->GetEntries()));
  pt->Draw();

  TPaveText *ptBin = new TPaveText(0.15, 0.9, 0.60, 0.95, "NDC");
  ptBin->SetFillColor(0);
  ptBin->SetBorderSize(0);
  ptBin->SetTextAlign(12);
  ptBin->SetTextSize(0.07);
  ptBin->AddText(Form("Bin %d: [%.0f-%.0f] MeV", bin, mass_min, mass_max));
  ptBin->Draw();

  TLegend *leg = new TLegend(0.15, 0.60, 0.5, 0.8);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(h, sample_type, "lep");
  leg->AddEntry(core, "Core Gaussian", "l");
  leg->Draw();
}

// ----------------------------------------------------------------------------
// Main macro
// ----------------------------------------------------------------------------
void trackmass_scan(const TString tree_type = "TDATA",
               const TString sample_type = "Data",
               const TString fit_type = "pull",
               bool draw_bins = true,
               const TString input_file_nm = "/home/bo/Desktop/bdt_tuning_TDATA_norm_true_temp/cut/tree_pre.root",
               const TString fit_model = "gausPoly",
	       const TString pull_type = "new")
{
  TString pdf_name = Form("../trackmass_scan/bin_histograms_%s_%s.pdf", tree_type.Data(), pull_type.Data());
  
  TString root_name = "";
  if (pull_type == "new") {
    root_name = Form("../trackmass_scan/pull_scan_%s_new.root", tree_type.Data());
  }
  else {
    root_name = Form("../trackmass_scan/pull_scan_%s.root", tree_type.Data());
  }
  
  gSystem->Exec("mkdir -p ../trackmass_scan");

  cout << "\n========================================" << endl;
  cout << "  PULL AND ppIM SCAN (invariant mass of pi+ pi-)" << endl;
  cout << "  Sample: " << sample_type << endl;
  cout << "  Tree:   " << tree_type << endl;
  cout << "  Fit type: " << fit_type << endl;
  cout << "  Fit model: " << fit_model << endl;
  cout << "========================================\n" << endl;

  gROOT->GetListOfCanvases()->Delete();
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(6);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();
  TH1::AddDirectory(kFALSE);

  // Input file
  TFile* tree_file = new TFile(input_file_nm);
  if (!tree_file || tree_file->IsZombie()) {
    cerr << "ERROR: Cannot open " << input_file_nm << endl;
    return;
  }

  TTree *INPUT_TREE = (TTree*)tree_file->Get(tree_type);
  if (!INPUT_TREE) {
    cerr << "ERROR: Cannot find " << tree_type << endl;
    tree_file->Close();
    return;
  }
  cout << "✓ Loaded " << tree_type << " with " << INPUT_TREE->GetEntries() << " entries" << endl;

  // Output file
  TFile *fout = new TFile(root_name, "RECREATE");
  fout->cd();


}
