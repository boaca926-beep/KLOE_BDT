// ============================================================================
// trackmass_scan.C
//
// Fits residual distributions for invariant mass of pi+p- system in mass bins.
// Options:
//   - MC signal: fits (ppIM_rec - ppIM_true) -> gives bias & resolution
//   - Data:      uses histogram mean & RMS of ppIM_rec -> gives mean mass & RMS
// Usage:
//   .x trackmass_scan.C("TISR3PI_SIG_PEAK", "Signal", "pull", true, "gausPoly")
//   .x trackmass_scan.C("TDATA", "Data", "pull", true, "gausPoly")
// If sample_type == "Data" (case-insensitive), data mode is forced.
// Options for fit_model:
//   "gausPoly"    : Gaussian + linear polynomial (default)
//   "doubleGaus"  : Double Gaussian (core + tail)
//   "crystalBall" : Crystal Ball (asymmetric tail)
//   "gausCheb2"   : Gaussian + 2nd-order polynomial (Chebyshev-like)

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
// Fit a single histogram with the chosen model (only for MC)
// ----------------------------------------------------------------------------
bool fitHist(TH1D *h, double &mean, double &sigma,
             double &mean_err, double &sigma_err,
             double &core_amp, double &fit_min, double &fit_max,
             double &chi2ndf, const TString &model = "crystalBall") {
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
    fit = new TF1("gausCheb2", "gaus(0)+pol2(3)", fit_min, fit_max);
    double amp_g = peak * 0.9;
    fit->SetParameters(amp_g, mean0, sigma_guess_core, 0.0, 0.0, 0.0);
    fit->SetParLimits(0, 0, peak * 3);
    fit->SetParLimits(1, mean0 - 1.0, mean0 + 1.0);
    fit->SetParLimits(2, 0.005, 5);
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
// Draw a single bin histogram into a given TPad
// ----------------------------------------------------------------------------
void drawBinHistInPad(TH1D *h, int bin, double mean, double sigma,
                      double core_amp, double fit_min, double fit_max,
                      double mass_min, double mass_max,
                      double chi2ndf, const TString &fit_type, const TString &sample_type,
                      TPad *pad, bool isMC) {
  if (!h || !pad) return;

  pad->cd();
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.02);

  h->SetLineWidth(1);
  h->SetLineColor(kBlue);

  h->GetXaxis()->SetNdivisions(505);
  h->GetXaxis()->SetTitle(isMC ? "M^{Data}_{#pi#pi} - M^{true}_{#pi#pi} (MeV/c^{2})" : "M_{#pi#pi} (MeV/c^{2})");
  h->GetXaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetXaxis()->SetLabelSize(0.07);
  h->GetYaxis()->SetLabelSize(0.07);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetRangeUser(0.01, 1.4 * h->GetMaximum());
  h->GetYaxis()->SetNdivisions(505);

  h->Draw("E0");

  TF1 *core = new TF1("core", "gaus", fit_min, fit_max);
  core->SetParameters(core_amp, mean, sigma);
  core->SetLineColor(kRed);
  core->SetLineWidth(1);

  TPaveText *pt = new TPaveText(0.67, 0.55, 0.95, 0.85, "NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.05);
  if (isMC) {
    pt->AddText(Form("Bias = %.3f", mean));
    pt->AddText(Form("Sigma = %.3f", sigma));
  } else {
    pt->AddText(Form("Mean = %.3f", mean));
    pt->AddText(Form("RMS = %.3f", sigma));
  }
  pt->AddText(Form("#chi^{2}/NDF = %.2f", chi2ndf));
  pt->AddText(Form("Entries = %d", (int)h->GetEntries()));
  pt->Draw();

  TPaveText *ptBin = new TPaveText(0.15, 0.9, 0.60, 0.95, "NDC");
  ptBin->SetFillColor(0);
  ptBin->SetBorderSize(0);
  ptBin->SetTextAlign(12);
  ptBin->SetTextSize(0.07);
  ptBin->AddText(Form("Bin %d: [%.0f-%.0f] MeV/c^{2}", bin, mass_min, mass_max));
  ptBin->Draw();

  TLegend *leg = new TLegend(0.15, 0.60, 0.5, 0.8);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(h, sample_type, "lep");
  if (isMC) leg->AddEntry(core, "Core Gaussian", "l");
  leg->Draw();
}

// ----------------------------------------------------------------------------
// Main macro
// ----------------------------------------------------------------------------
void trackmass_scan(const TString tree_type = "TISR3PI_SIG_PEAK",
                    const TString sample_type = "Signal",
                    const TString fit_type = "pull",
                    bool draw_bins = true,
                    //const TString input_file_nm = "/home/bo/Desktop/bdt_tuning_TDATA_chain_false/cut/tree_pre.root",
		    const TString input_file_nm = "/home/bo/Desktop/bdt_tuning_TDATA_norm_false/cut/tree_pre.root",
                    const TString fit_model = "gausPoly", //gausCheb2, gausPoly 
                    const TString pull_type = "old")
{

/*
void trackmass_scan(const TString tree_type = "TDATA",
                    const TString sample_type = "Data",
                    const TString fit_type = "pull",
                    bool draw_bins = true,
                    const TString input_file_nm = "/home/bo/Desktop/bdt_tuning_TDATA_chain_false/cut/tree_pre.root",
                    const TString fit_model = "gausPoly",
                    const TString pull_type = "")
{
*/

  TString pdf_name = Form("../trackmass_scan/bin_histograms_%s_%s.pdf", tree_type.Data(), pull_type.Data());

  TString root_name = "";
  if (pull_type == "new") {
    root_name = Form("../trackmass_scan/pull_scan_%s_new.root", tree_type.Data());
  } else {
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

  // ----------------------------------------------------------------------
  // Determine if we are in MC or data mode
  // ----------------------------------------------------------------------
  bool isMC = true;
  // If sample_type contains "Data" (case-insensitive), force data mode
  if (sample_type.Contains("Data", TString::kIgnoreCase)) {
    isMC = false;
    cout << "✓ Sample type is Data: using histogram mean & RMS of ppIM distribution" << endl;
  } else {
    // Otherwise, auto-detect based on branch existence
    if (INPUT_TREE->GetBranch("Br_ppIM_true") != nullptr) {
      isMC = true;
      cout << "✓ MC tree detected: will fit pull distribution (ppIM - ppIM_true)" << endl;
    } else {
      isMC = false;
      cout << "✓ Data tree (no truth branch): using histogram mean & RMS of ppIM distribution" << endl;
    }
  }

  // Output file
  TFile *fout = new TFile(root_name, "RECREATE");
  fout->cd();

  // Merged TTree (stores the quantity per event)
  TTree *ppIMTree = new TTree("ppIMTree", "merged branches");
  double ppIM_quantity = 0.;
  ppIMTree->Branch("ppIM_quantity", &ppIM_quantity, "ppIM_quantity/D");

  TH1D *h_overall = new TH1D("h_overall", "", 200, isMC ? -10. : 250., isMC ? 10. : 650.);
  if (isMC) {
    h_overall->SetTitle("Pull distribution (ppIM - ppIM_true)");
  } else {
    h_overall->SetTitle("M_{#pi#pi} distribution (data)");
  }

  // Binning in reconstructed mass
  const double mass_min = 260.;
  const double mass_max = 650.;
  const double binwidth = 2 * 5.6;   // 11.2 MeV/c²
  const int nbins = (int)((mass_max - mass_min) / binwidth);
  cout << "Binning in M_{#pi#pi}: [" << mass_min << ", " << mass_max << "] MeV/c^{2}, "
       << "bins: " << nbins << ", width: " << binwidth << " MeV/c^{2}" << endl;

  vector<TH1D*> hist_hists(nbins);
  vector<double> bin_center(nbins);
  vector<double> bin_mean(nbins), bin_sigma(nbins);
  vector<double> bin_mean_err(nbins), bin_sigma_err(nbins);
  vector<double> bin_core_amp(nbins);
  vector<double> bin_fit_min(nbins), bin_fit_max(nbins);
  vector<double> bin_chi2ndf(nbins);
  vector<int> bin_entries(nbins);

  // Set histogram range appropriately
  double hmin = isMC ? -5.0 : 250.0;
  double hmax = isMC ?  5.0 : 650.0;

  for (int b = 0; b < nbins; ++b) {
    double lo = mass_min + b * binwidth;
    double hi = lo + binwidth;
    bin_center[b] = (lo + hi) / 2.0;
    TString name = Form("hist_bin_%d", b);
    hist_hists[b] = new TH1D(name, "", 100, hmin, hmax);
    hist_hists[b]->Sumw2();
    bin_mean[b] = 0; bin_sigma[b] = 0;
    bin_mean_err[b] = 0; bin_sigma_err[b] = 0;
    bin_core_amp[b] = 0;
    bin_fit_min[b] = 0; bin_fit_max[b] = 0;
    bin_chi2ndf[b] = 0;
    bin_entries[b] = 0;
  }

  // Branch addresses
  double ppIM = 0., ppIM_true = 0.;
  INPUT_TREE->SetBranchAddress("Br_ppIM", &ppIM);
  if (isMC) {
    INPUT_TREE->SetBranchAddress("Br_ppIM_true", &ppIM_true);
  }

  Long64_t nentries = INPUT_TREE->GetEntries();
  for (Long64_t irow = 0; irow < nentries; ++irow) {
    INPUT_TREE->GetEntry(irow);

    // Determine quantity to histogram
    double value;
    if (isMC) {
      value = ppIM - ppIM_true;   // pull
    } else {
      value = ppIM;               // reconstructed mass
    }

    h_overall->Fill(value);

    // Fill per‑bin histograms based on ppIM (reconstructed mass)
    if (ppIM >= mass_min && ppIM < mass_max) {
      int b = (int)((ppIM - mass_min) / binwidth);
      if (b >= 0 && b < nbins) {
        hist_hists[b]->Fill(value);
        bin_entries[b]++;
      }
    }

    ppIM_quantity = value;
    ppIMTree->Fill();
  }

  // Per‑bin analysis
  cout << "\n=== Per‑bin results (" << (isMC ? "fit" : "mean & RMS") << ") ===" << endl;
  if (isMC) {
    cout << "Bin  M0 (MeV/c²) M1 (MeV/c²)  Entries  Mean ± err  Sigma ± err  χ²/NDF" << endl;
  } else {
    cout << "Bin  M0 (MeV/c²) M1 (MeV/c²)  Entries  Mean ± err  RMS ± err" << endl;
  }

  vector<double> M0_LIST(nbins), M1_LIST(nbins);

  for (int b = 0; b < nbins; ++b) {
    bin_entries[b] = hist_hists[b]->GetEntries();
    bool ok = false;

    if (isMC) {
      // Fit for MC (pull distribution)
      ok = fitHist(hist_hists[b], bin_mean[b], bin_sigma[b],
                   bin_mean_err[b], bin_sigma_err[b],
                   bin_core_amp[b], bin_fit_min[b], bin_fit_max[b],
                   bin_chi2ndf[b], fit_model);
      if (ok) {
        M0_LIST[b] = mass_min + b * binwidth;
        M1_LIST[b] = mass_min + (b+1) * binwidth;
        printf("%3d  %6.1f %6.1f  %6d  %7.3f±%-7.3f  %7.3f±%-7.3f  %7.3f\n",
               b, M0_LIST[b], M1_LIST[b],
               bin_entries[b], bin_mean[b], bin_mean_err[b],
               bin_sigma[b], bin_sigma_err[b], bin_chi2ndf[b]);
      }
    } else {
      // Data: use histogram statistics
      if (bin_entries[b] >= 500) {
        bin_mean[b] = hist_hists[b]->GetMean();
        bin_sigma[b] = hist_hists[b]->GetRMS();
        bin_mean_err[b] = hist_hists[b]->GetMeanError(); // standard error of mean
        // Approximate RMS error (for Gaussian): RMS / sqrt(2*N)
        bin_sigma_err[b] = (bin_entries[b] > 0) ? bin_sigma[b] / sqrt(2.0 * bin_entries[b]) : 0;
        bin_chi2ndf[b] = 0; // not applicable
        ok = true;
        M0_LIST[b] = mass_min + b * binwidth;
        M1_LIST[b] = mass_min + (b+1) * binwidth;
        printf("%3d  %6.1f %6.1f  %6d  %7.3f±%-7.3f  %7.3f±%-7.3f\n",
               b, M0_LIST[b], M1_LIST[b],
               bin_entries[b], bin_mean[b], bin_mean_err[b],
               bin_sigma[b], bin_sigma_err[b]);
      }
    }

    // If not ok, keep zeros and skip from TGraphs later.
  }

  // Build TGraphs
  TGraphErrors *g_bias = new TGraphErrors();
  TGraphErrors *g_sigma = new TGraphErrors();
  g_bias->SetName("g_bias_vs_M");
  g_sigma->SetName("g_resolution_vs_M");

  int point = 0;
  for (int b = 0; b < nbins; ++b) {
    // For MC, we require fit success (ok); for data we require entries>=500.
    bool valid = (isMC) ? (bin_entries[b] >= 500 && bin_sigma[b] > 0.001 && bin_sigma_err[b] > 0)
                         : (bin_entries[b] >= 500 && bin_sigma[b] > 0);
    if (valid) {
      g_bias->SetPoint(point, bin_center[b], bin_mean[b]);
      g_bias->SetPointError(point, 0, bin_mean_err[b]);
      g_sigma->SetPoint(point, bin_center[b], bin_sigma[b]);
      g_sigma->SetPointError(point, 0, bin_sigma_err[b]);
      point++;
    }
  }
  cout << "Added " << point << " points to TGraphs" << endl;

  // ----------------------------------------------------------------------
  // Draw bin histograms in PDF (6x6 grid)
  // ----------------------------------------------------------------------
  if (draw_bins) {
    cout << "\n=== Generating bin histogram PDF (6x6) ===" << endl;
    TCanvas *c_pdf = new TCanvas("c_pdf", "Bin Histograms", 1600, 1200);
    c_pdf->Print(pdf_name + "[");

    vector<int> bins_to_draw;
    for (int b = 0; b < nbins; ++b) {
      if (bin_entries[b] >= 100) bins_to_draw.push_back(b);
    }

    const int nCols = 6;
    const int nRows = 6;
    const int nPerPage = nCols * nRows;

    for (size_t i = 0; i < bins_to_draw.size(); i += nPerPage) {
      c_pdf->Clear();
      c_pdf->Divide(nCols, nRows, 1e-5, 1e-5);

      for (int j = 0; j < nPerPage && (i+j) < bins_to_draw.size(); ++j) {
        int b = bins_to_draw[i+j];
        TPad *pad = (TPad*)c_pdf->GetPad(j+1);
        pad->cd();
        drawBinHistInPad(hist_hists[b], b, bin_mean[b], bin_sigma[b],
                         bin_core_amp[b], bin_fit_min[b], bin_fit_max[b],
                         M0_LIST[b], M1_LIST[b],
                         bin_chi2ndf[b], fit_type, sample_type, pad, isMC);
      }
      c_pdf->Update();
      c_pdf->Print(pdf_name);
    }

    c_pdf->Print(pdf_name + "]");
    delete c_pdf;
    cout << "✅ Bin histograms saved to " << pdf_name << endl;
  }

  // ----------------------------------------------------------------------
  // Summary canvases
  // ----------------------------------------------------------------------
  // Overall distribution
  TCanvas *c2 = new TCanvas("c2", "Distribution", 900, 700);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.15);
  h_overall->SetLineWidth(2);
  h_overall->SetLineColor(4);
  h_overall->GetXaxis()->SetTitle(isMC ? "Pull (MeV)" : "M_{#pi#pi} (MeV)");
  h_overall->GetYaxis()->SetTitle("Entries");
  h_overall->GetXaxis()->CenterTitle();
  h_overall->GetYaxis()->CenterTitle();
  h_overall->GetYaxis()->SetRangeUser(0.01, 1.6 * h_overall->GetMaximum());
  h_overall->Draw("hist");

  // Resolution vs mass
  TCanvas *c3 = new TCanvas("c3", "Resolution vs mass", 900, 700);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.15);
  g_sigma->SetMarkerStyle(20);
  g_sigma->SetMarkerSize(1.2);
  g_sigma->GetXaxis()->SetTitle("M_{#pi#pi} (MeV/c^{2})");
  g_sigma->GetYaxis()->SetTitle(isMC ? "Resolution (sigma) [MeV]" : "RMS [MeV]");
  g_sigma->GetXaxis()->CenterTitle();
  g_sigma->GetYaxis()->CenterTitle();
  g_sigma->Draw("AP");

  // Bias/Mean vs mass
  TCanvas *c4 = new TCanvas("c4", "Mean vs mass", 900, 700);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.15);
  g_bias->SetMarkerStyle(20);
  g_bias->SetMarkerSize(1.2);
  g_bias->GetXaxis()->SetTitle("M_{#pi#pi} (MeV/c^{2})");
  g_bias->GetYaxis()->SetTitle(isMC ? "Bias (Mean residual) [MeV]" : "Mean M_{#pi#pi} [MeV]");
  g_bias->GetXaxis()->CenterTitle();
  g_bias->GetYaxis()->CenterTitle();
  g_bias->Draw("AP");

  // Write output
  fout->cd();
  g_bias->Write();
  g_sigma->Write();
  ppIMTree->Write();
  h_overall->Write();
  fout->Close();

  cout << "\n✅ Output written to " << root_name << endl;
}
