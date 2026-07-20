// ============================================================================
// pull_scan.C
//
// Fits residual distributions for π⁰ photons in energy bins with Gaussian + linear polynomial.
// Options for fit_type:
//   "pull"     : (E_raw - E_fit)/sigma   (simplified pull for π⁰ photons)
//   "rawdiff"  : E_nofit - E_fit  (absolute difference, MeV)
//   "ratio"    : (E_nofit - E_fit)/E_fit
//   "symdiff"  : (E_nofit - E_fit)/(E_nofit + E_fit)
// Usage:
//   .x pull_scan.C("TISR3PI_SIG_PEAK", "Signal", "pull", true)
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
#include <iostream>
#include <vector>

using namespace std;

// ----------------------------------------------------------------------------
// Fit a single histogram with Gaussian + linear polynomial, fallback to single Gaussian
// ----------------------------------------------------------------------------
bool fitHist(TH1D *h, double &mean, double &sigma,
             double &mean_err, double &sigma_err,
             double &core_amp, double &fit_min, double &fit_max,
             double &chi2ndf) {
  if (!h || h->GetEntries() < 500) return false;

  double xmin = h->GetXaxis()->GetXmin();
  double xmax = h->GetXaxis()->GetXmax();
  double rms = h->GetRMS();
  double peak = h->GetMaximum();
  double mean0 = h->GetMean();

  // Fit range: ±3.0 sigma to include tails
  const double fit_factor = 2.;
  fit_min = mean0 - fit_factor * rms;
  fit_max = mean0 + fit_factor * rms;
  if (fit_min < xmin) fit_min = xmin;
  if (fit_max > xmax) fit_max = xmax;

  // Estimate core sigma from FWHM
  int bin_peak = h->GetMaximumBin();
  double half_max = peak / 2.0;
  int bin_left = bin_peak, bin_right = bin_peak;
  while (bin_left > 1 && h->GetBinContent(bin_left) > half_max) bin_left--;
  while (bin_right < h->GetNbinsX() && h->GetBinContent(bin_right) > half_max) bin_right++;
  double fwhm = h->GetBinCenter(bin_right) - h->GetBinCenter(bin_left);
  double sigma_guess_core = fwhm / 2.3548;
  if (sigma_guess_core <= 0 || sigma_guess_core > rms) sigma_guess_core = rms * 0.7;

  // ------------------------------------------------
  // Try Gaussian + linear polynomial
  // ------------------------------------------------
  TF1 *fit = new TF1("gausPoly", "gaus(0)+pol1(3)", fit_min, fit_max);
  // Gaussian: [0]=amp, [1]=mean, [2]=sigma
  // Polynomial: [3]=p0, [4]=p1
  double amp_g = peak * 0.9;
  fit->SetParameters(amp_g, mean0, sigma_guess_core,
                     0.0, 0.0);   // initial polynomial coefficients
  fit->SetParLimits(0, 0, peak * 3);
  fit->SetParLimits(1, mean0 - 1.0, mean0 + 1.0);
  fit->SetParLimits(2, 0.005, 5);      // sigma up to 5 (fraction)
  // Polynomial coefficients – no strict limits (allow flexibility)
  fit->SetParLimits(3, -peak, peak);
  fit->SetParLimits(4, -peak, peak);

  Int_t status = h->Fit(fit, "RQS");
  double chi2 = fit->GetChisquare();
  int ndf = fit->GetNDF();
  chi2ndf = (ndf > 0) ? chi2 / ndf : 0;

  // Extract Gaussian parameters
  bool ok = (status == 0 && chi2ndf < 5.0);
  if (ok) {
    mean = fit->GetParameter(1);
    sigma = fit->GetParameter(2);
    mean_err = fit->GetParError(1);
    sigma_err = fit->GetParError(2);
    core_amp = fit->GetParameter(0);
    if (sigma > 0.005 && sigma_err > 0) {
      delete fit;
      return true;
    }
  }
  delete fit;

  // ------------------------------------------------
  // Fallback: single Gaussian (if polynomial fit fails)
  // ------------------------------------------------
  TF1 *sgaus = new TF1("singleGaus", "gaus", fit_min, fit_max);
  sgaus->SetParameters(peak, mean0, sigma_guess_core);
  sgaus->SetParLimits(0, 0, peak * 3);
  sgaus->SetParLimits(1, mean0 - 1.0, mean0 + 1.0);
  sgaus->SetParLimits(2, 0.005, 5);

  status = h->Fit(sgaus, "RQS");
  chi2 = sgaus->GetChisquare();
  ndf = sgaus->GetNDF();
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
// Draw a single bin histogram with the core Gaussian fit
// ----------------------------------------------------------------------------
void drawBinHist(TH1D *h, int bin, double mean, double sigma,
                 double core_amp, double fit_min, double fit_max,
		 double mass_min, double mass_max,
                 double chi2ndf, const TString &fit_type,
                 TCanvas *c, TPad *pad) {
  if (!h || h->GetEntries() < 500) return;

  pad->cd();
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.15);

  h->SetLineWidth(2);
  h->SetLineColor(4);
  h->GetXaxis()->SetTitle(Form("Residual (%s)", fit_type.Data()));
  h->GetYaxis()->SetTitle("Entries");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetRangeUser(0.01, 1.6 * h->GetMaximum());
  h->Draw("hist");

  // Core Gaussian
  TF1 *core = new TF1("core", "gaus", fit_min, fit_max);
  core->SetParameters(core_amp, mean, sigma);
  core->SetLineColor(kRed);
  core->SetLineWidth(2);
  core->Draw("same");

  // Text box with results
  TPaveText *pt = new TPaveText(0.6, 0.55, 0.85, 0.85, "NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.045);
  pt->AddText(Form("Mean = %.3f", mean));
  pt->AddText(Form("Sigma = %.3f", sigma));
  pt->AddText(Form("#chi^{2}/NDF = %.2f", chi2ndf));
  pt->AddText(Form("Entries = %d", (int)h->GetEntries()));
  pt->Draw();

  // Bin label
  //double lo = h->GetXaxis()->GetXmin();
  //double hi = h->GetXaxis()->GetXmax();
  TPaveText *ptBin = new TPaveText(0.2, 0.75, 0.45, 0.85, "NDC");
  ptBin->SetFillColor(0);
  ptBin->SetBorderSize(0);
  ptBin->SetTextAlign(12);
  ptBin->SetTextSize(0.05);
  ptBin->AddText(Form("Bin %d: %.0f-%.0f MeV", bin, mass_min, mass_max));
  ptBin->Draw();

  // Legend
  TLegend *leg = new TLegend(0.6, 0.4, 0.9, 0.5);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h, "Data", "l");
  leg->AddEntry(core, "Core Gaussian", "l");
  leg->Draw();
}

// .x pull_scan.C("TISR3PI_SIG_PEAK", "Signal", "pull", true)
//
// ----------------------------------------------------------------------------
// Main macro – only π⁰ photons (e1, e2)
// ----------------------------------------------------------------------------
/*
void pull_scan(const TString tree_type = "TISR3PI_SIG_PEAK",
               const TString sample_type = "Signal",
               const TString fit_type = "pull",
               bool draw_bins = true) {
*/

void pull_scan(const TString tree_type = "TDATA",
               const TString sample_type = "Data",
               const TString fit_type = "pull",
               bool draw_bins = true) {


  // Construct output filenames with tree_type
  TString pdf_name = Form("../pull_scan/bin_histograms_%s.pdf", tree_type.Data());
  TString root_name = Form("../pull_scan/pull_scan_%s.root", tree_type.Data());
  gSystem->Exec("mkdir -p ../pull_scan");
  
  cout << "\n========================================" << endl;
  cout << "  PULL SCAN (π⁰ photons only)" << endl;
  cout << "  Sample: " << sample_type << endl;
  cout << "  Tree:   " << tree_type << endl;
  cout << "  Fitting: " << fit_type << endl;
  if (fit_type == "pull")
    cout << "    Pull" << endl;
  else if (fit_type == "rawdiff")
    cout << "    (E_nofit - E_fit) in MeV" << endl;
  else if (fit_type == "ratio")
    cout << "    (E_nofit - E_fit) / E_fit" << endl;
  else if (fit_type == "symdiff")
    cout << "    (E_nofit - E_fit) / (E_nofit + E_fit)" << endl;
  else
    cout << "    Unknown fit_type!" << endl;
  cout << "========================================\n" << endl;

  gROOT->GetListOfCanvases()->Delete();
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();
  TH1::AddDirectory(kFALSE);
  //gSystem->Exec("mkdir -p ../pull_scan");

  // ----------------------------------------------------------------------
  // Input file
  // ----------------------------------------------------------------------
  const TString input_file_nm = "/home/bo/Desktop/bdt_tuning_TDATA_chain/cut/tree_pre.root";
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
  // Output file
  // ----------------------------------------------------------------------
  TFile *fout = new TFile(root_name, "RECREATE");
  fout->cd();

  // Merged TTree – store only two photons
  TTree *EPhoTree = new TTree("EPhoTree", "merged branches");
  double EPho_fit[2], EPho_nofit[2], EPho_value[2];
  EPhoTree->Branch("EPho_fit", EPho_fit, "EPho_fit[2]/D");
  EPhoTree->Branch("EPho_nofit", EPho_nofit, "EPho_nofit[2]/D");
  EPhoTree->Branch("EPho_value", EPho_value, "EPho_value[2]/D");

  // Global histograms
  TH1D *h_phoE      = new TH1D("h_phoE", "", 300, 0., 350.);
  TH1D *h_pull      = new TH1D("h_pull", "", 200, -10, 10);

  // ----------------------------------------------------------------------
  // Binning
  // ----------------------------------------------------------------------
  const double Emin = 0.0;
  const double Emax = 350.0;
  const double binwidth = 2 * 5.6;   // 28 MeV
  const int nbins = (int)((Emax - Emin) / binwidth);
  cout << "Using " << nbins << " bins of width " << binwidth << " MeV" << endl;

  vector<TH1D*> pull_hists(nbins);
  vector<double> bin_center(nbins);
  vector<double> bin_mean(nbins), bin_sigma(nbins);
  vector<double> bin_mean_err(nbins), bin_sigma_err(nbins);
  vector<double> bin_core_amp(nbins);
  vector<double> bin_fit_min(nbins), bin_fit_max(nbins);
  vector<double> bin_chi2ndf(nbins);
  vector<int> bin_entries(nbins);

  for (int b = 0; b < nbins; ++b) {
    double lo = Emin + b * binwidth;
    double hi = lo + binwidth;
    bin_center[b] = (lo + hi) / 2.0;
    TString name = Form("pull_bin_%d", b);
    pull_hists[b] = new TH1D(name, "", 200, -5, 5);
    pull_hists[b]->Sumw2();
    bin_mean[b] = 0; bin_sigma[b] = 0;
    bin_mean_err[b] = 0; bin_sigma_err[b] = 0;
    bin_core_amp[b] = 0;
    bin_fit_min[b] = 0; bin_fit_max[b] = 0;
    bin_chi2ndf[b] = 0;
    bin_entries[b] = 0;
  }

  // ----------------------------------------------------------------------
  // Branch addresses – only π⁰ photons (e1, e2)
  // ----------------------------------------------------------------------
  double e1_fit=0, e2_fit=0;
  double sigma_e1=0, sigma_e2=0;
  double e1_nofit=0, e2_nofit=0;
  double e1_pull=0, e2_pull=0;
  double m_gg_bdt = 0.;
  
  INPUT_TREE->SetBranchAddress("Br_e1_fit", &e1_fit);
  INPUT_TREE->SetBranchAddress("Br_e2_fit", &e2_fit);

  INPUT_TREE->SetBranchAddress("Br_sigma_e1_bdt", &sigma_e1);
  INPUT_TREE->SetBranchAddress("Br_sigma_e2_bdt", &sigma_e2);
  
  INPUT_TREE->SetBranchAddress("Br_e1_nofit_bdt", &e1_nofit);
  INPUT_TREE->SetBranchAddress("Br_e2_nofit_bdt", &e2_nofit);

  INPUT_TREE->SetBranchAddress("Br_e1_pull_bdt", &e1_pull);
  INPUT_TREE->SetBranchAddress("Br_e2_pull_bdt", &e2_pull);

  INPUT_TREE->SetBranchAddress("Br_m_gg_bdt", &m_gg_bdt);

  // ----------------------------------------------------------------------
  // Loop over events
  // ----------------------------------------------------------------------
  const double pi0_true = 134.977;
  
  Long64_t nentries = INPUT_TREE->GetEntries();
  for (Long64_t irow = 0; irow < nentries; ++irow) {
    INPUT_TREE->GetEntry(irow);

    if (m_gg_bdt > pi0_true + 5. || m_gg_bdt < pi0_true - 5.) continue;

    double sigma[2] = {sigma_e1, sigma_e2};
    double fit[2]  = {e1_fit, e2_fit};
    double nofit[2] = {e1_nofit, e2_nofit};
    double pull[2] = {e1_pull, e2_pull};
    
    
    for (int ip = 0; ip < 2; ++ip) {
      EPho_fit[ip]   = fit[ip];
      EPho_nofit[ip] = nofit[ip];
      if (fit_type == "pull") {
        EPho_value[ip] = pull[ip];
      } else if (fit_type == "rawdiff") {
        EPho_value[ip] = nofit[ip] - fit[ip];
      } else if (fit_type == "ratio") {
        EPho_value[ip] = (nofit[ip] - fit[ip]) / fit[ip];
      } else if (fit_type == "symdiff") {
        double denom = nofit[ip] + fit[ip];
        EPho_value[ip] = (denom != 0) ? (nofit[ip] - fit[ip]) / denom : 0;
      } else {
        EPho_value[ip] = 0;
      }
    }
    EPhoTree->Fill();

    // Fill histograms for both photons
    for (int ip = 0; ip < 2; ++ip) {
      double Ef = fit[ip];
      if (Ef <= 0) continue;

      double value = EPho_value[ip];
      h_pull->Fill(value);
      h_phoE->Fill(Ef);

      if (Ef < Emin || Ef >= Emax) continue;
      int b = (int)((Ef - Emin) / binwidth);
      if (b < 0 || b >= nbins) continue;
      pull_hists[b]->Fill(value);
    }
  }

  // ----------------------------------------------------------------------
  // Fit each bin histogram
  // ----------------------------------------------------------------------
  cout << "\n=== Per‑bin Gaussian + Linear Polynomial fit results (" << fit_type << ") ===" << endl;
  cout << "Bin  E0 (MeV) E1 (MeV)  Entries  Mean ± err  Sigma ± err  χ²/NDF" << endl;
  vector<double> E0_LIST(nbins);
  vector<double> E1_LIST(nbins);
 
  for (int b = 0; b < nbins; ++b) {
    bin_entries[b] = pull_hists[b]->GetEntries();
    bool ok = fitHist(pull_hists[b], bin_mean[b], bin_sigma[b],
                      bin_mean_err[b], bin_sigma_err[b],
                      bin_core_amp[b], bin_fit_min[b], bin_fit_max[b],
                      bin_chi2ndf[b]);
    if (ok) {
      E0_LIST[b] = Emin+b*binwidth;
      E1_LIST[b] = Emin+(b+1)*binwidth;
      
      printf("%3d  %6.1f %6.1f  %6d  %7.3f±%-7.3f  %7.3f±%-7.3f  %7.3f\n",
             b, Emin+b*binwidth, Emin+(b+1)*binwidth,
             bin_entries[b], bin_mean[b], bin_mean_err[b],
             bin_sigma[b], bin_sigma_err[b], bin_chi2ndf[b]);
    }
  }

  // ----------------------------------------------------------------------
  // Build TGraphs
  // ----------------------------------------------------------------------
  TGraphErrors *g_bias = new TGraphErrors();
  TGraphErrors *g_sigma = new TGraphErrors();
  g_bias->SetName("g_bias_vs_E");
  g_bias->SetTitle("Bias (mean residual) vs Energy");
  g_sigma->SetName("g_resolution_vs_E");
  g_sigma->SetTitle("Resolution (sigma) vs Energy");

  int point = 0;
  for (int b = 0; b < nbins; ++b) {
    if (bin_entries[b] >= 500 && bin_sigma[b] > 0.001 && bin_sigma_err[b] > 0) {
      g_bias->SetPoint(point, bin_center[b], bin_mean[b]);
      g_bias->SetPointError(point, 0, bin_mean_err[b]);
      g_sigma->SetPoint(point, bin_center[b], bin_sigma[b]);
      g_sigma->SetPointError(point, 0, bin_sigma_err[b]);
      point++;
    }
  }

  // ----------------------------------------------------------------------
  // Draw all bin histograms (if requested)
  // ----------------------------------------------------------------------
  if (draw_bins) {
    cout << "\n=== Generating bin histogram PDF ===" << endl;
    TCanvas *c_pdf = new TCanvas("c_pdf", "Bin Histograms", 900, 700);
    c_pdf->Print(pdf_name + "[");
    for (int b = 0; b < nbins; ++b) {
      if (bin_entries[b] < 500) continue;
      //cout << "mass range " << b + 1 << ": "<< E0_LIST[b] << ", " << E1_LIST[b] << endl;
      drawBinHist(pull_hists[b], b, bin_mean[b], bin_sigma[b],
                  bin_core_amp[b], bin_fit_min[b], bin_fit_max[b],
		  E0_LIST[b], E1_LIST[b],
                  bin_chi2ndf[b], fit_type, c_pdf, c_pdf);
      c_pdf->Print(pdf_name);
    }
    c_pdf->Print(pdf_name + "]");
    delete c_pdf;
    cout << "✅ Bin histograms saved to " << pdf_name << endl;
  }

  // ----------------------------------------------------------------------
  // Draw summary plots
  // ----------------------------------------------------------------------
  TCanvas *c1 = new TCanvas("c1", "Photon energy spectrum", 900, 700);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.15);
  h_phoE->SetLineWidth(2);
  h_phoE->SetLineColor(4);
  h_phoE->GetXaxis()->SetTitle("E_{fit} (MeV)");
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

  TCanvas *c2 = new TCanvas("c2", "Distribution", 900, 700);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.15);
  h_pull->SetLineWidth(2);
  h_pull->SetLineColor(4);
  h_pull->GetXaxis()->SetTitle(Form("Residual (%s)", fit_type.Data()));
  h_pull->GetYaxis()->SetTitle("Entries");
  h_pull->GetXaxis()->CenterTitle();
  h_pull->GetYaxis()->CenterTitle();
  h_pull->GetYaxis()->SetRangeUser(0.01, 1.6 * h_pull->GetMaximum());
  h_pull->Draw("hist");

  TCanvas *c3 = new TCanvas("c3", "Resolution vs energy", 900, 700);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.15);
  g_sigma->SetMarkerStyle(20);
  g_sigma->SetMarkerSize(1.2);
  g_sigma->GetXaxis()->SetTitle("E_{fit} (MeV)");
  g_sigma->GetYaxis()->SetTitle("#sigma (fraction)");
  g_sigma->GetXaxis()->CenterTitle();
  g_sigma->GetYaxis()->CenterTitle();
  g_sigma->GetXaxis()->SetRangeUser(0, 350);
  g_sigma->Draw("AP");

  TCanvas *c4 = new TCanvas("c4", "Bias vs energy", 900, 700);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.15);
  g_bias->SetMarkerStyle(20);
  g_bias->SetMarkerSize(1.2);
  g_bias->GetXaxis()->SetTitle("E_{fit} (MeV)");
  g_bias->GetYaxis()->SetTitle("Mean residual");
  g_bias->GetXaxis()->CenterTitle();
  g_bias->GetYaxis()->CenterTitle();
  g_bias->GetXaxis()->SetRangeUser(0, 350);
  g_bias->Draw("AP");

  // ----------------------------------------------------------------------
  // Write output
  // ----------------------------------------------------------------------
  fout->cd();
  EPhoTree->Write();
  h_phoE->Write();
  h_pull->Write();
  for (int b = 0; b < nbins; ++b) pull_hists[b]->Write();
  g_bias->Write();
  g_sigma->Write();
  fout->Close();

  tree_file->Close();

  cout << "\n✅ Output written to ../pull_scan/pull_scan.root" << endl;
}
