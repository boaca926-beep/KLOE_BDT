// ============================================================================
// omega_resolution.C
//
// Fits the ω peak in the 3π invariant mass (Br_m3pi_bdt) for data and MC.
// Extracts the Gaussian width (sigma) which reflects the detector resolution.
// Usage:
//   .x omega_resolution.C("TDATA", "Data")
//   .x omega_resolution.C("TISR3PI_SIG_PEAK", "Signal")
//
// The macro writes a PDF plot and a ROOT file with the histogram and fit result.
// ============================================================================

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>

using namespace std;

// ----------------------------------------------------------------------------
// Fit ω peak with Gaussian + 2nd-order polynomial background
// ----------------------------------------------------------------------------
bool fitOmega(TH1D *h, double &mean, double &sigma,
              double &mean_err, double &sigma_err,
              double &chi2ndf) {
  if (!h || h->GetEntries() < 100) return false;

  // Fit range: 740–810 MeV (adjust if needed)
  double xmin = 740.0;
  double xmax = 810.0;
  h->GetXaxis()->SetRangeUser(xmin, xmax);

  // Estimate initial parameters from histogram
  double peak = h->GetMaximum();
  double mean0 = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
  double rms = h->GetRMS();
  double sigma_guess = (rms > 0) ? rms * 0.6 : 5.0; // rough estimate

  // Define fit function: Gaussian + polynomial (order 2)
  TF1 *fit = new TF1("omega_fit", "gaus(0) + pol2(3)", xmin, xmax);
  // Parameters: [0] = amplitude, [1] = mean, [2] = sigma,
  //             [3] = p0, [4] = p1, [5] = p2
  fit->SetParameters(peak, mean0, sigma_guess, 0.0, 0.0, 0.0);
  fit->SetParLimits(0, 0, peak * 3);
  fit->SetParLimits(1, mean0 - 2.0, mean0 + 2.0);
  fit->SetParLimits(2, 0.5, 15.0);
  fit->SetParLimits(3, -peak/2, peak/2);
  fit->SetParLimits(4, -peak/10, peak/10);
  fit->SetParLimits(5, -peak/100, peak/100);

  // Perform fit
  Int_t status = h->Fit(fit, "RQS");
  double chi2 = fit->GetChisquare();
  int ndf = fit->GetNDF();
  chi2ndf = (ndf > 0) ? chi2 / ndf : 0;

  bool ok = false;
  if (status == 0 && chi2ndf < 3.0) {
    mean = fit->GetParameter(1);
    sigma = fit->GetParameter(2);
    mean_err = fit->GetParError(1);
    sigma_err = fit->GetParError(2);
    ok = (sigma > 0 && sigma_err > 0);
  }
  delete fit;
  return ok;
}

// ----------------------------------------------------------------------------
// Main macro
// ----------------------------------------------------------------------------
void omega_resolution(const TString tree_type = "TDATA",
                      const TString sample_type = "Data") {
  // ----------------------------------------------------------------------
  // Configuration
  // ----------------------------------------------------------------------
  // Adjust input file path as needed – this is the output of tree_cut_bdt_tuning.C
  TString input_file_nm = "/home/bo/Desktop/bdt_tuning_TDATA_chain_false/cut/tree_pre.root";
  // If you used a different folder, change the path above.

  TString output_file = Form("../omega_resolution/omega_%s.root", tree_type.Data());
  TString pdf_name   = Form("../omega_resolution/omega_fit_%s.pdf", tree_type.Data());

  gSystem->Exec("mkdir -p ../omega_resolution");

  cout << "\n========================================" << endl;
  cout << "  OMEGA RESOLUTION EXTRACTION" << endl;
  cout << "  Sample: " << sample_type << endl;
  cout << "  Tree:   " << tree_type << endl;
  cout << "========================================\n" << endl;

  gROOT->GetListOfCanvases()->Delete();
  gErrorIgnoreLevel = kError;
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  TH1::AddDirectory(kFALSE);

  // ----------------------------------------------------------------------
  // Open input file and get tree
  // ----------------------------------------------------------------------
  TFile *tree_file = new TFile(input_file_nm);
  if (!tree_file || tree_file->IsZombie()) {
    cerr << "ERROR: Cannot open " << input_file_nm << endl;
    return;
  }

  TTree *INPUT_TREE = (TTree*)tree_file->Get(tree_type);
  if (!INPUT_TREE) {
    cerr << "ERROR: Cannot find tree " << tree_type << endl;
    tree_file->Close();
    return;
  }
  cout << "✓ Loaded " << tree_type << " with " << INPUT_TREE->GetEntries() << " entries" << endl;

  // ----------------------------------------------------------------------
  // Histogram of M3π (Br_m3pi_bdt)
  // ----------------------------------------------------------------------
  TH1D *h_m3pi = new TH1D("h_m3pi", "", 200, 700, 850);
  h_m3pi->Sumw2();

  double m3pi = 0.;
  INPUT_TREE->SetBranchAddress("Br_m3pi_bdt", &m3pi);
  if (!INPUT_TREE->GetBranch("Br_m3pi_bdt")) {
    cerr << "ERROR: Branch Br_m3pi_bdt not found in tree " << tree_type << endl;
    tree_file->Close();
    return;
  }

  Long64_t nentries = INPUT_TREE->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    INPUT_TREE->GetEntry(i);
    h_m3pi->Fill(m3pi);
  }

  // ----------------------------------------------------------------------
  // Fit the ω peak
  // ----------------------------------------------------------------------
  double mean = 0., sigma = 0.;
  double mean_err = 0., sigma_err = 0.;
  double chi2ndf = 0.;

  bool ok = fitOmega(h_m3pi, mean, sigma, mean_err, sigma_err, chi2ndf);

  if (!ok) {
    cerr << "ERROR: Fit failed for " << tree_type << endl;
    tree_file->Close();
    return;
  }

  // ----------------------------------------------------------------------
  // Draw and save
  // ----------------------------------------------------------------------
  TCanvas *c = new TCanvas("c", "Omega fit", 900, 700);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.15);

  h_m3pi->SetLineWidth(2);
  h_m3pi->SetLineColor(kBlack);
  h_m3pi->GetXaxis()->SetTitle("M_{3#pi} (MeV)");
  h_m3pi->GetYaxis()->SetTitle("Entries / 0.75 MeV");
  h_m3pi->GetXaxis()->CenterTitle();
  h_m3pi->GetYaxis()->CenterTitle();
  h_m3pi->GetYaxis()->SetRangeUser(0.01, 1.4 * h_m3pi->GetMaximum());
  h_m3pi->Draw("E0");

  // Overlay the fit function
  TF1 *fitFunc = (TF1*)h_m3pi->GetFunction("omega_fit");
  if (fitFunc) {
    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineWidth(2);
    fitFunc->Draw("same");
  }

  // Draw the core Gaussian separately for visualization
  TF1 *core = new TF1("core_gaus", "gaus", mean - 3*sigma, mean + 3*sigma);
  core->SetParameters(fitFunc->GetParameter(0), mean, sigma);
  core->SetLineColor(kBlue);
  core->SetLineStyle(2);
  core->Draw("same");

  // Text box
  TPaveText *pt = new TPaveText(0.55, 0.70, 0.90, 0.90, "NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.04);
  pt->AddText(Form("Sample: %s", sample_type.Data()));
  pt->AddText(Form("Mean = %.3f ± %.3f MeV", mean, mean_err));
  pt->AddText(Form("#sigma = %.3f ± %.3f MeV", sigma, sigma_err));
  pt->AddText(Form("#chi^{2}/NDF = %.2f", chi2ndf));
  pt->AddText(Form("Entries = %d", (int)h_m3pi->GetEntries()));
  pt->Draw();

  TLegend *leg = new TLegend(0.15, 0.70, 0.45, 0.90);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_m3pi, sample_type, "lep");
  leg->AddEntry(fitFunc, "Total fit (Gaus + poly2)", "l");
  leg->AddEntry(core, "Core Gaussian", "l");
  leg->Draw();

  c->SaveAs(pdf_name);

  // ----------------------------------------------------------------------
  // Save results to ROOT file
  // ----------------------------------------------------------------------
  TFile *fout = new TFile(output_file, "RECREATE");
  fout->cd();
  h_m3pi->Write();
  if (fitFunc) fitFunc->Write("fitFunc");
  core->Write("coreGaus");
  // Store fit parameters as a TTree
  TTree *resultTree = new TTree("result", "Fit results");
  double mean_val = mean, sigma_val = sigma;
  double mean_err_val = mean_err, sigma_err_val = sigma_err;
  resultTree->Branch("mean", &mean_val, "mean/D");
  resultTree->Branch("mean_err", &mean_err_val, "mean_err/D");
  resultTree->Branch("sigma", &sigma_val, "sigma/D");
  resultTree->Branch("sigma_err", &sigma_err_val, "sigma_err/D");
  resultTree->Branch("chi2ndf", &chi2ndf, "chi2ndf/D");
  resultTree->Fill();
  resultTree->Write();
  fout->Close();

  cout << "\n✅ Fit results:" << endl;
  cout << "   Mean   = " << mean << " ± " << mean_err << " MeV" << endl;
  cout << "   Sigma  = " << sigma << " ± " << sigma_err << " MeV" << endl;
  cout << "   χ²/NDF = " << chi2ndf << endl;
  cout << "   Output saved to " << output_file << " and " << pdf_name << endl;

  tree_file->Close();
  delete c;
}
