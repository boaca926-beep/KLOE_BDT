// ============================================================================
// pull_tuning.C
//
// Fits pull distributions of E1+E2 (π⁰ photons) and E3 (ISR photon)
// with a pure Gaussian (signal) fit. Extracts bias (mean) and scale (sigma)
// for MC and writes them to a header file for use in tree_cut_bdt_tuning.C.
// ============================================================================

#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TKey.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TString.h>

// ----------------------------------------------------------------------------
// Helper to list contents of a ROOT file (for debugging)
// ----------------------------------------------------------------------------
void checkFile(TFile *f_input){
  TIter next_tree(f_input->GetListOfKeys());
  TString objnm_tree, classnm_tree;
  int i = 0;
  TKey *key;
  while ( (key = (TKey *) next_tree()) ) {
    i++;
    objnm_tree   = key->GetName();
    classnm_tree = key->GetClassName();
    key->GetSeekKey();
    cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
  }
}

// ----------------------------------------------------------------------------
// Structure to store fit results
// ----------------------------------------------------------------------------
struct FitResult {
    TString name;
    double mean, mean_err;
    double sigma, sigma_err;
    double chi2_ndf;
    int entries;
};

// ----------------------------------------------------------------------------
// Main tuning macro
// ----------------------------------------------------------------------------
void pull_tuning() {

  gROOT->GetListOfCanvases()->Delete();
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(5);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  gSystem->Exec("mkdir -p ../pull_tuning");

  // --- Path to the histogram files (produced by tree_cut_bdt_raw.C) ---
  const TString tuning_type = "tuning";
  const TString input_folder = "../output_bdt_" + tuning_type;   // <-- adjust to your hist folder

  TString fin_E1_nm = input_folder + "_pull_E1/hist_pull_E1.root";
  TFile *fin_E1 = new TFile(fin_E1_nm);
  if (!fin_E1 || fin_E1->IsZombie()) {
    std::cerr << "ERROR: cannot open " << fin_E1_nm << std::endl;
    return;
  }

  TString fin_E2_nm = input_folder + "_pull_E2/hist_pull_E2.root";
  TFile *fin_E2 = new TFile(fin_E2_nm);
  if (!fin_E2 || fin_E2->IsZombie()) {
    std::cerr << "ERROR: cannot open " << fin_E2_nm << std::endl;
    return;
  }

  TString fin_E3_nm = input_folder + "_pull_E3/hist_pull_E3.root";
  TFile *fin_E3 = new TFile(fin_E3_nm);
  if (!fin_E3 || fin_E3->IsZombie()) {
    std::cerr << "ERROR: cannot open " << fin_E3_nm << std::endl;
    return;
  }

  // --- Retrieve histograms ---
  // MC (signal) and Data histograms from each pull file
  TH1D* hE1_MC   = (TH1D*)fin_E1->Get("hist_isr3pi_sc");
  TH1D* hE1_DATA = (TH1D*)fin_E1->Get("hist_data");

  TH1D* hE2_MC   = (TH1D*)fin_E2->Get("hist_isr3pi_sc");
  TH1D* hE2_DATA = (TH1D*)fin_E2->Get("hist_data");

  TH1D* hE3_MC   = (TH1D*)fin_E3->Get("hist_isr3pi_sc");
  TH1D* hE3_DATA = (TH1D*)fin_E3->Get("hist_data");

  // --- Combine E1 and E2 for the two π⁰ photons ---
  TH1D *hE12_MC = (TH1D*)hE1_MC->Clone("hE12_MC");
  hE12_MC->Add(hE2_MC, 1.);

  TH1D *hE12_DATA = (TH1D*)hE1_DATA->Clone("hE12_DATA");
  hE12_DATA->Add(hE2_DATA, 1.);

  hE3_MC->SetName("hE3_MC");
  hE3_DATA->SetName("hE3_DATA");

  // --- Prepare arrays for fitting ---
  const int nb_hist = 4;
  TH1D *HLIST[nb_hist] = {hE12_MC, hE12_DATA, hE3_MC, hE3_DATA};
  FitResult results[nb_hist];
  TF1 *gaus_fits[nb_hist];

  // --- Loop over histograms and fit ---
  for (int i = 0; i < nb_hist; i++) {

    TH1D *h_tmp = (TH1D*)HLIST[i]->Clone(Form("%s", HLIST[i]->GetName()));
    double mean_est = h_tmp->GetMean();
    double rms_est  = h_tmp->GetRMS();
    double amp_est  = h_tmp->GetBinContent(h_tmp->FindBin(mean_est));

    // Fit range: ±1.5σ around mean to focus on the signal peak
    const double fit_width = 1.5;   // narrower than before
    const double fit_min = mean_est - fit_width * rms_est;
    const double fit_max = mean_est + fit_width * rms_est;

    std::cout << "========================================" << std::endl;
    std::cout << "Sample: " << h_tmp->GetName() << std::endl;
    std::cout << "Number of entries: " << h_tmp->GetEntries() << std::endl;
    std::cout << "Mean estimate: " << mean_est << ", RMS: " << rms_est << std::endl;
    std::cout << "Fit range: [" << fit_min << ", " << fit_max << "]" << std::endl;

    // ----------------------------------------------------------------
    // Fit function: PURE GAUSSIAN (signal only) – no background
    //   f(x) = A * exp(-0.5*((x-mean)/sigma)^2)
    // Parameters: [0]=A, [1]=mean, [2]=sigma
    // ----------------------------------------------------------------
    TF1 *fitFunc = new TF1("fitFunc", "gaus(0)", fit_min, fit_max);
    fitFunc->SetParameters(amp_est, mean_est, rms_est * 0.8);
    // Parameter limits for stability
    fitFunc->SetParLimits(1, mean_est - rms_est, mean_est + rms_est);   // constrain mean
    fitFunc->SetParLimits(2, 0.1, 10.0);                                // sigma > 0
    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineWidth(2);

    // Perform the fit (R=range, Q=quiet, S=store result)
    h_tmp->Fit(fitFunc, "RQS");

    double mean = fitFunc->GetParameter(1);
    double sigma = fitFunc->GetParameter(2);
    double mean_err = fitFunc->GetParError(1);
    double sigma_err = fitFunc->GetParError(2);
    double chi2ndf = fitFunc->GetChisquare() / fitFunc->GetNDF();

    std::cout << "Gaussian fit: mean = " << mean << " +/- " << mean_err << std::endl;
    std::cout << "Gaussian fit: sigma = " << sigma << " +/- " << sigma_err << std::endl;
    std::cout << "χ²/ndf = " << chi2ndf << std::endl;

    // Store results
    results[i].name = TString(h_tmp->GetName());
    results[i].mean = mean;
    results[i].mean_err = mean_err;
    results[i].sigma = sigma;
    results[i].sigma_err = sigma_err;
    results[i].chi2_ndf = chi2ndf;
    results[i].entries = h_tmp->GetEntries();
    gaus_fits[i] = fitFunc;
  }

  // --- Set colors for the fit lines (for plotting) ---
  gaus_fits[0]->SetLineColor(kRed);
  gaus_fits[0]->SetLineWidth(2);
  gaus_fits[1]->SetLineColor(kGreen+2);
  gaus_fits[1]->SetLineWidth(2);
  gaus_fits[2]->SetLineColor(kRed);
  gaus_fits[2]->SetLineWidth(2);
  gaus_fits[3]->SetLineColor(kGreen+2);
  gaus_fits[3]->SetLineWidth(2);

  // ========== Canvas 1: E1+E2 Pull ==========
  TCanvas *c_E12 = new TCanvas("c_E12", "E1+E2 Pull Distributions", 900, 900);
  c_E12->cd(1);
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);

  const double ymax_E12 = hE12_DATA->GetMaximum();
  hE12_DATA->GetYaxis()->SetTitle("Events");
  hE12_DATA->GetYaxis()->SetRangeUser(0.01, ymax_E12 * 1.6);
  hE12_DATA->GetYaxis()->CenterTitle();
  hE12_DATA->GetYaxis()->SetTitleSize(0.05);
  hE12_DATA->GetYaxis()->SetTitleOffset(1.4);
  hE12_DATA->GetYaxis()->SetLabelSize(0.04);
  hE12_DATA->GetXaxis()->SetTitle("E_{1}+E_{2} Pull [MeV]");
  hE12_DATA->GetXaxis()->SetTitleSize(0.05);
  hE12_DATA->GetXaxis()->SetTitleOffset(1.2);
  hE12_DATA->GetXaxis()->SetLabelSize(0.04);
  hE12_DATA->GetXaxis()->CenterTitle();

  hE12_DATA->Draw();
  gaus_fits[1]->Draw("same");
  hE12_MC->Draw("same hist");
  gaus_fits[0]->Draw("same");

  TLegend *leg_E12 = new TLegend(0.2, 0.7, 0.9, 0.9);
  leg_E12->SetFillStyle(0);
  leg_E12->SetBorderSize(0);
  leg_E12->SetNColumns(2);
  leg_E12->SetTextSize(0.03);
  leg_E12->AddEntry(hE12_MC, "MC", "f");
  leg_E12->AddEntry(gaus_fits[0], Form("#mu = %.3f, #sigma = %.3f", results[0].mean, results[0].sigma), "l");
  leg_E12->AddEntry(hE12_DATA, "Data", "f");
  leg_E12->AddEntry(gaus_fits[1], Form("#mu = %.3f, #sigma = %.3f", results[1].mean, results[1].sigma), "l");
  leg_E12->Draw();
  gPad->Update();

  // ========== Canvas 2: E3 Pull ==========
  TCanvas *c_E3 = new TCanvas("c_E3", "E3 Pull Distributions", 900, 900);
  c_E3->cd(1);
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);

  const double ymax_E3 = hE3_DATA->GetMaximum();
  hE3_DATA->GetYaxis()->SetTitle("Events");
  hE3_DATA->GetYaxis()->SetRangeUser(0.01, ymax_E3 * 1.6);
  hE3_DATA->GetYaxis()->CenterTitle();
  hE3_DATA->GetYaxis()->SetTitleSize(0.05);
  hE3_DATA->GetYaxis()->SetTitleOffset(1.4);
  hE3_DATA->GetYaxis()->SetLabelSize(0.04);
  hE3_DATA->GetXaxis()->SetTitle("E_{3} Pull [MeV]");
  hE3_DATA->GetXaxis()->SetTitleSize(0.05);
  hE3_DATA->GetXaxis()->SetTitleOffset(1.2);
  hE3_DATA->GetXaxis()->SetLabelSize(0.04);
  hE3_DATA->GetXaxis()->CenterTitle();

  hE3_DATA->Draw();
  gaus_fits[3]->Draw("same");
  hE3_MC->Draw("same hist");
  gaus_fits[2]->Draw("same");

  TLegend *leg_E3 = new TLegend(0.2, 0.7, 0.9, 0.9);
  leg_E3->SetFillStyle(0);
  leg_E3->SetBorderSize(0);
  leg_E3->SetNColumns(2);
  leg_E3->SetTextSize(0.03);
  leg_E3->AddEntry(hE3_MC, "MC", "f");
  leg_E3->AddEntry(gaus_fits[2], Form("#mu = %.3f, #sigma = %.3f", results[2].mean, results[2].sigma), "l");
  leg_E3->AddEntry(hE3_DATA, "Data", "f");
  leg_E3->AddEntry(gaus_fits[3], Form("#mu = %.3f, #sigma = %.3f", results[3].mean, results[3].sigma), "l");
  leg_E3->Draw();
  gPad->Update();

  // --- Save canvases ---
  c_E12->SaveAs("../pull_tuning/pull_E12.pdf");
  c_E3->SaveAs("../pull_tuning/pull_E3.pdf");

  // --- Print summary of fit results ---
  std::cout << "\n========================================" << std::endl;
  std::cout << "Summary of E Pull Fit Results (Pure Gaussian):" << std::endl;
  std::cout << "========================================" << std::endl;
  for (int i = 0; i < nb_hist; i++) {
    std::cout << Form("%-15s: mean = %6.3f +/- %6.3f, sigma = %6.3f +/- %6.3f, χ²/ndf = %.3f",
                      results[i].name.Data(),
                      results[i].mean, results[i].mean_err,
                      results[i].sigma, results[i].sigma_err,
                      results[i].chi2_ndf) << std::endl;
  }

  // --- Extract MC parameters (results[0] = MC E1+E2, results[2] = MC E3) ---
  double bias_E12         = results[0].mean;        // MC bias for π⁰ photons
  double sigma_scale_E12  = results[0].sigma;       // MC scale for π⁰ photons
  double bias_E3          = results[2].mean;        // MC bias for ISR photon
  double sigma_scale_E3   = results[2].sigma;       // MC scale for ISR photon

  // --- Write the MC parameters to the tuning header ---
  std::ofstream myfile;
  TString myfile_nm = "../header_bdt/tuning.h";
  myfile.open(myfile_nm.Data());
  myfile << "// MC pull parameters extracted from BDT-selected signal MC\n";
  myfile << "// Fitted with pure Gaussian (pull_tuning.C)\n";
  myfile << "const double bias_E12 = " << bias_E12 << ";\n";
  myfile << "const double sigma_scale_E12 = " << sigma_scale_E12 << ";\n";
  myfile << "const double bias_E3 = " << bias_E3 << ";\n";
  myfile << "const double sigma_scale_E3 = " << sigma_scale_E3 << ";\n\n";
  myfile.close();

  std::cout << "\nParameters written to " << myfile_nm << std::endl;
  std::cout << "  bias_E12         = " << bias_E12 << std::endl;
  std::cout << "  sigma_scale_E12  = " << sigma_scale_E12 << std::endl;
  std::cout << "  bias_E3          = " << bias_E3 << std::endl;
  std::cout << "  sigma_scale_E3   = " << sigma_scale_E3 << std::endl;
  std::cout << "Parameters are saved at " << myfile_nm << std::endl;
}
