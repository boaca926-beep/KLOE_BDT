// MC tuning on final states particles gamma1, 2 and 3. Energy, position and cluster time
// mc_tuning.C
// Fit combined pull distributions for data and signal MC to extract mean and sigma
#include "../header_bdt/plot_resol.h"
#include <TF1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>

// Define FitResult struct BEFORE using it
struct FitResult {
    TString name;
    double mean, mean_err;
    double sigma, sigma_err;
    double chi2_ndf;
    int entries;
};

void mc_tuning() {
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  gSystem->Exec("mkdir -p ../plots_tuning");

  TFile *fout = new TFile("../plots_tuning/pull_fits.root", "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: cannot create output file." << std::endl;
    return;
  }

  //TString treeFile = "/home/bo/Desktop/input_bdt_TDATA_norm/cut/tree_pre_bdt.root";
  TFile *ftree = TFile::Open(treeFile);
  if (!ftree || ftree->IsZombie()) {
    std::cerr << "ERROR: cannot open " << treeFile << std::endl;
    return;
  }

  const int nb_mc_type = 2;
  const TString MC_TYPE[nb_mc_type] = {"TISR3PI_SIG", "TDATA"};
  const TString MC_LABEL[nb_mc_type] = {"Signal MC", "Data"};
  const int nb_pull_type = 2;
  const TString pull_names[2] = {"E1", "E2"};
  
  const int bin_size = 100;
  const double XMIN = -5.0;
  const double XMAX = 5.0;
  
  // Different fit widths for MC and Data
  const double fit_width_mc = 1.5;
  const double fit_width_data = 2.0;

  // ---- SIDE BY SIDE PLOT ----
  TCanvas *c_side = new TCanvas("c_side", "Pull Distributions: MC vs Data", 1400, 700);
  c_side->Divide(nb_mc_type, 1);

  // Store histograms for overlay
  TH1D *h_mc = nullptr;
  TH1D *h_data = nullptr;
  TF1 *gaus_mc = nullptr;
  TF1 *gaus_data = nullptr;
  FitResult result_mc, result_data;

  // Loop over data/MC types
  for (int i = 0; i < nb_mc_type; i++) {
    TString mc_type = MC_TYPE[i];
    TString mc_label = MC_LABEL[i];
    TTree *ttree = (TTree*) ftree->Get(mc_type);
    if (!ttree) {
      std::cerr << "Tree " << mc_type << " not found, skipping." << std::endl;
      continue;
    }

    // Create histogram for COMBINED E1+E2 pulls
    TH1D *h_pull_combined = new TH1D(Form("h_pull_%s_E1plusE2", mc_type.Data()), 
                                      "", bin_size, XMIN, XMAX);
    h_pull_combined->SetMarkerStyle(20);
    h_pull_combined->SetMarkerSize(0.6);
    h_pull_combined->Sumw2();

    // Fill combined pulls from E1 and E2
    double pull_val = 0.0;
    for (int p = 0; p < nb_pull_type; p++) {
      TString pull_name = pull_names[p];
      TString branch_name = "Br_pull_" + pull_name;
      ttree->SetBranchAddress(branch_name, &pull_val);

      Long64_t nentries = ttree->GetEntries();
      for (Long64_t j = 0; j < nentries; j++) {
        ttree->GetEntry(j);
        h_pull_combined->Fill(pull_val);
      }
    }

    // Check if histogram has entries
    if (h_pull_combined->GetEntries() < 10) {
      std::cerr << "WARNING: Too few entries for " << mc_type 
                << " (" << h_pull_combined->GetEntries() << "), skipping." << std::endl;
      continue;
    }

    // Fit the combined distribution
    double mean_est = h_pull_combined->GetMean();
    double rms_est = h_pull_combined->GetRMS();
    double amp_est = h_pull_combined->GetBinContent(h_pull_combined->FindBin(mean_est));
    
    // Choose fit width based on sample type
    double fit_width = (mc_type == "TDATA") ? fit_width_data : fit_width_mc;
    const double fit_min = mean_est - fit_width;
    const double fit_max = mean_est + fit_width;
    
    std::cout << "========================================" << std::endl;
    std::cout << "Sample: " << mc_label << " (" << mc_type << ")" << std::endl;
    std::cout << "Combined Pull: E1 + E2" << std::endl;
    std::cout << "Number of entries: " << h_pull_combined->GetEntries() << std::endl;
    std::cout << "Mean estimate: " << mean_est << ", RMS: " << rms_est << std::endl;
    std::cout << "Fit range: [" << fit_min << ", " << fit_max << "]" << std::endl;
    
    TF1 *gaus = new TF1(Form("gaus_%s", mc_type.Data()), "gaus", fit_min, fit_max);
    gaus->SetParameters(amp_est, mean_est, rms_est * 0.7);
    gaus->SetParLimits(2, 0.2, 5.0);
    
    h_pull_combined->Fit(gaus, "RQS");
    
    double mean = gaus->GetParameter(1);
    double sigma = gaus->GetParameter(2);
    double mean_err = gaus->GetParError(1);
    double sigma_err = gaus->GetParError(2);
    double chi2_ndf = gaus->GetChisquare() / gaus->GetNDF();
    
    std::cout << "Mean = " << mean << " +/- " << mean_err << std::endl;
    std::cout << "Sigma = " << sigma << " +/- " << sigma_err << std::endl;
    std::cout << "χ²/ndf = " << chi2_ndf << std::endl;
    
    // Quality checks
    if (fabs(mean) > 0.1) {
      std::cout << "WARNING: Mean is significantly shifted from 0!" << std::endl;
    }
    if (fabs(sigma - 1.0) > 0.1) {
      std::cout << "WARNING: Sigma deviates from 1.0 by " << (sigma - 1.0) << std::endl;
    }
    if (chi2_ndf > 2.0) {
      std::cout << "WARNING: Poor fit quality (χ²/ndf = " << chi2_ndf << ")" << std::endl;
    }
    std::cout << "========================================" << std::endl;
    
    // Store for overlay
    if (mc_type == "TISR3PI_SIG") {
      h_mc = h_pull_combined;
      gaus_mc = gaus;
      result_mc.name = mc_type;
      result_mc.mean = mean;
      result_mc.mean_err = mean_err;
      result_mc.sigma = sigma;
      result_mc.sigma_err = sigma_err;
      result_mc.chi2_ndf = chi2_ndf;
      result_mc.entries = h_pull_combined->GetEntries();
    } else {
      h_data = h_pull_combined;
      gaus_data = gaus;
      result_data.name = mc_type;
      result_data.mean = mean;
      result_data.mean_err = mean_err;
      result_data.sigma = sigma;
      result_data.sigma_err = sigma_err;
      result_data.chi2_ndf = chi2_ndf;
      result_data.entries = h_pull_combined->GetEntries();
    }
    
    // ---- SIDE BY SIDE PLOT - FIXED ----
    // IMPORTANT: cd() to the correct pad BEFORE drawing
    c_side->cd(i + 1);  // Pad 1: MC, Pad 2: Data
    
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.12);
    
    // Normalize histograms to unit area for comparison
    TH1D *h_norm = (TH1D*) h_pull_combined->Clone(Form("h_norm_%s", mc_type.Data()));
    double integral = h_norm->Integral();
    if (integral > 0) h_norm->Scale(1.0 / integral);
    h_norm->SetMarkerStyle(20);
    h_norm->SetMarkerSize(0.6);
    h_norm->SetLineColor(kBlack);
    
    h_norm->GetYaxis()->SetTitle("Normalized Entries");
    h_norm->GetXaxis()->SetTitle("Pull (E1+E2)");
    h_norm->GetYaxis()->SetRangeUser(0, h_norm->GetMaximum() * 1.3);
    h_norm->Draw("E");
    
    // Draw fit (scaled to match normalized histogram)
    TF1 *gaus_norm = new TF1(Form("gaus_norm_%s", mc_type.Data()), "gaus", fit_min, fit_max);
    if (integral > 0) {
      gaus_norm->SetParameters(gaus->GetParameter(0) / integral, 
                                gaus->GetParameter(1), 
                                gaus->GetParameter(2));
    } else {
      gaus_norm->SetParameters(gaus->GetParameter(0), gaus->GetParameter(1), gaus->GetParameter(2));
    }
    gaus_norm->SetLineColor(kRed);
    gaus_norm->SetLineWidth(2);
    gaus_norm->Draw("same");
    
    // Draw vertical line at 0 for reference
    TLine *line0 = new TLine(0, 0, 0, h_norm->GetMaximum() * 1.15);
    line0->SetLineStyle(2);
    line0->SetLineColor(kGray);
    line0->Draw();
    
    // Draw vertical lines at ±1 for reference (ideal pull distribution)
    TLine *line1a = new TLine(-1, 0, -1, h_norm->GetMaximum() * 1.15);
    line1a->SetLineStyle(3);
    line1a->SetLineColor(kGray+1);
    line1a->Draw();
    TLine *line1b = new TLine(1, 0, 1, h_norm->GetMaximum() * 1.15);
    line1b->SetLineStyle(3);
    line1b->SetLineColor(kGray+1);
    line1b->Draw();
    
    // Legend - position adjusted for each pad
    TLegend *leg = new TLegend(0.55, 0.70, 0.92, 0.88);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(h_norm, Form("%s (E1+E2)", mc_label.Data()), "lep");
    leg->AddEntry(gaus_norm, Form("#mu = %.3f, #sigma = %.3f", mean, sigma), "l");
    leg->Draw();
    
    // Add text with χ²/ndf
    TPaveText *pt = new TPaveText(0.15, 0.78, 0.48, 0.88, "NDC");
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->SetTextSize(0.04);
    pt->AddText(Form("#chi^{2}/ndf = %.2f", chi2_ndf));
    pt->Draw();
    
    // Update the pad
    gPad->Update();
    
    // Write histogram to output file
    fout->cd();
    h_pull_combined->Write();
    gaus->Write();
  }
  
  // Important: Update canvas before saving
  c_side->Update();
  c_side->SaveAs("../plots_tuning/pull_fits_MC_vs_Data.pdf");
  
  // ---- OVERLAY PLOT (MC vs Data on same canvas) ----
  if (h_mc && h_data) {
    TCanvas *c_overlay = new TCanvas("c_overlay", "Pull Distributions: MC vs Data Overlay", 1200, 700);
    c_overlay->SetLeftMargin(0.12);
    c_overlay->SetRightMargin(0.05);
    c_overlay->SetTopMargin(0.05);
    c_overlay->SetBottomMargin(0.12);
    
    // Normalize both histograms
    TH1D *h_mc_norm = (TH1D*) h_mc->Clone("h_mc_norm");
    TH1D *h_data_norm = (TH1D*) h_data->Clone("h_data_norm");
    double mc_integral = h_mc_norm->Integral();
    double data_integral = h_data_norm->Integral();
    if (mc_integral > 0) h_mc_norm->Scale(1.0 / mc_integral);
    if (data_integral > 0) h_data_norm->Scale(1.0 / data_integral);
    
    h_mc_norm->SetLineColor(kBlue);
    h_mc_norm->SetLineWidth(2);
    h_mc_norm->SetMarkerStyle(20);
    h_mc_norm->SetMarkerSize(0.6);
    h_mc_norm->SetMarkerColor(kBlue);
    
    h_data_norm->SetLineColor(kRed);
    h_data_norm->SetLineWidth(2);
    h_data_norm->SetMarkerStyle(24);
    h_data_norm->SetMarkerSize(0.8);
    h_data_norm->SetMarkerColor(kRed);
    
    double max_y = std::max(h_mc_norm->GetMaximum(), h_data_norm->GetMaximum());
    h_mc_norm->GetYaxis()->SetTitle("Normalized Entries");
    h_mc_norm->GetXaxis()->SetTitle("Pull (E1+E2)");
    h_mc_norm->GetYaxis()->SetRangeUser(0, max_y * 1.25);
    h_mc_norm->Draw("E");
    h_data_norm->Draw("E same");
    
    // Draw vertical line at 0
    TLine *line0 = new TLine(0, 0, 0, max_y * 1.15);
    line0->SetLineStyle(2);
    line0->SetLineColor(kGray);
    line0->Draw();
    
    TLegend *leg_overlay = new TLegend(0.55, 0.72, 0.92, 0.92);
    leg_overlay->SetFillStyle(0);
    leg_overlay->SetBorderSize(0);
    leg_overlay->SetTextSize(0.04);
    leg_overlay->AddEntry(h_mc_norm, Form("Signal MC: #mu=%.3f, #sigma=%.3f", 
                           result_mc.mean, result_mc.sigma), "lep");
    leg_overlay->AddEntry(h_data_norm, Form("Data: #mu=%.3f, #sigma=%.3f", 
                           result_data.mean, result_data.sigma), "lep");
    leg_overlay->Draw();
    
    c_overlay->Update();
    c_overlay->SaveAs("../plots_tuning/pull_fits_MC_vs_Data_overlay.pdf");
    delete c_overlay;
  }
  
  // Print summary table
  std::cout << "\n========================================" << std::endl;
  std::cout << "SUMMARY OF PULL FIT RESULTS" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "Sample\t\tEntries\tMean\t\tSigma\t\tχ²/ndf" << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  if (h_mc) {
    std::cout << "Signal MC\t" << result_mc.entries << "\t" 
              << result_mc.mean << " ± " << result_mc.mean_err << "\t"
              << result_mc.sigma << " ± " << result_mc.sigma_err << "\t"
              << result_mc.chi2_ndf << std::endl;
  }
  if (h_data) {
    std::cout << "Data\t\t" << result_data.entries << "\t" 
              << result_data.mean << " ± " << result_data.mean_err << "\t"
              << result_data.sigma << " ± " << result_data.sigma_err << "\t"
              << result_data.chi2_ndf << std::endl;
  }
  std::cout << "========================================" << std::endl;
  
  // Clean up
  delete c_side;
  fout->Close();
  ftree->Close();
  
  std::cout << "\nDone. Results saved to:" << std::endl;
  std::cout << "  - ../plots_tuning/pull_fits.root" << std::endl;
  std::cout << "  - ../plots_tuning/pull_fits_MC_vs_Data.pdf (side by side)" << std::endl;
  std::cout << "  - ../plots_tuning/pull_fits_MC_vs_Data_overlay.pdf (overlay)" << std::endl;
}
