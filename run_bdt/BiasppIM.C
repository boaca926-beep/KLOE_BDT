// ============================================================================
// track_compr.C
//
// Compare track bias (median of ppIM)  between signal and data.
// Reads:
//   ../pull_scan/pull_scan_TISR3PI_SIG_PEAK.root
//   ../pull_scan/pull_scan_TDATA.root
// Plots:
//   - cv1: track bias vs ppIM for signal and data - signal (real points only).
// Saves:
//   - ../Bias_ppIM/data_mc_comparison_bkg_sub_ppIM.pdf
// ============================================================================

#include <TFile.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TGaxis.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TLine.h>
#include <TPad.h>
#include <iostream>
#include <vector>

#include "../header_method/method.h"
#include "../header_plot/plot.h"

using namespace std;

struct FitResult {
  TString name;
  double mean;
  double rms;
  double median;
  int entries;
};

int BiasppIM(const TString tuning_type = "raw_false", //raw_false, tuning_false
	     const TString var_nm = "ppIM",
	     const TString var_symb = "M_{2#pi} [MeV/c^{2}]"
	     ) {

  const TString tree_file_nm = "../" + tuning_type + "_" + var_nm + "/hist.root";

  const TString out_dir = "../BiasppIM_" + tuning_type;

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  // ---- Ensure output directory exists ----
  gSystem->mkdir(out_dir, kTRUE);

  TFile* tree_file = new TFile(tree_file_nm);
  if (!tree_file || tree_file->IsZombie()) {
    cerr << "ERROR: Cannot open " << tree_file_nm << endl;
    return 1;
  }

  checkFile(tree_file); // optional, shows contents

  // ---- Retrieve raw data and signal histograms ----
  TH1D *hist_eeg = (TH1D*)tree_file->Get("hist_eeg_sc");
  TH1D *hist_signal = (TH1D*)tree_file->Get("hist_isr3pi_sc");
  TH1D *hist_omegapi = (TH1D*)tree_file->Get("hist_omegapi_sc");
  TH1D *hist_nonreson = (TH1D*)tree_file->Get("hist_nonreson_sc");
  TH1D *hist_ksl = (TH1D*)tree_file->Get("hist_ksl_sc");
  TH1D *hist_mcrest = (TH1D*)tree_file->Get("hist_mcrest_sc");
  TH1D *hist_data = (TH1D*)tree_file->Get("hist_data");

  if (!hist_data) {
    cerr << "ERROR: Data histogram not found!" << endl;
    tree_file->Close();
    return 1;
  }
  if (!hist_signal || !hist_eeg || !hist_omegapi || !hist_nonreson || !hist_ksl || !hist_mcrest) {
    cerr << "ERROR: MC histogram not found!" << endl;
    tree_file->Close();
    return 1;
  }

  TH1D *hist_bkg_sum = (TH1D*)hist_data->Clone("hist_bkg_sum");
  hist_bkg_sum->Reset();

  // Add backgrounds with optional scaling (set to 1.0 for now)
  if (hist_eeg) hist_bkg_sum->Add(hist_eeg, 1.0);
  if (hist_omegapi) hist_bkg_sum->Add(hist_omegapi, 1.0);
  if (hist_nonreson) hist_bkg_sum->Add(hist_nonreson, 1.0);
  if (hist_ksl) hist_bkg_sum->Add(hist_ksl, 1.0);
  if (hist_mcrest) hist_bkg_sum->Add(hist_mcrest, 1.0);

  // ---- Subtract backgrounds from data ----
  TH1D *hist_data_sub = (TH1D*)hist_data->Clone("hist_data_sub");
  hist_data_sub->Add(hist_bkg_sum, -1.0);

  const int nb_mass = 2;
  TH1D *hMassList[nb_mass] = {hist_signal, hist_data_sub};
  TString massNameList[nb_mass] = {"MC", "Data - Background"};
  int massColor[nb_mass] = {kRed, kRed};
  FitResult massResults[nb_mass];

  for (int i = 0; i < nb_mass; i++) {
    TH1D *h_mass = hMassList[i];
    if (!h_mass) {
      std::cerr << "Skipping mass histogram " << massNameList[i] << " (null)." << std::endl;
      continue;
    }

    // Clone to avoid modifying original
    TH1D *h_mass_copy = (TH1D*)h_mass->Clone(Form("h_mass_%s", massNameList[i].Data()));
    h_mass_copy->SetDirectory(0);

    h_mass_copy->GetXaxis()->SetRangeUser(300, 650); // avoid tail effects
    double mass_mean = h_mass_copy->GetMean();
    double mass_rms = h_mass_copy->GetRMS();

    double quantile = 0.5;
    double val;
    
    double median = h_mass_copy->GetQuantiles(1, &val, &quantile);
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Mass histogram: " << massNameList[i] << std::endl;
    std::cout << "Mean estimate: " << mass_mean << ", RMS: " << mass_rms << std::endl;
    std::cout << "Median estimate: " << val << std::endl;

    massResults[i].name = massNameList[i];
    massResults[i].mean = mass_mean;
    massResults[i].rms = mass_rms;
    massResults[i].median = val;

    for (int i = 1; i <= h_mass_copy->GetNbinsX(); ++i) {
      double x = h_mass_copy->GetBinCenter(i);
      double y = h_mass_copy->GetBinContent(i);
      //std::cout << x << "\t" << y << std::endl;
    }
    
  }

  // ---- Summary ----
  std::cout << "\n========================================" << std::endl;
  std::cout << "Summary of trak mass parameters:" << std::endl;
  std::cout << "========================================" << std::endl;
  for (int i = 0; i < nb_mass; i++) {
    if (hMassList[i]) {
      std::cout << Form("%-10s: mean = %6.3f, rms = %6.3f, median = %6.3f",
                        massResults[i].name.Data(),
                        massResults[i].mean,
                        massResults[i].rms,
                        massResults[i].median) << std::endl;
    }
  }

  // ---- Write track scale ----
  //sigma²_data=sigma²_MC + (s * mu_MC)²
  double track_scale = massResults[1].median / massResults[0].median;
  double track_scale_err = 0.0;

  double track_smearing     = TMath::Sqrt(TMath::Power(massResults[1].rms, 2) - TMath::Power(massResults[0].rms, 2)) / massResults[0].median;
  double track_smearing_err = 0.0;

  //cout << track_smearing << endl;

  
  std::ofstream myfile;
  TString myfile_nm = "../pull_scan/track_scale.txt";
  myfile.open(myfile_nm.Data());
  myfile << "const double track_scale = " << track_scale << ";\n";
  myfile << "const double track_scale_err = " << track_scale_err << ";\n";
  myfile << "const double track_smearing = " << track_smearing << ";\n";
  myfile << "const double track_smearing_err = " << track_smearing_err << ";\n";
  
  myfile.close();

  cout << "Comparison plot saved to: ../pull_scan/track_sale.txt" << endl;
  
  // Plots
  TCanvas *c1 = new TCanvas("c1", "Data/MC Comparison (Background Subtracted)", 900, 900);
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.15);

  hist_data_sub->SetMarkerStyle(20);
  hist_data_sub->SetMarkerSize(0.6);
  hist_data_sub->GetYaxis()->SetTitle("Events");
  hist_data_sub->GetYaxis()->SetRangeUser(0.01, hist_data_sub->GetMaximum() * 1.6);
  hist_data_sub->GetYaxis()->CenterTitle();
  hist_data_sub->GetYaxis()->SetTitleSize(0.05);
  hist_data_sub->GetYaxis()->SetTitleOffset(1.4);
  hist_data_sub->GetYaxis()->SetLabelSize(0.04);
  hist_data_sub->GetXaxis()->SetTitle(var_symb);
  hist_data_sub->GetXaxis()->SetTitleSize(0.05);
  hist_data_sub->GetXaxis()->SetTitleOffset(1.2);
  hist_data_sub->GetXaxis()->SetLabelSize(0.04);
  hist_data_sub->GetXaxis()->CenterTitle();

  hist_data_sub->Draw("E1");
  hist_signal->Draw("HIST SAME");

  TPaveText *pt = new TPaveText(0.2, 0.8, 0.85, 0.89, "NDC");
  //TPaveText *pt = new TPaveText(0.18, 0.7, 0.85, 0.89, "NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.03);
  pt->SetTextFont(42);
  pt->AddText(Form("Median_{MC} = %.3f   Median_{Data} = %.3f [MeV/c^{2}]", massResults[0].median, massResults[1].median));
  pt->AddText(Form("RMS_{MC} = %.3f   RMS_{Data} = %.3f [MeV/c^{2}]", massResults[0].rms, massResults[1].rms));
  //pt->AddText(Form("#alpha^{trk}_{smear} = %.3f  #beta^{trk}_{smear} = %.3f", track_scale, track_smearing));
  
  //pt->AddText(Form("#Gamma^{Data}_{#omega}/#Gamma^{MC}_{#omega} = %.3f #pm %.3f", width_ratio, width_ratio_err));
  pt->Draw();

  TLegend *leg = new TLegend(0.2, 0.7, 0.7, 0.8);
  leg->SetTextFont(132);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hist_data_sub, "Data - Background", "lep");
  leg->AddEntry(hist_signal, "Signal MC", "l");
  leg->Draw();

  c1->SaveAs(out_dir + "/data_mc_comparison_bkg_sub_ppIM.pdf");
  cout << "Comparison plot saved to: " << out_dir << "/data_mc_comparison_bkg_sub_ppIM.pdf" << endl;
  

  delete c1;
  
  return 0;

}
