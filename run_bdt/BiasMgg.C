// massBias_bdt.C – full code with background subtraction
// Update compr_bdt.sh with corresponding tuning_type

#include "../header_method/method.h"
#include "../header_plot/plot.h"

struct FitResult {
    TString name;
    double mean, mean_err;
    double sigma, sigma_err;
    double chi2_ndf;
    int entries;
};

int BiasMgg(const TString tuning_type = "tuning_false",
            const TString var_nm = "m_gg_bdt",
            const TString var_symb = "M_{#gamma#gamma} [MeV/c^{2}]"
            ) {
  
  const TString tree_file_nm = "../" + tuning_type + "_" + var_nm + "/hist.root";
  const TString out_dir = "../BiasMgg_" + tuning_type;

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  gSystem->mkdir(out_dir, kTRUE);

  TFile* tree_file = new TFile(tree_file_nm);
  if (!tree_file || tree_file->IsZombie()) {
    cerr << "ERROR: Cannot open " << tree_file_nm << endl;
    return 1;
  }

  checkFile(tree_file);

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

  if (hist_eeg) hist_bkg_sum->Add(hist_eeg, 1.0);
  if (hist_omegapi) hist_bkg_sum->Add(hist_omegapi, 1.0);
  if (hist_nonreson) hist_bkg_sum->Add(hist_nonreson, 1.0);
  if (hist_ksl) hist_bkg_sum->Add(hist_ksl, 1.0);
  if (hist_mcrest) hist_bkg_sum->Add(hist_mcrest, 1.0);

  TH1D *hist_data_sub = (TH1D*)hist_data->Clone("hist_data_sub");
  hist_data_sub->Add(hist_bkg_sum, -1.0);

  // ------------------------------------------------------------------
  // Fit both MC and background‑subtracted data with Gaussian+pol2
  // ------------------------------------------------------------------
  const int nb_mass = 2;
  TH1D *hMassList[nb_mass] = {hist_signal, hist_data_sub};
  TString massNameList[nb_mass] = {"MC", "Data - Background "};
  int massColor[nb_mass] = {kRed, kRed};
  FitResult massResults[nb_mass];

  std::cout << "MC signal entries: " << hist_signal->GetEntries() << std::endl;
  std::cout << "MC signal mean: " << hist_signal->GetMean() << std::endl;
  std::cout << "Data entries: " << hist_data_sub->GetEntries() << std::endl;
  std::cout << "Data mean: " << hist_data_sub->GetMean() << std::endl;

  for (int i = 0; i < nb_mass; i++) {
    TH1D *h_mass = hMassList[i];
    if (!h_mass) {
      std::cerr << "Skipping mass histogram " << massNameList[i] << " (null)." << std::endl;
      continue;
    }

    TCanvas *c_mass = new TCanvas("c_mass_" + massNameList[i], "Mass Distributions (Gaussian+pol2 fits)", 700, 700);
    c_mass->Divide(nb_mass, 1);

    TH1D *h_mass_copy = (TH1D*)h_mass->Clone(Form("h_mass_%s", massNameList[i].Data()));
    h_mass_copy->SetDirectory(0);

    double mass_mean = h_mass_copy->GetMean();
    double mass_rms = h_mass_copy->GetRMS();
    double mass_peak = h_mass_copy->GetBinContent(h_mass_copy->GetMaximumBin());
    double mass_peak_pos = h_mass_copy->GetBinCenter(h_mass_copy->GetMaximumBin());

    double mass_min = h_mass_copy->GetXaxis()->GetXmin();
    double mass_max = h_mass_copy->GetXaxis()->GetXmax();

    // Fit range: ±1σ around mean, but constrained to the histogram range
    double fit_min = TMath::Max(mass_min, mass_mean - mass_rms);
    double fit_max = TMath::Min(mass_max, mass_mean + mass_rms);

    // If the range is too narrow (e.g., near edges), widen to ±2σ
    if (fit_max - fit_min < 0.5 * mass_rms) {
      fit_min = TMath::Max(mass_min, mass_mean - 2.0 * mass_rms);
      fit_max = TMath::Min(mass_max, mass_mean + 2.0 * mass_rms);
    }

    // Safety: if still invalid, fallback to full range
    if (fit_min >= fit_max) {
      fit_min = mass_min;
      fit_max = mass_max;
    }

    std::cout << "\n========================================" << std::endl;
    std::cout << "Mass histogram: " << massNameList[i] << std::endl;
    std::cout << "Number of entries: " << h_mass_copy->GetEntries() << std::endl;
    std::cout << "Mean estimate: " << mass_mean << ", RMS: " << mass_rms << std::endl;
    std::cout << "Fit range: [" << fit_min << ", " << fit_max << "] MeV/c^{2}" << std::endl;

    // ---- Gaussian + 2nd-order polynomial ----
    TF1 *fitFunc = new TF1(Form("fit_%s", massNameList[i].Data()), "gaus(0)+pol2(3)", fit_min, fit_max);
    double amp_g = mass_peak * 0.9;
    double sigma_guess = mass_rms * 0.7;
    fitFunc->SetParameters(amp_g, mass_peak_pos, sigma_guess, 0.0, 0.0, 0.0);
    fitFunc->SetParLimits(0, 0, mass_peak * 5);
    fitFunc->SetParLimits(1, fit_min, fit_max);
    fitFunc->SetParLimits(2, 0.1, 20.0);
    fitFunc->SetParLimits(3, -mass_peak*0.5, mass_peak*0.5);
    fitFunc->SetParLimits(4, -mass_peak*0.1, mass_peak*0.1);
    fitFunc->SetParLimits(5, -mass_peak*0.01, mass_peak*0.01);
    fitFunc->SetLineColor(massColor[i]);
    fitFunc->SetLineWidth(2);

    h_mass_copy->Fit(fitFunc, "RQS");

    // Check if fit succeeded (mean within range, positive sigma)
    if (fitFunc->GetParameter(1) <= fit_min || fitFunc->GetParameter(1) >= fit_max ||
        fitFunc->GetParameter(2) <= 0.1) {
      std::cerr << "Warning: fit for " << massNameList[i] << " may have failed. Using histogram mean as fallback." << std::endl;
      massResults[i].mean = mass_mean;
      massResults[i].mean_err = h_mass_copy->GetMeanError();
      massResults[i].sigma = mass_rms;
      massResults[i].sigma_err = 0.0;
      massResults[i].chi2_ndf = 0.0;
      massResults[i].entries = h_mass_copy->GetEntries();
    } else {
      massResults[i].name = massNameList[i];
      massResults[i].mean = fitFunc->GetParameter(1);
      massResults[i].mean_err = fitFunc->GetParError(1);
      massResults[i].sigma = fitFunc->GetParameter(2);
      massResults[i].sigma_err = fitFunc->GetParError(2);
      massResults[i].chi2_ndf = fitFunc->GetChisquare() / fitFunc->GetNDF();
      massResults[i].entries = h_mass_copy->GetEntries();
    }

    // ---- Draw ----
    c_mass->cd(i+1);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);

    h_mass_copy->SetMarkerStyle(20);
    h_mass_copy->SetMarkerSize(0.6);
    h_mass_copy->GetYaxis()->SetTitle("Events");
    h_mass_copy->GetYaxis()->SetRangeUser(0.01, h_mass_copy->GetMaximum() * 1.6);
    h_mass_copy->GetYaxis()->CenterTitle();
    h_mass_copy->GetYaxis()->SetTitleSize(0.05);
    h_mass_copy->GetYaxis()->SetTitleOffset(1.4);
    h_mass_copy->GetYaxis()->SetLabelSize(0.04);
    h_mass_copy->GetXaxis()->SetTitle(var_symb);
    h_mass_copy->GetXaxis()->SetTitleSize(0.05);
    h_mass_copy->GetXaxis()->SetTitleOffset(1.2);
    h_mass_copy->GetXaxis()->SetLabelSize(0.04);
    h_mass_copy->GetXaxis()->CenterTitle();

    if (massNameList[i] == "MC") {
      h_mass_copy->Draw("hist");
    } else {
      h_mass_copy->Draw("E");
    }

    // If the fit was successful, draw it; otherwise we already used fallback
    if (fitFunc->GetParameter(1) > fit_min && fitFunc->GetParameter(1) < fit_max &&
        fitFunc->GetParameter(2) > 0.1) {
      fitFunc->Draw("same");
    }

    TLegend *leg_mass = new TLegend(0.2, 0.75, 0.65, 0.9);
    leg_mass->SetFillStyle(0);
    leg_mass->SetBorderSize(0);
    leg_mass->SetTextSize(0.03);
    leg_mass->AddEntry(h_mass_copy, Form("%s ", massNameList[i].Data()), "lep");
    if (fitFunc->GetParameter(1) > fit_min && fitFunc->GetParameter(1) < fit_max) {
      /*
      leg_mass->AddEntry(fitFunc, Form("Gaus: #mu = %.2f#pm%.2f [MeV], #sigma = %.2f [MeV], #chi^{2}/ndf=%.2f",
                                       massResults[i].mean, massResults[i].mean_err,
                                       massResults[i].sigma, massResults[i].chi2_ndf), "l");*/
      leg_mass->AddEntry(fitFunc, Form("#mu = %.2f#pm%.2f [MeV], #sigma = %.2f [MeV]", massResults[i].mean, massResults[i].mean_err, massResults[i].sigma), "l");
    } else {
      leg_mass->AddEntry((TObject*)0, "No valid fit", "");
    }
    leg_mass->Draw();

    c_mass->Update();
    c_mass->SaveAs(out_dir + "/mass_fit_" + massNameList[i] + ".pdf");
    delete c_mass;
  }

  // ---- Summary ----
  std::cout << "\n========================================" << std::endl;
  std::cout << "Summary of Mass Fit Results (Gaussian+pol2):" << std::endl;
  std::cout << "========================================" << std::endl;
  for (int i = 0; i < nb_mass; i++) {
    if (hMassList[i]) {
      std::cout << Form("%-15s: mean = %6.3f +/-%6.3f, sigma = %6.3f +/-%6.3f, χ²/ndf = %.3f",
                        massResults[i].name.Data(),
                        massResults[i].mean, massResults[i].mean_err,
                        massResults[i].sigma, massResults[i].sigma_err,
                        massResults[i].chi2_ndf) << std::endl;
    }
  }

  double mass_bias = -(massResults[0].mean - massResults[1].mean);
  double mass_bias_err = TMath::Sqrt(TMath::Power(massResults[0].mean_err, 2) + TMath::Power(massResults[1].mean_err, 2));
  double mass_bias_Z = TMath::Abs(mass_bias) / mass_bias_err;

  // ---- correct error propagation for mass ratio ----
  double R0 = massResults[1].mean / massResults[0].mean;
  double R0_err = R0 * TMath::Sqrt(
      TMath::Power(massResults[1].mean_err / massResults[1].mean, 2) +
      TMath::Power(massResults[0].mean_err / massResults[0].mean, 2)
  );
  
  double width_data = massResults[1].sigma;
  double width_data_err = massResults[1].sigma_err;
  
  double width_mc = massResults[0].sigma;
  double width_mc_err = massResults[0].sigma_err;
  
  // ---- correct error propagation for width_ratio ----
  double width_ratio = massResults[1].sigma / massResults[0].sigma;
  double width_ratio_err = width_ratio * TMath::Sqrt(
      TMath::Power(massResults[1].sigma_err / massResults[1].sigma, 2) +
      TMath::Power(massResults[0].sigma_err / massResults[0].sigma, 2)
  );

  // ----------------------------------------------------
  cout << "mass bias = " << mass_bias << " +/- " << mass_bias_err << "\n"
       << "R0 = " << R0 << " +/- " << R0_err << "\n"
       << "width_ratio = " << width_ratio << " +/- " << width_ratio_err << "\n";

  // ---- Write residual bias ----
  std::ofstream myfile;
  TString myfile_nm = "../pull_scan/massbias_bdt.txt";
  myfile.open(myfile_nm.Data());
  myfile << "const double mpi0_data = " << massResults[1].mean << ";\n";
  myfile << "const double mpi0_data_err = " << massResults[1].mean_err << ";\n\n";
  myfile << "const double mpi0_mc = " << massResults[0].mean << ";\n";
  myfile << "const double mpi0_mc_err = " << massResults[0].mean_err << ";\n\n";
  myfile << "const double mass_bias = " << mass_bias << ";\n";
  myfile << "const double mass_bias_err = " << mass_bias_err << ";\n";
  myfile << "const double mass_bias_Z = " << mass_bias_Z << ";\n\n";

  myfile << "const double R0 = " << R0 << ";\n";
  myfile << "const double R0_err = " << R0_err << ";\n\n";

  myfile << "const double width_data = " << width_data << ";\n";
  myfile << "const double width_data_err = " << width_data_err << ";\n\n";

  myfile << "const double width_mc = " << width_mc << ";\n";
  myfile << "const double width_mc_err = " << width_mc_err << ";\n\n";
  
  myfile << "const double width_ratio = " << width_ratio << ";\n";
  myfile << "const double width_ratio_err = " << width_ratio_err << ";\n";
  
  myfile.close();

  // ---- Optional: Data/MC comparison plot with background‑subtracted data ----
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

  hist_data_sub->Draw("E0");
  hist_signal->Draw("HIST SAME");

  TPaveText *pt = new TPaveText(0.3, 0.8, 0.85, 0.89, "NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.03);
  pt->SetTextFont(42);
  pt->AddText(Form("m^{MC}_{0}-m^{Data}_{0} = %.3f #pm %.3f [MeV/c^{2}]", -1 * mass_bias, mass_bias_err));
  pt->AddText(Form("#Gamma^{Data}_{0}/#Gamma^{MC}_{0} = %.3f #pm %.3f", width_ratio, width_ratio_err));
  pt->Draw();

  TLegend *leg = new TLegend(0.2, 0.65, 0.6, 0.8);
  leg->SetTextFont(132);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(hist_data_sub, "Data - Background", "lep");
  leg->AddEntry(hist_signal, "Signal", "l");
  leg->Draw();

  c1->SaveAs(out_dir + "/data_mc_comparison_bkg_sub_mgg.pdf");
  cout << "Comparison plot saved to: " << out_dir << "/data_mc_comparison_bkg_sub_mgg.pdf" << endl;
  //cout << "Parameter mpi0 saved to: " << myfile_nm << "\n";

  delete c1;
  delete leg;
  delete pt;
  tree_file->Close();
  delete tree_file;

  return 0;
}
