// comparison for backgrounds including etagam
// Uses the same 1D histograms from compr_bdt, scaled with fit results
#include "../header_bdt/compr.h"
#include "../header_plot/plot.h"
#include "../header_method/method.h"
#include "../header_bdt/sfw2d_bdt.txt"   // provides eeg_sfw, isr3pi_sfw, nonReson_sfw, etc.

int plot_compr_roofit() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  // Open the file containing 1D histograms (same as used by compr_bdt)
  TString filename = "../output_" + var_nm + "/hist_" + var_nm + ".root";
  TFile* intree = new TFile(filename);
  if (!intree || intree->IsZombie()) {
    std::cerr << "ERROR: Cannot open input file " << filename << std::endl;
    return 1;
  }

  TList *Hlist = (TList*)intree->Get("Hlist");
  if (!Hlist) {
    std::cerr << "ERROR: Hlist not found in file" << std::endl;
    return 1;
  }

  auto getHist = [&](const char* name) -> TH1D* {
    TH1D* h = (TH1D*)Hlist->FindObject(name);
    if (!h) std::cerr << "WARNING: histogram " << name << " not found" << std::endl;
    return h;
  };

  // Get histograms - these are already 1D projections
  TH1D *hist_data              = getHist("hist_data");
  TH1D *hist_eeg               = getHist("hist_eeg_sc");        // EEG template
  TH1D *hist_isr3pi_peak       = getHist("hist_isr3pi_peak_sc"); // ISR3pi peak template
  TH1D *hist_isr3pi_nonReson   = getHist("hist_isr3pi_nonReson_sc"); // Non-resonant template
  TH1D *hist_omegapi           = getHist("hist_omegapi_sc");     // OmegaPi template
  TH1D *hist_etagam            = getHist("hist_etagam_sc");      // EtaGamma template
  TH1D *hist_ksl               = getHist("hist_ksl_sc");         // KSL template
  TH1D *hist_mcrest            = getHist("hist_mcrest_sc");      // MC Rest (kpm + rhopi + bkgrest) NO etagam!
  TH1D *hist_bkgsum            = getHist("hist_bkgsum_sc");      // Background sum (from compr_bdt)

  if (!hist_data) {
    std::cerr << "ERROR: hist_data not found!" << std::endl;
    return 1;
  }

  // Print raw integrals before scaling
  std::cout << "\n===== Raw integrals (before scaling) =====" << std::endl;
  std::cout << "Data                : " << hist_data->Integral() << std::endl;
  if (hist_eeg) std::cout << "EEG                 : " << hist_eeg->Integral() << std::endl;
  if (hist_isr3pi_peak) std::cout << "ISR3pi peak         : " << hist_isr3pi_peak->Integral() << std::endl;
  if (hist_isr3pi_nonReson) std::cout << "ISR3pi non-resonant : " << hist_isr3pi_nonReson->Integral() << std::endl;
  if (hist_omegapi) std::cout << "OmegaPi             : " << hist_omegapi->Integral() << std::endl;
  if (hist_etagam) std::cout << "EtaGamma            : " << hist_etagam->Integral() << std::endl;
  if (hist_ksl) std::cout << "KSL                 : " << hist_ksl->Integral() << std::endl;
  if (hist_mcrest) std::cout << "MC Rest             : " << hist_mcrest->Integral() << std::endl;

  // Force bin errors to sqrt(content) for data
  for (int bin = 1; bin <= hist_data->GetNbinsX(); ++bin) {
    double content = hist_data->GetBinContent(bin);
    hist_data->SetBinError(bin, TMath::Sqrt(content));
  }

  // ------------------------------------------------------------------
  // IMPORTANT: mcrest in the fit INCLUDES etagam!
  // So we need to add etagam to mcrest for consistency with the fit
  // ------------------------------------------------------------------
  if (hist_mcrest && hist_etagam) {
    hist_mcrest->Add(hist_etagam, 1.0);
    std::cout << "\nNOTE: Added etagam to mcrest (to match fit definition)" << std::endl;
  }

  // Apply scaling factors from the fit
  if (hist_eeg) hist_eeg->Scale(eeg_sfw);
  if (hist_isr3pi_peak) hist_isr3pi_peak->Scale(isr3pi_sfw);
  if (hist_isr3pi_nonReson) hist_isr3pi_nonReson->Scale(nonReson_sfw);
  if (hist_omegapi) hist_omegapi->Scale(omegapi_sfw);
  if (hist_ksl) hist_ksl->Scale(ksl_sfw);
  if (hist_mcrest) hist_mcrest->Scale(mcrest_sfw);  // mcrest now includes etagam

  // Print scaled integrals
  std::cout << "\n===== Scaled integrals (after applying fit scaling factors) =====" << std::endl;
  if (hist_eeg) std::cout << "EEG                 : " << hist_eeg->Integral() << std::endl;
  if (hist_isr3pi_peak) std::cout << "ISR3pi peak         : " << hist_isr3pi_peak->Integral() << std::endl;
  if (hist_isr3pi_nonReson) std::cout << "ISR3pi non-resonant : " << hist_isr3pi_nonReson->Integral() << std::endl;
  if (hist_omegapi) std::cout << "OmegaPi             : " << hist_omegapi->Integral() << std::endl;
  if (hist_ksl) std::cout << "KSL                 : " << hist_ksl->Integral() << std::endl;
  if (hist_mcrest) std::cout << "MC Rest (incl etagam): " << hist_mcrest->Integral() << std::endl;

  // Recompute background sum (all backgrounds, etagam now in mcrest)
  TH1D *hist_bkgsum_scaled = nullptr;
  if (hist_eeg) {
    hist_bkgsum_scaled = (TH1D*)hist_eeg->Clone("hist_bkgsum_scaled");
    if (hist_omegapi) hist_bkgsum_scaled->Add(hist_omegapi, 1.);
    if (hist_ksl) hist_bkgsum_scaled->Add(hist_ksl, 1.);
    if (hist_mcrest) hist_bkgsum_scaled->Add(hist_mcrest, 1.);
  }

  // Total MC sum (signal + backgrounds)
  TH1D *hist_mcsum_sc = nullptr;
  if (hist_isr3pi_peak) {
    hist_mcsum_sc = (TH1D*)hist_isr3pi_peak->Clone("hist_mcsum_sc");
    if (hist_isr3pi_nonReson) hist_mcsum_sc->Add(hist_isr3pi_nonReson, 1.);
    if (hist_bkgsum_scaled) hist_mcsum_sc->Add(hist_bkgsum_scaled, 1.);
  }

  if (!hist_mcsum_sc) {
    std::cerr << "ERROR: Could not create total MC sum histogram" << std::endl;
    return 1;
  }

  std::cout << "\nTotal MC sum integral: " << hist_mcsum_sc->Integral() << std::endl;
  std::cout << "Data integral (unchanged): " << hist_data->Integral() << std::endl;

  if (hist_mcsum_sc->Integral() == 0) {
    std::cerr << "ERROR: Total MC sum is zero. Nothing to plot." << std::endl;
    return 1;
  }

  // Style histograms
  hist_mcsum_sc->SetLineColor(kBlack);
  hist_mcsum_sc->SetLineWidth(3);
  hist_mcsum_sc->SetLineStyle(kSolid);
  hist_mcsum_sc->SetFillStyle(0);

  if (hist_isr3pi_peak) {
    hist_isr3pi_peak->SetLineColor(kRed);
    hist_isr3pi_peak->SetLineWidth(2);
    hist_isr3pi_peak->SetLineStyle(kSolid);
  }
  
  if (hist_isr3pi_nonReson) {
    hist_isr3pi_nonReson->SetLineColor(kOrange);
    hist_isr3pi_nonReson->SetLineWidth(2);
    hist_isr3pi_nonReson->SetLineStyle(kDashed);
  }
  
  if (hist_omegapi) {
    hist_omegapi->SetLineColor(kGreen);
    hist_omegapi->SetLineWidth(2);
    hist_omegapi->SetLineStyle(kDashed);
  }
  
  if (hist_ksl) {
    hist_ksl->SetLineColor(kCyan);
    hist_ksl->SetLineWidth(2);
    hist_ksl->SetLineStyle(kDashDotted);
  }
  
  if (hist_eeg) {
    hist_eeg->SetLineColor(kBlue);
    hist_eeg->SetLineWidth(2);
    hist_eeg->SetLineStyle(kDashed);
  }
  
  if (hist_mcrest) {
    hist_mcrest->SetLineColor(kGray+2);
    hist_mcrest->SetLineWidth(2);
    hist_mcrest->SetLineStyle(kDotted);
  }

  hist_data->SetMarkerStyle(20);
  hist_data->SetMarkerSize(0.8);
  hist_data->SetLineWidth(1);
  hist_data->SetMarkerColor(kBlack);

  // --- Residuals ---
  const double ymax = hist_data->GetMaximum();
  TH1D *hresidul = new TH1D("hresidul", "", binsize, var_min, var_max);
  TH1D *hresidul_distr = new TH1D("hresidul_distr", "", 200, -10, 10);

  for (int j = 1; j <= binsize; ++j) {
    double nb_data = hist_data->GetBinContent(j);
    double nb_mcsum = hist_mcsum_sc->GetBinContent(j);
    double nb_data_err = hist_data->GetBinError(j);
    double nb_mcsum_err = hist_mcsum_sc->GetBinError(j);
    double evnt_err = TMath::Sqrt(nb_data_err*nb_data_err + nb_mcsum_err*nb_mcsum_err);
    if (evnt_err > 0) {
      double residul = (nb_data - nb_mcsum) / evnt_err;
      hresidul->SetBinContent(j, residul);
      hresidul_distr->Fill(residul);
    }
  }

  // --- Canvas ---
  TCanvas *cv = new TCanvas("cv", "Comparison", 900, 900);
  cv->SetBottomMargin(0.12);
  cv->SetLeftMargin(0.12);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
  pad1->SetBottomMargin(0.01);
  pad1->SetLeftMargin(0.12);
  pad1->Draw();
  pad1->cd();

  // Draw components (backgrounds first, then signal, then data on top)
  if (hist_mcrest) hist_mcrest->Draw("hist");
  if (hist_eeg) hist_eeg->Draw("hist same");
  if (hist_omegapi) hist_omegapi->Draw("hist same");
  if (hist_ksl) hist_ksl->Draw("hist same");
  if (hist_isr3pi_nonReson) hist_isr3pi_nonReson->Draw("hist same");
  if (hist_isr3pi_peak) hist_isr3pi_peak->Draw("hist same");
  hist_mcsum_sc->Draw("hist same");
  hist_data->Draw("E1 same");

  // Axis formatting
  const double limit_factor = 0.5;
  hist_data->GetXaxis()->SetRangeUser(0, var_max * limit_factor);
  hist_data->GetXaxis()->SetTitle("");
  hist_data->GetYaxis()->SetTitle("Events / (0.005 GeV^{2})");
  hist_data->GetYaxis()->CenterTitle();
  hist_data->GetYaxis()->SetTitleSize(0.05);
  hist_data->GetYaxis()->SetTitleOffset(1.2);
  hist_data->GetYaxis()->SetLabelSize(0.04);
  
  if (ymax > 0)
    hist_data->GetYaxis()->SetRangeUser(0.01, ymax * 1.6);
  else
    hist_data->GetYaxis()->SetRangeUser(0.01, 1);

  // Legend
  TLegend *legd_cv = new TLegend(0.65, 0.35, 0.9, 0.9);
  legd_cv->SetTextFont(132);
  legd_cv->SetFillStyle(0);
  legd_cv->SetBorderSize(0);
  legd_cv->SetNColumns(1);
  legd_cv->AddEntry(hist_data, "Data", "lep");
  legd_cv->AddEntry(hist_mcsum_sc, "Total MC", "l");
  legd_cv->AddEntry(hist_isr3pi_peak, "#pi^{+}#pi^{-}#pi^{0}#gamma (peak)", "l");
  legd_cv->AddEntry(hist_isr3pi_nonReson, "#pi^{+}#pi^{-}#pi^{0}#gamma (non-reson)", "l");
  legd_cv->AddEntry(hist_omegapi, "#omega#pi^{0}", "l");
  legd_cv->AddEntry(hist_ksl, "K_{L}K_{S}", "l");
  legd_cv->AddEntry(hist_eeg, "e^{+}e^{-}#gamma", "l");
  legd_cv->AddEntry(hist_mcrest, "Others (incl. #eta#gamma)", "l");
  legd_cv->Draw("Same");

  // Lower pad (residuals)
  cv->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.3);
  pad2->SetLeftMargin(0.12);
  pad2->Draw();
  pad2->cd();

  hresidul->SetMarkerStyle(20);
  hresidul->SetMarkerSize(0.6);
  hresidul->GetXaxis()->SetRangeUser(0, var_max * limit_factor);
  hresidul->GetXaxis()->SetTitle(var_symb + " " + unit);
  hresidul->GetXaxis()->SetTitleSize(0.12);
  hresidul->GetXaxis()->SetTitleOffset(1.0);
  hresidul->GetXaxis()->SetLabelSize(0.1);
  hresidul->GetXaxis()->CenterTitle();
  hresidul->GetYaxis()->SetTitle("Pull");
  hresidul->GetYaxis()->SetTitleSize(0.12);
  hresidul->GetYaxis()->SetTitleOffset(0.5);
  hresidul->GetYaxis()->SetLabelSize(0.08);
  hresidul->GetYaxis()->SetRangeUser(-5, 5);
  hresidul->GetYaxis()->SetNdivisions(505);
  hresidul->GetYaxis()->CenterTitle();
  hresidul->Draw("P");

  pad2->cd();
  gPad->SetGrid();

  TLine *line = new TLine(var_min, 0, var_max * limit_factor, 0);
  line->SetLineStyle(2);
  line->Draw();

  cv->Update();
  cv->SaveAs("../output_" + var_nm + "/cv_compr_rootfit_" + var_nm + ".pdf");

  // Residual distribution canvas
  TCanvas *cv_res = new TCanvas("cv_res", "Residual Distribution", 700, 500);
  hresidul_distr->Draw();
  gStyle->SetOptFit(1);
  hresidul_distr->Fit("gaus", "Q");
  cv_res->SaveAs("../output_" + var_nm + "/residuals_distr_rootfit_" + var_nm + ".pdf");

  // Print summary
  std::cout << "\n===== FIT SUMMARY =====" << std::endl;
  std::cout << "Data events: " << hist_data->Integral() << std::endl;
  std::cout << "Total MC events: " << hist_mcsum_sc->Integral() << std::endl;
  std::cout << "Pull mean: " << hresidul_distr->GetMean() << std::endl;
  std::cout << "Pull sigma: " << hresidul_distr->GetRMS() << std::endl;
  std::cout << "======================\n" << std::endl;

  delete cv;
  delete cv_res;
  intree->Close();
  delete intree;

  return 0;
}
