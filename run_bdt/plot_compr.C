#include "../header_bdt/compr.h"
#include "../header_plot/plot.h"
#include "../header_method/method.h"
#include "../header_bdt/sfw2d_bdt.txt"   // provides eeg_sfw, isr3pi_sfw, etc.

int plot_compr() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  // Open the file created by compr_bdt
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

  TH1D *hist_data        = getHist("hist_data");
  TH1D *hist_eeg_sc      = getHist("hist_eeg_sc");
  TH1D *hist_isr3pi_sc   = getHist("hist_isr3pi_sc");
  TH1D *hist_omegapi_sc  = getHist("hist_omegapi_sc");
  TH1D *hist_etagam_sc   = getHist("hist_etagam_sc");
  TH1D *hist_ksl_sc      = getHist("hist_ksl_sc");
  TH1D *hist_mcrest_sc   = getHist("hist_mcrest_sc");
  TH1D *hist_bkgsum_sc   = getHist("hist_bkgsum_sc");

  if (!hist_data) {
    std::cerr << "ERROR: hist_data not found!" << std::endl;
    return 1;
  }

  // --- DEBUG: Print raw integrals before scaling ---
  std::cout << "\n===== Raw integrals (before scaling) =====" << std::endl;
  std::cout << "Data     : " << hist_data->Integral() << std::endl;
  std::cout << "EEG      : " << hist_eeg_sc->Integral() << std::endl;
  std::cout << "ISR3pi   : " << hist_isr3pi_sc->Integral() << std::endl;
  std::cout << "OmegaPi  : " << hist_omegapi_sc->Integral() << std::endl;
  std::cout << "EtaGamma : " << hist_etagam_sc->Integral() << std::endl;
  std::cout << "KSL      : " << hist_ksl_sc->Integral() << std::endl;
  std::cout << "MC Rest  : " << hist_mcrest_sc->Integral() << std::endl;

  // Force bin errors to sqrt(content) for data
  for (int bin = 1; bin <= hist_data->GetNbinsX(); ++bin) {
    double content = hist_data->GetBinContent(bin);
    hist_data->SetBinError(bin, TMath::Sqrt(content));
  }

  // Apply scaling factors
  hist_eeg_sc->Scale(eeg_sfw);
  hist_isr3pi_sc->Scale(isr3pi_sfw);
  hist_omegapi_sc->Scale(omegapi_sfw);
  hist_etagam_sc->Scale(etagam_sfw);
  hist_ksl_sc->Scale(ksl_sfw);
  hist_mcrest_sc->Scale(mcrest_sfw);

  // --- DEBUG: Print scaled integrals ---
  std::cout << "\n===== Scaled integrals =====" << std::endl;
  std::cout << "EEG      : " << hist_eeg_sc->Integral() << std::endl;
  std::cout << "ISR3pi   : " << hist_isr3pi_sc->Integral() << std::endl;
  std::cout << "OmegaPi  : " << hist_omegapi_sc->Integral() << std::endl;
  std::cout << "EtaGamma : " << hist_etagam_sc->Integral() << std::endl;
  std::cout << "KSL      : " << hist_ksl_sc->Integral() << std::endl;
  std::cout << "MC Rest  : " << hist_mcrest_sc->Integral() << std::endl;

  // Recompute background sum
  TH1D *hist_bkgsum_scaled = (TH1D*)hist_eeg_sc->Clone();
  hist_bkgsum_scaled->Add(hist_omegapi_sc, 1.);
  hist_bkgsum_scaled->Add(hist_ksl_sc, 1.);
  hist_bkgsum_scaled->Add(hist_etagam_sc, 1.);
  hist_bkgsum_scaled->Add(hist_mcrest_sc, 1.);
  hist_bkgsum_scaled->SetName("hist_bkgsum_scaled");

  // Total MC sum
  TH1D *hist_mcsum_sc = (TH1D*)hist_bkgsum_scaled->Clone();
  hist_mcsum_sc->Add(hist_isr3pi_sc, 1.);
  hist_mcsum_sc->SetName("hist_mcsum_sc");

  std::cout << "\nTotal MC sum integral: " << hist_mcsum_sc->Integral() << std::endl;
  std::cout << "Data integral after scaling (unchanged): " << hist_data->Integral() << std::endl;

  // If total MC is zero, abort
  if (hist_mcsum_sc->Integral() == 0) {
    std::cerr << "ERROR: Total MC sum is zero. Nothing to plot." << std::endl;
    return 1;
  }

  // Style histograms
  hist_mcsum_sc->SetLineColor(kRed);
  hist_mcsum_sc->SetLineWidth(2);
  hist_mcsum_sc->SetLineStyle(kSolid);

  hist_isr3pi_sc->SetLineStyle(kDashed);
  hist_etagam_sc->SetLineStyle(kDotted);
  hist_ksl_sc->SetLineStyle(kDashDotted);
  hist_omegapi_sc->SetLineStyle(kDashed);
  hist_mcrest_sc->SetLineStyle(kDotted);
  hist_eeg_sc->SetLineStyle(kDashDotted);

  hist_data->SetMarkerStyle(21);
  hist_data->SetMarkerSize(0.7);
  hist_data->SetLineWidth(1);

  // --- Residuals ---
  const double ymax = hist_data->GetMaximum();
  TH1D *hresidul = new TH1D("hresidul", "", binsize, var_min, var_max);
  TH1D *hresidul_distr = new TH1D("hresidul_distr", "", 200, -10, 10);

  for (int j = 1; j <= binsize; ++j) {
    double nb_data = hist_data->GetBinContent(j);
    double nb_mcsum = hist_mcsum_sc->GetBinContent(j);
    double evnt_err = TMath::Sqrt(nb_data + nb_mcsum);
    if (evnt_err > 0) {
      double residul = (nb_data - nb_mcsum) / evnt_err;
      hresidul->SetBinContent(j, residul);
      hresidul_distr->Fill(residul);
    }
  }

  // --- Canvas ---
  TCanvas *cv = new TCanvas("cv", "Comparison", 800, 800);
  cv->SetBottomMargin(0.12);
  cv->SetLeftMargin(0.12);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
  pad1->SetBottomMargin(0.01);
  pad1->SetLeftMargin(0.12);
  pad1->Draw();
  pad1->cd();

  // Draw components (corrected draw options)
  //hist_data->SetMinimum(0);
  hist_data->Draw("E1");
  hist_isr3pi_sc->Draw("hist same");
  hist_omegapi_sc->Draw("hist same");
  hist_etagam_sc->Draw("hist same");
  hist_ksl_sc->Draw("hist same");
  hist_mcrest_sc->Draw("hist same");
  hist_eeg_sc->Draw("hist same");
  //hist_mcsum_sc->Draw("hist same");
  //gPad->SetLogy();
  
  // Axis formatting
  hist_data->GetXaxis()->SetTitle("");
  hist_data->GetYaxis()->SetTitle("Events");
  hist_data->GetYaxis()->CenterTitle();
  hist_data->GetYaxis()->SetTitleSize(0.05);
  hist_data->GetYaxis()->SetTitleOffset(1.2);
  hist_data->GetYaxis()->SetLabelSize(0.04);
  
  // Set y-axis range if ymax > 0, otherwise use auto
  if (ymax > 0)
    hist_data->GetYaxis()->SetRangeUser(0.01, ymax * 1.6);
  else
    hist_data->GetYaxis()->SetRangeUser(0.01, 1);  // fallback, but will show empty

  // Legend
  TLegend *legd_cv = new TLegend(0.65, 0.35, 0.9, 0.9);
  //TLegend *legd_cv = new TLegend(0.15, 0.35, 0.6, 0.9);
  
  legd_cv->SetTextFont(132);
  legd_cv->SetFillStyle(0);
  legd_cv->SetBorderSize(0);
  legd_cv->SetNColumns(1);
  legd_cv->AddEntry(hist_data, "Data", "lep");
  //legd_cv->AddEntry(hist_mcsum_sc, "MC sum", "l");
  legd_cv->AddEntry(hist_isr3pi_sc, "#pi^{+}#pi^{-}#pi^{0}#gamma", "l");
  legd_cv->AddEntry(hist_omegapi_sc, "#omega#pi^{0}", "l");
  legd_cv->AddEntry(hist_etagam_sc, "#eta#gamma", "l");
  legd_cv->AddEntry(hist_ksl_sc, "K_{L}K_{S}", "l");
  legd_cv->AddEntry(hist_eeg_sc, "e^{+}e^{-}#gamma", "l");
  legd_cv->AddEntry(hist_mcrest_sc, "Others", "l");
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
  hresidul->GetXaxis()->SetTitle(var_symb + " " + unit);
  hresidul->GetXaxis()->SetTitleSize(0.12);
  hresidul->GetXaxis()->SetTitleOffset(1.0);
  hresidul->GetXaxis()->SetLabelSize(0.1);
  hresidul->GetXaxis()->CenterTitle();
  hresidul->GetYaxis()->SetTitle("Pull");
  hresidul->GetYaxis()->SetTitleSize(0.12);
  hresidul->GetYaxis()->SetTitleOffset(0.5);
  hresidul->GetYaxis()->SetLabelSize(0.08);
  hresidul->GetYaxis()->SetRangeUser(-20, 20);
  hresidul->GetYaxis()->SetNdivisions(505);
  hresidul->GetYaxis()->CenterTitle();
  hresidul->Draw("P");

  pad2->cd();
  gPad->SetGrid();        // add this line
  

  TLine *line = new TLine(var_min, 0, var_max, 0);
  line->SetLineStyle(2);
  line->Draw();

  cv->Update();
  cv->SaveAs("../output_" + var_nm + "/cv_compr_" + var_nm + ".pdf");

  TCanvas *cv_res = new TCanvas("cv_res", "Residual Distribution", 700, 500);
  hresidul_distr->Draw();
  cv_res->SaveAs("../output_" + var_nm + "/residuals_distr_" + var_nm + ".pdf");

  delete cv;
  delete cv_res;
  intree->Close();
  delete intree;

  return 0;
}
