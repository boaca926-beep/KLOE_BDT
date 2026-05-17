#include "../header_bdt/compr.h"
#include "../header_plot/plot.h"
#include "../header_method/method.h"
#include "../header_bdt/sfw2d.txt"   // provides eeg_sfw, isr3pi_sfw, etc.

int plot_compr() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.5);            // makes error bar caps visible
  TH1::SetDefaultSumw2();

  // Open the file created by compr_bdt
  TString filename = "../output_" + var_nm + "/hist_" + var_nm + ".root";
  TFile* intree = new TFile(filename);
  if (!intree || intree->IsZombie()) {
    std::cerr << "ERROR: Cannot open input file " << filename << std::endl;
    return 1;
  }

  // Retrieve the TList (written by compr_bdt)
  TList *Hlist = (TList*)intree->Get("Hlist");
  if (!Hlist) {
    std::cerr << "ERROR: Hlist not found in file" << std::endl;
    return 1;
  }

  // Helper to get histograms from the list
  auto getHist = [&](const char* name) -> TH1D* {
    TH1D* h = (TH1D*)Hlist->FindObject(name);
    if (!h) std::cerr << "WARNING: histogram " << name << " not found" << std::endl;
    return h;
  };

  // Get the histograms (raw and scaled versions)
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

  // Force bin errors to sqrt(content) for data (already set, but safe)
  for (int bin = 1; bin <= hist_data->GetNbinsX(); ++bin) {
    double content = hist_data->GetBinContent(bin);
    hist_data->SetBinError(bin, TMath::Sqrt(content));
  }

  // Apply scaling factors (from sfw2d.txt)
  hist_eeg_sc->Scale(eeg_sfw);
  hist_isr3pi_sc->Scale(isr3pi_sfw);
  hist_omegapi_sc->Scale(omegapi_sfw);
  hist_etagam_sc->Scale(etagam_sfw);
  hist_ksl_sc->Scale(ksl_sfw);
  hist_mcrest_sc->Scale(mcrest_sfw);

  // Recompute background sum (eeg + all backgrounds except isr3pi)
  TH1D *hist_bkgsum_scaled = (TH1D*)hist_eeg_sc->Clone();
  hist_bkgsum_scaled->Add(hist_omegapi_sc, 1.);
  hist_bkgsum_scaled->Add(hist_ksl_sc, 1.);
  hist_bkgsum_scaled->Add(hist_etagam_sc, 1.);
  hist_bkgsum_scaled->Add(hist_mcrest_sc, 1.);
  hist_bkgsum_scaled->SetName("hist_bkgsum_scaled");

  // Total MC sum (background + signal)
  TH1D *hist_mcsum_sc = (TH1D*)hist_bkgsum_scaled->Clone();
  hist_mcsum_sc->Add(hist_isr3pi_sc, 1.);
  hist_mcsum_sc->SetName("hist_mcsum_sc");
  format_h(hist_mcsum_sc, 1, 2);          // black solid line

  // Style the component histograms
  hist_isr3pi_sc->SetLineStyle(2);   // dashed
  hist_etagam_sc->SetLineStyle(3);   // dotted
  hist_ksl_sc->SetLineStyle(4);      // dash-dotted
  hist_omegapi_sc->SetLineStyle(5);  // long-dashed
  hist_mcrest_sc->SetLineStyle(6);   // double-dashed
  hist_eeg_sc->SetLineStyle(7);      // dash-double-dotted

  // Prepare data histogram style
  hist_data->SetMarkerStyle(21);
  hist_data->SetMarkerSize(0.7);
  //hist_data->SetLineWidth(0);

  // --- Residuals (pulls per bin) ---
  const double ymax = hist_data->GetMaximum();
  TH1D *hresidul = new TH1D("hresidul", "", binsize, var_min, var_max);
  TH1D *hresidul_distr = new TH1D("hresidul_distr", "", 200, -10, 10);

  double nb_data = 0., nb_mcsum = 0., residul = 0.;

  for (int j = 1; j <= binsize; ++j) {
    nb_data = hist_data->GetBinContent(j);
    nb_mcsum = hist_mcsum_sc->GetBinContent(j);
    double evnt_err = TMath::Sqrt(nb_data + nb_mcsum);   // Poisson approximation
    if (evnt_err > 0) {
      residul = (nb_data - nb_mcsum) / evnt_err;
      hresidul->SetBinContent(j, residul);
      hresidul_distr->Fill(residul);
    }
  }

  // --- Create canvas with two pads (main plot and residuals) ---
  TCanvas *cv = new TCanvas("cv", "Comparison", 800, 800);
  cv->SetBottomMargin(0.12);
  cv->SetLeftMargin(0.12);

  // Upper pad (main plot)
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
  pad1->SetBottomMargin(0.01);
  pad1->SetLeftMargin(0.12);
  pad1->Draw();
  pad1->cd();

  // Draw data with error bars (E1 ensures markers + error bars)
  hist_data->Draw("E");
  hist_mcsum_sc->Draw("SameHist");
  // Draw components
  hist_isr3pi_sc->Draw("SameHist");
  hist_omegapi_sc->Draw("SameHist");
  hist_etagam_sc->Draw("SameHist");
  hist_ksl_sc->Draw("SameHist");
  hist_mcrest_sc->Draw("SameHist");
  hist_eeg_sc->Draw("SameHist");

  // Axis formatting for upper pad
  hist_data->GetXaxis()->SetTitle("");
  hist_data->GetYaxis()->SetTitle("Events");
  hist_data->GetYaxis()->CenterTitle();
  hist_data->GetYaxis()->SetTitleSize(0.05);
  hist_data->GetYaxis()->SetTitleOffset(1.2);
  hist_data->GetYaxis()->SetLabelSize(0.04);
  // Allow extra headroom for error bars
  hist_data->GetYaxis()->SetRangeUser(0.01, ymax * 1.6);

  // TPaveText label
  TPaveText *pt34 = new TPaveText(0.2, 0.75, 0.4, 0.8, "NDC");
  PteAttr(pt34);
  pt34->SetTextSize(0.1);
  pt34->AddText("(a)");
  pt34->Draw("Same");

  // Legend
  TLegend *legd_cv = new TLegend(0.65, 0.45, 0.9, 0.9);
  legd_cv->SetTextFont(132);
  legd_cv->SetFillStyle(0);
  legd_cv->SetBorderSize(0);
  legd_cv->SetNColumns(1);
  legd_cv->AddEntry(hist_data, "Data", "lep");
  legd_cv->AddEntry(hist_mcsum_sc, "MC sum", "l");
  legd_cv->AddEntry(hist_isr3pi_sc, "#pi^{+}#pi^{-}#pi^{0}#gamma", "l");
  legd_cv->AddEntry(hist_omegapi_sc, "#omega#pi^{0}", "l");
  legd_cv->AddEntry(hist_etagam_sc, "#eta#gamma", "l");
  legd_cv->AddEntry(hist_ksl_sc, "K_{L}K_{S}", "l");
  legd_cv->AddEntry(hist_eeg_sc, "e^{+}e^{-}#gamma", "l");
  legd_cv->AddEntry(hist_mcrest_sc, "Others", "l");
  legd_cv->Draw("Same");
  legtextsize(legd_cv, 0.04);

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
  hresidul->SetLineWidth(0);
  hresidul->GetXaxis()->SetTitle(var_symb + " " + unit);
  hresidul->GetXaxis()->SetTitleSize(0.12);
  hresidul->GetXaxis()->SetTitleOffset(1.0);
  hresidul->GetXaxis()->SetLabelSize(0.1);
  hresidul->GetYaxis()->SetTitle("Pull");
  hresidul->GetYaxis()->SetTitleSize(0.12);
  hresidul->GetYaxis()->SetTitleOffset(0.5);
  hresidul->GetYaxis()->SetLabelSize(0.1);
  hresidul->GetYaxis()->SetRangeUser(-5, 5);
  hresidul->Draw("P");

  // Line at zero
  TLine *line = new TLine(var_min, 0, var_max, 0);
  line->SetLineStyle(2);
  line->Draw();

  // Save the main plot
  cv->SaveAs("../output_" + var_nm + "/cv_compr_" + var_nm + ".pdf");

  // Save residual distribution (optional)
  TCanvas *cv_res = new TCanvas("cv_res", "Residual Distribution", 700, 500);
  hresidul_distr->Draw();
  cv_res->SaveAs("../output_" + var_nm + "/residuals_distr_" + var_nm + ".pdf");

  // Cleanup
  delete cv;
  delete cv_res;
  intree->Close();
  delete intree;

  return 0;
}
