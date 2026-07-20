// Resolution – fit to m3pi_diff distribution with inner Gaussian values (normalized)
#include "../header_bdt/plot_resol.h"

void plot_resol() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  // Create output directory if it doesn't exist
  gSystem->Exec("mkdir -p ../plots_resol");

  // Open tree file
  TFile *ftree = TFile::Open(treeFile);
  if (!ftree || ftree->IsZombie()) {
    std::cerr << "ERROR: cannot open " << treeFile << std::endl;
    return;
  }

  TTree *ttree = (TTree*) ftree->Get("TISR3PI_SIG_PEAK");
  
  if (!ttree) { std::cerr << "No such tree.\n"; return; }

  // Histogram for difference: (reco – true) mass
  TH1D *h_diff = new TH1D("h_diff", "", bin_size, XMIN, XMAX);
  h_diff->SetMarkerStyle(21);
  h_diff->SetMarkerSize(0.7);
  h_diff->SetLineColor(1);
  h_diff->Sumw2();

  double var = 0.0, var_true = 0.0;
  double var_diff = 0.0;
  //int recon_indx = -1, bkg_indx = -1; //No need, purity conditions are applied at the level of selection
  
  //ttree->SetBranchAddress("Br_recon_indx_bdt", &recon_indx);
  //ttree->SetBranchAddress("Br_bkg_indx", &bkg_indx);
  ttree->SetBranchAddress(var_type, &var);
  ttree->SetBranchAddress(var_type_true, &var_true);
  
  for (Long64_t i = 0; i < ttree->GetEntries(); ++i) {
    ttree->GetEntry(i);
    //if (recon_indx == 2 && bkg_indx == 1) {
    var_diff = var - var_true;
    h_diff->Fill(var_diff);
      //}
  }

  // --- Normalize the histogram (integral = 1) ---
  double integral = h_diff->Integral();
  if (integral > 0) {
    h_diff->Scale(1.0 / integral);
  } else {
    std::cerr << "ERROR: histogram has zero entries." << std::endl;
    return;
  }

  // Fit range: mean ± factor * RMS
  double mean = 0.; //h_diff->GetMean();
  double rms  = h_diff->GetRMS();
  double fit_min = mean - fit_factor * rms;
  double fit_max = mean + fit_factor * rms;
  double hist_min = h_diff->GetXaxis()->GetXmin();
  double hist_max = h_diff->GetXaxis()->GetXmax();
  if (fit_min < hist_min) fit_min = hist_min;
  if (fit_max > hist_max) fit_max = hist_max;
  std::cout << "Fit range for normalized var_diff: [" << fit_min << ", " << fit_max << "] " << unit << std::endl;

  double peak = h_diff->GetMaximum();

  // ------------------------------------------------------------------
  // 1. Fit single Gaussian to get baseline χ²
  // ------------------------------------------------------------------
  TF1 *singleGaus = new TF1("singleGaus", "gaus", fit_min, fit_max);
  singleGaus->SetParameters(peak, mean, rms * 0.7);
  singleGaus->SetParLimits(2, 0.1, 30.0);
  singleGaus->SetLineColor(kRed);
  singleGaus->SetLineWidth(2);
  h_diff->Fit(singleGaus, "RQS");
  double chi2_single = singleGaus->GetChisquare();
  int ndf_single = singleGaus->GetNDF();

  // ------------------------------------------------------------------
  // 2. Fit double Gaussian
  // ------------------------------------------------------------------
  TF1 *doubleGaus = new TF1("doubleGaus", "gaus(0)+gaus(3)", fit_min, fit_max);
  doubleGaus->SetParameter(0, peak * 0.8);
  doubleGaus->SetParameter(1, mean);
  doubleGaus->SetParameter(2, rms * 0.7);
  doubleGaus->SetParameter(3, peak * 0.2);
  doubleGaus->SetParameter(4, mean + 2.0);
  doubleGaus->SetParameter(5, rms * 1.2);
  doubleGaus->SetParLimits(0, 0.0, peak * 3.0);
  doubleGaus->SetParLimits(1, mean - 3.0, mean + 3.0);
  doubleGaus->SetParLimits(2, 0.1, 30.0);
  doubleGaus->SetParLimits(3, 0.0, peak * 2.0);
  doubleGaus->SetParLimits(4, mean - 5.0, mean + 5.0);
  doubleGaus->SetParLimits(5, 0.1, 30.0);
  doubleGaus->SetLineColor(kRed);
  doubleGaus->SetLineWidth(2);
  h_diff->Fit(doubleGaus, "RQS");

  double chi2_double = doubleGaus->GetChisquare();
  int ndf_double = doubleGaus->GetNDF();

  double amp1 = doubleGaus->GetParameter(0);
  double amp2 = doubleGaus->GetParameter(3);
  double sigma1 = doubleGaus->GetParameter(2);
  double sigma2 = doubleGaus->GetParameter(5);

  std::cout << "Double Gaussian parameters:" << std::endl;
  std::cout << "  amp1 = " << amp1 << " +/- " << doubleGaus->GetParError(0) << std::endl;
  std::cout << "  mean1 = " << doubleGaus->GetParameter(1) << " +/- " << doubleGaus->GetParError(1) << std::endl;
  std::cout << "  sigma1 = " << sigma1 << " +/- " << doubleGaus->GetParError(2) << std::endl;
  std::cout << "  amp2 = " << amp2 << " +/- " << doubleGaus->GetParError(3) << std::endl;
  std::cout << "  mean2 = " << doubleGaus->GetParameter(4) << " +/- " << doubleGaus->GetParError(4) << std::endl;
  std::cout << "  sigma2 = " << sigma2 << " +/- " << doubleGaus->GetParError(5) << std::endl;
  std::cout << "  χ²/ndf = " << doubleGaus->GetChisquare() / doubleGaus->GetNDF() << std::endl;
  
  
  // ------------------------------------------------------------------
  // 3. Decide which fit to use
  // ------------------------------------------------------------------
  double chi2_improvement = chi2_single - chi2_double;
  double amp_ratio = TMath::Max(amp1, amp2) / TMath::Min(amp1, amp2);
  double sigma_ratio = TMath::Max(sigma1, sigma2) / TMath::Min(sigma1, sigma2);
  bool distinct = (amp_ratio > 1.2) && (sigma_ratio > 1.1);

  cout << "chi2_improvement = " << chi2_improvement << ", amp_ratio = " << amp_ratio << ", sigma_ratio = " << sigma_ratio << ", chi2_improvement / chi2_single = " << chi2_improvement / chi2_single << endl;
  
  bool double_better = (chi2_improvement > 10.0) && (chi2_improvement / chi2_single > 0.05);

  TF1 *finalFit;
  TString fitType;
  bool doubleUsed;
  double chi2ndf;

  if (double_better && distinct) {
  //if (distinct) { // only double gaussian
    std::cout << "Double Gaussian is significantly better, keeping it." << std::endl;
    finalFit = doubleGaus;
    fitType = "Double Gaussian";
    doubleUsed = true;
    chi2ndf = chi2_double / ndf_double;
    delete singleGaus;  // clean up
  } else {
    std::cout << "Double Gaussian not significantly better, using single Gaussian." << std::endl;
    finalFit = singleGaus;
    fitType = "Single Gaussian";
    doubleUsed = false;
    chi2ndf = chi2_single / ndf_single;
    delete doubleGaus;
  }

  // Draw histogram and fit
  TCanvas *c1 = new TCanvas("c1", "var_diff resolution (normalized)", 700, 700);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);

  double ymax = h_diff -> GetMaximum();
  
  h_diff->SetLineWidth(2);
  h_diff->SetLineColor(4);
  h_diff->GetXaxis()->SetTitle(x_title);
  h_diff->GetYaxis()->SetTitle("Normalized Entries");
  h_diff->GetXaxis()->CenterTitle();
  h_diff->GetYaxis()->CenterTitle();
  h_diff->GetXaxis()->SetTitleSize(0.05);
  h_diff->GetYaxis()->SetTitleSize(0.06);
  h_diff->GetXaxis()->SetLabelSize(0.05);
  h_diff->GetYaxis()->SetLabelSize(0.05);
  h_diff->GetYaxis()->SetTitleOffset(1.3);
  h_diff->GetYaxis()->SetRangeUser(0., 1.6 * ymax);
  // ============================================================

  h_diff->GetYaxis()->SetNdivisions(505);
  h_diff->GetXaxis()->SetNdivisions(505);

  h_diff->Draw("hist");
  finalFit->Draw("same");

  // ----- Extract inner (narrower) Gaussian values if double Gaussian -----
  TString line1, line2, line3;
  if (doubleUsed) {
    double sigma1 = finalFit->GetParameter(2);
    double sigma2 = finalFit->GetParameter(5);
    if (sigma1 < sigma2) {
      double mean_inner = finalFit->GetParameter(1);
      double sigma_inner = sigma1;
      double err_mean = finalFit->GetParError(1);
      double err_sigma = finalFit->GetParError(2);
      line1 = Form("Inner Gaussian (core):");
      line2 = Form("#mu = %.2f#pm%.3f " + unit, mean_inner, err_mean);
      line3 = Form("#sigma = %.2f#pm%.3f " + unit, sigma_inner, err_sigma);
    } else {
      double mean_inner = finalFit->GetParameter(4);
      double sigma_inner = sigma2;
      double err_mean = finalFit->GetParError(4);
      double err_sigma = finalFit->GetParError(5);
      line1 = Form("Inner Gaussian:");
      line2 = Form("#mu = %.2f#pm%.3f " + unit, mean_inner, err_mean);
      line3 = Form("#sigma = %.2f#pm%.3f " + unit, sigma_inner, err_sigma);
    }
  } else {
    // Single Gaussian
    double mean_sg = finalFit->GetParameter(1);
    double sigma_sg = finalFit->GetParameter(2);
    double err_mean = finalFit->GetParError(1);
    double err_sigma = finalFit->GetParError(2);
    line1 = Form("Gaussian fit:");
    line2 = Form("#mu = %.2e#pm%.2e %s ", mean_sg, err_mean, unit.Data());
    line3 = Form("#sigma = %.2e#pm%.2e %s ", sigma_sg, err_sigma, unit.Data());
  }
  TString line4 = Form("#chi^{2}/NDF = %.2f", chi2ndf);

  TPaveText *pt = new TPaveText(0.3, 0.72, 0.9, 0.88, "NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.04);
  pt->SetTextFont(42);
  if (doubleUsed) pt->AddText(line1);
  pt->AddText(line2);
  pt->AddText(line3);
  pt->Draw();

  TLegend *leg = new TLegend(0.55, 0.60, 0.9, 0.68);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(finalFit, fitType, "l");
  leg->Draw();

  // Save canvas
  c1->SaveAs("../plots_resol/" + var_type + "_diff_fit_normalized.png");
  std::cout << "\nPlot saved to ../plots_resol/" + var_type + "_diff_fit_normalized.png" << std::endl;
}
