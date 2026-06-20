// At the beginning of the function, set output formatting
std::cout << std::fixed;  // Fixed notation (not scientific)
std::cout << std::setprecision(6);  // 6 decimal places

// Resolution – fit to m3pi_diff distribution with inner Gaussian values (normalized)
#include "../header_bdt/plot_resol.h"


void plot_resol_better() {

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

  TTree *ttree = (TTree*) ftree->Get("TISR3PI_SIG");
  if (!ttree) { std::cerr << "No such tree.\n"; return; }

  // Histogram for difference: (reco – true) mass
  TH1D *h_diff = new TH1D("h_diff", "", bin_size, XMIN, XMAX);
  h_diff->SetMarkerStyle(21);
  h_diff->SetMarkerSize(0.7);
  h_diff->SetLineColor(1);
  h_diff->Sumw2();

  double var = 0.0, var_true = 0.0;
  double var_diff = 0.0;
  int recon_indx_bdt = -1, bkg_indx = -1;
  
  ttree->SetBranchAddress("Br_recon_indx_bdt", &recon_indx_bdt);
  ttree->SetBranchAddress("Br_bkg_indx", &bkg_indx);
  ttree->SetBranchAddress(var_type, &var);
  ttree->SetBranchAddress(var_type_true, &var_true);
  
  Long64_t nEntries = ttree->GetEntries();
  for (Long64_t i = 0; i < nEntries; ++i) {
    ttree->GetEntry(i);
    if (recon_indx_bdt == 2 && bkg_indx == 1) {
      var_diff = var - var_true;
      h_diff->Fill(var_diff);
    }
  }

  // Check if histogram has entries
  double integral = h_diff->Integral();
  if (integral < 10) {
    std::cerr << "ERROR: histogram has too few entries (" << integral << ")" << std::endl;
    return;
  }
  
  // --- Normalize the histogram (integral = 1) ---
  h_diff->Scale(1.0 / integral);

  // --- Calculate robust statistics ---
  double mean = h_diff->GetMean();
  double rms = h_diff->GetRMS();
  
  // Use more robust initial estimates from median and quantiles
  double median = 0.0;
  double q16 = 0.0, q84 = 0.0;
  
  // Get quantiles for robust sigma estimate
  Double_t quantiles[3] = {0.16, 0.50, 0.84};
  Double_t values[3];
  h_diff->GetQuantiles(3, values, quantiles);
  q16 = values[0];
  median = values[1];
  q84 = values[2];
  double robust_sigma = (q84 - q16) / 2.0; // ~1 sigma for Gaussian
  
  // Use median as better initial mean estimate if distribution is symmetric
  if (TMath::Abs(mean - median) < 2.0 * robust_sigma) {
    mean = median; // Use median if close to mean
  }
  
  // Fit range: mean ± N * sigma
  double fit_factor_robust = fit_factor;
  double fit_min = mean - fit_factor_robust * robust_sigma;
  double fit_max = mean + fit_factor_robust * robust_sigma;
  
  // Clamp to histogram range
  double hist_min = h_diff->GetXaxis()->GetXmin();
  double hist_max = h_diff->GetXaxis()->GetXmax();
  fit_min = TMath::Max(fit_min, hist_min);
  fit_max = TMath::Min(fit_max, hist_max);
  
  std::cout << "\n=== Fit Information ===" << std::endl;
  std::cout << "Mean: " << mean << " " << unit << std::endl;
  std::cout << "RMS: " << rms << " " << unit << std::endl;
  std::cout << "Robust sigma: " << robust_sigma << " " << unit << std::endl;
  std::cout << "Fit range: [" << fit_min << ", " << fit_max << "] " << unit << std::endl;
  std::cout << "Entries in fit range: " << h_diff->Integral(h_diff->FindBin(fit_min), h_diff->FindBin(fit_max)) << std::endl;

  // --- Fitting strategy with stability checks ---
  double peak = h_diff->GetMaximum();
  
  // Try double Gaussian fit with better initial parameters
  TF1 *doubleGaus = new TF1("doubleGaus", "gaus(0)+gaus(3)", fit_min, fit_max);
  
  // Initial parameters: core Gaussian (70%) + tail Gaussian (30%)
  double core_frac = 0.7;
  double tail_frac = 0.3;
  double core_sigma = robust_sigma * 0.6;
  double tail_sigma = robust_sigma * 1.5;
  
  doubleGaus->SetParameter(0, peak * core_frac);
  doubleGaus->SetParameter(1, mean);
  doubleGaus->SetParameter(2, core_sigma);
  doubleGaus->SetParameter(3, peak * tail_frac);
  doubleGaus->SetParameter(4, mean + 2.0);
  doubleGaus->SetParameter(5, tail_sigma);
  
  // Set parameter limits
  doubleGaus->SetParLimits(0, 0.0, peak * 2.0);
  doubleGaus->SetParLimits(1, mean - 5.0 * robust_sigma, mean + 5.0 * robust_sigma);
  doubleGaus->SetParLimits(2, 0.5, 30.0);
  doubleGaus->SetParLimits(3, 0.0, peak * 1.5);
  doubleGaus->SetParLimits(4, mean - 20.0, mean + 20.0);
  doubleGaus->SetParLimits(5, 1.0, 50.0);
  
  doubleGaus->SetLineColor(kRed);
  doubleGaus->SetLineWidth(2);

  // Perform double Gaussian fit
  Int_t fitStatus = h_diff->Fit(doubleGaus, "RQSN"); // R=range, Q=quiet, S=sumw2, N=no drawing
  
  // --- Evaluate fit stability ---
  bool doubleGausStable = false;
  double chi2ndf_double = 0.0;
  
  if (fitStatus == 0) { // Fit converged
    chi2ndf_double = doubleGaus->GetChisquare() / doubleGaus->GetNDF();
    
    // Get fit parameters and errors
    double amp1 = doubleGaus->GetParameter(0);
    double amp2 = doubleGaus->GetParameter(3);
    double err_amp1 = doubleGaus->GetParError(0);
    double err_amp2 = doubleGaus->GetParError(3);
    double sigma1 = doubleGaus->GetParameter(2);
    double sigma2 = doubleGaus->GetParameter(5);
    double err_sigma1 = doubleGaus->GetParError(2);
    double err_sigma2 = doubleGaus->GetParError(5);
    
    // Stability criteria
    bool ampStable = (err_amp1 / (amp1 + 1e-6) < 1.5) && (err_amp2 / (amp2 + 1e-6) < 2.0);
    bool sigmaStable = (err_sigma1 / (sigma1 + 1e-6) < 1.0) && (err_sigma2 / (sigma2 + 1e-6) < 1.0);
    bool positiveAmps = (amp1 > 0) && (amp2 > 0);
    bool positiveSigmas = (sigma1 > 0) && (sigma2 > 0);
    bool reasonableChi2 = (chi2ndf_double < 10.0) && (chi2ndf_double > 0.1);
    bool separatedGaussians = (TMath::Abs(sigma1 - sigma2) / TMath::Min(sigma1, sigma2) > 0.2);
    bool noSaturation = (amp1 + amp2 < peak * 2.0);
    
    std::cout << "\n=== Double Gaussian Stability Checks ===" << std::endl;
    std::cout << "  chi2/NDF: " << chi2ndf_double << std::endl;
    std::cout << "  amp1: " << amp1 << " ± " << err_amp1 << " (err: " << err_amp1/(amp1+1e-6)*100 << "%)" << std::endl;
    std::cout << "  amp2: " << amp2 << " ± " << err_amp2 << " (err: " << err_amp2/(amp2+1e-6)*100 << "%)" << std::endl;
    std::cout << "  sigma1: " << sigma1 << " ± " << err_sigma1 << std::endl;
    std::cout << "  sigma2: " << sigma2 << " ± " << err_sigma2 << std::endl;
    std::cout << "  amp ratio: " << amp1/(amp1+amp2) << std::endl;
    std::cout << "  sigma ratio: " << sigma2/sigma1 << std::endl;
    
    doubleGausStable = ampStable && sigmaStable && positiveAmps && positiveSigmas && 
                        reasonableChi2 && separatedGaussians && noSaturation;
  }
  
  // --- Final fit selection ---
  TF1 *finalFit = nullptr;
  TString fitType = "";
  bool doubleUsed = false;
  double chi2ndf = 0.0;
  
  if (doubleGausStable) {
    // Use double Gaussian
    finalFit = doubleGaus;
    fitType = "Double Gaussian";
    doubleUsed = true;
    chi2ndf = chi2ndf_double;
    std::cout << "✅ Double Gaussian fit is stable - using it." << std::endl;
  } else {
    // Fallback to single Gaussian
    std::cout << "⚠️ Double Gaussian unstable, falling back to single Gaussian." << std::endl;
    
    // Try to get better initial parameters for single Gaussian
    TF1 *singleGaus = new TF1("singleGaus", "gaus", fit_min, fit_max);
    singleGaus->SetParameters(peak, mean, robust_sigma);
    singleGaus->SetParLimits(0, 0.0, peak * 2.0);
    singleGaus->SetParLimits(1, mean - 5.0 * robust_sigma, mean + 5.0 * robust_sigma);
    singleGaus->SetParLimits(2, 0.5, 30.0);
    singleGaus->SetLineColor(kRed);
    singleGaus->SetLineWidth(2);
    
    Int_t singleStatus = h_diff->Fit(singleGaus, "RQSN");
    if (singleStatus == 0) {
      finalFit = singleGaus;
      fitType = "Single Gaussian";
      doubleUsed = false;
      chi2ndf = singleGaus->GetChisquare() / singleGaus->GetNDF();
      std::cout << "✅ Single Gaussian fit succeeded." << std::endl;
    } else {
      // Last resort: use unbinned or simple fit
      std::cout << "⚠️ Single Gaussian also failed. Using simplified fit." << std::endl;
      TF1 *simpleGaus = new TF1("simpleGaus", "gaus", fit_min, fit_max);
      simpleGaus->SetParameters(peak, mean, robust_sigma);
      simpleGaus->SetLineColor(kRed);
      simpleGaus->SetLineWidth(2);
      h_diff->Fit(simpleGaus, "RQSN");
      finalFit = simpleGaus;
      fitType = "Simple Gaussian";
      doubleUsed = false;
      chi2ndf = simpleGaus->GetChisquare() / simpleGaus->GetNDF();
    }
  }

  // --- Draw histogram and fit ---
  TCanvas *c1 = new TCanvas("c1", "var_diff resolution (normalized)", 700, 700);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);

  double ymax = h_diff->GetMaximum();
  
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
  h_diff->GetXaxis()->SetRangeUser(range_factor * fit_min, range_factor * fit_max);
  h_diff->GetYaxis()->SetNdivisions(505);
  h_diff->GetXaxis()->SetNdivisions(505);

  h_diff->Draw("hist");
  if (finalFit) finalFit->Draw("same");

  // --- Extract and display results ---
  TString line1, line2, line3, line4;
  
  if (doubleUsed && finalFit) {
    // Extract inner (narrower) Gaussian values
    double sigma1 = finalFit->GetParameter(2);
    double sigma2 = finalFit->GetParameter(5);
    double err_sigma1 = finalFit->GetParError(2);
    double err_sigma2 = finalFit->GetParError(5);
    
    // Choose the narrower sigma as the core resolution
    if (sigma1 < sigma2) {
      double mean_inner = finalFit->GetParameter(1);
      double sigma_inner = sigma1;
      double err_mean = finalFit->GetParError(1);
      double err_sigma = err_sigma1;
      line1 = Form("Inner Gaussian (core):");
      line2 = Form("#mu = %.2f#pm%.3f %s", mean_inner, err_mean, unit.Data());
      line3 = Form("#sigma = %.2f#pm%.3f %s", sigma_inner, err_sigma, unit.Data());
      // Add fraction info
      double amp1 = finalFit->GetParameter(0);
      double amp2 = finalFit->GetParameter(3);
      double fraction = amp1 / (amp1 + amp2);
      line4 = Form("Core fraction: %.1f%%", fraction * 100);
    } else {
      double mean_inner = finalFit->GetParameter(4);
      double sigma_inner = sigma2;
      double err_mean = finalFit->GetParError(4);
      double err_sigma = err_sigma2;
      line1 = Form("Inner Gaussian (core):");
      line2 = Form("#mu = %.2f#pm%.3f %s", mean_inner, err_mean, unit.Data());
      line3 = Form("#sigma = %.2f#pm%.3f %s", sigma_inner, err_sigma, unit.Data());
      double amp1 = finalFit->GetParameter(0);
      double amp2 = finalFit->GetParameter(3);
      double fraction = amp2 / (amp1 + amp2);
      line4 = Form("Core fraction: %.1f%%", fraction * 100);
    }
  } else if (finalFit) {
    // Single Gaussian
    double mean_sg = finalFit->GetParameter(1);
    double sigma_sg = finalFit->GetParameter(2);
    double err_mean = finalFit->GetParError(1);
    double err_sigma = finalFit->GetParError(2);
    line1 = Form("Gaussian fit:");
    // Format mean
    line2 = Form("#mu = %.3e#pm%.3e %s", mean_sg, err_mean, unit.Data());
    line3 = Form("#sigma = %.3e#pm%.3e %s", sigma_sg, err_sigma, unit.Data());
    
    //line2 = Form("#mu = %.2f#pm%.3f %s", mean_sg, err_mean, unit.Data());
    //line3 = Form("#sigma = %.2f#pm%.3f %s", sigma_sg, err_sigma, unit.Data());
    //line4 = "";
    //cout << mean_sg << endl;
  }
  
  TString line5 = Form("#chi^{2}/NDF = %.2f", chi2ndf);

  // Create transparent info box
  TPaveText *pt = new TPaveText(0.25, 0.70, 0.90, 0.88, "NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.04);
  pt->SetTextFont(42);
  pt->AddText(line1);
  pt->AddText(line2);
  pt->AddText(line3);
  //if (line4.Length() > 0) pt->AddText(line4);
  //pt->AddText(line5);
  pt->Draw();

  // Legend
  TLegend *leg = new TLegend(0.55, 0.65, 0.90, 0.7);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  //leg->AddEntry(h_diff, "Data", "l");
  leg->AddEntry(finalFit, fitType, "l");
  leg->Draw();

  // Add information about fit stability
  TString stabilityText = doubleUsed ? "Fit: Double Gaussian (stable)" : "Fit: Single Gaussian";
  TLatex *lat = new TLatex(0.15, 0.92, stabilityText);
  lat->SetNDC();
  lat->SetTextSize(0.035);
  lat->SetTextFont(42);
  lat->Draw();

  // Save canvas
  c1->SaveAs("../plots_resol/" + var_type + "_diff_fit_normalized.png");
  std::cout << "\n✅ Plot saved to ../plots_resol/" + var_type + "_diff_fit_normalized.png" << std::endl;
  
  // Close file
  ftree->Close();
}
