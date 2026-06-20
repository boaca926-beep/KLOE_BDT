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
  //TString treeFile = "/home/bo/Desktop/input_bdt_TDATA_norm/cut/tree_pre_bdt.root";
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
  
  for (Long64_t i = 0; i < ttree->GetEntries(); ++i) {
    ttree->GetEntry(i);
    if (recon_indx_bdt == 2 && bkg_indx == 1) {
      var_diff = var - var_true;
      h_diff->Fill(var_diff);
      //cout << var_diff << endl;
    }
  }

  // --- Normalize the histogram (integral = 1) ---
  double integral = h_diff->Integral();
  if (integral > 0) {
    h_diff->Scale(1.0 / integral);
  } else {
    std::cerr << "ERROR: histogram has zero entries." << std::endl;
    return;
  }

  // Fit range: mean ± 3σ
  //double mean = h_diff->GetMean();
  double mean = 0.0;
  double rms  = h_diff->GetRMS();
  double fit_min = mean - fit_factor * rms;
  double fit_max = mean + fit_factor * rms;
  double hist_min = h_diff->GetXaxis()->GetXmin();
  double hist_max = h_diff->GetXaxis()->GetXmax();
  if (fit_min < hist_min) fit_min = hist_min;
  if (fit_max > hist_max) fit_max = hist_max;
  std::cout << "Fit range for normalized var_diff: [" << fit_min << ", " << fit_max << "] " << unit << std::endl;

  double peak = h_diff->GetMaximum();
  
  // Double Gaussian fit
  TF1 *doubleGaus = new TF1("doubleGaus", "gaus(0)+gaus(3)", fit_min, fit_max);
  doubleGaus->SetParameter(0, peak * 0.8);
  doubleGaus->SetParameter(1, mean);
  doubleGaus->SetParameter(2, rms * 0.7);
  doubleGaus->SetParameter(3, peak * 0.2);
  doubleGaus->SetParameter(4, mean + 2.0);
  doubleGaus->SetParameter(5, rms * 1.2);
  doubleGaus->SetParLimits(2, 0.5, 30.0);
  doubleGaus->SetParLimits(5, 1.0, 50.0);
  doubleGaus->SetParLimits(4, mean - 10.0, mean + 20.0);
  doubleGaus->SetParLimits(0, 0.0, peak * 2.0);
  doubleGaus->SetParLimits(3, 0.0, peak * 1.5);
  doubleGaus->SetLineColor(kRed);
  doubleGaus->SetLineWidth(2);

  h_diff->Fit(doubleGaus, "R");
  
  double chi2ndf = doubleGaus->GetChisquare() / doubleGaus->GetNDF();
  double err_amp2 = doubleGaus->GetParError(3);
  double amp2 = doubleGaus->GetParameter(3);
  bool stable = (chi2ndf < 10.0) && (err_amp2 / (amp2 + 1e-6) < 2.0);
  
  TF1 *finalFit = doubleGaus;
  TString fitType = "Double Gaussian";
  bool doubleUsed = true;

  if (!stable) {
    std::cout << "Double Gaussian unstable, falling back to single Gaussian." << std::endl;
    TF1 *singleGaus = new TF1("singleGaus", "gaus", fit_min, fit_max);
    singleGaus->SetParameters(peak, mean, rms * 0.7);
    singleGaus->SetParLimits(2, 0.5, 30.0);
    singleGaus->SetLineColor(kRed);
    singleGaus->SetLineWidth(2);
    h_diff->Fit(singleGaus, "R");
    finalFit = singleGaus;
    fitType = "Single Gaussian";
    doubleUsed = false;
    chi2ndf = singleGaus->GetChisquare() / singleGaus->GetNDF();
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
  h_diff->GetXaxis()->SetRangeUser(range_factor * fit_min, range_factor * fit_max); // or -50,50
  //h_diff->GetXaxis()->SetRangeUser(XMIN, XMAX);   
  h_diff->GetYaxis()->SetNdivisions(505);
  h_diff->GetXaxis()->SetNdivisions(505);

  h_diff->Draw("hist");
  finalFit->Draw("same");

  //cout << "x range: " << range_factor * fit_min << ", " << range_factor * fit_max << endl;

  // ----- Extract inner (narrower) Gaussian values if double Gaussian -----
  TString line1, line2, line3;
  if (doubleUsed) {
    double sigma1 = finalFit->GetParameter(2);
    double sigma2 = finalFit->GetParameter(5);
    // Choose the narrower sigma (smaller value) as the inner resolution
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

  // Create a transparent TPaveText
  TPaveText *pt = new TPaveText(0.3, 0.72, 0.9, 0.88, "NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.04);
  pt->SetTextFont(42);
  if (doubleUsed) pt->AddText(line1);
  pt->AddText(line2);
  pt->AddText(line3);
  //pt->AddText(line4);
  pt->Draw();

  // Legend for fit line
  TLegend *leg = new TLegend(0.55, 0.60, 0.9, 0.68);
  //leg->SetTextFont(132);
  leg->SetTextSize(0.035);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(finalFit, fitType, "l");
  leg->Draw();

  // Save canvas
  c1->SaveAs("../plots_resol/" + var_type + "_diff_fit_normalized.png");
  std::cout << "\nPlot saved to ../plots_resol/" + var_type + "diff_fit_normalized.png" << std::endl;
}
