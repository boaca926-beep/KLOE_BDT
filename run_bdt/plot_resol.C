// Resolution – fit to m3pi_diff distribution with inner Gaussian values (normalized)

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
  TString treeFile = "/home/kloe/Desktop/input_bdt_TDATA_chain/cut/tree_pre_bdt.root";
  TFile *ftree = TFile::Open(treeFile);
  if (!ftree || ftree->IsZombie()) {
    std::cerr << "ERROR: cannot open " << treeFile << std::endl;
    return;
  }

  TTree *ttree = (TTree*) ftree->Get("TISR3PI_SIG");
  if (!ttree) { std::cerr << "No such tree.\n"; return; }

  // Histogram for difference: (reco – true) mass
  TH1D *h_m3pi_diff = new TH1D("h_m3pi_diff", "", 200, -50, 50);
  h_m3pi_diff->SetMarkerStyle(21);
  h_m3pi_diff->SetMarkerSize(0.7);
  h_m3pi_diff->SetLineColor(1);
  h_m3pi_diff->Sumw2();

  double m3pi = 0.0, m3pi_true = 0.0;
  double m3pi_diff = 0.0;
  int recon_indx_bdt = -1, bkg_indx = -1;
  
  ttree->SetBranchAddress("Br_recon_indx_bdt", &recon_indx_bdt);
  ttree->SetBranchAddress("Br_bkg_indx", &bkg_indx);
  ttree->SetBranchAddress("Br_m3pi_bdt", &m3pi);
  ttree->SetBranchAddress("Br_m3pi_true_bdt", &m3pi_true);
  
  for (Long64_t i = 0; i < ttree->GetEntries(); ++i) {
    ttree->GetEntry(i);
    if (recon_indx_bdt == 2 && bkg_indx == 1) {
      m3pi_diff = m3pi - m3pi_true;
      h_m3pi_diff->Fill(m3pi_diff);
    }
  }

  // --- Normalize the histogram (integral = 1) ---
  double integral = h_m3pi_diff->Integral();
  if (integral > 0) {
    h_m3pi_diff->Scale(1.0 / integral);
  } else {
    std::cerr << "ERROR: histogram has zero entries." << std::endl;
    return;
  }

  // Fit range: mean ± 3σ
  double mean = h_m3pi_diff->GetMean();
  double rms  = h_m3pi_diff->GetRMS();
  double fit_min = mean - 1.5 * rms;
  double fit_max = mean + 1.5 * rms;
  double hist_min = h_m3pi_diff->GetXaxis()->GetXmin();
  double hist_max = h_m3pi_diff->GetXaxis()->GetXmax();
  if (fit_min < hist_min) fit_min = hist_min;
  if (fit_max > hist_max) fit_max = hist_max;
  std::cout << "Fit range for normalized m3pi_diff: [" << fit_min << ", " << fit_max << "] MeV/c^2" << std::endl;

  double peak = h_m3pi_diff->GetMaximum();
  
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

  h_m3pi_diff->Fit(doubleGaus, "R");
  
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
    singleGaus->SetLineColor(kBlue);
    singleGaus->SetLineWidth(2);
    h_m3pi_diff->Fit(singleGaus, "R");
    finalFit = singleGaus;
    fitType = "Single Gaussian";
    doubleUsed = false;
    chi2ndf = singleGaus->GetChisquare() / singleGaus->GetNDF();
  }

  // Draw histogram and fit
  TCanvas *c1 = new TCanvas("c1", "m3pi_diff resolution (normalized)", 700, 700);
  c1->SetLeftMargin(0.15);
  c1->SetBottomMargin(0.15);

  h_m3pi_diff->SetLineWidth(2);
  h_m3pi_diff->GetXaxis()->SetTitle("M^{true}_{3#pi}-M^{rec}_{3#pi} [MeV/c^{2}]");
  h_m3pi_diff->GetYaxis()->SetTitle("Normalized Entries");
  h_m3pi_diff->GetXaxis()->CenterTitle();
  h_m3pi_diff->GetYaxis()->CenterTitle();
  h_m3pi_diff->GetXaxis()->SetTitleSize(0.05);
  h_m3pi_diff->GetYaxis()->SetTitleSize(0.06);
  h_m3pi_diff->GetXaxis()->SetLabelSize(0.05);
  h_m3pi_diff->GetYaxis()->SetLabelSize(0.05);
  h_m3pi_diff->GetYaxis()->SetTitleOffset(1.3);
  h_m3pi_diff->GetYaxis()->SetRangeUser(0., 0.1);
  h_m3pi_diff->GetXaxis()->SetRangeUser(3 * fit_min, 3 * fit_max); // or -50,50
  h_m3pi_diff->GetYaxis()->SetNdivisions(505);

  h_m3pi_diff->Draw("hist");
  finalFit->Draw("same");

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
      line2 = Form("#mu = %.2f #pm %.2f MeV/c^{2}", mean_inner, err_mean);
      line3 = Form("#deltaM_{3#pi} = %.2f #pm %.2f MeV/c^{2}", sigma_inner, err_sigma);
    } else {
      double mean_inner = finalFit->GetParameter(4);
      double sigma_inner = sigma2;
      double err_mean = finalFit->GetParError(4);
      double err_sigma = finalFit->GetParError(5);
      line1 = Form("Inner Gaussian:");
      line2 = Form("#mu = %.2f #pm %.2f MeV/c^{2}", mean_inner, err_mean);
      line3 = Form("#deltaM_{3#pi} = %.2f #pm %.2f MeV/c^{2}", sigma_inner, err_sigma);
    }
  } else {
    // Single Gaussian
    double mean_sg = finalFit->GetParameter(1);
    double sigma_sg = finalFit->GetParameter(2);
    double err_mean = finalFit->GetParError(1);
    double err_sigma = finalFit->GetParError(2);
    line1 = Form("Gaussian fit:");
    line2 = Form("#mu = %.2f #pm %.2f MeV/c^{2}", mean_sg, err_mean);
    line3 = Form("#sigma = %.2f #pm %.2f MeV/c^{2}", sigma_sg, err_sigma);
  }
  TString line4 = Form("#chi^{2}/NDF = %.2f", chi2ndf);

  // Create a transparent TPaveText
  TPaveText *pt = new TPaveText(0.45, 0.72, 0.9, 0.88, "NDC");
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
  c1->SaveAs("../plots_resol/m3pi_diff_fit_normalized.png");
  std::cout << "\nPlot saved to ../plots_resol/m3pi_diff_fit_normalized.png" << std::endl;
}
