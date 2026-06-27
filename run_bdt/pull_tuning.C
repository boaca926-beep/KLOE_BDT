void checkFile(TFile *f_input){

  TIter next_tree(f_input -> GetListOfKeys());

  TString objnm_tree, classnm_tree;

  int i = 0;
  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

}

// Define FitResult struct BEFORE using it
struct FitResult {
    TString name;
    double mean, mean_err;
    double sigma, sigma_err;
    double chi2_ndf;
    int entries;
};

// Breit-Wigner function: p0 / ((x-p1)^2 + p2^2) with p0 = normalization, p1 = mean, p2 = gamma/2
Double_t breitwigner(Double_t *x, Double_t *par) {
    return par[0] / ((x[0] - par[1]) * (x[0] - par[1]) + par[2] * par[2]);
}

void pull_tuning() {

  gROOT->GetListOfCanvases()->Delete();
  
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(5);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  gSystem->Exec("mkdir -p ../pull_tuning");

  TFile *fin_E1 = new TFile("../output_pull_E1/hist_pull_E1.root");
  if (!fin_E1 || fin_E1->IsZombie()) {
    std::cerr << "ERROR: cannot open output_pull_E1/hist_pull_E1.root" << std::endl;
    return;
  }

  TFile *fin_E2 = new TFile("../output_pull_E2/hist_pull_E2.root");
  if (!fin_E2 || fin_E2->IsZombie()) {
    std::cerr << "ERROR: cannot open output_pull_E2/hist_pull_E2.root" << std::endl;
    return;
  }

  TFile *fin_E3 = new TFile("../output_pull_E3/hist_pull_E3.root");
  if (!fin_E3 || fin_E3->IsZombie()) {
    std::cerr << "ERROR: cannot open output_pull_E3/hist_pull_E3.root" << std::endl;
    return;
  }

  TFile *fin_m3pi = new TFile("../output_m3pi_bdt/hist_m3pi_bdt.root");
  if (!fin_m3pi || fin_m3pi->IsZombie()) {
    std::cerr << "ERROR: cannot open output_m3pi_bdt/hist_m3pi_bdt.root" << std::endl;
    return;
  }

  // Get mass histograms
  TH1D *hm3pi_MC = (TH1D*)fin_m3pi->Get("hist_isr3pi_sc"); 
  TH1D *hm3pi_Data = (TH1D*)fin_m3pi->Get("hist_data"); 
  if (!hm3pi_MC) std::cerr << "WARNING: hm3pi_MC not found." << std::endl;
  if (!hm3pi_Data) std::cerr << "WARNING: hm3pi_Data not found." << std::endl;

  // Print file contents for debugging
  // checkFile(fin_m3pi);
  
  TH1D* hE1_MC = (TH1D*)fin_E1->Get("hist_isr3pi_sc");
  TH1D* hE1_DATA = (TH1D*)fin_E1->Get("hist_data");

  TH1D* hE2_MC = (TH1D*)fin_E2->Get("hist_isr3pi_sc");
  TH1D* hE2_DATA = (TH1D*)fin_E2->Get("hist_data");

  TH1D* hE3_MC = (TH1D*)fin_E3->Get("hist_isr3pi_sc");
  TH1D* hE3_DATA = (TH1D*)fin_E3->Get("hist_data");

  // histos E1+E2, E3
  TH1D *hE12_MC = (TH1D*)hE1_MC->Clone("hE12_MC");
  hE12_MC->Add(hE2_MC, 1.);

  TH1D *hE12_DATA = (TH1D*)hE1_DATA->Clone("hE12_DATA");
  hE12_DATA->Add(hE2_DATA, 1.);

  hE3_MC->SetName("hE3_MC");
  hE3_DATA->SetName("hE3_DATA");
  
  // Fit the combined distribution
  const int nb_hist = 4;
  TH1D *HLIST[nb_hist] = {hE12_MC, hE12_DATA, hE3_MC, hE3_DATA};
  FitResult results[nb_hist];
  TF1 *gaus_fits[nb_hist];
  
  for (int i = 0; i < nb_hist; i ++) {

    TH1D *h_tmp = (TH1D*)HLIST[i]->Clone(Form("%s", HLIST[i]->GetName()));
    double mean_est = h_tmp->GetMean();
    double rms_est = h_tmp->GetRMS();
    double amp_est = h_tmp->GetBinContent(h_tmp->FindBin(mean_est));

    // ============================================================
    // Fit range: ±1.5σ around mean (symmetric core fit)
    // ============================================================
    const double fit_width = 1.5;
    const double fit_min = mean_est - fit_width;
    const double fit_max = mean_est + fit_width;
    
    std::cout << "========================================" << std::endl;
    std::cout << "Sample: " << h_tmp->GetName() << std::endl;
    std::cout << "Number of entries: " << h_tmp->GetEntries() << std::endl;
    std::cout << "Mean estimate: " << mean_est << ", RMS: " << rms_est << std::endl;
    std::cout << "Fit range: [" << fit_min << ", " << fit_max << "]" << std::endl;

    // ============================================================
    // SINGLE GAUSSIAN FIT (for both MC and Data)
    // High signal purity → no need for double Gaussian
    // ============================================================
    TF1 *singleGaus = new TF1("singleGaus", "gaus", fit_min, fit_max);
    singleGaus->SetParameters(amp_est, mean_est, rms_est * 0.7);
    singleGaus->SetParLimits(2, 0.1, 10.0);
    singleGaus->SetLineColor(kRed);
    singleGaus->SetLineWidth(2);
    
    h_tmp->Fit(singleGaus, "RQS");
    
    double mean = singleGaus->GetParameter(1);
    double sigma = singleGaus->GetParameter(2);
    double mean_err = singleGaus->GetParError(1);
    double sigma_err = singleGaus->GetParError(2);
    double chi2ndf = singleGaus->GetChisquare() / singleGaus->GetNDF();
    
    std::cout << "Single Gaussian: mean = " << mean << " +/- " << mean_err << std::endl;
    std::cout << "Single Gaussian: sigma = " << sigma << " +/- " << sigma_err << std::endl;
    std::cout << "χ²/ndf = " << chi2ndf << std::endl;

    // Store results
    results[i].name = TString(h_tmp->GetName());
    results[i].mean = mean;
    results[i].mean_err = mean_err;
    results[i].sigma = sigma;
    results[i].sigma_err = sigma_err;
    results[i].chi2_ndf = chi2ndf;
    results[i].entries = h_tmp->GetEntries();
    gaus_fits[i] = singleGaus;
    
  }
  
  // Set colors for fits (unchanged)
  gaus_fits[0]->SetLineColor(kRed);      // hE12_MC
  gaus_fits[0]->SetLineWidth(2);
  gaus_fits[1]->SetLineColor(kGreen+2);  // hE12_DATA
  gaus_fits[1]->SetLineWidth(2);
  gaus_fits[2]->SetLineColor(kRed);      // hE3_MC
  gaus_fits[2]->SetLineWidth(2);
  gaus_fits[3]->SetLineColor(kGreen+2);  // hE3_DATA
  gaus_fits[3]->SetLineWidth(2);

  // ========== Canvas 1: E1+E2 Pull ==========
  TCanvas *c_E12 = new TCanvas("c_E12", "E1+E2 Pull Distributions", 900, 900);

  c_E12->cd(1);
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);

  const double ymax_E12 = hE12_DATA->GetMaximum();
  hE12_DATA->GetYaxis()->SetTitle("Events");
  hE12_DATA->GetYaxis()->SetRangeUser(0.01, ymax_E12 * 1.6);
  hE12_DATA->GetYaxis()->CenterTitle();
  hE12_DATA->GetYaxis()->SetTitleSize(0.05);
  hE12_DATA->GetYaxis()->SetTitleOffset(1.4);
  hE12_DATA->GetYaxis()->SetLabelSize(0.04);
  hE12_DATA->GetXaxis()->SetTitle("E_{1}+E_{2} Pull [MeV]");
  hE12_DATA->GetXaxis()->SetTitleSize(0.05);
  hE12_DATA->GetXaxis()->SetTitleOffset(1.2);
  hE12_DATA->GetXaxis()->SetLabelSize(0.04);
  hE12_DATA->GetXaxis()->CenterTitle();
  
  hE12_DATA->Draw();
  gaus_fits[1]->Draw("same");
  hE12_MC->Draw("same hist");
  gaus_fits[0]->Draw("same");

  TLegend *leg_E12 = new TLegend(0.2, 0.7, 0.9, 0.9);
  leg_E12->SetFillStyle(0);
  leg_E12->SetBorderSize(0);
  leg_E12->SetNColumns(2);
  leg_E12->SetTextSize(0.03);
  leg_E12->AddEntry(hE12_MC, "MC", "f");
  leg_E12->AddEntry(gaus_fits[0], Form("#mu = %.3f, #sigma = %.3f", results[0].mean, results[0].sigma), "l");
  leg_E12->AddEntry(hE12_DATA, "Data", "f");
  leg_E12->AddEntry(gaus_fits[1], Form("#mu = %.3f, #sigma = %.3f", results[1].mean, results[1].sigma), "l");
  leg_E12->Draw();
    
  gPad->Update();

  // ========== Canvas 2: E3 Pull ==========
  TCanvas *c_E3 = new TCanvas("c_E3", "E3 Pull Distributions", 900, 900);

  c_E3->cd(1);
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);

  const double ymax_E3 = hE3_DATA->GetMaximum();
  hE3_DATA->GetYaxis()->SetTitle("Events");
  hE3_DATA->GetYaxis()->SetRangeUser(0.01, ymax_E3 * 1.6);
  hE3_DATA->GetYaxis()->CenterTitle();
  hE3_DATA->GetYaxis()->SetTitleSize(0.05);
  hE3_DATA->GetYaxis()->SetTitleOffset(1.4);
  hE3_DATA->GetYaxis()->SetLabelSize(0.04);
  hE3_DATA->GetXaxis()->SetTitle("E_{3} Pull [MeV]");
  hE3_DATA->GetXaxis()->SetTitleSize(0.05);
  hE3_DATA->GetXaxis()->SetTitleOffset(1.2);
  hE3_DATA->GetXaxis()->SetLabelSize(0.04);
  hE3_DATA->GetXaxis()->CenterTitle();
  
  hE3_DATA->Draw();
  gaus_fits[3]->Draw("same");
  hE3_MC->Draw("same hist");
  gaus_fits[2]->Draw("same");

  TLegend *leg_E3 = new TLegend(0.2, 0.7, 0.9, 0.9);
  leg_E3->SetFillStyle(0);
  leg_E3->SetBorderSize(0);
  leg_E3->SetNColumns(2);
  leg_E3->SetTextSize(0.03);
  leg_E3->AddEntry(hE3_MC, "MC", "f");
  leg_E3->AddEntry(gaus_fits[2], Form("#mu = %.3f, #sigma = %.3f", results[2].mean, results[2].sigma), "l");
  leg_E3->AddEntry(hE3_DATA, "Data", "f");
  leg_E3->AddEntry(gaus_fits[3], Form("#mu = %.3f, #sigma = %.3f", results[3].mean, results[3].sigma), "l");
  leg_E3->Draw();
    
  gPad->Update();

  // Save canvases
  c_E12->SaveAs("../pull_tuning/pull_E12.pdf");
  c_E3->SaveAs("../pull_tuning/pull_E3.pdf");

  // ============================================================
  // MASS FITS: MC and Data side-by-side (FIXED)
  // ============================================================
  const int nb_mass = 2;
  TH1D *hMassList[nb_mass] = {hm3pi_MC, hm3pi_Data};
  TString massNameList[nb_mass] = {"MC", "Data"};
  int massColor[nb_mass] = {kBlue, kRed};
  FitResult massResults[nb_mass];
  TF1 *bw_fits[nb_mass];

  // Create canvas for mass fits
  for (int i = 0; i < nb_mass; i++) {
    TH1D *h_mass = hMassList[i];
    if (!h_mass) {
      std::cerr << "Skipping mass histogram " << massNameList[i] << " (null)." << std::endl;
      continue;
    }

    TCanvas *c_mass = new TCanvas("c_mass_" + massNameList[i], "3π Mass Distributions (Breit-Wigner fits)", 700, 700);
    c_mass->Divide(nb_mass, 1);

  // Clone to avoid modifying original
    TH1D *h_mass_copy = (TH1D*)h_mass->Clone(Form("h_mass_%s", massNameList[i].Data()));
    h_mass_copy->SetDirectory(0);
    h_mass_copy->SetLineColor(massColor[i]);

    double mass_mean = h_mass_copy->GetMean();
    double mass_rms = h_mass_copy->GetRMS();
    double mass_peak = h_mass_copy->GetBinContent(h_mass_copy->GetMaximumBin());
    double mass_peak_pos = h_mass_copy->GetBinCenter(h_mass_copy->GetMaximumBin());

    // Fit range: ±1.5σ around mean (like pulls) but constrained
    double fit_min_mass = mass_mean - .5 * mass_rms;
    double fit_max_mass = mass_mean + .5 * mass_rms;
    if (fit_min_mass < 760) fit_min_mass = 760;
    if (fit_max_mass > 810) fit_max_mass = 810;

    std::cout << "\n========================================" << std::endl;
    std::cout << "Mass histogram: " << massNameList[i] << std::endl;
    std::cout << "Number of entries: " << h_mass_copy->GetEntries() << std::endl;
    std::cout << "Mean estimate: " << mass_mean << ", RMS: " << mass_rms << std::endl;
    std::cout << "Fit range: [" << fit_min_mass << ", " << fit_max_mass << "] MeV/c^{2}" << std::endl;

    // Breit-Wigner fit
    TF1 *bw = new TF1(Form("bw_%s", massNameList[i].Data()), breitwigner, fit_min_mass, fit_max_mass, 3);
    bw->SetParameters(mass_peak * 4.0, mass_peak_pos, 4.0);
    bw->SetParLimits(1, 780, 786);
    bw->SetParLimits(2, 0.5, 10.0);
    bw->SetLineColor(massColor[i]);
    bw->SetLineWidth(2);

    h_mass_copy->Fit(bw, "RQS");
    double mass_mean_fit = bw->GetParameter(1);
    double mass_gamma_half = bw->GetParameter(2);
    double mass_mean_err = bw->GetParError(1);
    double mass_gamma_err = bw->GetParError(2);
    double chi2ndf_mass = bw->GetChisquare() / bw->GetNDF();

    std::cout << "Breit-Wigner fit: mean = " << mass_mean_fit << " +/- " << mass_mean_err << " MeV/c^{2}" << std::endl;
    std::cout << "Breit-Wigner gamma/2 = " << mass_gamma_half << " +/- " << mass_gamma_err << " MeV" << std::endl;
    std::cout << "χ²/ndf = " << chi2ndf_mass << std::endl;

    // Store results
    massResults[i].name = massNameList[i];
    massResults[i].mean = mass_mean_fit;
    massResults[i].mean_err = mass_mean_err;
    massResults[i].sigma = mass_gamma_half;
    massResults[i].sigma_err = mass_gamma_err;
    massResults[i].chi2_ndf = chi2ndf_mass;
    massResults[i].entries = h_mass_copy->GetEntries();
    bw_fits[i] = bw;

    // Draw in pad
    c_mass->cd(i+1);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);

    double ymax_mass = h_mass_copy->GetMaximum();
    h_mass_copy->SetMarkerStyle(20);
    h_mass_copy->SetMarkerSize(0.6);
    h_mass_copy->SetLineColor(kBlack);
    h_mass_copy->GetYaxis()->SetTitle("Events");
    h_mass_copy->GetYaxis()->SetRangeUser(0.01, ymax_mass * 1.6);
    h_mass_copy->GetYaxis()->CenterTitle();
    h_mass_copy->GetYaxis()->SetTitleSize(0.05);
    h_mass_copy->GetYaxis()->SetTitleOffset(1.4);
    h_mass_copy->GetYaxis()->SetLabelSize(0.04);
    h_mass_copy->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
    h_mass_copy->GetXaxis()->SetTitleSize(0.05);
    h_mass_copy->GetXaxis()->SetTitleOffset(1.2);
    h_mass_copy->GetXaxis()->SetLabelSize(0.04);
    h_mass_copy->GetXaxis()->CenterTitle();

    if (massNameList[i] == "MC") {
      h_mass_copy->Draw("hist");
    }
    else {
        h_mass_copy->Draw("E");
    }

    bw->Draw("same");

    TLegend *leg_mass = new TLegend(0.2, 0.7, 0.65, 0.9);
    leg_mass->SetFillStyle(0);
    leg_mass->SetBorderSize(0);
    leg_mass->SetTextSize(0.035);
    leg_mass->AddEntry(h_mass_copy, Form("%s 3#pi mass", massNameList[i].Data()), "lep");
    leg_mass->AddEntry(bw, Form("BW: M = %.2f, #Gamma/2 = %.2f", mass_mean_fit, mass_gamma_half), "l");
    leg_mass->Draw();

    c_mass->Update();
    c_mass->SaveAs("../pull_tuning/mass_fit_" + massNameList[i] + ".pdf");
    //delete c_mass;
  }

  
  // Print summary (only for pull fits)
  std::cout << "\n========================================" << std::endl;
  std::cout << "Summary of Pull Fit Results (Single Gaussian):" << std::endl;
  std::cout << "========================================" << std::endl;
  for (int i = 0; i < nb_hist; i++) {
    std::cout << Form("%-15s: mean = %6.3f +/- %6.3f, sigma = %6.3f +/- %6.3f, χ²/ndf = %.3f", 
                      results[i].name.Data(), 
                      results[i].mean, results[i].mean_err,
                      results[i].sigma, results[i].sigma_err,
                      results[i].chi2_ndf) << std::endl;
  }

  std::cout << "\n========================================" << std::endl;
  std::cout << "Summary of Mass Fit Results (Breit-Wigner):" << std::endl;
  std::cout << "========================================" << std::endl;
  for (int i = 0; i < nb_mass; i++) {
    if (hMassList[i]) {
      std::cout << Form("%-10s: mean = %6.3f +/- %6.3f, gamma/2 = %6.3f +/- %6.3f, χ²/ndf = %.3f", 
                        massResults[i].name.Data(), 
                        massResults[i].mean, massResults[i].mean_err,
                        massResults[i].sigma, massResults[i].sigma_err,
                        massResults[i].chi2_ndf) << std::endl;
    }
  }

  // Use the single Gaussian parameters for MC only
  double bias_E12 = results[0].mean;     // hE12_MC single mean
  double sigma_scale_E12 = results[0].sigma; // hE12_MC single sigma

  double bias_E3 = results[2].mean;      // hE3_MC single mean
  double sigma_scale_E3 = results[2].sigma; // hE3_MC single sigma

  std::ofstream myfile;
  TString myfile_nm = "../header_bdt/tuning.txt";
  myfile.open(myfile_nm.Data());
  
  myfile << "const double bias_E12 = " << bias_E12 << ";\n"
	 << "const double sigma_scale_E12 = " << sigma_scale_E12 << ";\n"
	 << "const double bias_E3 = " << bias_E3 << ";\n"
    	 << "const double sigma_scale_E3 = " << sigma_scale_E3 << ";\n\n";
  //<< "const double energy_shift = " << massResults[1].mean - massResults[0].mean << ";\n";
	 
  myfile.close();
  
}
