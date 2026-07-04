// Normalized pulls of the photon final states distributions using KLOE raw data. E1 and E2 are the first and the second paried photon. E3 is the the unpaired photon. M(E1+E2) gives the reconstructed pi0 invaraint mass.
// Determine 

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

void pull_tuning() {

  gROOT->GetListOfCanvases()->Delete();
  
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(5);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  gSystem->Exec("mkdir -p ../pull_tuning");

  const TString tuning_type = "raw";
  
  const TString fin_E1_nm = "../output_kloe_" + tuning_type + "_pull_E1/hist_pull_E1.root";

  TFile *fin_E1 = new TFile(fin_E1_nm);
  if (!fin_E1 || fin_E1->IsZombie()) {
    std::cerr << "ERROR: cannot open " + fin_E1_nm << std::endl;
    return;
  }

  const TString fin_E2_nm = "../output_kloe_" + tuning_type + "_pull_E2/hist_pull_E2.root";

  TFile *fin_E2 = new TFile(fin_E2_nm);
  if (!fin_E2 || fin_E2->IsZombie()) {
    std::cerr << "ERROR: cannot open output_pull_E2_kloe/hist_pull_E2.root" << std::endl;
    return;
  }

  const TString fin_E3_nm = "../output_kloe_" + tuning_type + "_pull_E3/hist_pull_E3.root";

  TFile *fin_E3 = new TFile(fin_E3_nm);
  if (!fin_E3 || fin_E3->IsZombie()) {
    std::cerr << "ERROR: cannot open output_pull_E3_kloe/hist_pull_E3.root" << std::endl;
    return;
  }

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

  // Print summary (only for pull fits)
  std::cout << "\n========================================" << std::endl;
  std::cout << "Summary of E Pull Fit Results (Single Gaussian):" << std::endl;
  std::cout << "========================================" << std::endl;
  for (int i = 0; i < nb_hist; i++) {
    std::cout << Form("%-15s: mean = %6.3f +/- %6.3f, sigma = %6.3f +/- %6.3f, χ²/ndf = %.3f", 
                      results[i].name.Data(), 
                      results[i].mean, results[i].mean_err,
                      results[i].sigma, results[i].sigma_err,
                      results[i].chi2_ndf) << std::endl;
  }

  // Use Data parameters (not MC)
  double bias_E12 = results[1].mean;        // hE12_DATA mean
  double sigma_scale_E12 = results[1].sigma; // hE12_DATA sigma
  double bias_E3 = results[3].mean;         // hE3_DATA mean
  double sigma_scale_E3 = results[3].sigma; // hE3_DATA sigma
 
  std::ofstream myfile;
  TString myfile_nm = "../header/tuning_" + tuning_type + ".h";
  myfile.open(myfile_nm.Data());
  myfile << "const double bias_E12 = " << bias_E12 << ";\n"
	 << "const double sigma_scale_E12 = " << sigma_scale_E12 << ";\n"
	 << "const double bias_E3 = " << bias_E3 << ";\n"
    	 << "const double sigma_scale_E3 = " << sigma_scale_E3 << ";\n\n";
  //<< "const double energy_shift = " << massResults[1].mean - massResults[0].mean << ";\n";
	 
  myfile.close();
  
}
