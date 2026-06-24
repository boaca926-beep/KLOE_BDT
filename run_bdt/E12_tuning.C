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
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  gSystem->Exec("mkdir -p ../pull_tuning");

  TFile *fin_E1 = new TFile("../output_pull_E1/hist_pull_E1.root");
  if (!fin_E1 || fin_E1->IsZombie()) {
    std::cerr << "ERROR: cannot create output file." << std::endl;
    return;
  }

  TFile *fin_E2 = new TFile("../output_pull_E2/hist_pull_E2.root");
  if (!fin_E2 || fin_E2->IsZombie()) {
    std::cerr << "ERROR: cannot create output file." << std::endl;
    return;
  }

  // Check file content
  //checkFile(fin_E1);
  //checkFile(fin_E2);
  
  TH1D* hE1_MC = (TH1D*)fin_E1->Get("hist_isr3pi_sc");
  TH1D* hE1_DATA = (TH1D*)fin_E1->Get("hist_data");

  TH1D* hE2_MC = (TH1D*)fin_E2->Get("hist_isr3pi_sc");
  TH1D* hE2_DATA = (TH1D*)fin_E2->Get("hist_data");

  // histos E1+E2
  TH1D *hE12_MC = (TH1D*)hE1_MC->Clone("hE12_MC");
  hE12_MC->Add(hE2_MC, 1.);

  TH1D *hE12_DATA = (TH1D*)hE1_DATA->Clone("hE12_DATA");
  hE12_DATA->Add(hE2_DATA, 1.);

  // Fit the combined distribution
  const int nb_hist = 2;
  TH1D *HLIST[nb_hist] = {hE12_MC, hE12_DATA};
  FitResult result_mc, result_data;
  TF1 *gaus_mc = nullptr;
  TF1 *gaus_data = nullptr;
  
  for (int i = 0; i < nb_hist; i ++) {

    TH1D *h_tmp = (TH1D*)HLIST[i]->Clone();
    double mean_est = h_tmp->GetMean();
    double rms_est = h_tmp->GetRMS();
    double amp_est = h_tmp->GetBinContent(hE12_MC->FindBin(mean_est));

    const double fit_width = 1.5;
    const double fit_min = mean_est - fit_width;
    const double fit_max = mean_est + fit_width;
    
    std::cout << "========================================" << std::endl;
    std::cout << "Sample: " << h_tmp->GetName() << std::endl;
    std::cout << "Combined Pull: E1 + E2" << std::endl;
    std::cout << "Number of entries: " << h_tmp->GetEntries() << std::endl;
    std::cout << "Mean estimate: " << mean_est << ", RMS: " << rms_est << std::endl;
    std::cout << "Fit range: [" << fit_min << ", " << fit_max << "]" << std::endl;

    TF1 *gaus = new TF1("gaus", "gaus", fit_min, fit_max);
    gaus->SetParameters(amp_est, mean_est, rms_est * 0.7);
    gaus->SetParLimits(2, 0.2, 5.0);
    h_tmp->Fit(gaus, "RQS");

    double mean = gaus->GetParameter(1);
    double sigma = gaus->GetParameter(2);
    double mean_err = gaus->GetParError(1);
    double sigma_err = gaus->GetParError(2);
    double chi2_ndf = gaus->GetChisquare() / gaus->GetNDF();
    
    std::cout << "Mean = " << mean << " +/- " << mean_err << std::endl;
    std::cout << "Sigma = " << sigma << " +/- " << sigma_err << std::endl;
    std::cout << "χ²/ndf = " << chi2_ndf << std::endl;

    if (TString(h_tmp->GetName()) == "hE12_MC") {
      cout << h_tmp->GetName() << endl;
      gaus_mc = gaus;
      result_mc.name = TString(h_tmp->GetName());
      result_mc.mean = mean;
      result_mc.mean_err = mean_err;
      result_mc.sigma = sigma;
      result_mc.sigma_err = sigma_err;
      result_mc.chi2_ndf = chi2_ndf;
      result_mc.entries = h_tmp->GetEntries();
    } else if (TString(h_tmp->GetName()) == "hE12_DATA") {
      cout << h_tmp->GetName() << endl;
      gaus_data = gaus;
      result_data.name = TString(h_tmp->GetName());
      result_data.mean = mean;
      result_data.mean_err = mean_err;
      result_data.sigma = sigma;
      result_data.sigma_err = sigma_err;
      result_data.chi2_ndf = chi2_ndf;
      result_data.entries = h_tmp->GetEntries();
    }
    
  }
  
  //hE1_DATA->Draw();
  //hE1_MC->Draw("same");

  //hE2_DATA->Draw();
  //hE2_MC->Draw("same");

  // Add text with χ²/ndf
  //TPaveText *pt = new TPaveText(0.15, 0.78, 0.48, 0.88, "NDC");
  //pt->SetFillStyle(0);
  //pt->SetBorderSize(0);
  //pt->SetTextSize(0.04);
  //pt->AddText(Form("#chi^{2}/ndf = %.2f", chi2_ndf));
  //pt->Draw();
  
 
  gaus_mc->SetLineColor(kRed);
  gaus_mc->SetLineWidth(2);

  gaus_data->SetLineColor(kBlack);
  gaus_data->SetLineWidth(2);

  TCanvas *c_side = new TCanvas("c_side", "Pull Distributions", 1400, 700);

  c_side->Divide(2, 1);

  c_side->cd(1);
  gPad->SetBottomMargin(0.12);
  gPad->SetLeftMargin(0.12);

  
  const double ymax = hE12_DATA->GetMaximum();
  hE12_DATA->GetYaxis()->SetTitle("Events");
  hE12_DATA->GetYaxis()->SetRangeUser(0.01, ymax * 1.6);
  hE12_DATA->GetYaxis()->CenterTitle();
  hE12_DATA->GetYaxis()->SetTitleSize(0.05);
  hE12_DATA->GetYaxis()->SetTitleOffset(1.2);
  hE12_DATA->GetYaxis()->SetLabelSize(0.04);
  hE12_DATA->GetXaxis()->SetTitle("E_{1}+E_{2} Pull [MeV]");
  hE12_DATA->GetXaxis()->SetTitleSize(0.05);
  hE12_DATA->GetXaxis()->SetTitleOffset(1.0);
  hE12_DATA->GetXaxis()->SetLabelSize(0.04);
  hE12_DATA->GetXaxis()->CenterTitle();
  
  hE12_DATA->Draw();
  gaus_data->Draw("same");
  hE12_MC->Draw("same");
  gaus_mc->Draw("same");

  // Legend - position adjusted for each pad
  TLegend *leg = new TLegend(0.15, 0.7, 0.9, 0.9);
  //TLegend *leg = new TLegend(0.6, 0.50, 0.92, 0.88);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetNColumns(2);
  leg->SetTextSize(0.03);
  leg->AddEntry(hE12_MC, Form("%s (E1+E2)", result_mc.name.Data()), "lep");
  leg->AddEntry(gaus_mc, Form("#mu = %.3f, #sigma = %.3f", result_mc.mean, result_mc.sigma), "l");
  leg->AddEntry(hE12_DATA, Form("%s (E1+E2)", result_data.name.Data()), "lep");
  leg->AddEntry(gaus_data, Form("#mu = %.3f, #sigma = %.3f", result_data.mean, result_data.sigma), "l");
  leg->Draw();
    
  // Update the pad
  gPad->Update();
    
    
  
}
