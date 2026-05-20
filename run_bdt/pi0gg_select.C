// pi0gg_select.C – compare χ² and BDT π⁰ photon identification
void pi0gg_select() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFitFormat("6.4g");

  const char* data_filename = "/home/kloe/Desktop/input_bdt_TDATA_chain/cut/tree_pre_bdt.root";
  const char* tree_name = "TISR3PI_SIG";
    
  TFile* file = TFile::Open(data_filename);
  if (!file || file->IsZombie()) return;
  TTree* tree = (TTree*)file->Get(tree_name);
  if (!tree) return;
  std::cout << "Tree " << tree_name << " has " << tree->GetEntries() << " entries." << std::endl;

  int recon_indx = -1, recon_indx_bdt = -1;
  tree->SetBranchAddress("Br_recon_indx", &recon_indx);
  tree->SetBranchAddress("Br_recon_indx_bdt", &recon_indx_bdt);
    
  TH1D* h_chi2 = new TH1D("h_chi2", "; Number of correct-selected #pi^{0} photons;Events", 3, 0, 3);
  TH1D* h_bdt  = new TH1D("h_bdt",  ";  Number of correct-selected #pi^{0} photons (BDT pairing);Events", 3, 0, 3);
  h_chi2->SetLineColor(kBlue);
  h_bdt->SetLineColor(kRed);
  h_chi2->SetFillStyle(3001);
  h_bdt->SetFillStyle(3001);

  // Integer bin labels
  for (int i = 1; i <= h_chi2->GetXaxis()->GetNbins(); ++i) {
    h_chi2->GetXaxis()->SetBinLabel(i, Form("%d", i-1));
  }

  Long64_t nentries = tree->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    tree->GetEntry(i);
    h_chi2->Fill(recon_indx);
    h_bdt->Fill(recon_indx_bdt);
  }

  std::cout << "\n=== Summary ===\n";
  std::cout << "χ² pairing:  " << h_chi2->GetBinContent(3) << " events with 2 correct photons ("
            << 100.*h_chi2->GetBinContent(3)/nentries << "%)\n";
  std::cout << "BDT pairing: " << h_bdt->GetBinContent(3)  << " events with 2 correct photons ("
            << 100.*h_bdt->GetBinContent(3)/nentries << "%)\n";

  // 1D comparison
  TCanvas* c1 = new TCanvas("c1","Comparison",900,900);
  c1->SetBottomMargin(0.12);
  c1->SetLeftMargin(0.12);

  const double ymax = h_bdt->GetMaximum();

  h_chi2->GetXaxis()->CenterTitle();
  h_chi2->GetYaxis()->CenterTitle();
  h_chi2->GetXaxis()->SetLabelSize(0.05);
  h_chi2->GetYaxis()->SetLabelSize(0.04);
  h_chi2->GetXaxis()->SetTitleSize(0.05);
  h_chi2->GetYaxis()->SetTitleSize(0.05);
  h_chi2->GetXaxis()->SetTitleOffset(1.0);
  h_chi2->GetYaxis()->SetTitleOffset(1.2);
  
  h_chi2->SetLineWidth(2);
  h_chi2->GetYaxis()->SetRangeUser(0, ymax * 1.1);
  h_chi2->Draw("hist");
  h_bdt->SetLineWidth(2);
  h_bdt->Draw("hist same");
  TLegend* leg = new TLegend(0.2,0.7,0.6,0.88);
  leg->AddEntry(h_chi2,"#chi^{2} pairing","f");
  leg->AddEntry(h_bdt,"BDT pairing","f");
  leg->Draw();
  c1->SaveAs("../plots_select/pi0gg_recon_compare.pdf");

  // 2D correlation (absolute counts)
  TH2D* h_corr = new TH2D("h_corr",";BDT-selected correct #pi^{0} photons;#chi^{2}-selected correct #pi^{0} photons",3,0,3,3,0,3);
  for (Long64_t i = 0; i < nentries; ++i) {
    tree->GetEntry(i);
    h_corr->Fill(recon_indx_bdt, recon_indx);
  }
  // Integer bin labels
  for (int i = 1; i <= 3; ++i) {
    h_corr->GetXaxis()->SetBinLabel(i, Form("%d", i-1));
    h_corr->GetYaxis()->SetBinLabel(i, Form("%d", i-1));
  }

  // ========== ABSOLUTE COUNTS (manual text to avoid .000) ==========
  TCanvas* c2 = new TCanvas("c2","Absolute counts",900,900);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.18);
  h_corr->SetStats(0);
  h_corr->GetXaxis()->CenterTitle();
  h_corr->GetYaxis()->CenterTitle();
  h_corr->Draw("colz");                 // draw colour map only (no text)
  // Manually draw integer bin contents
  gStyle->SetTextSize(0.03);
  for (int i = 1; i <= h_corr->GetNbinsX(); ++i) {
    for (int j = 1; j <= h_corr->GetNbinsY(); ++j) {
      double val = h_corr->GetBinContent(i, j);
      if (val == 0) continue;
      double x = h_corr->GetXaxis()->GetBinCenter(i);
      double y = h_corr->GetYaxis()->GetBinCenter(j);
      TText *t = new TText(x, y, Form("%.0f", val));
      t->SetTextAlign(22);
      t->SetTextSize(0.05);
      t->SetTextColor(kBlack);
      t->Draw();
    }
  }
  c2->SaveAs("../plots_select/pi0gg_recon_correlation_absolute.png");

  // ========== FRACTION BY ROW ==========
  TH2D* h_frac = (TH2D*)h_corr->Clone("h_frac");
  for (int ix = 1; ix <= h_frac->GetNbinsX(); ++ix) {
    double row_sum = 0;
    for (int iy = 1; iy <= h_frac->GetNbinsY(); ++iy)
      row_sum += h_frac->GetBinContent(ix, iy);
    if (row_sum > 0) {
      for (int iy = 1; iy <= h_frac->GetNbinsY(); ++iy) {
	double val = h_frac->GetBinContent(ix, iy);
	h_frac->SetBinContent(ix, iy, val / row_sum);
      }
    }
  }

  
  TCanvas* c3 = new TCanvas("c3", "Fraction by row", 900, 900);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.18);
  h_frac->SetStats(0);
  h_frac->GetXaxis()->CenterTitle();
  h_frac->GetYaxis()->CenterTitle();
  // Draw colour map only (no automatic text)
  h_frac->Draw("colz");
  
  // Manually draw bin contents as text with desired size
  double textSize = 0.05;   // adjust as needed
  for (int i = 1; i <= h_frac->GetNbinsX(); ++i) {
    for (int j = 1; j <= h_frac->GetNbinsY(); ++j) {
      double val = h_frac->GetBinContent(i, j);
      if (val == 0) continue;
      double x = h_frac->GetXaxis()->GetBinCenter(i);
        double y = h_frac->GetYaxis()->GetBinCenter(j);
        TText *t = new TText(x, y, Form("%.2f", val));  // two decimal places
        t->SetTextAlign(22);
        t->SetTextSize(textSize);
        t->SetTextColor(kBlack);
        t->Draw();
    }
  }
  
  c3->SaveAs("../plots_select/pi0gg_recon_correlation_fraction.png");
  
}
