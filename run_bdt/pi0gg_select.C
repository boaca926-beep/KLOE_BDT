#include "../header_bdt/binning.h"

/*
Why BDT correct is not suitable for background

    In signal MC, when the BDT correctly identifies both π⁰ photons (recon_indx_bdt == 2), the event is almost always a true signal event with the correct pairing.

    The mass distribution m3pi_bdt from these events shows a clear ω peak around 780 MeV, not the smooth, falling shape expected from combinatorial background.

    Using the BDT‑correct sample as a background template would artificially inject the signal peak into the background model, ruining any fit that attempts to separate signal and background.
*/

// pi0gg_select.C – compare χ² and BDT π⁰ photon identification
void pt_style(TPaveText *pt, TString text) {
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.04);
  pt->SetTextFont(42);
  pt->AddText(text);
}

void pi0gg_select() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFitFormat("6.4g");

  const char* data_filename = "/home/kloe/Desktop/input_bdt_TDATA_norm/cut/tree_pre_bdt.root";
  const char* tree_name = "TISR3PI_SIG";
    
  TFile* file = TFile::Open(data_filename);
  if (!file || file->IsZombie()) return;
  TTree* tree = (TTree*)file->Get(tree_name);
  if (!tree) return;
  std::cout << "Tree " << tree_name << " has " << tree->GetEntries() << " entries." << std::endl;

  int recon_indx = -1, recon_indx_bdt = -1;
  int bkg_indx = -1;
  double m3pi_bdt = 0.0, IM3pi_7C = 0.0;
  tree->SetBranchAddress("Br_bkg_indx", &bkg_indx);
  tree->SetBranchAddress("Br_recon_indx", &recon_indx);
  tree->SetBranchAddress("Br_recon_indx_bdt", &recon_indx_bdt);
  tree->SetBranchAddress("Br_m3pi_bdt", &m3pi_bdt);
  tree->SetBranchAddress("Br_IM3pi_7C", &IM3pi_7C);
      
  TH1D* h_chi2 = new TH1D("h_chi2", "; Number of correct-selected #pi^{0} photons;Events", 3, 0, 3);
  TH1D* h_bdt  = new TH1D("h_bdt",  ";  Number of correct-selected #pi^{0} photons (BDT pairing);Events", 3, 0, 3);
  h_chi2->SetLineColor(kBlue);
  h_bdt->SetLineColor(kRed);
  h_chi2->SetFillStyle(3001);
  h_bdt->SetFillStyle(3001);

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
  double purity_chi2 = 100.*h_chi2->GetBinContent(3)/nentries;
  double purity_bdt = 100.*h_bdt->GetBinContent(3)/nentries; 
    
  std::cout << "χ² pairing:  " << h_chi2->GetBinContent(3) << " events with 2 correct photons ("
            << purity_chi2 << "%)\n";
  std::cout << "BDT pairing: " << h_bdt->GetBinContent(3)  << " events with 2 correct photons ("
            << purity_bdt << "%)\n";

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

  TPaveText *pt = new TPaveText(0.15, 0.6, 0.5, 0.68, "NDC");
  TPaveText *pt1 = new TPaveText(0.15, 0.67, 0.5, 0.72, "NDC");
  
  TString line = Form("Purity #chi^{2}-selection = %.1f%%", purity_chi2);
  TString line1 = Form("Purity BDT-selection = %.1f%%", purity_bdt);

  pt_style(pt, line);
  pt->SetTextColor(kBlue);
        
  pt_style(pt1, line1);
  pt1->SetTextColor(kRed);
  
  pt->Draw();
  pt1->Draw("same");

  TLegend* leg = new TLegend(0.2,0.8,0.5,0.88);
  leg->SetTextSize(0.03);
  leg->SetBorderSize(0);
  leg->AddEntry(h_chi2,"#chi^{2} pairing","l");
  leg->AddEntry(h_bdt,"BDT pairing","l");
  leg->Draw();
  c1->SaveAs(Form("../plots_select/pi0gg_recon_compare_%s.pdf", tree_name));

  // Diagonal correct and wrong
  TH1D* h_diag_correct = new TH1D("h_diag_correct","BDT-selected purely correct #pi^{0} photons (M_{3#pi})",NBINS,MASS_MIN,MASS_MAX);
  TH1D* h_diag_wrong   = new TH1D("h_diag_wrong",  "BDT-selected purely wrong #pi^{0} photons (M_{3#pi})",NBINS,MASS_MIN,MASS_MAX);

  // Off-diagonal histograms (renamed for clarity)
  TH1D* h_off_diag_bdt2 = new TH1D("h_off_diag_bdt2",    "BDT = 2 correct (off-diagonal)", NBINS, MASS_MIN, MASS_MAX); // events where BDT is correct (recon_indx_bdt == 2) but χ² disagrees (recon_indx != 2).
  TH1D* h_off_diag_bdt_lt2 = new TH1D("h_off_diag_bdt_lt2", "BDT < 2 correct (off-diagonal)", NBINS, MASS_MIN, MASS_MAX);
  
  // off-diagonal BDT-wrong only (recon_indx_bdt < 2 and recon_indx > recon_indx_bdt)
  TH1D* h_offdiag_bdt_wrong = new TH1D("h_offdiag_bdt_wrong", "Off-diag BDT wrong (recon<2 & recon>recon_bdt)", NBINS, MASS_MIN, MASS_MAX); // events where BDT is wrong (recon_indx_bdt < 2) and χ² correctness is higher (recon_indx > recon_indx_bdt).
  
  // 2D correlation (3pi invariant mass)
  TH2D* h_corr_m3pi_bdt_chi2 = new TH2D("h_corr_m3pi_bdt_chi2",";BDT M_{3#pi} [MeV]; #chi^{2} M_{3#pi} [MeV]",NBINS,MASS_MIN,MASS_MAX,NBINS,MASS_MIN,MASS_MAX);
  
  // 2D correlation (absolute counts)
  TH2D* h_corr = new TH2D("h_corr",";BDT-selected correct #pi^{0} photons;#chi^{2}-selected correct #pi^{0} photons",3,0,3,3,0,3);
  for (Long64_t i = 0; i < nentries; ++i) {
    tree->GetEntry(i);
    h_corr_m3pi_bdt_chi2->Fill(m3pi_bdt, IM3pi_7C);

    // Diagonal vs off-diagonal
    if (recon_indx_bdt == recon_indx) {
      if (recon_indx_bdt == 2)
        h_diag_correct->Fill(m3pi_bdt);
      else
        h_diag_wrong->Fill(m3pi_bdt);
    } else {
      if (recon_indx_bdt == 2) {
        h_off_diag_bdt2->Fill(m3pi_bdt);
      } else if (recon_indx_bdt < 2) {
        h_off_diag_bdt_lt2->Fill(m3pi_bdt);
        // NEW: fill off-diagonal BDT-wrong
        if (recon_indx > recon_indx_bdt)
          h_offdiag_bdt_wrong->Fill(m3pi_bdt);
      }
    }
    
    h_corr->Fill(recon_indx_bdt, recon_indx);
  }
  // Integer bin labels
  for (int i = 1; i <= 3; ++i) {
    h_corr->GetXaxis()->SetBinLabel(i, Form("%d", i-1));
    h_corr->GetYaxis()->SetBinLabel(i, Form("%d", i-1));
  }
  
  // ========== ABSOLUTE COUNTS ==========
  TCanvas* c2 = new TCanvas("c2","Absolute counts",900,900);
  c2->cd();
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.18);
  h_corr->SetStats(0);
  h_corr->GetXaxis()->CenterTitle();
  h_corr->GetYaxis()->CenterTitle();
  h_corr->Draw("colz");
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
  c2->SaveAs(Form("../plots_select/pi0gg_recon_correlation_absolute_%s.png", tree_name));

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
  c3->cd();
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.18);
  h_frac->SetStats(0);
  h_frac->GetXaxis()->CenterTitle();
  h_frac->GetYaxis()->CenterTitle();
  h_frac->Draw("colz");
  
  double textSize = 0.05;
  for (int i = 1; i <= h_frac->GetNbinsX(); ++i) {
    for (int j = 1; j <= h_frac->GetNbinsY(); ++j) {
      double val = h_frac->GetBinContent(i, j);
      if (val == 0) continue;
      double x = h_frac->GetXaxis()->GetBinCenter(i);
      double y = h_frac->GetYaxis()->GetBinCenter(j);
      TText *t = new TText(x, y, Form("%.2f", val));
      t->SetTextAlign(22);
      t->SetTextSize(textSize);
      t->SetTextColor(kBlack);
      t->Draw();
    }
  }
  c3->SaveAs(Form("../plots_select/pi0gg_recon_correlation_fraction_%s.pdf", tree_name));

  // ========== CORRELATION PLOTS ==========
  TCanvas* c4 = new TCanvas("c4","correlation matrix; BDT vs chi2 selection (omega peak, correct paired pi0 photons)",900,900);
  c4->Divide(2,2);
  
  c4->cd(1);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.12);
  h_corr_m3pi_bdt_chi2->SetStats(0);
  h_corr_m3pi_bdt_chi2->GetXaxis()->CenterTitle();
  h_corr_m3pi_bdt_chi2->GetYaxis()->CenterTitle();
  h_corr_m3pi_bdt_chi2->Draw("colz");
  gPad->SetLogz();
  
  c4->cd(2);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.12);
  // Off-diagonal shapes (original)
  h_off_diag_bdt2->SetLineColor(kBlue);
  h_off_diag_bdt_lt2->SetLineColor(kRed);
  h_off_diag_bdt2->Draw("hist");
  h_off_diag_bdt_lt2->Draw("same hist");
  TLegend *leg_off = new TLegend(0.65,0.7,0.9,0.85);
  leg_off->AddEntry(h_off_diag_bdt2,"BDT = 2 correct (off-diag)","l");
  leg_off->AddEntry(h_off_diag_bdt_lt2,"BDT < 2 correct (off-diag)","l");
  leg_off->Draw();

  gPad->Update();

  std::cout << "h_diag_correct entries: " << h_diag_correct->Integral() << std::endl;
  std::cout << "h_diag_wrong entries: " << h_diag_wrong->Integral() << std::endl;
  std::cout << "h_off_diag_bdt2 entries: " << h_off_diag_bdt2->Integral() << std::endl;
  std::cout << "h_off_diag_bdt_lt2 entries: " << h_off_diag_bdt_lt2->Integral() << std::endl;
  std::cout << "h_offdiag_bdt_wrong entries: " << h_offdiag_bdt_wrong->Integral() << std::endl;

  c4->cd(3);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.12);
  h_diag_correct->SetStats(0);
  h_diag_correct->GetXaxis()->CenterTitle();
  h_diag_correct->GetYaxis()->CenterTitle();
  h_diag_correct->Draw("hist");
  TLegend *leg3 = new TLegend(0.65,0.7,0.9,0.85);
  leg3->AddEntry(h_diag_correct,"Purely correct (diag)","l");
  leg3->Draw();
 
  c4->cd(4);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.12);
  h_diag_wrong->SetStats(0);
  h_diag_wrong->GetXaxis()->CenterTitle();
  h_diag_wrong->GetYaxis()->CenterTitle();
  h_diag_wrong->Draw("hist");
  h_off_diag_bdt_lt2->SetLineColor(kRed);
  h_off_diag_bdt_lt2->Draw("hist same");
  
  TLegend *leg4 = new TLegend(0.65,0.7,0.9,0.85);
  leg4->AddEntry(h_diag_wrong,"Purely wrong (diag)","l");
  leg4->AddEntry(h_diag_correct,"BDT < 2 correct (off-diag)","l");
  leg4->Draw();

  c4->SaveAs("../plots_select/pi0gg_correlation_matrix.pdf");

  // Optional: off-diagonal normalized (original)
  TCanvas* c5 = new TCanvas("c5","Off-diagonal M_{3#pi} (normalized)",900,600);
  TH1D* h_off_diag_bdt2_norm = (TH1D*)h_off_diag_bdt2->Clone("h_off_diag_bdt2_norm");
  TH1D* h_off_diag_bdt_lt2_norm = (TH1D*)h_off_diag_bdt_lt2->Clone("h_off_diag_bdt_lt2_norm");
  h_off_diag_bdt2_norm->Scale(1.0/h_off_diag_bdt2_norm->Integral());
  h_off_diag_bdt_lt2_norm->Scale(1.0/h_off_diag_bdt_lt2_norm->Integral());
  //h_off_diag_bdt2_norm->SetLineColor(kBlue);
  h_off_diag_bdt_lt2_norm->SetLineColor(kRed);
  h_off_diag_bdt2_norm->GetYaxis()->SetTitle("Normalized events");
  h_off_diag_bdt2_norm->Draw("hist");
  h_off_diag_bdt_lt2_norm->Draw("same hist");
  TLegend *leg2 = new TLegend(0.65,0.7,0.85,0.85);
  leg2->AddEntry(h_off_diag_bdt2_norm,"BDT = 2 correct (off-diag)","l");
  leg2->AddEntry(h_off_diag_bdt_lt2_norm,"BDT < 2 correct (off-diag)","l");
  leg2->Draw();
  c5->SaveAs("../plots_select/pi0gg_off_diagonal_m3pi_normalized.pdf");

  // ========== NEW: Off-diagonal BDT-wrong shape for template ==========
  TCanvas* c6 = new TCanvas("c6","Off-diagonal BDT-wrong combinatorial template",900,600);
  h_offdiag_bdt_wrong->SetLineColor(kMagenta);
  h_offdiag_bdt_wrong->SetFillStyle(3001);
  h_offdiag_bdt_wrong->SetFillColor(kMagenta);
  h_offdiag_bdt_wrong->Draw("hist");
  c6->SaveAs("../plots_select/pi0gg_offdiag_bdt_wrong_template.pdf");

  // Save the new histogram to a ROOT file for template fits
  TFile* fout = new TFile("../plots_select/combinatorial_template_offdiag_wrong.root", "RECREATE");
  h_offdiag_bdt_wrong->Write();
  fout->Close();

  // Print yield for information
  std::cout << "\n=== Off-diagonal BDT-wrong yield ===\n";
  std::cout << "Events in h_offdiag_bdt_wrong: " << h_offdiag_bdt_wrong->Integral() << std::endl;
}
