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
  pt->SetTextSize(0.035);
  pt->SetTextFont(42);
  pt->AddText(text);
}

void pi0gg_select() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFitFormat("6.4g");

  const char* data_filename = "/home/bo/Desktop/input_bdt_TDATA_chain/cut/tree_pre_bdt.root";
  const char* tree_name = "TISR3PI_SIG";
    
  TFile* file = TFile::Open(data_filename);
  if (!file || file->IsZombie()) return;
  TTree* tree = (TTree*)file->Get(tree_name);
  if (!tree) return;
  std::cout << "Tree " << tree_name << " has " << tree->GetEntries() << " entries." << std::endl;

  int recon_indx = -1, recon_indx_bdt = -1, total_recon_quality = -1;
  int bkg_indx = -1;
  double m3pi_bdt = 0.0, IM3pi_7C = 0.0;
  double ppIM = 0.;
  tree->SetBranchAddress("Br_bkg_indx", &bkg_indx);
  tree->SetBranchAddress("Br_recon_indx", &recon_indx);
  tree->SetBranchAddress("Br_recon_indx_bdt", &recon_indx_bdt);
  tree->SetBranchAddress("Br_total_recon_quality", &total_recon_quality);
  tree->SetBranchAddress("Br_m3pi_bdt", &m3pi_bdt);
  tree->SetBranchAddress("Br_IM3pi_7C", &IM3pi_7C);
  tree->SetBranchAddress("Br_ppIM", &ppIM);
      
  // ===== PPIM BINNING (for background templates) =====
  const int PPIM_BINS = 100;
  const double PPIM_MIN = 300;
  const double PPIM_MAX = 650;
  
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

  // Diagonal correct and wrong (using ppIM)
  TH1D* h_diag_correct = new TH1D("h_diag_correct","BDT-selected purely correct #pi^{0} photons (ppIM)", PPIM_BINS, PPIM_MIN, PPIM_MAX);
  TH1D* h_diag_wrong   = new TH1D("h_diag_wrong",  "BDT-selected purely wrong #pi^{0} photons (ppIM)", PPIM_BINS, PPIM_MIN, PPIM_MAX);

  // Off-diagonal histograms (using ppIM)
  TH1D* h_off_diag_bdt2 = new TH1D("h_off_diag_bdt2",    "BDT = 2 correct (off-diagonal)", PPIM_BINS, PPIM_MIN, PPIM_MAX);
  TH1D* h_off_diag_bdt_lt2 = new TH1D("h_off_diag_bdt_lt2", "BDT < 2 correct (off-diagonal)", PPIM_BINS, PPIM_MIN, PPIM_MAX);
  TH1D* h_offdiag_bdt_wrong = new TH1D("h_offdiag_bdt_wrong", "Off-diag BDT wrong (recon<2 & recon>recon_bdt)", PPIM_BINS, PPIM_MIN, PPIM_MAX);

  // ===== BACKGROUND TEMPLATES (using ppIM) =====
  // PURE: Both methods fail - USE THIS FOR FITS
  TH1D* h_background_pure = new TH1D("h_background_pure", 
      "Pure background (both methods fail: recon<2 && recon_bdt<2) [ppIM]", 
      PPIM_BINS, PPIM_MIN, PPIM_MAX);
  h_background_pure->SetLineColor(kBlack);
  h_background_pure->SetLineWidth(2);
  h_background_pure->SetFillStyle(3001);
  h_background_pure->SetFillColor(kGray);

  // CONTAMINATED: BDT < 2 off-diagonal - DO NOT USE FOR FITS
  TH1D* h_background_contaminated = new TH1D("h_background_contaminated", 
      "Contaminated (BDT<2 & off-diagonal) - DO NOT USE [ppIM]", 
      PPIM_BINS, PPIM_MIN, PPIM_MAX);
  h_background_contaminated->SetLineColor(kRed);
  h_background_contaminated->SetLineWidth(2);
  h_background_contaminated->SetLineStyle(2);

  // ===== PEAK IDENTIFICATION HISTOGRAMS =====
  // Diagonal components (both fail and agree)
  TH1D* h_diag_00 = new TH1D("h_diag_00", "(0,0): both see 0 correct [ppIM]", PPIM_BINS, PPIM_MIN, PPIM_MAX);
  h_diag_00->SetLineColor(kBlue);
  h_diag_00->SetLineWidth(2);

  TH1D* h_diag_11 = new TH1D("h_diag_11", "(1,1): both see 1 correct [ppIM]", PPIM_BINS, PPIM_MIN, PPIM_MAX);
  h_diag_11->SetLineColor(kGreen);
  h_diag_11->SetLineWidth(2);

  // Off-diagonal components (both fail but disagree)
  TH1D* h_offdiag_01 = new TH1D("h_offdiag_01", "(0,1): BDT=0, χ²=1 [ppIM]", PPIM_BINS, PPIM_MIN, PPIM_MAX);
  h_offdiag_01->SetLineColor(kRed);
  h_offdiag_01->SetLineWidth(2);

  TH1D* h_offdiag_10 = new TH1D("h_offdiag_10", "(1,0): BDT=1, χ²=0 [ppIM]", PPIM_BINS, PPIM_MIN, PPIM_MAX);
  h_offdiag_10->SetLineColor(kOrange);
  h_offdiag_10->SetLineWidth(2);

  // χ² correct reference (to confirm no peak)
  TH1D* h_chi2_correct = new TH1D("h_chi2_correct", "χ²=2 correct (any BDT) [ppIM]", PPIM_BINS, PPIM_MIN, PPIM_MAX);
  h_chi2_correct->SetLineColor(kMagenta);
  h_chi2_correct->SetLineWidth(2);
  h_chi2_correct->SetLineStyle(2);
  
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
        h_diag_correct->Fill(ppIM);
      else
        h_diag_wrong->Fill(ppIM);
    } else {
      if (recon_indx_bdt == 2) {
        h_off_diag_bdt2->Fill(ppIM);
      } else if (recon_indx_bdt < 2) {
        h_off_diag_bdt_lt2->Fill(ppIM);
        if (recon_indx > recon_indx_bdt)
          h_offdiag_bdt_wrong->Fill(ppIM);
      }
    }

    // ===== FILL BACKGROUND TEMPLATES =====
    //if (!(total_recon_quality == 3)) {
    if (!(recon_indx_bdt == 2 && bkg_indx == 1)) {
      // PURE: Both methods fail - USE THIS FOR FITS
      if (recon_indx_bdt < 2 && recon_indx < 2) {
        //h_background_pure->Fill(ppIM);

        // Identify which component contributes to the peak
        if (recon_indx_bdt == recon_indx) {
          // DIAGONAL: both agree they're wrong
          if (recon_indx_bdt == 0 && recon_indx == 0) {
            h_diag_00->Fill(ppIM);
          } else if (recon_indx_bdt == 1 && recon_indx == 1) {
            h_diag_11->Fill(ppIM);
          }
        } else {
          // OFF-DIAGONAL: both fail but disagree
          if (recon_indx_bdt == 0 && recon_indx == 1) {
            h_offdiag_01->Fill(ppIM);
          } else if (recon_indx_bdt == 1 && recon_indx == 0) {
            h_offdiag_10->Fill(ppIM);
          }
        }
      }

      // CONTAMINATED: BDT < 2 off-diagonal - DO NOT USE
      if (recon_indx_bdt < 2 && recon_indx != recon_indx_bdt) {
	h_background_contaminated->Fill(ppIM);
      }
    }

    // Fill χ² correct for reference
    if (recon_indx == 2) {
      h_chi2_correct->Fill(ppIM);
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
  h_off_diag_bdt2->SetLineColor(kBlue);
  h_off_diag_bdt_lt2->SetLineColor(kRed);
  h_off_diag_bdt2->Draw("hist");
  h_off_diag_bdt_lt2->Draw("same hist");
  TLegend *leg1 = new TLegend(0.65,0.7,0.9,0.85);
  leg1->AddEntry(h_off_diag_bdt2,"BDT = 2 correct (off-diag)","l");
  leg1->AddEntry(h_off_diag_bdt_lt2,"BDT < 2 correct (off-diag)","l");
  leg1->Draw();

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
  h_off_diag_bdt_lt2_norm->SetLineColor(kRed);
  h_off_diag_bdt2_norm->GetYaxis()->SetTitle("Normalized events");
  h_off_diag_bdt2_norm->Draw("hist");
  h_off_diag_bdt_lt2_norm->Draw("same hist");
  TLegend *leg2 = new TLegend(0.65,0.7,0.85,0.85);
  leg2->AddEntry(h_off_diag_bdt2_norm,"BDT = 2 correct (off-diag)","l");
  leg2->AddEntry(h_off_diag_bdt_lt2_norm,"BDT < 2 correct (off-diag)","l");
  leg2->Draw();
  c5->SaveAs("../plots_select/pi0gg_off_diagonal_m3pi_normalized.pdf");

  // ========== Off-diagonal BDT-wrong shape for template ==========
  TCanvas* c6 = new TCanvas("c6","Off-diagonal BDT-wrong combinatorial template",900,600);
  h_offdiag_bdt_wrong->SetLineColor(kMagenta);
  h_offdiag_bdt_wrong->SetFillStyle(3001);
  h_offdiag_bdt_wrong->SetFillColor(kMagenta);
  h_offdiag_bdt_wrong->Draw("hist");
  c6->SaveAs("../plots_select/pi0gg_offdiag_bdt_wrong_template.pdf");

  // ========== Pure background template (ppIM) ==========
  TCanvas* c7 = new TCanvas("c7","Pure Background Template (both methods fail) [ppIM]",900,600);
  h_background_pure->SetLineColor(kBlack);
  h_background_pure->SetLineWidth(2);
  h_background_pure->SetFillStyle(3001);
  h_background_pure->SetFillColor(kGray);
  h_background_pure->GetXaxis()->SetTitle("ppIM [MeV]");
  h_background_pure->Draw("hist");
  c7->SaveAs("../plots_select/background_pure_ppIM.pdf");

  // ========== Comparison: Pure vs Contaminated ==========
  TCanvas* c8 = new TCanvas("c8","Background Template Comparison: PURE vs CONTAMINATED (ppIM)",900,600);
  
  // Draw highest statistics first
  if (h_background_contaminated->Integral() > h_background_pure->Integral()) {
    h_background_contaminated->Draw("hist");
    h_background_pure->Draw("hist same");
  } else {
    h_background_pure->Draw("hist");
    h_background_contaminated->Draw("hist same");
  }
  
  h_background_pure->GetXaxis()->SetTitle("ppIM [MeV]");
  
  TLegend* leg_compare = new TLegend(0.55, 0.7, 0.9, 0.85);
  leg_compare->AddEntry(h_background_pure, "PURE (both fail) - USE THIS", "l");
  leg_compare->AddEntry(h_background_contaminated, "CONTAMINATED (has peak)", "l");
  leg_compare->Draw();
  
  c8->SaveAs("../plots_select/background_pure_vs_contaminated_ppIM.pdf");

  // ========== PEAK IDENTIFICATION (NO NORMALIZATION) ==========
  TCanvas* c_peak = new TCanvas("c_peak", "Peak Source Identification", 1200, 900);
  c_peak->Divide(2,2);

  // Panel 1: All components (raw counts)
  c_peak->cd(1);
  gPad->SetLeftMargin(0.12);
  
  // Draw highest statistics first
  double int_00 = h_diag_00->Integral();
  double int_11 = h_diag_11->Integral();
  double int_01 = h_offdiag_01->Integral();
  double int_10 = h_offdiag_10->Integral();
  
  // Find max and draw first
  double max_int = TMath::Max(TMath::Max(int_00, int_11), TMath::Max(int_01, int_10));
  
  if (int_00 == max_int) {
    h_diag_00->Draw("hist");
    h_diag_11->Draw("hist same");
    h_offdiag_01->Draw("hist same");
    h_offdiag_10->Draw("hist same");
  } else if (int_11 == max_int) {
    h_diag_11->Draw("hist");
    h_diag_00->Draw("hist same");
    h_offdiag_01->Draw("hist same");
    h_offdiag_10->Draw("hist same");
  } else if (int_01 == max_int) {
    h_offdiag_01->Draw("hist");
    h_diag_00->Draw("hist same");
    h_diag_11->Draw("hist same");
    h_offdiag_10->Draw("hist same");
  } else {
    h_offdiag_10->Draw("hist");
    h_diag_00->Draw("hist same");
    h_diag_11->Draw("hist same");
    h_offdiag_01->Draw("hist same");
  }
  
  h_diag_00->SetLineColor(kBlue);
  h_diag_11->SetLineColor(kGreen);
  h_offdiag_01->SetLineColor(kRed);
  h_offdiag_10->SetLineColor(kOrange);
  
  h_diag_00->GetXaxis()->SetTitle("ppIM [MeV]");
  h_diag_00->GetYaxis()->SetTitle("Events");
  
  TLegend* leg_peak = new TLegend(0.55, 0.6, 0.9, 0.85);
  leg_peak->AddEntry(h_diag_00, Form("(0,0): both see 0 (%.0f)", int_00), "l");
  leg_peak->AddEntry(h_diag_11, Form("(1,1): both see 1 (%.0f)", int_11), "l");
  leg_peak->AddEntry(h_offdiag_01, Form("(0,1): BDT=0, χ²=1 (%.0f)", int_01), "l");
  leg_peak->AddEntry(h_offdiag_10, Form("(1,0): BDT=1, χ²=0 (%.0f)", int_10), "l");
  leg_peak->Draw();

  // Panel 2: Diagonal (0,0) vs (1,1)
  c_peak->cd(2);
  gPad->SetLeftMargin(0.12);
  
  if (int_00 > int_11) {
    h_diag_00->Draw("hist");
    h_diag_11->Draw("hist same");
  } else {
    h_diag_11->Draw("hist");
    h_diag_00->Draw("hist same");
  }
  
  h_diag_00->GetXaxis()->SetTitle("ppIM [MeV]");
  h_diag_00->GetYaxis()->SetTitle("Events");

  TLegend* leg_diag = new TLegend(0.55, 0.7, 0.9, 0.85);
  leg_diag->AddEntry(h_diag_00, Form("(0,0): both see 0 (%.0f)", int_00), "l");
  leg_diag->AddEntry(h_diag_11, Form("(1,1): both see 1 (%.0f)", int_11), "l");
  leg_diag->Draw();

  // Panel 3: Off-diagonal (0,1) vs (1,0)
  c_peak->cd(3);
  gPad->SetLeftMargin(0.12);
  
  if (int_01 > int_10) {
    h_offdiag_01->Draw("hist");
    h_offdiag_10->Draw("hist same");
  } else {
    h_offdiag_10->Draw("hist");
    h_offdiag_01->Draw("hist same");
  }
  
  h_offdiag_01->GetXaxis()->SetTitle("ppIM [MeV]");
  h_offdiag_01->GetYaxis()->SetTitle("Events");

  TLegend* leg_off = new TLegend(0.55, 0.7, 0.9, 0.85);
  leg_off->AddEntry(h_offdiag_01, Form("(0,1): BDT=0, χ²=1 (%.0f)", int_01), "l");
  leg_off->AddEntry(h_offdiag_10, Form("(1,0): BDT=1, χ²=0 (%.0f)", int_10), "l");
  leg_off->Draw();

  // Panel 4: Pure vs Contaminated with χ²=2 reference (raw counts)
  c_peak->cd(4);
  gPad->SetLeftMargin(0.12);
  
  double int_pure = h_background_pure->Integral();
  double int_contam = h_background_contaminated->Integral();
  double int_chi2 = h_chi2_correct->Integral();
  
  // Draw highest statistics first
  if (int_pure >= int_contam && int_pure >= int_chi2) {
    h_background_pure->Draw("hist");
    h_background_contaminated->Draw("hist same");
    //h_chi2_correct->Draw("hist same");
  } else if (int_contam >= int_pure && int_contam >= int_chi2) {
    h_background_contaminated->Draw("hist");
    h_background_pure->Draw("hist same");
    //h_chi2_correct->Draw("hist same");
  } else {
    //h_chi2_correct->Draw("hist");
    h_background_pure->Draw("hist same");
    h_background_contaminated->Draw("hist same");
  }
  
  h_background_pure->SetLineColor(kBlack);
  h_background_contaminated->SetLineColor(kRed);
  h_chi2_correct->SetLineColor(kMagenta);
  h_chi2_correct->SetLineStyle(2);
  
  h_background_pure->GetXaxis()->SetTitle("ppIM [MeV]");
  h_background_pure->GetYaxis()->SetTitle("Events");

  TLegend* leg_ref = new TLegend(0.55, 0.6, 0.9, 0.85);
  leg_ref->AddEntry(h_background_pure, Form("Pure (both fail) (%.0f)", int_pure), "l");
  leg_ref->AddEntry(h_background_contaminated, Form("Contaminated (%.0f)", int_contam), "l");
  leg_ref->AddEntry(h_chi2_correct, Form("#chi^{2} correct (NO peak) (%.0f)", int_chi2), "l");
  leg_ref->Draw();

  c_peak->SaveAs("../plots_select/peak_source_identification.pdf");

  // ========== SAVE TEMPLATES ==========
  TFile* fout = new TFile("../plots_select/background_templates.root", "RECREATE");
  h_background_pure->Write("background_pure");              // USE THIS for final fits
  h_background_contaminated->Write("background_contaminated"); // Shows the peak (comparison only)
  h_offdiag_bdt_wrong->Write("background_bdt_wrong");       // For comparison only
  fout->Close();

  // ========== PRINT STATISTICS ==========
  std::cout << "\n=== Background Template Statistics (ppIM) ===\n";
  std::cout << "PURE (both methods fail): " 
            << h_background_pure->Integral() << " events" << std::endl;
  std::cout << "CONTAMINATED (has peak): " 
            << h_background_contaminated->Integral() << " events" << std::endl;
  std::cout << "Difference (peak events to exclude): " 
            << h_background_contaminated->Integral() - h_background_pure->Integral() 
            << " events" << std::endl;
  
  std::cout << "\n=== Peak Source Identification ===\n";
  std::cout << "Diagonal (0,0): " << h_diag_00->Integral() << " events" << std::endl;
  std::cout << "Diagonal (1,1): " << h_diag_11->Integral() << " events" << std::endl;
  std::cout << "Off-diagonal (0,1): " << h_offdiag_01->Integral() << " events" << std::endl;
  std::cout << "Off-diagonal (1,0): " << h_offdiag_10->Integral() << " events" << std::endl;
  std::cout << "χ²=2 correct: " << h_chi2_correct->Integral() << " events (NO peak)" << std::endl;
  
  std::cout << "\n=== Saved Templates ===\n";
  std::cout << "File: ../plots_select/background_templates.root" << std::endl;
  std::cout << "  - background_pure        : USE THIS (no peak)" << std::endl;
  std::cout << "  - background_contaminated: Shows the peak (for comparison)" << std::endl;
  std::cout << "  - background_bdt_wrong   : For comparison only" << std::endl;
}
