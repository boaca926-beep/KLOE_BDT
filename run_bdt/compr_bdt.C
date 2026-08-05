// comparsion for backgrounds including etagam
#include "../header_method/method.h"
#include "../header_bdt/compr.h"   // defines var_nm, binsize, var_min, var_max
#include "../header_plot/plot.h"
//#include "../header_bdt/path.h"    // for outputSfw2D

/*
int compr_bdt(const TString tree_file_nm = "/home/bo/Desktop/bdt_tuning_TDATA_chain_false/cut/tree_pre.root",
	      const TString out_dir = "../tuning_false_m_gg_bdt",
	      const TString outputSfw2D = "/home/bo/Desktop/bdt_tuning_TDATA_chain_false/sfw2d/",
	      const TString var_nm = "m_gg_bdt",
	      const TString var_symb = "M_{#gamma#gamma}",
	      const TString unit = "",
	      
	      const int binsize = 120,
	      const double var_min = 120,
	      const double var_max = 150
) {
*/

/*
int compr_bdt(const TString tree_file_nm = "/home/bo/Desktop/bdt_tuning_TDATA_chain_false/cut/tree_pre.root",
	      const TString out_dir = "../tuning_false_m3pi_bdt",
	      const TString outputSfw2D = "/home/bo/Desktop/bdt_tuning_TDATA_chain_false/sfw2d/",
	      const TString var_nm = "m3pi_bdt",
	      const TString var_symb = "M_{3#pi}",
	      const TString unit = "",
	      
	      const int binsize = 120,
	      const double var_min = 760,
	      const double var_max = 800
) {
*/

int compr_bdt() {
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  TFile* tree_file = new TFile(tree_file_nm);
  if (!tree_file || tree_file->IsZombie()) {
    cerr << "ERROR: Cannot open " << tree_file_nm << endl;
    return 1;
  }

  const int TLSize = 10;
  TTree *TDATA      = static_cast<TTree*>(tree_file->Get("TDATA"));
  TTree *TEEG       = static_cast<TTree*>(tree_file->Get("TEEG"));
  TTree *TOMEGAPI   = static_cast<TTree*>(tree_file->Get("TOMEGAPI"));
  TTree *TKSL       = static_cast<TTree*>(tree_file->Get("TKSL"));
  TTree *TKPM       = static_cast<TTree*>(tree_file->Get("TKPM"));
  TTree *TRHOPI     = static_cast<TTree*>(tree_file->Get("TRHOPI"));
  TTree *TETAGAM    = static_cast<TTree*>(tree_file->Get("TETAGAM"));
  TTree *TBKGREST   = static_cast<TTree*>(tree_file->Get("TBKGREST"));
  TTree *TISR3PI_SIG_PEAK = static_cast<TTree*>(tree_file->Get("TISR3PI_SIG_PEAK"));
  TTree *TISR3PI_SIG_NON_RESON = static_cast<TTree*>(tree_file->Get("TISR3PI_SIG_NON_RESON"));
  
  TTree *TrList[TLSize] = {TDATA, TEEG, TOMEGAPI, TKSL, TKPM, TRHOPI, TETAGAM, TBKGREST, TISR3PI_SIG_PEAK, TISR3PI_SIG_NON_RESON};
  const TString TrNm[TLSize] = {"data", "eeg", "omegapi", "ksl", "kpm", "rhopi", "etagam", "bkgrest", "isr3pi", "nonreson"};
  int color_list[TLSize] = {1, 6, 7, 28, 46, 42, 3, 37, 4, 9};

  TList *Hlist = new TList();
  double var_value = 0.;

  for (int i = 0; i < TLSize; ++i) {
    if (!TrList[i]) {
      cerr << "WARNING: Tree " << TrNm[i] << " not found!" << endl;
      continue;
    }

    TH1D *h = new TH1D("hist_" + TrNm[i], "", binsize, var_min, var_max);
    h->Sumw2();
    TrList[i]->SetBranchAddress("Br_" + var_nm, &var_value);

    Long64_t nentries = TrList[i]->GetEntries();
    for (Long64_t irow = 0; irow < nentries; ++irow) {
      TrList[i]->GetEntry(irow);
      h->Fill(var_value);
    }

    format_h(h, color_list[i], 2);
    Hlist->Add(h);
  }
  
  auto getHist = [&](const char* name) -> TH1D* {
    TH1D* h = (TH1D*)Hlist->FindObject(name);
    if (!h) cerr << "WARNING: histogram " << name << " not found in list" << endl;
    return h;
  };

  TH1D *hist_data    = getHist("hist_data");
  TH1D *hist_eeg     = getHist("hist_eeg");
  TH1D *hist_omegapi = getHist("hist_omegapi");
  TH1D *hist_ksl     = getHist("hist_ksl");
  TH1D *hist_kpm     = getHist("hist_kpm");
  TH1D *hist_rhopi   = getHist("hist_rhopi");
  TH1D *hist_etagam  = getHist("hist_etagam");
  TH1D *hist_bkgrest = getHist("hist_bkgrest");
  TH1D *hist_isr3pi  = getHist("hist_isr3pi");
  TH1D *hist_nonreson = getHist("hist_nonreson");

  // MC rest
  TH1D *hist_mcrest = (TH1D*)hist_bkgrest->Clone("hist_mcrest");
  hist_mcrest->Add(hist_kpm, 1.);
  hist_mcrest->Add(hist_rhopi, 1.);
  hist_mcrest->Add(hist_etagam, 1.);
  Hlist->Add(hist_mcrest);

  // ===== LOAD SCALING FACTORS FROM SFW2D =====
  TFile *f_sfw2d = TFile::Open(outputSfw2D + "sfw2d.root");
  double scale_eeg = 1.0, scale_isr3pi = 1.0, scale_omegapi = 1.0;
  double scale_nonReson = 1.0, scale_ksl = 1.0, scale_mcrest = 1.0;
  
  if (f_sfw2d && !f_sfw2d->IsZombie()) {
    // FIXED: Use "TRESULT" instead of "fit_results"
    TTree *fitTree = (TTree*)f_sfw2d->Get("TRESULT");
    if (fitTree) {
      double feeg, fisr3pi, fomegapi, fnonReson, fksl, fmcrest;
      double eeg_sfw, isr3pi_sfw, omegapi_sfw, nonReson_sfw, ksl_sfw, mcrest_sfw;
      // FIXED: Add "Br_" prefix to all branch names
      fitTree->SetBranchAddress("Br_feeg", &feeg);
      fitTree->SetBranchAddress("Br_eeg_sfw", &eeg_sfw);
      fitTree->SetBranchAddress("Br_fisr3pi", &fisr3pi);
      fitTree->SetBranchAddress("Br_isr3pi_sfw", &isr3pi_sfw);
      fitTree->SetBranchAddress("Br_fomegapi", &fomegapi);
      fitTree->SetBranchAddress("Br_omegapi_sfw", &omegapi_sfw);
      fitTree->SetBranchAddress("Br_fnonReson", &fnonReson);
      fitTree->SetBranchAddress("Br_nonReson_sfw", &nonReson_sfw);
      fitTree->SetBranchAddress("Br_fksl", &fksl);
      fitTree->SetBranchAddress("Br_ksl_sfw", &ksl_sfw);
      fitTree->SetBranchAddress("Br_fmcrest", &fmcrest);
      fitTree->SetBranchAddress("Br_mcrest_sfw", &mcrest_sfw);

      fitTree->GetEntry(0);
      
      double nb_data_sum = hist_data->GetSumOfWeights();
      double nb_eeg_sum = hist_eeg->GetSumOfWeights();
      double nb_isr3pi_sum = hist_isr3pi->GetSumOfWeights();
      double nb_omegapi_sum = hist_omegapi->GetSumOfWeights();
      double nb_nonReson_sum = hist_nonreson->GetSumOfWeights();
      double nb_ksl_sum = hist_ksl->GetSumOfWeights();
      double nb_mcrest_sum = hist_mcrest->GetSumOfWeights();
      
      auto getscale = [](double Nd, double fra, double N) -> double {
        return (N == 0.0) ? 0.0 : Nd * fra / N;
      };
      
      scale_eeg     = eeg_sfw; //getscale(nb_data_sum, feeg,     nb_eeg_sum);
      scale_isr3pi  = isr3pi_sfw; //getscale(nb_data_sum, fisr3pi,  nb_isr3pi_sum);
      scale_omegapi = omegapi_sfw; //getscale(nb_data_sum, fomegapi, nb_omegapi_sum);
      scale_nonReson= nonReson_sfw; //getscale(nb_data_sum, fnonReson, nb_nonReson_sum);
      scale_ksl     = ksl_sfw; //getscale(nb_data_sum, fksl,     nb_ksl_sum);
      scale_mcrest  = mcrest_sfw; //getscale(nb_data_sum, fmcrest,  nb_mcrest_sum);
      
      cout << "\n=== Scaling factors from SFW2D ===" << endl;
      cout << "scale_eeg = " << scale_eeg << ", eeg_sfw = " << eeg_sfw << endl;
      cout << "scale_isr3pi = " << scale_isr3pi << ", isr3pi_sfw = " << isr3pi_sfw << endl;
      cout << "scale_omegapi = " << scale_omegapi << ", omegapi_sfw = " << omegapi_sfw << endl;
      cout << "scale_nonReson = " << scale_nonReson << ", nonReson_sfw = " << nonReson_sfw << endl;
      cout << "scale_ksl = " << scale_ksl << ", ksl_sfw = " << ksl_sfw << endl;
      cout << "scale_mcrest = " << scale_mcrest << ", mcrest_sfw = " << mcrest_sfw << endl;
    } else {
      cout << "WARNING: TRESULT tree not found in sfw2d.root" << endl;
    }
    f_sfw2d->Close();
  } else {
    cout << "WARNING: sfw2d.root not found, using unscaled MC" << endl;
  }
  
  // Create scaled histograms
  TH1D *hist_eeg_sc      = (TH1D*)hist_eeg->Clone("hist_eeg_sc");
  TH1D *hist_isr3pi_sc   = (TH1D*)hist_isr3pi->Clone("hist_isr3pi_sc");
  TH1D *hist_omegapi_sc  = (TH1D*)hist_omegapi->Clone("hist_omegapi_sc");
  TH1D *hist_nonreson_sc = (TH1D*)hist_nonreson->Clone("hist_nonreson_sc");
  TH1D *hist_ksl_sc      = (TH1D*)hist_ksl->Clone("hist_ksl_sc");
  TH1D *hist_mcrest_sc   = (TH1D*)hist_mcrest->Clone("hist_mcrest_sc");
  
  hist_eeg_sc->Scale(scale_eeg);
  hist_isr3pi_sc->Scale(scale_isr3pi);
  hist_omegapi_sc->Scale(scale_omegapi);
  hist_nonreson_sc->Scale(scale_nonReson);
  hist_ksl_sc->Scale(scale_ksl);
  hist_mcrest_sc->Scale(scale_mcrest);

  cout << "\n=== Summary ===" << endl;
  double peak_nb = hist_isr3pi_sc->Integral();
  double distorted_nb = hist_nonreson_sc->Integral();
  double signal_sum = peak_nb + distorted_nb;
  double purity = peak_nb / signal_sum;
  double hist_low = hist_data->GetXaxis()->GetXmin();
  double hist_max = hist_data->GetXaxis()->GetXmax();
  
  std::cout << "Mass range [" << hist_low << ", " << hist_max << "] MeV/c²\n";
  std::cout << "DATA: " << hist_data->Integral() << "\n"
	    << "SIGNAL: " << signal_sum << ", purity = " << purity * 100. << "%\n"
	    << "\tpeak: " << peak_nb << ", sfw = "<< scale_isr3pi << "\n"
	    << "\tdistorted: " << distorted_nb << "\n"
	    << "EEG: " << hist_eeg_sc->Integral() << "\n"
	    << "OMEGAPI: " << hist_omegapi_sc->Integral() << "\n"
	    << "KSL: " << hist_ksl_sc->Integral() << "\n"
	    << "MCREST: " << hist_mcrest_sc->Integral() << "\n"
    	    << std::endl;

  // Total scaled MC
  TH1D *hist_mcsum = (TH1D*)hist_eeg_sc->Clone("hist_mcsum");
  hist_mcsum->Add(hist_isr3pi_sc, 1.);
  hist_mcsum->Add(hist_omegapi_sc, 1.);
  hist_mcsum->Add(hist_nonreson_sc, 1.);
  hist_mcsum->Add(hist_ksl_sc, 1.);
  hist_mcsum->Add(hist_mcrest_sc, 1.);
  hist_mcsum->SetLineColor(kRed);

  hist_nonreson_sc->SetLineColor(kBlue);
  hist_nonreson_sc->SetLineStyle(3);

  // Bkg sum
  TH1D *hist_bkgsum = (TH1D*)hist_mcsum->Clone("hist_bkgsum");
  hist_bkgsum->Add(hist_isr3pi_sc, -1.);
  hist_bkgsum->Add(hist_nonreson_sc, -1.);
  hist_bkgsum->SetLineColor(37);

  // Create output directory
  //TString out_dir = "../output_bdt_" + TString(var_nm);
  //gSystem->mkdir(out_dir, kTRUE);

  // Save histograms (including scaled ones)
  TString outfile_name = out_dir + "/hist.root";
  cout << outfile_name << endl;
  TFile *f_out = new TFile(outfile_name, "recreate");
  Hlist->Write("Hlist", TObject::kSingleKey);
  hist_eeg_sc->Write();
  hist_isr3pi_sc->Write();
  hist_omegapi_sc->Write();
  hist_nonreson_sc->Write();
  hist_ksl_sc->Write();
  hist_mcrest_sc->Write();
  hist_mcsum->Write();
  hist_data->Write();
  
  f_out->Close();

  // ===== SCALED Data/MC Comparison Plot =====
  TCanvas *c1 = new TCanvas("c1", "Data/MC Comparison (Scaled)", 900, 900);
  c1->SetBottomMargin(0.12);
  c1->SetLeftMargin(0.12);

  c1->Divide(1, 2);

  // --- Residuals ---
  TH1D *hist_ratio = new TH1D("hist_ratio", "", binsize, var_min, var_max);
  TH1D *hist_ratio_distr = new TH1D("hist_ratio_distr", "", 200, -10, 10);

  for (int j = 1; j <= binsize; ++j) {
    double nb_data = hist_data->GetBinContent(j);
    double nb_mcsum = hist_mcsum->GetBinContent(j);
    double evnt_err = TMath::Sqrt(nb_data + nb_mcsum);
    if (evnt_err > 0) {
      double residul = (nb_data - nb_mcsum) / evnt_err;
      hist_ratio->SetBinContent(j, residul);
      hist_ratio_distr->Fill(residul);
    }
  }
  
  // Upper pad
  TPad *pad1 = (TPad*)c1->cd(1);
  pad1->SetPad(0, 0.3, 1, 1);
  pad1->SetBottomMargin(0.01);
  pad1->SetLeftMargin(0.12);

  hist_data->SetMarkerStyle(20);
  hist_data->SetMarkerSize(0.8);
  
  hist_data->Draw("E1");
  hist_mcsum->Draw("HIST SAME");
  hist_mcrest_sc->Draw("HIST SAME");
  hist_ksl_sc->Draw("HIST SAME");
  hist_omegapi_sc->Draw("HIST SAME");
  hist_eeg_sc->Draw("HIST SAME");
  hist_isr3pi_sc->Draw("HIST SAME");
  hist_nonreson_sc->Draw("HIST SAME");

  //hist_bkgsum->Draw("HIST SAME");

  const double mpi0_data = 135.118;
  Double_t maxVal = hist_data->GetMaximum();
  cout << "maxVal = " << maxVal << endl;
  
  TLine *line = new TLine(var_min, 0, var_max, 0);
  //TLine *line = new TLine(0.28, 0, 0.28, 5e3);
  //TLine *line = new TLine(mpi0_data, 0, mpi0_data, maxVal);
  
  line->SetLineColor(2);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  //line->Draw();
  
  const double ymax = hist_data->GetMaximum();
  hist_data->GetYaxis()->SetTitle("Events");
  hist_data->GetYaxis()->SetRangeUser(0.01, ymax * 1.8);
  //hist_data->GetYaxis()->SetRangeUser(3, ymax * 1.6);
  hist_data->GetYaxis()->CenterTitle();
  hist_data->GetYaxis()->SetTitleSize(0.05);
  hist_data->GetYaxis()->SetTitleOffset(1.2);
  hist_data->GetYaxis()->SetLabelSize(0.04);
  hist_data->GetYaxis()->SetLabelOffset(10.);

  TPaveText *pt0 = new TPaveText(0.65, 0.7, 0.85, 0.85, "NDC");
  pt0->SetFillColor(0);
  pt0->SetBorderSize(0);
  pt0->SetTextAlign(12);
  pt0->SetTextSize(0.04);
  pt0->SetTextFont(42);
  pt0->AddText(Form("Purity = %.1f%%", purity * 100.));
  //pt0->Draw();

  //gPad->SetLogy();

  //line->Draw();
  
  //TLegend *leg = new TLegend(0.7, 0.35, 0.9, 0.9);
  TLegend *leg = new TLegend(0.15, 0.6, 0.9, 0.9);
  leg->SetTextFont(132);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetNColumns(2);
  leg->AddEntry(hist_data, "Data", "lep");
  leg->AddEntry(hist_bkgsum, "Background", "l");
  leg->AddEntry(hist_isr3pi_sc, "#pi^{+}#pi^{-}#pi^{0}#gamma (signal)", "l");
  leg->AddEntry(hist_nonreson_sc, "Mis-recon. signal", "l");
  leg->AddEntry(hist_eeg_sc, "e^{+}e^{-}#gamma", "l");
  leg->AddEntry(hist_omegapi_sc, "#omega#pi^{0}", "l");
  leg->AddEntry(hist_ksl_sc, "K_{L}K_{S}", "l");
  leg->AddEntry(hist_mcrest_sc, "MC Rest", "l");
  leg->AddEntry(hist_mcsum, "MC sum", "l");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw();

  
  // Lower pad - ratio
  TPad *pad2 = (TPad*)c1->cd(2);
  pad2->SetPad(0, 0, 1, 0.3);
  pad2->SetTopMargin(0.05);
  pad2->SetBottomMargin(0.3);
  pad2->SetLeftMargin(0.12);
  pad2->SetGridy(1);
  
  hist_ratio->SetMarkerStyle(20);
  hist_ratio->SetMarkerSize(0.6);
  hist_ratio->GetXaxis()->SetTitle(var_symb + " " + unit);
  hist_ratio->GetXaxis()->SetTitleSize(0.12);
  hist_ratio->GetXaxis()->SetTitleOffset(1.0);
  hist_ratio->GetXaxis()->SetLabelSize(0.1);
  hist_ratio->GetXaxis()->CenterTitle();
  hist_ratio->GetYaxis()->SetTitle("Pull");
  hist_ratio->GetYaxis()->SetTitleSize(0.12);
  hist_ratio->GetYaxis()->SetTitleOffset(0.5);
  hist_ratio->GetYaxis()->SetLabelSize(0.08);
  hist_ratio->GetYaxis()->SetRangeUser(-5., 5.);
  hist_ratio->GetYaxis()->SetNdivisions(505);
  hist_ratio->GetYaxis()->CenterTitle();
  hist_ratio->Draw("EP");

  line->Draw();
  
  /*
  // Prompt plot
  TCanvas *c2 = new TCanvas("c2", "Data/MC Comparison (Scaled)", 1000, 700);
  c2->SetBottomMargin(0.12);
  c2->SetLeftMargin(0.12);

  hist_data->GetXaxis()->SetTitle(var_symb + " " + unit);
  hist_data->GetXaxis()->SetTitleSize(0.05);
  hist_data->GetXaxis()->SetTitleOffset(1.1);
  hist_data->GetXaxis()->SetLabelSize(0.05);
  hist_data->GetXaxis()->CenterTitle();
  hist_data->GetYaxis()->SetLabelSize(0.05);
  hist_data->GetYaxis()->SetTitleSize(0.06);
  hist_data->GetYaxis()->SetTitleOffset(.8);
  
  hist_data->Draw("E1");
  hist_isr3pi_sc->Draw("HIST SAME");
  //hist_nonreson_sc->Draw("HIST SAME");
  hist_bkgsum->Draw("HIST SAME");
  hist_mcsum->Draw("HIST SAME");
  line->Draw();
  
  TPaveText *pt = new TPaveText(0.15, 0.65, 0.85, 0.75, "NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.045);
  pt->SetTextFont(42);
  pt->AddText(Form("m_{#pi^{0}} = %.2f MeV/c^{2}   BDT-Threshold ('any') = %.2f", mpi0_data, 0.35));
  
  //gPad->SetLogy();
  TLegend *leg1 = new TLegend(0.15, 0.7, 0.9, 0.9);
  leg1->SetTextFont(132);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.05);
  leg1->SetNColumns(3);
  
  leg1->AddEntry(hist_data, "Data", "lep");
  leg1->AddEntry(hist_bkgsum, "Background", "l");
  leg1->AddEntry(hist_isr3pi_sc, "#pi^{+}#pi^{-}#pi^{0}#gamma (signal)", "l");
  //leg1->AddEntry(line, "Threshold", "l");
  
  leg1->Draw();
  pt->Draw();
  */
  
  // Save plots
  c1->SaveAs(out_dir + "/comparison_" + var_nm + ".pdf");
  //c2->SaveAs(out_dir + "/prompt_comparison_" + var_nm + ".pdf");
  
  cout << "Scaled plot saved to: " << out_dir << "/data_mc_comparison_scaled_" << var_nm << ".pdf" << endl;
  
  delete c1;
  //delete c2;
  delete leg;
  delete line;
  delete f_out;
  delete Hlist;
  tree_file->Close();
  delete tree_file;
  
  return 0;
}
