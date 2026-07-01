// comparsion for backgrounds including etagam
#include "../header_method/method.h"
#include "../header/compr.h"   // defines var_nm, binsize, var_min, var_max
#include "../header_plot/plot.h"

int compr_kloe() {

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

  checkFile(tree_file);

  const int TLSize = 10;
  TTree *TDATA      = static_cast<TTree*>(tree_file->Get("TDATA"));
  TTree *TEEG       = static_cast<TTree*>(tree_file->Get("TEEG"));
  TTree *TOMEGAPI   = static_cast<TTree*>(tree_file->Get("TOMEGAPI"));
  TTree *TKSL       = static_cast<TTree*>(tree_file->Get("TKSL"));
  TTree *TKPM       = static_cast<TTree*>(tree_file->Get("TKPM"));
  TTree *TRHOPI     = static_cast<TTree*>(tree_file->Get("TRHOPI"));
  TTree *TETAGAM    = static_cast<TTree*>(tree_file->Get("TETAGAM"));
  TTree *TBKGREST   = static_cast<TTree*>(tree_file->Get("TBKGREST"));
  TTree *TISR3PI_SIG = static_cast<TTree*>(tree_file->Get("TISR3PI_SIG"));
  
  TTree *TrList[TLSize] = {TDATA, TEEG, TOMEGAPI, TKSL, TKPM, TRHOPI, TETAGAM, TBKGREST, TISR3PI_SIG};
  const TString TrNm[TLSize] = {"data", "eeg", "omegapi", "ksl", "kpm", "rhopi", "etagam", "bkgrest", "isr3pi"};
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
  
  // MC rest
  TH1D *hist_mcrest = (TH1D*)hist_bkgrest->Clone("hist_mcrest");
  hist_mcrest->Add(hist_kpm, 1.);
  hist_mcrest->Add(hist_rhopi, 1.);
  //hist_mcrest->Add(hist_etagam, 1.);
  Hlist->Add(hist_mcrest);

  // ===== LOAD SCALING FACTORS FROM SFW2D =====
  TFile *f_sfw2d = TFile::Open(outputSfw2D + "sfw2d.root");
  double scale_eeg = 1.0, scale_isr3pi = 1.0, scale_omegapi = 1.0;
  double scale_etagam = 1.0, scale_ksl = 1.0, scale_mcrest = 1.0;
  
  if (f_sfw2d && !f_sfw2d->IsZombie()) {
    // FIXED: Use "TRESULT" instead of "fit_results"
    TTree *fitTree = (TTree*)f_sfw2d->Get("TRESULT");
    if (fitTree) {
      double feeg, fisr3pi, fomegapi, fetagam, fksl, fmcrest;
      double eeg_sfw, isr3pi_sfw, omegapi_sfw, etagam_sfw, ksl_sfw, mcrest_sfw;
      // FIXED: Add "Br_" prefix to all branch names
      fitTree->SetBranchAddress("Br_feeg", &feeg);
      fitTree->SetBranchAddress("Br_eeg_sfw", &eeg_sfw);
      fitTree->SetBranchAddress("Br_fisr3pi", &fisr3pi);
      fitTree->SetBranchAddress("Br_isr3pi_sfw", &isr3pi_sfw);
      fitTree->SetBranchAddress("Br_fomegapi", &fomegapi);
      fitTree->SetBranchAddress("Br_omegapi_sfw", &omegapi_sfw);
      fitTree->SetBranchAddress("Br_fetagam", &fetagam);
      fitTree->SetBranchAddress("Br_etagam_sfw", &etagam_sfw);
      fitTree->SetBranchAddress("Br_fksl", &fksl);
      fitTree->SetBranchAddress("Br_ksl_sfw", &ksl_sfw);
      fitTree->SetBranchAddress("Br_fmcrest", &fmcrest);
      fitTree->SetBranchAddress("Br_mcrest_sfw", &mcrest_sfw);

      fitTree->GetEntry(0);
      
      double nb_data_sum = hist_data->GetSumOfWeights();
      double nb_eeg_sum = hist_eeg->GetSumOfWeights();
      double nb_isr3pi_sum = hist_isr3pi->GetSumOfWeights();
      double nb_omegapi_sum = hist_omegapi->GetSumOfWeights();
      double nb_etagam_sum = hist_etagam->GetSumOfWeights();
      double nb_ksl_sum = hist_ksl->GetSumOfWeights();
      double nb_mcrest_sum = hist_mcrest->GetSumOfWeights();
      
      auto getscale = [](double Nd, double fra, double N) -> double {
        return (N == 0.0) ? 0.0 : Nd * fra / N;
      };
      
      scale_eeg     = eeg_sfw * 2.; //getscale(nb_data_sum, feeg,     nb_eeg_sum);
      scale_isr3pi  = isr3pi_sfw; //getscale(nb_data_sum, fisr3pi,  nb_isr3pi_sum);
      scale_omegapi = omegapi_sfw; //getscale(nb_data_sum, fomegapi, nb_omegapi_sum);
      scale_etagam = etagam_sfw; //getscale(nb_data_sum, fetagam, nb_etagam_sum);
      scale_ksl     = ksl_sfw; //getscale(nb_data_sum, fksl,     nb_ksl_sum);
      scale_mcrest  = mcrest_sfw; //getscale(nb_data_sum, fmcrest,  nb_mcrest_sum);
      
      cout << "\n=== Scaling factors from SFW2D ===" << endl;
      cout << "scale_eeg = " << scale_eeg << ", eeg_sfw = " << eeg_sfw * 2. << endl;
      cout << "scale_isr3pi = " << scale_isr3pi << ", isr3pi_sfw = " << isr3pi_sfw << endl;
      cout << "scale_omegapi = " << scale_omegapi << ", omegapi_sfw = " << omegapi_sfw << endl;
      cout << "scale_etagam = " << scale_etagam << ", etagam_sfw = " << etagam_sfw << endl;
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
  TH1D *hist_etagam_sc = (TH1D*)hist_etagam->Clone("hist_etagam_sc");
  TH1D *hist_ksl_sc      = (TH1D*)hist_ksl->Clone("hist_ksl_sc");
  TH1D *hist_mcrest_sc   = (TH1D*)hist_mcrest->Clone("hist_mcrest_sc");

  hist_eeg_sc->Scale(scale_eeg);
  hist_isr3pi_sc->Scale(scale_isr3pi);
  hist_omegapi_sc->Scale(scale_omegapi);
  hist_etagam_sc->Scale(scale_etagam);
  hist_ksl_sc->Scale(scale_ksl);
  hist_mcrest_sc->Scale(scale_mcrest);

  hist_mcrest_sc->Add(hist_etagam_sc, 1.);
  
  cout << "\n=== Summary ===" << endl;
  double hist_low = hist_data->GetXaxis()->GetXmin();
  double hist_max = hist_data->GetXaxis()->GetXmax();
  
  std::cout << "Mass range [" << hist_low << ", " << hist_max << "] MeV/c²\n";
  std::cout << "DATA: " << hist_data->Integral() << "\n"
	    << "SIGNAL: " << hist_isr3pi_sc->Integral() << ", sfw = "<< scale_isr3pi << "\n"
	    << "ETAGAM (negligiable): " << hist_etagam_sc->Integral() << "\n"
	    << "EEG: " << hist_eeg_sc->Integral() << "\n"
	    << "OMEGAPI: " << hist_omegapi_sc->Integral() << "\n"
	    << "KSL: " << hist_ksl_sc->Integral() << "\n"
	    << "MCREST: " << hist_mcrest_sc->Integral() << "\n"
    	    << std::endl;

  // Total scaled MC
  TH1D *hist_mcsum = (TH1D*)hist_eeg_sc->Clone("hist_mcsum");
  hist_mcsum->Add(hist_isr3pi_sc, 1.);
  hist_mcsum->Add(hist_omegapi_sc, 1.);
  hist_mcsum->Add(hist_etagam_sc, 1.);
  hist_mcsum->Add(hist_ksl_sc, 1.);
  hist_mcsum->Add(hist_mcrest_sc, 1.);
  hist_mcsum->SetLineColor(kRed);

  // Create output directory
  //TString out_dir = "../output_" + TString(var_nm);
  gSystem->mkdir(out_dir, kTRUE);

  // Save histograms (including scaled ones)
  TString outfile_name = out_dir + "/hist_" + var_nm + ".root";
  cout << outfile_name << endl;
  TFile *f_out = new TFile(outfile_name, "recreate");
  Hlist->Write("Hlist", TObject::kSingleKey);
  hist_eeg_sc->Write();
  hist_isr3pi_sc->Write();
  hist_omegapi_sc->Write();
  hist_etagam_sc->Write();
  hist_ksl_sc->Write();
  hist_mcrest_sc->Write();
  hist_mcsum->Write();
  hist_data->Write();
  
  f_out->Close();

  // ===== SCALED Data/MC Comparison Plot =====
  //TCanvas *c1 = new TCanvas("c1", "Data/MC Comparison (Scaled)", 1200, 700);
  TCanvas *c1 = new TCanvas("c1", "Data/MC Comparison (Scaled)", 700, 700);
  
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

  hist_data->Draw("E1");
  
  hist_data->SetMarkerStyle(20);
  hist_data->SetMarkerSize(0.8);
  hist_mcsum->Draw("HIST SAME");
  hist_mcrest_sc->Draw("HIST SAME");
  hist_ksl_sc->Draw("HIST SAME");
  hist_omegapi_sc->Draw("HIST SAME");
  hist_etagam_sc->Draw("HIST SAME");
  hist_eeg_sc->Draw("HIST SAME");
  hist_isr3pi_sc->Draw("HIST SAME");

  //TLine *line = new TLine(var_min, 0, var_max, 0);
  TLine *line = new TLine(0.28, 0, 0.28, 5e3);
  line->SetLineColor(2);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  //line->Draw();
  
  const double ymax = hist_data->GetMaximum();
  hist_data->GetYaxis()->SetTitle("Events");
  hist_data->GetYaxis()->SetRangeUser(0.01, ymax * 1.6);
  hist_data->GetYaxis()->CenterTitle();
  hist_data->GetYaxis()->SetTitleSize(0.05);
  hist_data->GetYaxis()->SetTitleOffset(1.2);
  hist_data->GetYaxis()->SetLabelSize(0.04);
  
  //TLegend *leg = new TLegend(0.5, 0.35, 0.88, 0.9);
  TLegend *leg = new TLegend(0.15, 0.35, 0.6, 0.9);
  leg->SetTextFont(132);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetNColumns(1);
  leg->AddEntry(hist_data, "Data", "lep");
  leg->AddEntry(hist_mcsum, "MC sum", "l");
  leg->AddEntry(hist_isr3pi_sc, "#pi^{+}#pi^{-}#pi^{0}#gamma (signal)", "l");
  leg->AddEntry(hist_eeg_sc, "e^{+}e^{-}#gamma", "l");
  leg->AddEntry(hist_omegapi_sc, "#omega#pi^{0}", "l");
  leg->AddEntry(hist_ksl_sc, "K_{L}K_{S}", "l");
  leg->AddEntry(hist_mcrest_sc, "MC Rest", "l");
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
  
  c1->SaveAs(out_dir + "/data_mc_comparison_" + var_nm + "_kloe.pdf");
  cout << "Scaled plot saved to: " << out_dir << "/data_mc_comparison_scaled_" << var_nm << "_kloe.pdf" << endl;

  
  delete c1;
  delete leg;
  delete line;
  delete f_out;
  delete Hlist;
  tree_file->Close();
  delete tree_file;
  
  return 0;
}
