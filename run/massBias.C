// comparsion for backgrounds including etagam
#include "../header_method/method.h"
#include "../header_plot/plot.h"

const TString tree_file_nm = "/home/bo/Desktop/input_kloe_TDATA_chain/cut/tree_pre.root";
const TString outputSfw2D = "/home/bo/Desktop/input_kloe_TDATA_chain/sfw2d/";
const TString var_nm = "IM3pi_7C";
const TString unit = "[MeV/c^{2}]";
const TString var_symb = "M_{3#pi}";
const TString out_dir = "../output_IM3pi_7C_kloe";

const int binsize = 100;
const double var_min = 760;
const double var_max = 800;

const double IM3pi_min = 720; //760 720
const double IM3pi_max = 820; //800 620

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

int massBias() {

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
      fitTree->SetBranchAddress("Br_etagam", &fetagam);
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
  TString out_dir = "../massBias";
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

  // ------------------------------------------------------------------
  // * BW fit to determine 3pi mass peak position, mass bias [MeV/c^{2}]
  // ------------------------------------------------------------------
    
  const int nb_mass = 2;
  TH1D *hMassList[nb_mass] = {hist_isr3pi_sc, hist_data};
  TString massNameList[nb_mass] = {"MC", "Data"};
  int massColor[nb_mass] = {kRed, kRed};
  FitResult massResults[nb_mass];
  TF1 *bw_fits[nb_mass];

  // Create canvas for mass fits
  for (int i = 0; i < nb_mass; i++) {
    TH1D *h_mass = hMassList[i];
    if (!h_mass) {
      std::cerr << "Skipping mass histogram " << massNameList[i] << " (null)." << std::endl;
      continue;
    }

    TCanvas *c_mass = new TCanvas("c_mass_" + massNameList[i], "3#pi Mass Distributions (Breit-Wigner fits)", 700, 700);
    c_mass->Divide(nb_mass, 1);

    // Clone to avoid modifying original
    TH1D *h_mass_copy = (TH1D*)h_mass->Clone(Form("h_mass_%s", massNameList[i].Data()));
    h_mass_copy->SetDirectory(0);
    //h_mass_copy->SetLineColor(massColor[i]);

    double mass_mean = h_mass_copy->GetMean();
    double mass_rms = h_mass_copy->GetRMS();
    double mass_peak = h_mass_copy->GetBinContent(h_mass_copy->GetMaximumBin());
    double mass_peak_pos = h_mass_copy->GetBinCenter(h_mass_copy->GetMaximumBin());

    // Fit range: ±1.5σ around mean (like pulls) but constrained
    double fit_min_mass = mass_mean - 1. * mass_rms;
    double fit_max_mass = mass_mean + 1. * mass_rms;
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
    c_mass->SaveAs(out_dir + "/mass_fit_" + massNameList[i] + ".pdf");
    delete c_mass;
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
  double mass_bias = -(massResults[0].mean - massResults[1].mean);
  cout << "mass bias = " << mass_bias << endl;

  std::ofstream myfile;
  TString myfile_nm = "../header/massbias.h";
  myfile.open(myfile_nm.Data());
  myfile << "const double energy_shift = " << mass_bias << ";\n";
  myfile.close();

  // ===== SCALED Data/MC Comparison Plot =====
  //TCanvas *c1 = new TCanvas("c1", "Data/MC Comparison (Scaled)", 1200, 700);
  TCanvas *c1 = new TCanvas("c1", "Data/MC Comparison (Scaled)", 1400, 900);
  
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

  TPaveText *pt0 = new TPaveText(0.65, 0.72, 0.85, 0.82, "NDC");
  pt0->SetFillColor(0);
  pt0->SetBorderSize(0);
  pt0->SetTextAlign(12);
  pt0->SetTextSize(0.04);
  pt0->SetTextFont(42);
  pt0->AddText(Form("Mass bias = %.2f [MeV/c^{2}]", TMath::Abs(mass_bias)));
  pt0->Draw();
 
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
  
  c1->SaveAs(out_dir + "/data_mc_comparison_scaled_" + var_nm + "_kloe.pdf");
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
