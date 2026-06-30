// Normalized 3pi invariant mass distributions using KLOE raw data.
// omega region m3pi = [760, 800] MeV/c², where the signal is dominant
// Determine signal and data 3pi mass peak, and calculate the mass shift
// The 3pi mass shift gives energy correction for the pi0 decay photons, the mass bias value is stored in file ../header/massbias.h

#include "../header_method/method.h"
#include "../header_plot/plot.h"

const TString out_dir = "../output_IM3pi_7C_kloe_raw";
const TString tree_file_nm = out_dir + "/hist_IM3pi_7C.root";
const TString var_nm = "IM3pi_7C";
const TString unit = "[MeV/c^{2}]";
const TString var_symb = "M_{3#pi}";

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

  auto getHist = [&](const char* name) -> TH1D* {
    TH1D* h = (TH1D*)tree_file->Get(name);
    if (!h) cerr << "WARNING: histogram " << name << " not found in list" << endl;
    return h;
  };
  
  TH1D *hist_data    = getHist("hist_data");
  TH1D *hist_isr3pi_sc  = getHist("hist_isr3pi_sc");

  //hist_data->Draw();
  //hist_isr3pi_sc->Draw("same hist");
  
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

    h_mass_copy->SetMarkerStyle(20);
    h_mass_copy->SetMarkerSize(0.6);
    h_mass_copy->GetYaxis()->SetTitle("Events");
    h_mass_copy->GetYaxis()->SetRangeUser(0.01, h_mass_copy->GetMaximum() * 1.6);
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
  
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.15);

  hist_data->SetMarkerStyle(20);
  hist_data->SetMarkerSize(0.8);
  
  hist_data->SetMarkerStyle(20);
  hist_data->SetMarkerSize(0.6);
  hist_data->GetYaxis()->SetTitle("Events");
  hist_data->GetYaxis()->SetRangeUser(0.01, hist_data->GetMaximum() * 1.2);
  hist_data->GetYaxis()->CenterTitle();
  hist_data->GetYaxis()->SetTitleSize(0.05);
  hist_data->GetYaxis()->SetTitleOffset(1.4);
  hist_data->GetYaxis()->SetLabelSize(0.04);
  hist_data->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
  hist_data->GetXaxis()->SetTitleSize(0.05);
  hist_data->GetXaxis()->SetTitleOffset(1.2);
  hist_data->GetXaxis()->SetLabelSize(0.04);
  hist_data->GetXaxis()->CenterTitle();
  
  hist_data->Draw("E1");
  
  hist_isr3pi_sc->Draw("HIST SAME");

  //TLine *line = new TLine(var_min, 0, var_max, 0);
  TLine *line = new TLine(0.28, 0, 0.28, 5e3);
  line->SetLineColor(2);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  //line->Draw();

  TPaveText *pt0 = new TPaveText(0.55, 0.75, 0.85, 0.82, "NDC");
  pt0->SetFillColor(0);
  pt0->SetBorderSize(0);
  pt0->SetTextAlign(12);
  pt0->SetTextSize(0.04);
  pt0->SetTextFont(42);
  pt0->AddText(Form("Mass bias = %.2f [MeV/c^{2}]", TMath::Abs(mass_bias)));
  pt0->Draw();
 
  const double ymax = hist_data->GetMaximum();
  hist_data->GetYaxis()->SetTitle("Events");
  hist_data->GetYaxis()->SetRangeUser(0.01, ymax * 1.2);
  hist_data->GetYaxis()->CenterTitle();
  hist_data->GetYaxis()->SetTitleSize(0.05);
  hist_data->GetYaxis()->SetTitleOffset(1.2);
  hist_data->GetYaxis()->SetLabelSize(0.04);
  
  //TLegend *leg = new TLegend(0.5, 0.35, 0.88, 0.9);
  TLegend *leg = new TLegend(0.15, 0.6, 0.6, 0.9);
  leg->SetTextFont(132);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetNColumns(1);
  leg->AddEntry(hist_data, "Data", "lep");
  leg->AddEntry(hist_isr3pi_sc, "#pi^{+}#pi^{-}#pi^{0}#gamma (signal)", "l");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw();
  
  c1->SaveAs(out_dir + "/data_mc_comparison_scaled_" + var_nm + "_kloe.pdf");
  cout << "Scaled plot saved to: " << out_dir << "/data_mc_comparison_scaled_" << var_nm << "_kloe.pdf" << endl;

  
  delete c1;
  delete leg;
  delete line;
  tree_file->Close();
  delete tree_file;
  
  return 0;
}
