// massBias.C – full code with background subtraction
#include "../header_method/method.h"
#include "../header_plot/plot.h"
#include "../header/path.h"   // for outputHist, tuning_type

const TString tree_file_nm = outputHist + "hist.root";
const TString out_dir = "../massBias_" + tuning_type;

const TString var_nm = "IM3pi_7C";
const TString unit = "[MeV/c^{2}]";
const TString var_symb = "M_{3#pi}";

struct FitResult {
    TString name;
    double mean, mean_err;
    double sigma, sigma_err;
    double chi2_ndf;
    int entries;
};

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

  // ---- Ensure output directory exists ----
  gSystem->mkdir(out_dir, kTRUE);

  TFile* tree_file = new TFile(tree_file_nm);
  if (!tree_file || tree_file->IsZombie()) {
    cerr << "ERROR: Cannot open " << tree_file_nm << endl;
    return 1;
  }

  checkFile(tree_file); // optional, shows contents

  // ---- Get the TList with cross-section histograms ----
  TList *list_crx = (TList*)tree_file->Get("HIM3pi_crx");
  if (!list_crx) {
    cerr << "ERROR: TList HIM3pi_crx not found!" << endl;
    tree_file->Close();
    return 1;
  }

  // ---- Retrieve raw data and signal histograms ----
  TH1D *hist_data = (TH1D*)list_crx->FindObject("h1d_IM3pi_TDATA_CRX");
  TH1D *hist_signal   = (TH1D*)list_crx->FindObject("h1d_IM3pi_TISR3PI_SIG_CRX");

  // ---- Ready for background subtraction
  
  if (!hist_data) {
    cerr << "ERROR: Data histogram h1d_IM3pi_TDATA_CRX not found!" << endl;
    tree_file->Close();
    return 1;
  }
  if (!hist_signal) {
    cerr << "ERROR: Signal histogram h1d_IM3pi_TISR3PI_SIG_CRX not found!" << endl;
    tree_file->Close();
    return 1;
  }
  
  // ------------------------------------------------------------------
  // * BW fit to determine 3pi mass peak position, mass bias [MeV/c^{2}]
  // ------------------------------------------------------------------
    
  const int nb_mass = 2;
  TH1D *hMassList[nb_mass] = {hist_signal, hist_data};
  TString massNameList[nb_mass] = {"MC", "Data"};
  int massColor[nb_mass] = {kRed, kRed};
  FitResult massResults[nb_mass];

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

    double mass_mean = h_mass_copy->GetMean();
    double mass_rms = h_mass_copy->GetRMS();
    double mass_peak = h_mass_copy->GetBinContent(h_mass_copy->GetMaximumBin());
    double mass_peak_pos = h_mass_copy->GetBinCenter(h_mass_copy->GetMaximumBin());

    // Fit range: ±1σ around mean, constrained to [760,810]
    double fit_min = TMath::Max(760., mass_mean - mass_rms);
    double fit_max = TMath::Min(810., mass_mean + mass_rms);

    std::cout << "\n========================================" << std::endl;
    std::cout << "Mass histogram: " << massNameList[i] << std::endl;
    std::cout << "Number of entries: " << h_mass_copy->GetEntries() << std::endl;
    std::cout << "Mean estimate: " << mass_mean << ", RMS: " << mass_rms << std::endl;
    std::cout << "Fit range: [" << fit_min << ", " << fit_max << "] MeV/c^{2}" << std::endl;

    // Breit-Wigner fit
    TF1 *bw = new TF1(Form("bw_%s", massNameList[i].Data()), breitwigner, fit_min, fit_max, 3);
    bw->SetParameters(mass_peak * 4.0, mass_peak_pos, 4.0);
    bw->SetParLimits(1, 780, 786);
    bw->SetParLimits(2, 0.5, 10.0);
    bw->SetLineColor(massColor[i]);
    bw->SetLineWidth(2);

    h_mass_copy->Fit(bw, "RQS");
    massResults[i].name = massNameList[i];
    massResults[i].mean = bw->GetParameter(1);
    massResults[i].mean_err = bw->GetParError(1);
    massResults[i].sigma = bw->GetParameter(2);
    massResults[i].sigma_err = bw->GetParError(2);
    massResults[i].chi2_ndf = bw->GetChisquare() / bw->GetNDF();
    massResults[i].entries = h_mass_copy->GetEntries();

    // ---- Draw ----
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
    leg_mass->AddEntry(bw, Form("BW: M = %.2f, #Gamma/2 = %.2f", massResults[i].mean, massResults[i].sigma), "l");
    leg_mass->Draw();

    c_mass->Update();
    c_mass->SaveAs(out_dir + "/mass_" + massNameList[i] + ".pdf");
    delete c_mass;
  }

  // ---- Summary ----
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

  // ---- Write residual bias ----
  std::ofstream myfile;
  TString myfile_nm = "../header/massbias_" + tuning_type + ".h";
  myfile.open(myfile_nm.Data());
  myfile << "const double energy_shift = " << mass_bias << ";\n";
  myfile.close();

  tree_file->Close();
  delete tree_file;

  return 0;
}
