// MC tuning on final states particles gamma1, ,2 and 3. Energy, position and cluster time
// mc_tuning.C
// Fit pull distributions for data and signal MC to extract mean and sigma
// mc_tuning.C
// Fit pull distributions for data and signal MC to extract mean and sigma
#include "../header_bdt/plot_resol.h"
#include <TF1.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>

void mc_tuning() {
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  // Create output directory for plots (if needed)
  gSystem->Exec("mkdir -p ../plots_tuning");

  // Open output ROOT file in RECREATE mode (overwrites previous file)
  TFile *fout = new TFile("../plots_tuning/pull_fits.root", "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: cannot create output file." << std::endl;
    return;
  }

  // Open input tree file
  TString treeFile = "/home/kloe/Desktop/input_bdt_TDATA_norm_tuning/cut/tree_pre_bdt.root";
  TFile *ftree = TFile::Open(treeFile);
  if (!ftree || ftree->IsZombie()) {
    std::cerr << "ERROR: cannot open " << treeFile << std::endl;
    return;
  }

  // Types: data and signal MC
  const int nb_mc_type = 2;
  const TString MC_TYPE[nb_mc_type] = {"TDATA", "TISR3PI_SIG"};
  const int nb_pull_type = 3;
  const TString pull_names[3] = {"E1", "E2", "E3"};
  //const TString pull_names[nb_pull_type] = {"E1"};

  const int bin_size = 100;
  const double XMIN = -5.0;
  const double XMAX = 5.0;
  const double fit_min = -1.5;
  const double fit_max = 1.5;   // fit core only

  TCanvas *c = new TCanvas("c", "Pull distributions", 1200, 800);
  c->Divide(3, 2);

  // Loop over data/MC types and pull variables
  for (int i = 0; i < nb_mc_type; i++) {
    TString mc_type = MC_TYPE[i];
    TTree *ttree = (TTree*) ftree->Get(mc_type);
    if (!ttree) {
      std::cerr << "Tree " << mc_type << " not found, skipping." << std::endl;
      continue;
    }

    for (int p = 0; p < nb_pull_type; p++) {
      TString pull_name = pull_names[p];
      TH1D *h_pull = new TH1D(Form("h_pull_%s_%s", mc_type.Data(), pull_name.Data()),
                              "", bin_size, XMIN, XMAX);
      h_pull->SetMarkerStyle(20);
      h_pull->SetMarkerSize(0.6);
      h_pull->SetLineColor(kBlack);
      h_pull->Sumw2();

      double pull_val = 0.0;
      TString branch_name = "Br_pull_" + pull_name;
      ttree->SetBranchAddress(branch_name, &pull_val);

      Long64_t nentries = ttree->GetEntries();
      for (Long64_t j = 0; j < nentries; j++) {
        ttree->GetEntry(j);
        h_pull->Fill(pull_val);
      }

      // Fit only the core 
      TF1 *gaus = new TF1("gaus", "gaus", fit_min, fit_max);
      double mean_est = h_pull->GetMean();
      double rms_est = h_pull->GetRMS();
      double amp_est = h_pull->GetBinContent(h_pull->FindBin(mean_est));
      gaus->SetParameters(amp_est, mean_est, rms_est * 0.7);
      gaus->SetParLimits(2, 0.2, 5.0);
      h_pull->Fit(gaus, "RQS");
      
      double mean = gaus->GetParameter(1);
      double sigma = gaus->GetParameter(2);
      double mean_err = gaus->GetParError(1);
      double sigma_err = gaus->GetParError(2);
      
      std::cout << "----------------------------------------" << std::endl;
      std::cout << "Sample: " << mc_type << ", Pull: " << pull_name << std::endl;
      std::cout << "Mean = " << mean << " +/- " << mean_err << std::endl;
      std::cout << "Sigma = " << sigma << " +/- " << sigma_err << std::endl;
      
      // Draw
      c->cd(p + 1 + (i*3));
      h_pull->Draw("E");
      gaus->SetLineColor(kRed);
      gaus->Draw("same");
      TLegend *leg = new TLegend(0.65, 0.7, 0.95, 0.88);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(h_pull, Form("%s, pull %s", mc_type.Data(), pull_name.Data()), "lep");
      leg->AddEntry(gaus, Form("Gaussian: #mu = %.3f, #sigma = %.3f", mean, sigma), "l");
      leg->Draw();
      
      // Write histogram to the output ROOT file (now in RECREATE mode)
      fout->cd();
      h_pull->Write();
    }
  }
  
  c->SaveAs("../plots_tuning/pull_fits.pdf");
  delete c;
  fout->Close();
  ftree->Close();
  std::cout << "Done. Results saved to ../plots_tuning/pull_fits.root and .pdf" << std::endl;
}
