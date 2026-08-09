// Compare weight of signal MC 3pi invariant mass correction as a function M3pi from plots_m3pi_corr/*.root

#include "../header_method/method.h"

struct Results {
  TString name;
  double max_glob, x_glob;
};

int weight_compr(const TString input_folder = "/home/bo/Desktop/KLOE_BDT/plots_m3pi_corr/") {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  const TString out_dir = "../weight_compr/";

  // ---- Ensure output directory exists ----
  gSystem->mkdir(out_dir, kTRUE);

  // Get trees
  const int nb_files = 3;
  
  TH1D *hHistList[nb_files];
  TString fileNameList[nb_files] = {"raw", "track", "pull"};
  int ColorList[nb_files] = {kBlack, kRed, kGreen};
  Results WeightResults[nb_files];

  for (int i = 0; i < nb_files; i++) {

    TString file_nm = input_folder + "omega_fit_results_" + fileNameList[i] + ".root"; 
    cout << file_nm << endl;

    TFile* tree_file = new TFile(file_nm);
    if (!tree_file || tree_file->IsZombie()) {
      cerr << "ERROR: Cannot open " << file_nm << endl;
      return 1;
    }

    checkFile(tree_file); // optional, shows contents

    TH1D *hist = (TH1D*)tree_file->Get("h_weight_smooth");
    hist->SetName("h_weight_" + fileNameList[i]);
    hist->SetLineColor(ColorList[i]);
    hist->SetLineWidth(2);

    //cout << hist->GetName() << endl;
    hHistList[i] = hist;

    double max_val = hist->GetMaximum();
    int bin_max = hist->GetMaximumBin();
    double x_at_max = hist->GetBinCenter(bin_max);
    
    double min_val = hist->GetMinimum();
    int bin_min = hist->GetMinimumBin();
    double x_at_min = hist->GetBinCenter(bin_min);

    double max_glob = std::max({TMath::Abs(max_val), TMath::Abs(min_val)});
    double x_glob = (TMath::Abs(max_val) >= TMath::Abs(min_val)) ? x_at_max : x_at_min;
    
    cout << fileNameList[i] << ", max = " << max_val << " at " << x_at_max << ", min = " << min_val << " at " << x_at_min << ", global max = " << max_glob << " at " << x_glob << endl;

    WeightResults[i].name = fileNameList[i];
    WeightResults[i].max_glob = max_glob;
    WeightResults[i].x_glob = x_glob;
    
  }

  // Draw histos
  TH1D *h_weight_raw = hHistList[0];
  h_weight_raw->SetLineStyle(3);

  TH1D *h_weight_track = hHistList[1];
  h_weight_track->SetLineStyle(2);

  TH1D *h_weight_pull = hHistList[2];
  h_weight_pull->SetLineStyle(5);

  TAxis* xAxis = h_weight_raw->GetXaxis();
  double xMin = xAxis->GetXmin();
  double xMax = xAxis->GetXmax();
 
  TLine *line_weight = new TLine(xMin, 1, xMax, 1.);
  line_weight->SetLineColor(kGray + 2);
  line_weight->SetLineStyle(2);
  line_weight->SetLineWidth(2);
  
  TCanvas *c = new TCanvas("c_weight", "Signal MC 3#pi Mass Weight Distribution", 1200, 700);

  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);
  
  h_weight_raw->SetMarkerStyle(20);
  h_weight_raw->SetMarkerSize(0.6);
  h_weight_raw->GetXaxis()->SetNdivisions(505);
  h_weight_raw->GetYaxis()->SetTitle("Weight");
  h_weight_raw->GetYaxis()->SetRangeUser(0.01, h_weight_raw->GetMaximum() * 1.6);
  h_weight_raw->GetYaxis()->CenterTitle();
  h_weight_raw->GetYaxis()->SetTitleSize(0.05);
  h_weight_raw->GetYaxis()->SetTitleOffset(1.4);
  h_weight_raw->GetYaxis()->SetLabelSize(0.04);
  h_weight_raw->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
  h_weight_raw->GetXaxis()->SetTitleSize(0.05);
  h_weight_raw->GetXaxis()->SetTitleOffset(1.2);
  h_weight_raw->GetXaxis()->SetLabelSize(0.04);
  h_weight_raw->GetXaxis()->CenterTitle();
  
  h_weight_raw->Draw();
  h_weight_track->Draw("same hist");
  h_weight_pull->Draw("same hist");
  line_weight->Draw();

  TLegend *leg = new TLegend(0.6, 0.7, 0.85, 0.9);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  
  leg->AddEntry(h_weight_raw, "Kinematic fit", "l");
  leg->AddEntry(h_weight_track, "Track correction", "l");
  leg->AddEntry(h_weight_pull, "Pull correction", "l");
  leg->Draw();

  // Find overall maximum weight and its x-position
  double overall_max = -1e9;
  double overall_x = 0;
  for (int i = 0; i < nb_files; i++) {
    if (WeightResults[i].max_glob > overall_max) {
      overall_max = WeightResults[i].max_glob;
      overall_x = WeightResults[i].x_glob;
    }
  }

  double reduce_weight = TMath::Abs(WeightResults[0].max_glob - WeightResults[2].max_glob) / WeightResults[0].max_glob * 100.;
  
  cout << "raw max weight = " << WeightResults[0].max_glob << "\n"
       << "pull max weight = " << WeightResults[2].max_glob << "\n";
  
  TPaveText *pt = new TPaveText(0.2, 0.8, 0.35, 0.89, "NDC");
  pt->SetFillColor(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.04);
  pt->SetTextFont(42);
  //pt->AddText(Form("Maximum weight = %.3f at #pm %.3f", overall_max, overall_x));
  pt->AddText(Form("Maximum weight reduced by  %.1f%s", reduce_weight, "%"));
  
  pt->Draw();

  c->SaveAs(out_dir + "weight_compr.pdf");

  cout << "plots are save at " << out_dir << "weight_compr.pdf" << endl;
  return 0;
}
