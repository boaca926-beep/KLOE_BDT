#include "../header_plot/plot.h"
#include "../header_method/method.h"
#include "../header_bdt/plot.h"

// plot_type can be "sfw2d" (default) or "ppIM_vs_beta0"
int getplothist(const TString plot_type = "") {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  cout << "Creating 2D histos for " << plot_type << "... " << hist_root << "\n" << endl;

  TFile* infile = new TFile(hist_root);
  if (!infile || infile->IsZombie()) {
    cerr << "ERROR: Cannot open " << hist_root << endl;
    return 1;
  }

  // ------------------------------------------------------------------
  // Select the appropriate TList based on plot_type
  // ------------------------------------------------------------------
  TList *sourceList = nullptr;
  TString listName;
  if (plot_type == "sfw2d") {
    listName = "HSFW2D";
    sourceList = (TList*)infile->Get(listName);
  } else if (plot_type == "ppIM_vs_betapi0") {
    listName = "HppIM_vs_betapi0";
    sourceList = (TList*)infile->Get(listName);
  } else {
    cerr << "ERROR: Unknown plot_type '" << plot_type << "'. Use 'sfw2d' or 'ppIM_vs_betapi0'." << endl;
    return 1;
  }

  if (!sourceList) {
    cerr << "ERROR: " << listName << " not found in " << hist_root << endl;
    infile->Close();
    return 1;
  }

  // Helper function to get histogram from list with error checking
  auto getHist = [&](const char* name) -> TH2D* {
    TH2D* h = (TH2D*)sourceList->FindObject(name);
    if (!h) cerr << "WARNING: histogram " << name << " not found in " << listName << endl;
    return h;
  };

  // Determine histogram name prefix based on plot_type
  TString prefix = (plot_type == "sfw2d") ? "sfw" : "ppIM_vs_betapi0";
    
  // Extract individual histograms
  TH2D *h_sig     = getHist(Form("h2d_%s_TISR3PI_SIG", prefix.Data()));
  TH2D *h_omega   = getHist(Form("h2d_%s_TOMEGAPI",   prefix.Data()));
  TH2D *h_kpm     = getHist(Form("h2d_%s_TKPM",       prefix.Data()));
  TH2D *h_ksl     = getHist(Form("h2d_%s_TKSL",       prefix.Data()));
  TH2D *h_rhopi   = getHist(Form("h2d_%s_TRHOPI",     prefix.Data()));
  TH2D *h_etagam  = getHist(Form("h2d_%s_TETAGAM",    prefix.Data()));
  TH2D *h_bkgrest = getHist(Form("h2d_%s_TBKGREST",   prefix.Data()));
  TH2D *h_data    = getHist(Form("h2d_%s_TDATA",      prefix.Data()));
  TH2D *h_eeg     = getHist(Form("h2d_%s_TEEG",       prefix.Data()));
  
  if (!h_data) {
    cerr << "ERROR: Data histogram not found!" << endl;
    infile->Close();
    return 1;
  }

  // ------------------------------------------------------------------
  // Build combined histograms (MC rest, background sums, MC sum)
  // ------------------------------------------------------------------
  // MC rest (background without signal)
  TH2D *hmcrest = (TH2D*)h_bkgrest->Clone();
  hmcrest->Add(h_kpm, 1.);
  hmcrest->Add(h_rhopi, 1.);
  
  // Background sum without η→γγ
  TH2D *hbkgsum_noeta = (TH2D*)h_eeg->Clone();
  hbkgsum_noeta->Add(h_ksl, 1.);
  hbkgsum_noeta->Add(h_omega, 1.);
  hbkgsum_noeta->Add(hmcrest, 1.);
  
  // Background sum (including η→γγ)
  TH2D *hbkgsum = (TH2D*)hbkgsum_noeta->Clone();
  hbkgsum->Add(h_etagam, 1.);
  
  // MC sum (background + signal)
  TH2D *hmcsum = (TH2D*)hbkgsum->Clone();
  hmcsum->Add(h_sig, 1.);
  
  hmcrest->SetName("hmcrest");
  hbkgsum_noeta->SetName("hbkgsum_noeta");
  hbkgsum->SetName("hbkgsum");
  hmcsum->SetName("hmcsum");
 
  // ------------------------------------------------------------------
  // Write to output file
  // ------------------------------------------------------------------
  TString outfile_name = output_path + (plot_type == "sfw2d" ? "sfw2d_output.root" : "ppIM_vs_betapi0_output.root");
  TFile *f_out = new TFile(outfile_name, "recreate");
  h_data->Write();
  h_sig->Write();
  h_etagam->Write();
  hmcrest->Write();
  hbkgsum_noeta->Write();
  hbkgsum->Write();
  hmcsum->Write();
  f_out->Close();

  infile->Close();
  cout << "2D histograms for " << plot_type << " saved to " << outfile_name << endl;
  return 0;
}
