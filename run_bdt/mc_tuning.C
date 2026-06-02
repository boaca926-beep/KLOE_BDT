// MC tuning on final states particles gamma1, ,2 and 3. Energy, position and cluster time
#include "../header_bdt/plot_resol.h"

void mc_tuning() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetErrorX(0.8);
  TH1::SetDefaultSumw2();

  // Create output directory if it doesn't exist
  gSystem->Exec("mkdir -p ../plots_tuning");

  // Open tree file
  TString treeFile = "/home/kloe/Desktop/input_bdt_TDATA_norm_tuning/cut/tree_pre_bdt.root";
  TFile *ftree = TFile::Open(treeFile);
  if (!ftree || ftree->IsZombie()) {
    std::cerr << "ERROR: cannot open " << treeFile << std::endl;
    return;
  }
  
}
