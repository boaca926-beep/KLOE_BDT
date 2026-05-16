#include "../header_method/method.h"

TRandom *rnd = new TRandom(0); // seed 0 = use system time

int compr_bdt(){

  //gROOT->SetBatch(kTRUE);  
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);

  // switch on histogram errors
  TH1::SetDefaultSumw2();

  // random generator
  rnd=new TRandom3();

  // Inspect the input tree
  //cout << infile_nm << endl;
  const TString gen_file_nm = "/home/kloe/Desktop/input_bdt_TDATA_chain/gen/gen.root";
  const TString tree_file_nm = "/home/kloe/Desktop/input_bdt_TDATA_chain/cut/tree_pre_bdt.root";

  TFile* gen_file = new TFile(gen_file_nm);
  TFile* tree_file = new TFile(tree_file_nm);
  
  if (!gen_file || gen_file->IsZombie()) {
    cerr << "ERROR: Cannot open " << gen_file << endl;
    return 1;
  }

  if (!tree_file || tree_file->IsZombie()) {
    cerr << "ERROR: Cannot open " << tree_file << endl;
    return 1;
  }

  getObj(gen_file);
  getObj(tree_file);

  return 0;

}
