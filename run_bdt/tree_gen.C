#include "../header/path.h"

int tree_gen(){

  
  TFile *f_in = new TFile(sig_path + "sig.root");

  //checkFile(gen_in);

  TTree *ALLCHAIN_GEN = (TTree*)f_in -> Get("ALLCHAIN_GEN");

  double IM3pi_gen = 0.;

  TFile *f_output = new TFile(outputGen + "tree_gen.root", "recreate");

  TTree * TISR3PI_SIG_GEN = new TTree("TISR3PI_SIG_GEN", "recreate");
  TISR3PI_SIG_GEN -> SetAutoSave(0);

  TISR3PI_SIG_GEN -> Branch("Br_IM3pi_gen", &IM3pi_gen, "Br_IM3pi_gen/D");

  for (Int_t irow = 0; irow < ALLCHAIN_GEN -> GetEntries(); irow ++) {// loop trees

    //if (irow > 1e3) break;
    
    ALLCHAIN_GEN -> GetEntry(irow);

    IM3pi_gen = ALLCHAIN_GEN -> GetLeaf("Br_IM_3pi") -> GetValue(0);

    TISR3PI_SIG_GEN -> Fill();
	
    //cout << irow << endl;
    
  }

  TISR3PI_SIG_GEN -> Write();
	
  f_output -> Close();
  
  return 0;
  
}

