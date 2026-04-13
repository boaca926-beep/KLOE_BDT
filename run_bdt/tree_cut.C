#include "../header_bdt/sm_para.h"
#include "../header_bdt/path.h"
#include "../header_bdt/cut_para.h"
#include "../header_bdt/method.h"
#include <TStopwatch.h>

TStopwatch timer;
timer.Start();
  
int tree_cut(){

  //const TString f_path = folder_in + "/" + data_type + ".root";
  cout << "input path: " << sampleFile << endl;

  TFile *f_input = new TFile(sampleFile + ".root");

  //checkFile(f_input);

  TTree *ALLCHAIN_CUT = (TTree*)f_input -> Get("ALLCHAIN_CUT");
  
  // define variables
  double lagvalue_min_7C = 0.;
  double deltaE = 0.;
  double angle_pi0gam12 = 0.;  
  double betapi0 = 0.;
  double ppIM = 0.;
  double IM3pi_7C = 0., IM3pi_true = 0.;
  double IM_pi0_7C = 0.;
  double Eisr = 0., Epi0_pho1 = 0., Epi0_pho2 = 0.;
  double pull_E1 = 0.;
  double pull_x1 = 0.;
  double pull_y1 = 0.;
  double pull_z1 = 0.;
  double pull_t1 = 0.;
  double m02 = 0., mplus2 = 0.;
  
  int phid = 0, sig_type = 0;
  //int trigger_indx = 0;
  //int filfo_indx = 0;
  //int evtcls_indx = 0;
  //int fstate_indx = 0;
  int bkg_indx = 0, recon_indx = 0;

  double evnt_tot = 0; // total number of events
  double evnt_trigger = 0; // number of events after trigger
  double evnt_filfo = 0; // number of events after fiflo
  double evnt_evtcls = 0; // number of events after event classification
  double evnt_fstate = 0; // number of events after 2 tracks and 3 prompt photons
  double Eprompt_max = 0.;
  
  TFile *f_output = new TFile(outputCut + "tree_pre.root", "update");

  // ksl stream
  const int list_size = 11;
  const TString TNM[list_size] = {"TDATA", "TOMEGAPI", "TKPM", "TKSL", "T3PIGAM", "TRHOPI", "TETAGAM", "TBKGREST", "TUFO", "TEEG", "TISR3PI_SIG"};
  
  TTree *TTList[list_size];
  
  TCollection* tree_list = new TList;

  for (int i = 0; i < list_size; i ++) {// loop over the Tree list

    TTList[i] = new TTree(TNM[i], "recreate");
    TTList[i] -> SetAutoSave(0);
    
    tree_list -> Add(TTList[i]);
    
  } // end the Tree list
  
  // define branches
  TObject* treeout=0;
  TIter treeliter(tree_list);
  
  while((treeout=treeliter.Next()) != 0) {
    
    //treeout->Print();
    TTree* tree_tmp=dynamic_cast<TTree*>(treeout);

    tree_tmp -> Branch("Br_pull_E1", &pull_E1, "Br_pull_E1/D");
    tree_tmp -> Branch("Br_pull_x1", &pull_x1, "Br_pull_x1/D");
    tree_tmp -> Branch("Br_pull_y1", &pull_y1, "Br_pull_y1/D");
    tree_tmp -> Branch("Br_pull_z1", &pull_z1, "Br_pull_z1/D");
    tree_tmp -> Branch("Br_pull_t1", &pull_t1, "Br_pull_t1/D");

    tree_tmp -> Branch("Br_sig_type", &sig_type, "Br_sig_type/I");
    tree_tmp -> Branch("Br_bkg_indx", &bkg_indx, "Br_bkg_indx/I");
    tree_tmp -> Branch("Br_recon_indx", &recon_indx, "Br_recon_indx/I");

    tree_tmp -> Branch("Br_IM3pi_7C", &IM3pi_7C, "Br_IM3pi_7C/D");
    tree_tmp -> Branch("Br_IM3pi_true", &IM3pi_true, "Br_IM3pi_true/D");
    tree_tmp -> Branch("Br_IM_pi0_7C", &IM_pi0_7C, "Br_IM_pi0_7C/D");
    tree_tmp -> Branch("Br_mplus2", &mplus2, "Br_mplus2/D");
    tree_tmp -> Branch("Br_m02", &m02, "Br_m02/D");
    tree_tmp -> Branch("Br_ppIM", &ppIM, "Br_ppIM/D");

    tree_tmp -> Branch("Br_Eisr", &Eisr, "Br_Eisr/D");
    tree_tmp -> Branch("Br_Epi0_pho1", &Epi0_pho1, "Br_Epi0_pho1/D");
    tree_tmp -> Branch("Br_Epi0_pho2", &Epi0_pho2, "Br_Epi0_pho2/D");
    tree_tmp -> Branch("Br_angle_pi0gam12", &angle_pi0gam12, "Br_angle_pi0gam12/D");
    tree_tmp -> Branch("Br_betapi0", &betapi0, "Br_betapi0/D");
    tree_tmp -> Branch("Br_Eprompt_max", &Eprompt_max, "Br_Eprompt_max/D");

    tree_tmp -> Branch("Br_lagvalue_min_7C", &lagvalue_min_7C, "Br_lagvalue_min_7C/D");
    tree_tmp -> Branch("Br_deltaE", &deltaE, "Br_deltaE/D");
    //tree_tmp -> Branch("Br_angle_pi0gam12", &angle_pi0gam12, "Br_angle_pi0gam12/D");
    //tree_tmp -> Branch("Br_betapi0", &betapi0, "Br_betapi0/D");
    
  }
  

  for (Int_t irow = 0; irow < ALLCHAIN_CUT -> GetEntries(); irow ++) {// loop trees
	  
    ALLCHAIN_CUT -> GetEntry(irow);

    pull_E1 = ALLCHAIN_CUT -> GetLeaf("Br_PULLIST") -> GetValue(0);
    pull_x1 = ALLCHAIN_CUT -> GetLeaf("Br_PULLIST") -> GetValue(1);
    pull_y1 = ALLCHAIN_CUT -> GetLeaf("Br_PULLIST") -> GetValue(2);
    pull_z1 = ALLCHAIN_CUT -> GetLeaf("Br_PULLIST") -> GetValue(3);
    pull_t1 = ALLCHAIN_CUT -> GetLeaf("Br_PULLIST") -> GetValue(4);

    //trigger_indx = ALLCHAIN_CUT -> GetLeaf("Br_trigger_indx") -> GetValue(0);
    //filfo_indx = ALLCHAIN_CUT -> GetLeaf("Br_filfo_indx") -> GetValue(0);
    //evtcls_indx = ALLCHAIN_CUT -> GetLeaf("Br_evtcls_indx") -> GetValue(0);
    //fstate_indx = ALLCHAIN_CUT -> GetLeaf("Br_fstate_indx") -> GetValue(0);

    bkg_indx = ALLCHAIN_CUT -> GetLeaf("Br_bkg_indx") -> GetValue(0);
    recon_indx = ALLCHAIN_CUT -> GetLeaf("Br_recon_indx") -> GetValue(0);
    
    phid = ALLCHAIN_CUT -> GetLeaf("Br_phid") -> GetValue(0);
    sig_type = ALLCHAIN_CUT -> GetLeaf("Br_sig_type") -> GetValue(0);
    lagvalue_min_7C = ALLCHAIN_CUT -> GetLeaf("Br_lagvalue_min_7C") -> GetValue(0);
    deltaE = ALLCHAIN_CUT -> GetLeaf("Br_ENERGYLIST") -> GetValue(2); 
    angle_pi0gam12 = ALLCHAIN_CUT -> GetLeaf("Br_ANGLELIST") -> GetValue(0);
    betapi0 = ALLCHAIN_CUT -> GetLeaf("Br_betapi0") -> GetValue(0);
    ppIM = ALLCHAIN_CUT -> GetLeaf("Br_MASSLIST") -> GetValue(5);
    m02 = ALLCHAIN_CUT -> GetLeaf("Br_MASSLIST") -> GetValue(10);
    mplus2 = ALLCHAIN_CUT -> GetLeaf("Br_MASSLIST") -> GetValue(11);
    
    IM3pi_7C = ALLCHAIN_CUT -> GetLeaf("Br_IM3pi_7C") -> GetValue(0);
    IM_pi0_7C = ALLCHAIN_CUT -> GetLeaf("Br_IM_pi0_7C") -> GetValue(0);
    IM3pi_true = ALLCHAIN_CUT -> GetLeaf("Br_IM3pi_true") -> GetValue(0);
    Eisr = ALLCHAIN_CUT -> GetLeaf("Br_ENERGYLIST") -> GetValue(0);
    Epi0_pho1 = ALLCHAIN_CUT -> GetLeaf("Br_ENERGYLIST") -> GetValue(1);
    Epi0_pho2 = ALLCHAIN_CUT -> GetLeaf("Br_ENERGYLIST") -> GetValue(3);
    
    evnt_tot ++;

    Eprompt_max = 0.;

    if (Eisr > Eprompt_max) Eprompt_max = Eisr;
    if (Epi0_pho1 > Eprompt_max) Eprompt_max = Epi0_pho1;
    if (Epi0_pho2 > Eprompt_max) Eprompt_max = Epi0_pho2;

    //if (evnt_tot > 1e5) break;

    //// trigger cut
    //if (trigger_indx == 1) continue;
    //evnt_trigger ++;
    //cout << trigger_indx << endl;

    //// filfo cut
    //if (filfo_indx == 1) continue;
    //evnt_filfo ++;
    //cout << filfo_indx << endl;

    //// event classification cut
    //if (evtcls_indx == 1) continue;
    //evnt_evtcls ++;
    //cout << evtcls_indx << endl;

    //// 2 tracks and 3 prompt photons cut
    //if (fstate_indx == 1) continue;
    //evnt_fstate ++;
    //cout << fstate_indx << endl;

    //// background rejection (nominal condtions)
    if (lagvalue_min_7C > chi2_cut) continue;
    else if (deltaE > deltaE_cut) continue;
    else if (angle_pi0gam12 > angle_cut) continue;
    else if (betapi0 > GetFBeta(beta_cut, c0, c1, ppIM)) continue;
    //else if (betapi0 < GetFBeta(beta_cut, c0, c1, ppIM)) continue; // radiative background region

    //// background rejection (preserve etagam bkg)
    //else if (Eprompt_max < 320.) continue;
    
    //else if (IM3pi_7C > 850. || IM3pi_7C < 650.) continue; // reduce etagam background
    
    //cout << pull_E1 << endl;
    
    /*
    cout << chi2_cut << "\n"
	 << deltaE_cut << "\n"
	 << angle_cut << "\n"
	 << betapi0 << endl;
    */

    //cout << bkg_indx << ", " << recon_indx << endl;
    
    //cout << deltaE << endl;

    TTList[0]-> Fill(); // data
    TTList[8]-> Fill(); // ufo
    TTList[9]-> Fill(); // eeg
    TTList[10]-> Fill(); // sig
    

    if (phid == 0) {// omegapi
      TTList[1] -> Fill();
    }
    else if (phid == 1) {// kpm
      TTList[2] -> Fill();
    }
    else if (phid == 2) {// ksl
      TTList[3] -> Fill();
    }
    else if (phid == 3) {
      if (sig_type == 1) {// 3pi
	TTList[4] -> Fill();
      }
      else {// rho pi
	TTList[5] -> Fill();
      }
    }
    else if (phid == 5) {
      if (sig_type == 1) {// etagam 3pi
	TTList[6] -> Fill();
      }
      else {// etagam rest
	TTList[7] -> Fill();
      }
    }
    else {// bkg rest
      TTList[7] -> Fill();
    }
  }

  // save
  if (data_type == "exp") {
    TTList[0]-> Write();
    cout << "TDATA saved" << endl;
  }
  else if (data_type == "ufo") {
    TTList[8]-> Write();
    cout << "TUFO saved" << endl;

    // Fill and save cut parameters
    //TTree* TPARA_CUT = new TTree("TPARA_CUT", "recreate");
    //TPARA_CUT -> Branch("Br_egammamin", &egammamin, "Br_egammamin/D");
    //TPARA_CUT -> Fill();
    //TPARA_CUT -> Write();
  }
  else if (data_type == "eeg") {
    TTList[9]-> Write("TEEG");
    cout << "TEEG saved" << endl;
  }
  else if (data_type == "sig") {
    TTList[10]-> Write("TISR3PI_SIG");
    cout << "TISR3PI_SIG saved" << endl;
  }
  else if (data_type == "ksl") {
    TTList[1] -> Write();
    TTList[2] -> Write();
    TTList[3] -> Write();
    TTList[4] -> Write();
    TTList[5] -> Write();
    TTList[6] -> Write();
    TTList[7] -> Write();
    cout << "All trees in KSL saved" << endl;
  }

    
  
  // Summary
  cout << "=========================================\n"
       << f_output -> GetName() << endl;
  
  cout << "evnt_tot = " << evnt_tot << "\n"
       << "2 tacks and 3 prompt photons: " << evnt_fstate << "\n";

  cout << "Cut parameters: \n"
       << "chi2_cut = " << chi2_cut << "\n"
       << "angle_cut = " << angle_cut << "\n"
       << "deltaE_cut = " << deltaE_cut << "\n"
       << "beta_cut = " << beta_cut << "\n";
       
  /// histos

  /// save
  f_output -> Close();

  double realTime = timer.RealTime();  // Wall clock time in seconds
  double expect_time = realTime / 60; // Expected time used in mins
  //cout << "Real time = " << realTime << endl;
  cout << "Expected time: " << expect_time << " mins" << endl;

  timer.Stop();
  timer.Print();
  
  return 0;
  
}

  
  
