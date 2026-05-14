#include "../header/sm_para.h"
#include "../header/path.h"
#include "../header/cut_para.h"
#include "../header/method.h"
#include <TStopwatch.h>


int tree_cut(){

  TStopwatch timer;
  timer.Start();

  cout << "input path: " << sampleFile << endl;

  TFile *f_input = new TFile(sampleFile + ".root");
  TTree *ALLCHAIN_CUT = (TTree*)f_input -> Get("ALLCHAIN_CUT");
  
  // define variables
  double lagvalue_min_7C = 0.;
  double deltaE = 0.;
  double angle_pi0gam12 = 0.;  
  double betapi0 = 0.;
  double m02 = 0., mplus2 = 0.;
  double m3pi = 0.;
  double ppIM = 0.;
  double IM3pi_7C = 0., IM3pi_true = 0.;
  double IM_pi0_7C = 0.;
  double Eisr = 0., Epi0_pho1 = 0., Epi0_pho2 = 0.;
  // pulls
  double pull_E1 = 0.;
  double pull_x1 = 0.;
  double pull_y1 = 0.;
  double pull_z1 = 0.;
  double pull_t1 = 0.;
  // 4-momentum
  double ppl_E = 0., ppl_px = 0., ppl_py = 0., ppl_pz = 0.;
  double pmi_E = 0., pmi_px = 0., pmi_py = 0., pmi_pz = 0.;
  double pho_E1 = 0., pho_px1 = 0., pho_py1 = 0., pho_pz1 = 0.;
  double pho_E2 = 0., pho_px2 = 0., pho_py2 = 0., pho_pz2 = 0.;
  double pho_E3 = 0., pho_px3 = 0., pho_py3 = 0., pho_pz3 = 0.;
  // true 4-momentum
  double ppl_E_true = 0., ppl_px_true = 0., ppl_py_true = 0., ppl_pz_true = 0.;
  double pmi_E_true = 0., pmi_px_true = 0., pmi_py_true = 0., pmi_pz_true = 0.;
  double pho_E1_true = 0., pho_px1_true = 0., pho_py1_true = 0., pho_pz1_true = 0.;
  double pho_E2_true = 0., pho_px2_true = 0., pho_py2_true = 0., pho_pz2_true = 0.;
  double pho_E3_true = 0., pho_px3_true = 0., pho_py3_true = 0., pho_pz3_true = 0.;
    
  int phid = 0, sig_type = 0;
  int bkg_indx = 0, recon_indx = 0;
  
  double evnt_tot = 0;
  double Eprompt_max = 0.;
  
  TFile *f_output = new TFile(outputCut + "tree_pre.root", "update");

  const int list_size = 11;
  const TString TNM[list_size] = {"TDATA", "TOMEGAPI", "TKPM", "TKSL", "T3PIGAM", "TRHOPI", "TETAGAM", "TBKGREST", "TUFO", "TEEG", "TISR3PI_SIG"};
  
  TTree *TTList[list_size];
  TCollection* tree_list = new TList;

  for (int i = 0; i < list_size; i ++) {
    TTList[i] = new TTree(TNM[i], "recreate");
    TTList[i] -> SetAutoSave(0);
    tree_list -> Add(TTList[i]);
  }
  
  // define branches for all trees (same for all)
  TObject* treeout=0;
  TIter treeliter(tree_list);
  while((treeout=treeliter.Next()) != 0) {
    TTree* tree_tmp=dynamic_cast<TTree*>(treeout);
    tree_tmp -> Branch("Br_ppl_E", &ppl_E, "Br_ppl_E/D");
    tree_tmp -> Branch("Br_ppl_px", &ppl_px, "Br_ppl_px/D");
    tree_tmp -> Branch("Br_ppl_py", &ppl_py, "Br_ppl_py/D");
    tree_tmp -> Branch("Br_ppl_pz", &ppl_pz, "Br_ppl_pz/D");
    tree_tmp -> Branch("Br_pmi_E", &pmi_E, "Br_pmi_E/D");
    tree_tmp -> Branch("Br_pmi_px", &pmi_px, "Br_pmi_px/D");
    tree_tmp -> Branch("Br_pmi_py", &pmi_py, "Br_pmi_py/D");
    tree_tmp -> Branch("Br_pmi_pz", &pmi_pz, "Br_pmi_pz/D");

    tree_tmp -> Branch("Br_ppl_E_true", &ppl_E_true, "Br_ppl_E_true/D");
    tree_tmp -> Branch("Br_ppl_px_true", &ppl_px_true, "Br_ppl_px_true/D");
    tree_tmp -> Branch("Br_ppl_py_true", &ppl_py_true, "Br_ppl_py_true/D");
    tree_tmp -> Branch("Br_ppl_pz_true", &ppl_pz_true, "Br_ppl_pz_true/D");
    tree_tmp -> Branch("Br_pmi_E_true", &pmi_E_true, "Br_pmi_E_true/D");
    tree_tmp -> Branch("Br_pmi_px_true", &pmi_px_true, "Br_pmi_px_true/D");
    tree_tmp -> Branch("Br_pmi_py_true", &pmi_py_true, "Br_pmi_py_true/D");
    tree_tmp -> Branch("Br_pmi_pz_true", &pmi_pz_true, "Br_pmi_pz_true/D");
    
    tree_tmp -> Branch("Br_E1", &pho_E1, "Br_pho_E1/D");
    tree_tmp -> Branch("Br_px1", &pho_px1, "Br_pho_px1/D");
    tree_tmp -> Branch("Br_py1", &pho_py1, "Br_pho_py1/D");
    tree_tmp -> Branch("Br_pz1", &pho_pz1, "Br_pho_pz1/D");
    tree_tmp -> Branch("Br_E2", &pho_E2, "Br_pho_E2/D");
    tree_tmp -> Branch("Br_px2", &pho_px2, "Br_pho_px2/D");
    tree_tmp -> Branch("Br_py2", &pho_py2, "Br_pho_py2/D");
    tree_tmp -> Branch("Br_pz2", &pho_pz2, "Br_pho_pz2/D");
    tree_tmp -> Branch("Br_E3", &pho_E3, "Br_pho_E3/D");
    tree_tmp -> Branch("Br_E1_true", &pho_E1_true, "Br_pho_E1_true/D");
    tree_tmp -> Branch("Br_px1_true", &pho_px1_true, "Br_pho_px1_true/D");
    tree_tmp -> Branch("Br_py1_true", &pho_py1_true, "Br_pho_py1_true/D");
    tree_tmp -> Branch("Br_pz1_true", &pho_pz1_true, "Br_pho_pz1_true/D");
    tree_tmp -> Branch("Br_E2_true", &pho_E2_true, "Br_pho_E2_true/D");
    tree_tmp -> Branch("Br_px2_true", &pho_px2_true, "Br_pho_px2_true/D");
    tree_tmp -> Branch("Br_py2_true", &pho_py2_true, "Br_pho_py2_true/D");
    tree_tmp -> Branch("Br_pz2_true", &pho_pz2_true, "Br_pho_pz2_true/D");
    tree_tmp -> Branch("Br_E3_true", &pho_E3_true, "Br_pho_E3_true/D");
    tree_tmp -> Branch("Br_px3_true", &pho_px3_true, "Br_pho_px3_true/D");
    tree_tmp -> Branch("Br_py3_true", &pho_py3_true, "Br_pho_py3_true/D");
    tree_tmp -> Branch("Br_pz3_true", &pho_pz3_true, "Br_pho_pz3_true/D");
    tree_tmp -> Branch("Br_px3", &pho_px3, "Br_pho_px3/D");
    tree_tmp -> Branch("Br_py3", &pho_py3, "Br_pho_py3/D");
    tree_tmp -> Branch("Br_pz3", &pho_pz3, "Br_pho_pz3/D");
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
    tree_tmp -> Branch("Br_m3pi", &m3pi, "Br_m3pi/D");
    
  }

  TLorentzVector pi0gam1, pi0gam2, isrgam, trkplus, trkmin;
  
  // Event loop
  for (Int_t irow = 0; irow < ALLCHAIN_CUT -> GetEntries(); irow ++) {
    ALLCHAIN_CUT -> GetEntry(irow);

    ppl_E = ALLCHAIN_CUT -> GetLeaf("Br_ppl_E") -> GetValue(0);
    ppl_px = ALLCHAIN_CUT -> GetLeaf("Br_ppl_px") -> GetValue(0);
    ppl_py = ALLCHAIN_CUT -> GetLeaf("Br_ppl_py") -> GetValue(0);
    ppl_pz = ALLCHAIN_CUT -> GetLeaf("Br_ppl_pz") -> GetValue(0);

    pmi_E = ALLCHAIN_CUT -> GetLeaf("Br_pmi_E") -> GetValue(0);
    pmi_px = ALLCHAIN_CUT -> GetLeaf("Br_pmi_px") -> GetValue(0);
    pmi_py = ALLCHAIN_CUT -> GetLeaf("Br_pmi_py") -> GetValue(0);
    pmi_pz = ALLCHAIN_CUT -> GetLeaf("Br_pmi_pz") -> GetValue(0);

    ppl_E_true = ALLCHAIN_CUT -> GetLeaf("Br_ppl_E_true") -> GetValue(0);
    ppl_px_true = ALLCHAIN_CUT -> GetLeaf("Br_ppl_px_true") -> GetValue(0);
    ppl_py_true = ALLCHAIN_CUT -> GetLeaf("Br_ppl_py_true") -> GetValue(0);
    ppl_pz_true = ALLCHAIN_CUT -> GetLeaf("Br_ppl_pz_true") -> GetValue(0);

    pmi_E_true = ALLCHAIN_CUT -> GetLeaf("Br_pmi_E_true") -> GetValue(0);
    pmi_px_true = ALLCHAIN_CUT -> GetLeaf("Br_pmi_px_true") -> GetValue(0);
    pmi_py_true = ALLCHAIN_CUT -> GetLeaf("Br_pmi_py_true") -> GetValue(0);
    pmi_pz_true = ALLCHAIN_CUT -> GetLeaf("Br_pmi_pz_true") -> GetValue(0);

    pho_E1 = ALLCHAIN_CUT -> GetLeaf("Br_E1") -> GetValue(0);
    pho_px1 = ALLCHAIN_CUT -> GetLeaf("Br_px1") -> GetValue(0);
    pho_py1 = ALLCHAIN_CUT -> GetLeaf("Br_py1") -> GetValue(0);
    pho_pz1 = ALLCHAIN_CUT -> GetLeaf("Br_pz1") -> GetValue(0);

    pho_E2 = ALLCHAIN_CUT->GetLeaf("Br_E2")->GetValue(0);
    pho_px2 = ALLCHAIN_CUT->GetLeaf("Br_px2")->GetValue(0);
    pho_py2 = ALLCHAIN_CUT->GetLeaf("Br_py2")->GetValue(0);
    pho_pz2 = ALLCHAIN_CUT->GetLeaf("Br_pz2")->GetValue(0);

    pho_E3 = ALLCHAIN_CUT->GetLeaf("Br_E3")->GetValue(0);
    pho_px3 = ALLCHAIN_CUT->GetLeaf("Br_px3")->GetValue(0);
    pho_py3 = ALLCHAIN_CUT->GetLeaf("Br_py3")->GetValue(0);
    pho_pz3 = ALLCHAIN_CUT->GetLeaf("Br_pz3")->GetValue(0);

    pho_E1_true = ALLCHAIN_CUT -> GetLeaf("Br_E1_true") -> GetValue(0);
    pho_px1_true = ALLCHAIN_CUT -> GetLeaf("Br_px1_true") -> GetValue(0);
    pho_py1_true = ALLCHAIN_CUT -> GetLeaf("Br_py1_true") -> GetValue(0);
    pho_pz1_true = ALLCHAIN_CUT -> GetLeaf("Br_pz1_true") -> GetValue(0);

    pho_E2_true = ALLCHAIN_CUT -> GetLeaf("Br_E2_true") -> GetValue(0);
    pho_px2_true = ALLCHAIN_CUT -> GetLeaf("Br_px2_true") -> GetValue(0);
    pho_py2_true = ALLCHAIN_CUT -> GetLeaf("Br_py2_true") -> GetValue(0);
    pho_pz2_true = ALLCHAIN_CUT -> GetLeaf("Br_pz2_true") -> GetValue(0);

    pho_E3_true = ALLCHAIN_CUT -> GetLeaf("Br_E3_true") -> GetValue(0);
    pho_px3_true = ALLCHAIN_CUT -> GetLeaf("Br_px3_true") -> GetValue(0);
    pho_py3_true = ALLCHAIN_CUT -> GetLeaf("Br_py3_true") -> GetValue(0);
    pho_pz3_true = ALLCHAIN_CUT -> GetLeaf("Br_pz3_true") -> GetValue(0);
    
    pull_E1 = ALLCHAIN_CUT -> GetLeaf("Br_PULLIST") -> GetValue(0);
    pull_x1 = ALLCHAIN_CUT -> GetLeaf("Br_PULLIST") -> GetValue(1);
    pull_y1 = ALLCHAIN_CUT -> GetLeaf("Br_PULLIST") -> GetValue(2);
    pull_z1 = ALLCHAIN_CUT -> GetLeaf("Br_PULLIST") -> GetValue(3);
    pull_t1 = ALLCHAIN_CUT -> GetLeaf("Br_PULLIST") -> GetValue(4);

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

    pi0gam1.SetPxPyPzE(pho_px1, pho_py1, pho_pz1, pho_E1);
    pi0gam2.SetPxPyPzE(pho_px2, pho_py2, pho_pz2, pho_E2);
    isrgam.SetPxPyPzE(pho_px3, pho_py3, pho_pz3, pho_E3);
    trkplus.SetPxPyPzE(ppl_px, ppl_py, ppl_pz, ppl_E);
    trkmin.SetPxPyPzE(pmi_px, pmi_py, pmi_pz, pmi_E);
  
    m3pi = (pi0gam1 + pi0gam2 + trkplus + trkmin).M();

    evnt_tot ++;

    Eprompt_max = 0.;
    if (Eisr > Eprompt_max) Eprompt_max = Eisr;
    if (Epi0_pho1 > Eprompt_max) Eprompt_max = Epi0_pho1;
    if (Epi0_pho2 > Eprompt_max) Eprompt_max = Epi0_pho2;

    // Selection cuts
    if (lagvalue_min_7C > chi2_cut) continue;
    else if (deltaE > deltaE_cut) continue;
    else if (angle_pi0gam12 > angle_cut) continue;
    else if (betapi0 > GetFBeta(beta_cut, c0, c1, ppIM)) continue;

    // Fill only the tree(s) corresponding to the current sample type
    if (data_type == "exp") {
      TTList[0]->Fill();   // data tree
    } else if (data_type == "ufo") {
      TTList[8]->Fill();   // UFO tree
    } else if (data_type == "eeg") {
      TTList[9]->Fill();   // EEG tree
    } else if (data_type == "sig") {
      TTList[10]->Fill();  // signal tree
    } else if (data_type == "ksl") {
      // For KSL (full MC), classify events into the appropriate background/signal trees
      if (phid == 0) {
        TTList[1]->Fill();   // omegapi
      } else if (phid == 1) {
        TTList[2]->Fill();   // kpm
      } else if (phid == 2) {
        TTList[3]->Fill();   // ksl
      } else if (phid == 3) {
        if (sig_type == 1) TTList[4]->Fill();  // 3pi
        else TTList[5]->Fill();                // rho pi
      } else if (phid == 5) {
        if (sig_type == 1) TTList[6]->Fill();  // etagam 3pi
        else TTList[7]->Fill();                // etagam rest
      } else {
        TTList[7]->Fill();   // bkg rest
      }
    }
  }

  // Write only the trees that were filled
  if (data_type == "exp") {
    TTList[0]->Write();
    cout << "TDATA saved" << endl;
  } else if (data_type == "ufo") {
    TTList[8]->Write();
    cout << "TUFO saved" << endl;
  } else if (data_type == "eeg") {
    TTList[9]->Write("TEEG");
    cout << "TEEG saved" << endl;
  } else if (data_type == "sig") {
    TTList[10]->Write("TISR3PI_SIG");
    cout << "TISR3PI_SIG saved" << endl;
  } else if (data_type == "ksl") {
    TTList[1]->Write();
    TTList[2]->Write();
    TTList[3]->Write();
    TTList[4]->Write();
    TTList[5]->Write();
    TTList[6]->Write();
    TTList[7]->Write();
    cout << "All trees in KSL saved" << endl;
  }

  cout << "=========================================\n"
       << f_output->GetName() << endl;
  cout << "evnt_tot = " << evnt_tot << "\n";

  f_output->Close();
  f_input->Close();

  double realTime = timer.RealTime();
  cout << "Expected time: " << realTime/60 << " mins" << endl;
  timer.Stop();
  timer.Print();
  
  return 0;
}
