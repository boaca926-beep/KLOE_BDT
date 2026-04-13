#include "../header_bdt/sm_para.h"
#include "../header_bdt/path_sample.h"
#include "../header_bdt/cut_para.h"
#include "../header_bdt/method.h"
#include <TStopwatch.h>

TStopwatch timer;
timer.Start();

TLorentzVector Get4vector(double E, double px, double py, double pz) {

  //given a cluster index returns the 4-mom of a photon
  TVector3 gamma(px, py, pz);
  
  Double_t scale1;
  scale1=E/gamma.Mag();
  TLorentzVector gamma4mom(scale1*gamma, E);
  //cout << gamma4mom.M() << endl;
  return gamma4mom;

}

int tree_sample(){

  //const TString f_path = folder_in + "/" + data_type + ".root";
  cout << "input path: " << sampleFile << endl;

  TFile *f_input = new TFile(sampleFile + ".root");

  //checkFile(f_input);

  TTree *ALLCHAIN_CUT = (TTree*)f_input -> Get("ALLCHAIN_CUT");
  
  // define variables
  double pho_E1 = 0., pho_px1 = 0., pho_py1 = 0., pho_pz1 = 0.;
  double pho_E2 = 0., pho_px2 = 0., pho_py2 = 0., pho_pz2 = 0.;
  double pho_E3 = 0., pho_px3 = 0., pho_py3 = 0., pho_pz3 = 0.;
  double ppl_E = 0., ppl_px = 0., ppl_py = 0., ppl_pz = 0.;
  double pmi_E = 0., pmi_px = 0., pmi_py = 0., pmi_pz = 0.;
  double mpi0 = 0., m3pi = 0.;
  
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
  
  int phid = 0, sig_type = -1;
  int trigger_indx = -1;
  int filfo_indx = -1;
  //int filfo28_indx = -1;
  int evtcls_indx = -1;
  int sel_indx = -1;

  //int fstate_indx = -1;
  int bkg_indx = 0, recon_indx = 0;

  double evnt_tot = 0; // total number of events
  double evnt_sel = 0.; // number of selected events
  
  double evnt_trigger = 0; // number of events after trigger
  double evnt_filfo = 0; // number of events after fiflo
  double evnt_evtcls = 0; // number of events after event classification
  
  double Eprompt_max = 0.;
  
  TFile *f_output = new TFile(outputCut, "update");

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

    tree_tmp -> Branch("Br_E1", &pho_E1, "Br_E1/D");
    tree_tmp -> Branch("Br_px1", &pho_px1, "Br_px1/D");
    tree_tmp -> Branch("Br_py1", &pho_py1, "Br_py1/D");
    tree_tmp -> Branch("Br_pz1", &pho_pz1, "Br_pz1/D");

    tree_tmp -> Branch("Br_E2", &pho_E2, "Br_E2/D");
    tree_tmp -> Branch("Br_px2", &pho_px2, "Br_px2/D");
    tree_tmp -> Branch("Br_py2", &pho_py2, "Br_py2/D");
    tree_tmp -> Branch("Br_pz2", &pho_pz2, "Br_pz2/D");

    tree_tmp -> Branch("Br_E3", &pho_E3, "Br_E3/D");
    tree_tmp -> Branch("Br_px3", &pho_px3, "Br_px3/D");
    tree_tmp -> Branch("Br_py3", &pho_py3, "Br_py3/D");
    tree_tmp -> Branch("Br_pz3", &pho_pz3, "Br_pz3/D");

    tree_tmp -> Branch("Br_ppl_E", &ppl_E, "Br_ppl_E/D");
    tree_tmp -> Branch("Br_ppl_px", &ppl_px, "Br_ppl_px/D");
    tree_tmp -> Branch("Br_ppl_py", &ppl_py, "Br_ppl_py/D");
    tree_tmp -> Branch("Br_ppl_pz", &ppl_pz, "Br_ppl_pz/D");

    tree_tmp -> Branch("Br_pmi_E", &pmi_E, "Br_pmi_E/D");
    tree_tmp -> Branch("Br_pmi_px", &pmi_px, "Br_pmi_px/D");
    tree_tmp -> Branch("Br_pmi_py", &pmi_py, "Br_pmi_py/D");
    tree_tmp -> Branch("Br_pmi_pz", &pmi_pz, "Br_pmi_pz/D");
    
    tree_tmp -> Branch("Br_pull_E1", &pull_E1, "Br_pull_E1/D");
    tree_tmp -> Branch("Br_pull_x1", &pull_x1, "Br_pull_x1/D");
    tree_tmp -> Branch("Br_pull_y1", &pull_y1, "Br_pull_y1/D");
    tree_tmp -> Branch("Br_pull_z1", &pull_z1, "Br_pull_z1/D");
    tree_tmp -> Branch("Br_pull_t1", &pull_t1, "Br_pull_t1/D");

    tree_tmp -> Branch("Br_sig_type", &sig_type, "Br_sig_type/I");
    tree_tmp -> Branch("Br_bkg_indx", &bkg_indx, "Br_bkg_indx/I");
    tree_tmp -> Branch("Br_recon_indx", &recon_indx, "Br_recon_indx/I");
    tree_tmp -> Branch("Br_sel_indx", &sel_indx, "Br_sel_indx/I");

    tree_tmp -> Branch("Br_trigger_indx", &trigger_indx, "Br_trigger_indx/I");
    tree_tmp -> Branch("Br_filfo_indx", &filfo_indx, "Br_filfo_indx/I");
    //tree_tmp -> Branch("Br_filfo28_indx", &filfo28_indx, "Br_filfo28_indx/I");
    tree_tmp -> Branch("Br_evtcls_indx", &evtcls_indx, "Br_evtcls_indx/I");
    
    tree_tmp -> Branch("Br_IM3pi_7C", &IM3pi_7C, "Br_IM3pi_7C/D");
    tree_tmp -> Branch("Br_IM3pi_true", &IM3pi_true, "Br_IM3pi_true/D");
    tree_tmp -> Branch("Br_IM_pi0_7C", &IM_pi0_7C, "Br_IM_pi0_7C/D");
    tree_tmp -> Branch("Br_mplus2", &mplus2, "Br_mplus2/D");
    tree_tmp -> Branch("Br_m02", &m02, "Br_m02/D");
    tree_tmp -> Branch("Br_ppIM", &ppIM, "Br_ppIM/D");

    tree_tmp -> Branch("Br_Eisr", &Eisr, "Br_Eisr/D");
    tree_tmp -> Branch("Br_Epi0_pho1", &Epi0_pho1, "Br_Epi0_pho1/D");
    tree_tmp -> Branch("Br_Epi0_pho2", &Epi0_pho2, "Br_Epi0_pho2/D");
    tree_tmp -> Branch("Br_mpi0", &mpi0, "Br_mpi0/D");
    tree_tmp -> Branch("Br_m3pi", &m3pi, "Br_m3pi/D");
    
    tree_tmp -> Branch("Br_angle_pi0gam12", &angle_pi0gam12, "Br_angle_pi0gam12/D");
    tree_tmp -> Branch("Br_betapi0", &betapi0, "Br_betapi0/D");

    tree_tmp -> Branch("Br_lagvalue_min_7C", &lagvalue_min_7C, "Br_lagvalue_min_7C/D");
    tree_tmp -> Branch("Br_deltaE", &deltaE, "Br_deltaE/D");
    tree_tmp -> Branch("Br_angle_pi0gam12", &angle_pi0gam12, "Br_angle_pi0gam12/D");
    tree_tmp -> Branch("Br_betapi0", &betapi0, "Br_betapi0/D");
    tree_tmp -> Branch("Br_Eprompt_max", &Eprompt_max, "Br_Eprompt_max/D");
	  
  }
  

  TLorentzVector pi0gam1, pi0gam2, isrgam, trkplus, trkmin;
  
  for (Int_t irow = 0; irow < ALLCHAIN_CUT -> GetEntries(); irow ++) {// loop trees
	  
    ALLCHAIN_CUT -> GetEntry(irow);

    pho_E1 = ALLCHAIN_CUT -> GetLeaf("Br_E1") -> GetValue(0);
    pho_px1 = ALLCHAIN_CUT -> GetLeaf("Br_px1") -> GetValue(0);
    pho_py1 = ALLCHAIN_CUT -> GetLeaf("Br_py1") -> GetValue(0);
    pho_pz1 = ALLCHAIN_CUT -> GetLeaf("Br_pz1") -> GetValue(0);

    pho_E2 = ALLCHAIN_CUT -> GetLeaf("Br_E2") -> GetValue(0);
    pho_px2 = ALLCHAIN_CUT -> GetLeaf("Br_px2") -> GetValue(0);
    pho_py2 = ALLCHAIN_CUT -> GetLeaf("Br_py2") -> GetValue(0);
    pho_pz2 = ALLCHAIN_CUT -> GetLeaf("Br_pz2") -> GetValue(0);

    pho_E3 = ALLCHAIN_CUT -> GetLeaf("Br_E3") -> GetValue(0);
    pho_px3 = ALLCHAIN_CUT -> GetLeaf("Br_px3") -> GetValue(0);
    pho_py3 = ALLCHAIN_CUT -> GetLeaf("Br_py3") -> GetValue(0);
    pho_pz3 = ALLCHAIN_CUT -> GetLeaf("Br_pz3") -> GetValue(0);

    ppl_E = ALLCHAIN_CUT -> GetLeaf("Br_ppl_E") -> GetValue(0);
    ppl_px = ALLCHAIN_CUT -> GetLeaf("Br_ppl_px") -> GetValue(0);
    ppl_py = ALLCHAIN_CUT -> GetLeaf("Br_ppl_py") -> GetValue(0);
    ppl_pz = ALLCHAIN_CUT -> GetLeaf("Br_ppl_pz") -> GetValue(0);

    pmi_E = ALLCHAIN_CUT -> GetLeaf("Br_pmi_E") -> GetValue(0);
    pmi_px = ALLCHAIN_CUT -> GetLeaf("Br_pmi_px") -> GetValue(0);
    pmi_py = ALLCHAIN_CUT -> GetLeaf("Br_pmi_py") -> GetValue(0);
    pmi_pz = ALLCHAIN_CUT -> GetLeaf("Br_pmi_pz") -> GetValue(0);

    pull_E1 = ALLCHAIN_CUT -> GetLeaf("Br_PULLIST") -> GetValue(0);
    pull_x1 = ALLCHAIN_CUT -> GetLeaf("Br_PULLIST") -> GetValue(1);
    pull_y1 = ALLCHAIN_CUT -> GetLeaf("Br_PULLIST") -> GetValue(2);
    pull_z1 = ALLCHAIN_CUT -> GetLeaf("Br_PULLIST") -> GetValue(3);
    pull_t1 = ALLCHAIN_CUT -> GetLeaf("Br_PULLIST") -> GetValue(4);

    trigger_indx = ALLCHAIN_CUT -> GetLeaf("Br_trigger_indx") -> GetValue(0);
    filfo_indx = ALLCHAIN_CUT -> GetLeaf("Br_filfo_indx") -> GetValue(0);
    //filfo28_indx = ALLCHAIN_CUT -> GetLeaf("Br_filfo28_indx") -> GetValue(0);
    evtcls_indx = ALLCHAIN_CUT -> GetLeaf("Br_evtcls_indx") -> GetValue(0);
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

    Eprompt_max = 0.;
    if (Eisr > Eprompt_max) Eprompt_max = Eisr;
    if (Epi0_pho1 > Eprompt_max) Eprompt_max = Epi0_pho1;
    if (Epi0_pho2 > Eprompt_max) Eprompt_max = Epi0_pho2;

    //cout << Epi0_pho1 << ", " << Epi0_pho2 << ", " << Eisr << ", Eprompt_max = " << Eprompt_max << endl;

    pi0gam1 = Get4vector(pho_E1, pho_px1, pho_py1, pho_pz1);
    pi0gam2 = Get4vector(pho_E2, pho_px2, pho_py2, pho_pz2);
    isrgam = Get4vector(pho_E3, pho_px3, pho_py3, pho_pz3);
    trkplus = Get4vector(ppl_E, ppl_px, ppl_py, ppl_pz);
    trkmin = Get4vector(pmi_E, pmi_px, pmi_py, pmi_pz);
  
    mpi0 = (pi0gam1 + pi0gam2).M();
    m3pi = (pi0gam1 + pi0gam2 + trkplus + trkmin).M();

    //cout << mpi0 << ", " << m3pi << endl;
    
    evnt_tot ++;

    //if (evnt_tot > 1e5) break;

    //// preselection cut
    //if (trigger_indx == 0) continue; // CUT I
    evnt_trigger ++;
    //cout << trigger_indx << endl;

    //if (filfo_indx == 0) continue; // CUT II
    evnt_filfo ++;
    //cout << filfo_indx << endl;

    //if (evtcls_indx == 0) continue; // CUT III
    evnt_evtcls ++;
    
    //// background rejection
    //if (Eprompt_max > 320) continue;
    
    //if (lagvalue_min_7C > chi2_cut || deltaE > deltaE_cut || angle_pi0gam12 > angle_cut || betapi0 > GetFBeta(beta_cut, c0, c1, ppIM) || Eprompt_max > 320 || IM_pi0_7C > 142. || IM_pi0_7C < 126.) {
      
    if (lagvalue_min_7C > chi2_cut || deltaE > deltaE_cut || angle_pi0gam12 > angle_cut || betapi0 > GetFBeta(beta_cut, c0, c1, ppIM) || Eprompt_max > 320) {
    
      sel_indx = 0;
      //cout << sel_indx << endl;
      //cout << lagvalue_min_7C << endl;
    }
    else {
      sel_indx = 1;
    }

    /*
    if (lagvalue_min_7C > chi2_cut) continue;
    else if (deltaE > deltaE_cut) continue;
    else if (angle_pi0gam12 > angle_cut) continue;
    else if (betapi0 > GetFBeta(beta_cut, c0, c1, ppIM)) continue;
    */
    
    evnt_sel ++;
    
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
       << "evnt_trigger = " << evnt_trigger << "\n"
       << "evnt_filfo = " << evnt_filfo << "\n"
       << "evnt_evtcls = " << evnt_evtcls << "\n"
       << "evnt_sel = " << evnt_sel << "\n";
  
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

  
  
