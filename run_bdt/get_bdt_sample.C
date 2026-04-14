#include "../header_bdt/cut_para.h"

TLorentzVector Get4vector(double E, double px, double py, double pz) {

  //given a cluster index returns the 4-mom of a photon
  TVector3 gamma(px, py, pz);
  
  Double_t scale1;
  scale1=E/gamma.Mag();
  TLorentzVector gamma4mom(scale1*gamma, E);
  //cout << gamma4mom.M() << endl;
  return gamma4mom;

}

int get_bdt_sample(const TString data_type = "sig", const TString br_type = "TISR3PI_SIG", const TString inputFolder = "/home/kloe/Desktop/input_bdt_TDATA_chain/input/", const TString outputFolder = "/home/kloe/Desktop/KLOE_BDT/dataset/") {

  // from e+ e- -> omega gamma get 3 photon final state 4-momentum

  // Proceed.C
  TString sampleFile = inputFolder + data_type + ".root";
  //TString sampleFile = "/home/kloe/Desktop/" + data_type + ".root";
  
  TFile *f_input = new TFile(sampleFile);

  TTree *ALLCHAIN_CUT = (TTree*)f_input -> Get("ALLCHAIN_CUT");

  TFile *f_output = new TFile(outputFolder + data_type + "_bdt.root", "recreate");

  cout << "Creating BDT data ../n"
       << "Input root file: " << sampleFile << "\n"
       << "Output root file: " << outputFolder << f_output->GetName() << endl;
  
  
  const int list_size = 1;
  const TString TNM[list_size] = {br_type};

  TTree *TTList[list_size];
  
  TCollection* tree_list = new TList;

  for (int i = 0; i < list_size; i ++) {// loop over the Tree list

    TTList[i] = new TTree(TNM[i], "recreate");
    TTList[i] -> SetAutoSave(0);
    
    tree_list -> Add(TTList[i]);
    
  } // end the Tree list

  // define variables
  double pho_E1 = 0., pho_px1 = 0., pho_py1 = 0., pho_pz1 = 0.;
  double pho_E2 = 0., pho_px2 = 0., pho_py2 = 0., pho_pz2 = 0.;
  double pho_E3 = 0., pho_px3 = 0., pho_py3 = 0., pho_pz3 = 0.;
  double ppl_E = 0., ppl_px = 0., ppl_py = 0., ppl_pz = 0.;
  double pmi_E = 0., pmi_px = 0., pmi_py = 0., pmi_pz = 0.;
  double mpi0 = 0., mpi0_good = 0., mpi0_bad = 0.;
  double m3pi = 0., m3pi_good = 0., m3pi_bad = 0.;
  
  double lagvalue_min_7C = 0.;
  double deltaE = 0.;
  double angle_pi0gam12 = 0.;
  double betapi0 = 0.;
    
  int phid = -1, sig_type = -1;
  int recon_indx = -1, bkg_indx = -1;
  
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
    
    tree_tmp -> Branch("Br_mpi0", &mpi0, "Br_mpi0/D");
    tree_tmp -> Branch("Br_mpi0_good", &mpi0_good, "Br_mpi0_good/D");
    tree_tmp -> Branch("Br_mpi0_bad", &mpi0_bad, "Br_mpi0_bad/D");

    tree_tmp -> Branch("Br_m3pi", &m3pi, "Br_m3pi/D");
    tree_tmp -> Branch("Br_m3pi_good", &m3pi_good, "Br_m3pi_good/D");
    tree_tmp -> Branch("Br_m3pi_bad", &m3pi_bad, "Br_m3pi_bad/D");

    tree_tmp -> Branch("Br_phid", &phid, "Br_phid/I");
    tree_tmp -> Branch("Br_sig_type", &sig_type, "Br_sig_type/I");
    tree_tmp -> Branch("Br_recon_indx", &recon_indx, "Br_recon_indx/I");
    tree_tmp -> Branch("Br_bkg_indx", &bkg_indx, "Br_bkg_indx/I");
    
    tree_tmp -> Branch("Br_lagvalue_min_7C", &lagvalue_min_7C, "Br_lagvalue_min_7C/D");
    tree_tmp -> Branch("Br_betapi0", &betapi0, "Br_betapi0/D");
    tree_tmp -> Branch("Br_deltaE", &deltaE, "Br_deltaE/D");
    tree_tmp -> Branch("Br_angle_pi0gam12", &angle_pi0gam12, "Br_angle_pi0gam12/D");
    
  }

  TLorentzVector pi0gam1, pi0gam2, isrgam, trkplus, trkmin;

  const double m3pi_min = 0., m3pi_max = 1000.;
  const int bins = 300;
  
  TH1D *hmpi0 = new TH1D("hmpi0", "hmpi0", bins, m3pi_min, m3pi_max);
  TH1D *hmpi0_good = new TH1D("hmpi0_good", "hmpi0_good", bins, m3pi_min, m3pi_max);
  TH1D *hmpi0_bad = new TH1D("hmpi0_bad", "hmpi0_bad", bins, m3pi_min, m3pi_max);

  TH1D *hm3pi = new TH1D("hm3pi", "hm3pi", bins, m3pi_min, m3pi_max);
  TH1D *hm3pi_good = new TH1D("hm3pi_good", "hm3pi_good", bins, m3pi_min, m3pi_max);
  TH1D *hm3pi_bad = new TH1D("hm3pi_bad", "hm3pi_bad", bins, m3pi_min, m3pi_max);
  
  for (Int_t irow = 0; irow < ALLCHAIN_CUT -> GetEntries(); irow ++) {// loop trees
	  
    ALLCHAIN_CUT -> GetEntry(irow);

    phid = ALLCHAIN_CUT -> GetLeaf("Br_phid") -> GetValue(0);
    sig_type = ALLCHAIN_CUT -> GetLeaf("Br_sig_type") -> GetValue(0);
    recon_indx = ALLCHAIN_CUT -> GetLeaf("Br_recon_indx") -> GetValue(0);
    bkg_indx = ALLCHAIN_CUT -> GetLeaf("Br_bkg_indx") -> GetValue(0);

    lagvalue_min_7C = ALLCHAIN_CUT -> GetLeaf("Br_lagvalue_min_7C") -> GetValue(0);
    betapi0 = ALLCHAIN_CUT -> GetLeaf("Br_betapi0") -> GetValue(0);
    deltaE = ALLCHAIN_CUT -> GetLeaf("Br_ENERGYLIST") -> GetValue(2); 
    angle_pi0gam12 = ALLCHAIN_CUT -> GetLeaf("Br_ANGLELIST") -> GetValue(0);
    
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

    pi0gam1 = Get4vector(pho_E1, pho_px1, pho_py1, pho_pz1);
    pi0gam2 = Get4vector(pho_E2, pho_px2, pho_py2, pho_pz2);
    isrgam = Get4vector(pho_E3, pho_px3, pho_py3, pho_pz3);
    trkplus = Get4vector(ppl_E, ppl_px, ppl_py, ppl_pz);
    trkmin = Get4vector(pmi_E, pmi_px, pmi_py, pmi_pz);
  
    //cout << pho_E3 << ", " << pho_px3 << ", " << pho_py3 << ", " << pho_pz3 << endl;
    //cout << ppl_E << ", " << ppl_px << ", " << ppl_py << ", " << ppl_pz << endl;
    //cout << pmi_E << ", " << pmi_px << ", " << pmi_py << ", " << pmi_pz << endl;
    //cout << lagvalue_min_7C << endl;
    //cout << bkg_indx << ", " << recon_indx << endl;
      
    // pi0 invariant mass, identified pi0 photons
    mpi0 = (pi0gam1 + pi0gam2).M();
    m3pi = (pi0gam1 + pi0gam2 + trkplus + trkmin).M();

    // cuts
    //if (lagvalue_min_7C > chi2_cut) continue;
    
    hmpi0 -> Fill(mpi0);
    
    //cout << "pho_E1 = " << pho_E1 << endl;
      
    // fill trees for eta gamma
    if (!TMath::IsNaN(mpi0) && !TMath::IsNaN(m3pi)) {

      hm3pi -> Fill(m3pi);
    
      TTList[0]-> Fill();

      if (recon_indx == 2 && bkg_indx == 1) {// good events
	mpi0_good = mpi0;
	m3pi_good = m3pi;

	hmpi0_good -> Fill(mpi0);
	hm3pi_good -> Fill(m3pi);
    
	//cout << mpi0_good << endl;
      }
      else {// bad
	mpi0_bad = mpi0;
	m3pi_bad = m3pi;

	hmpi0_bad -> Fill(mpi0);
	hm3pi_bad -> Fill(m3pi);
    
      }
    
    }

    
  }
  
  //hm3pi -> Draw();
  //hm3pi_good -> Draw("HistSame");
  //hm3pi_bad -> Draw("HistSame");
  
  // save
  TTList[0]-> Write();
  
  return 0;

}
