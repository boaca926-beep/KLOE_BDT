#define MyClass_cxx
#include "../header_bdt/MyClass.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TDecompSVD.h"
#include <sstream>
#include "TError.h"

void MyClass::Main()
{

  gErrorIgnoreLevel = kFatal; //kFatal;

  cout << "PROCESSING ..." << endl;
  int sig_type = -1;
  int bkg_indx = -1;
  int recon_indx = -1;


  int fiduial_indx = -1;
  int trigger_indx = -1;
  int filfo_indx = -1;
  int evtcls_indx = -1;

  int test_indx = 0;     

  double evnt_sum = 0, evnt_pre_filtred = 0, evnt_pre_pass = 0;
  double evnt_unstr = 0, evnt_cls = 0;
  double evnt_sig = 0, evnt_bkg = 0;
  double evnt_pre = 0, evnt_str = 0, evnt_trig = 0, evnt_filfo = 0, evnt_trk = 0, evnt_class = 0, evnt_final = 0.;
  double evnt_phi5 = 0;
  double evnt_recon_type0 = 0., evnt_recon_type1 = 0., evnt_recon_type2 = 0.;
  // reconstruction
  int trknb_prompt = 0, trkindx1 = 0, trkindx2 = 0;

  int promptnb = 0;
  int evt_indx = -1;
  
  // vertex
  double Zv = 0;
  double Rhov = 0;
   
  // mass
  double IM_3pi = 0., IM_3pi_nofit = 0., IM3pi_7C = 0., IM3pi_true = 0., IM3pi = 0., IM3pi_beta = 0.;
  double IM_pi0 = 0., IM_pi0_nofit = 0., IM_pi0_7C = 0.;
  double trkmass_true = 0., trkmass = 0.;
  double ppIM_true = 0., ppIM = 0., ppIM_beta = 0.;
  double IMisrpho_miss = 0.;
  double IM3pi_pi12 = 0., IM3pi_pi13 = 0., IM3pi_pi23 = 0.; 
  double MASSLIST[100];
  double ANGLELIST[100];
  double PULLIST[100];
  double RESOLIST[100];
  double MOM4PHO1[4];
  double MOM4PHO2[4];
  double MOM4PHO3[4];
  double MOM4TRKPLUS[4];
  double MOM4TRKMINS[4];
  
  // energy
  double ENERGYLIST[100];
  double PI0PHORESD[100];
  double E_pho_isr = 0., Emax_pho = 0.;
  double EPI0GAM[2], EPI0NTMC[2], E_pi0gam1 = 0., E_pi0gam2 = 0.;
  double deltaE_true = 0.;
  double Emax_clust = 0., Esum_clust = 0.;
  double Esum = 0., E_radiv1 = 0., E_radiv2 = 0.;
  double deltaE = 0.;
  // angle
  double Angle_pho_isr = 0.;
  double Angle_clust = 0.;
  double angle_pi0gam12_true = 0.;
  double angle_ppl_3piboost = 0.;
  double minangle = 23.;
  
  // time
  double Tof_clust = 0.;
  double Sigma2_T_clust = 0., Sigma_T_clust = 0.;
  // kinfit
  const double xyres = 1.2, zres = 1.4, tres = 0.057, tcst = 0.140, eres = 0.057;
  // statistics
  double lagvalue_min_7C = 0., pvalue = 0.;
  double chi2mgg_min_7C = 0.;
  // others
  double tclusdiff = 0.;
  double beta_clu = 0.;
  double betapi0 = 0., betapi0_true = 0.;
  //
  double test_value = 0.;
  
  // define tree
  TTree ALLCHAIN_GEN("ALLCHAIN_GEN", "recreate"); ALLCHAIN_GEN.SetAutoSave(0);
  
  ALLCHAIN_GEN.Branch("Br_sig_type", &sig_type, "Br_sig_type/I");
  ALLCHAIN_GEN.Branch("Br_sel_type", &bit_select, "Br_sel_type/I");
  ALLCHAIN_GEN.Branch("Br_IM_3pi", &IM_3pi, "Br_IM_3pi/D");
  ALLCHAIN_GEN.Branch("Br_E_pho_isr", &E_pho_isr, "Br_E_pho_isr/D");
  ALLCHAIN_GEN.Branch("Br_Angle_pho_isr", &Angle_pho_isr, "Br_Angle_pho_isr/D");
  ALLCHAIN_GEN.Branch("Br_bpx", &bpx, "Br_bpx/F");
  ALLCHAIN_GEN.Branch("Br_Esum", &Esum, "Br_Esum/D");

  TTree ALLCHAIN_STR2("ALLCHAIN_STR2", "recreate"); ALLCHAIN_STR2.SetAutoSave(0);
  
  ALLCHAIN_STR2.Branch("Br_sig_type", &sig_type, "Br_sig_type/I");
  ALLCHAIN_STR2.Branch("Br_sel_type", &bit_select, "Br_sel_type/I");
  ALLCHAIN_STR2.Branch("Br_IM_3pi", &IM_3pi, "Br_IM_3pi/D");
  
  TTree TrSample ("TrSample", "recreate");
  TrSample.SetAutoSave(0);
  TrSample.Branch("Br_IM_3pi", &IM_3pi, "Br_IM_3pi/D");
  TrSample.Branch("Br_E_pho_isr", &E_pho_isr, "Br_E_pho_isr/D");
  TrSample.Branch("Br_Angle_pho_isr", &Angle_pho_isr, "Br_Angle_pho_isr/D");
  TrSample.Branch("Br_bpx", &bpx, "Br_bpx/F");
  TrSample.Branch("Br_Esum", &Esum, "Br_Esum/D");

  double pho_E1 = 0., pho_px1 = 0., pho_py1 = 0., pho_pz1 = 0.;
  double pho_E2 = 0., pho_px2 = 0., pho_py2 = 0., pho_pz2 = 0.;
  double pho_E3 = 0., pho_px3 = 0., pho_py3 = 0., pho_pz3 = 0.;
  double ppl_E = 0., ppl_px = 0., ppl_py = 0., ppl_pz = 0.;
  double pmi_E = 0., pmi_px = 0., pmi_py = 0., pmi_pz = 0.;

  double pho_E1_true = 0., pho_px1_true = 0., pho_py1_true = 0., pho_pz1_true = 0.;
  double pho_E2_true = 0., pho_px2_true = 0., pho_py2_true = 0., pho_pz2_true = 0.;
  double pho_E3_true = 0., pho_px3_true = 0., pho_py3_true = 0., pho_pz3_true = 0.;
  
  //
  TTree ALLCHAIN_CUT ("ALLCHAIN_CUT", "recreate");
  ALLCHAIN_CUT.SetAutoSave(0);
  ALLCHAIN_CUT.Branch("Br_sig_type", &sig_type, "Br_sig_type/I");
  ALLCHAIN_CUT.Branch("Br_sel_type", &bit_select, "Br_sel_type/I");
  ALLCHAIN_CUT.Branch("Br_phid", &phid, "Br_phid/I");
  ALLCHAIN_CUT.Branch("Br_kineid", &kineid, "Br_kineid/I");
  ALLCHAIN_CUT.Branch("Br_bkg_indx", &bkg_indx, "Br_bkg_indx/I");
  ALLCHAIN_CUT.Branch("Br_recon_indx", &recon_indx, "Br_recon_indx/I");
  ALLCHAIN_CUT.Branch("Br_trigger_indx", &trigger_indx, "Br_trigger_indx/I");
  ALLCHAIN_CUT.Branch("Br_evtcls_indx", &evtcls_indx, "Br_evtcls_indx/I");
  ALLCHAIN_CUT.Branch("Br_filfo_indx", &filfo_indx, "Br_filfo_indx/I");
  
  //
  //ALLCHAIN_CUT.Branch("Br_Rhov", &Rhov, "Br_Rhov/D");  
  //ALLCHAIN_CUT.Branch("Br_Zv", &Zv, "Br_Zv/D");  
  //ALLCHAIN_CUT.Branch("Br_Tof_clust", &Tof_clust, "Br_Tof_clust/D");
  ALLCHAIN_CUT.Branch("Br_egammamin", &egammamin, "Br_egammamin/D");  
  ALLCHAIN_CUT.Branch("Br_IM_3pi", &IM_3pi, "Br_IM_3pi/D");
  ALLCHAIN_CUT.Branch("Br_IM3pi_7C", &IM3pi_7C, "Br_IM3pi_7C/D");
  ALLCHAIN_CUT.Branch("Br_IM_3pi_nofit", &IM_3pi_nofit, "Br_IM_3pi_nofit/D");  
  ALLCHAIN_CUT.Branch("Br_IM3pi_pi12", &IM3pi_pi12, "Br_IM3pi_pi12/D"); 
  ALLCHAIN_CUT.Branch("Br_IM3pi_pi13", &IM3pi_pi13, "Br_IM3pi_pi13/D"); 
  ALLCHAIN_CUT.Branch("Br_IM3pi_pi23", &IM3pi_pi23, "Br_IM3pi_pi23/D"); 
  ALLCHAIN_CUT.Branch("Br_IM3pi_true", &IM3pi_true, "Br_IM3pi_true/D"); 
  ALLCHAIN_CUT.Branch("Br_IM_pi0", &IM_pi0, "Br_IM_pi0/D");
  ALLCHAIN_CUT.Branch("Br_IM_pi0_nofit", &IM_pi0_nofit, "Br_IM_pi0_nofit/D");  
  ALLCHAIN_CUT.Branch("Br_IM_pi0_7C", &IM_pi0_7C, "Br_IM_pi0_7C/D");
  ALLCHAIN_CUT.Branch("Br_trkmass_true", &trkmass_true, "Br_trkmass_true/D");
  ALLCHAIN_CUT.Branch("Br_trkmass", &trkmass, "Br_trkmass/D");
  ALLCHAIN_CUT.Branch("Br_ppIM_true", &ppIM_true, "Br_ppIM_true/D");
  ALLCHAIN_CUT.Branch("Br_ppIM", &ppIM, "Br_ppIM/D");
  ALLCHAIN_CUT.Branch("Br_ppIM_beta", &ppIM_beta, "Br_ppIM_beta/D");
  ALLCHAIN_CUT.Branch("Br_E_pho_isr", &E_pho_isr, "Br_E_pho_isr/D");
  ALLCHAIN_CUT.Branch("Br_deltaE_true", &deltaE_true, "Br_deltaE_true/D");
  ALLCHAIN_CUT.Branch("Br_Angle_pho_isr", &Angle_pho_isr, "Br_Angle_pho_isr/D"); 
  ALLCHAIN_CUT.Branch("Br_angle_pi0gam12_true", &angle_pi0gam12_true, "Br_angle_pi0gam12_true/D");
  ALLCHAIN_CUT.Branch("Br_Emax_clust", &Emax_clust, "Br_Emax_clust/D");
  ALLCHAIN_CUT.Branch("Br_Esum_clust", &Esum_clust, "Br_Esum_clust/D");

  //
  ALLCHAIN_CUT.Branch("Br_lagvalue_min_7C", &lagvalue_min_7C, "Br_lagvalue_min_7C/D");
  ALLCHAIN_CUT.Branch("Br_pvalue", &pvalue, "Br_pvalue/D");
  ALLCHAIN_CUT.Branch("Br_betapi0", &betapi0, "Br_betapi0/D");
  ALLCHAIN_CUT.Branch("Br_betapi0_true", &betapi0_true, "Br_betapi0_true/D");
  ALLCHAIN_CUT.Branch("Br_test", &test_value, "Br_test/D");
  ALLCHAIN_CUT.Branch("Br_chi2mgg_min_7C", &chi2mgg_min_7C, "Br_chi2mgg_min_7C/D");
  ALLCHAIN_CUT.Branch("Br_beta_clu", &beta_clu, "Br_beta_clu/D");
  ALLCHAIN_CUT.Branch("Br_tclusdiff", &tclusdiff, "Br_tclusdiff/D");
  
  //
  ALLCHAIN_CUT.Branch("Br_ENERGYLIST", &ENERGYLIST, "Br_ENERGYLIST[100]/D");
  //ENERGYLIST[0]: ISR photon energy, 7C kin. fit
  
  ALLCHAIN_CUT.Branch("Br_PI0PHORESD", &PI0PHORESD, "Br_PI0PHORESD[100]/D");
  
  ALLCHAIN_CUT.Branch("Br_MASSLIST", &MASSLIST, "Br_MASSLIST[100]/D");
  //MASSLIST[0]: pi0 mass no fit
  //MASSLIST[1]: 3pi mass no fit
  //MASSLIST[2]: pi0 mass 7C kin. fit
  //MASSLIST[3]: 3pi mass 7C kin. fit
  //MASSLIST[4]: track mass
  //MASSLIST[5]: pi+ pi- invariant mass

  ALLCHAIN_CUT.Branch("Br_ANGLELIST", &ANGLELIST, "Br_ANGLELIST[100]/D");
  //ANGLELIST[0]: pi0 gamma12
  
  ALLCHAIN_CUT.Branch("Br_RESOLIST", &RESOLIST, "Br_RESOLIST[100]/D");
  //RESOLIST[0]: pi0 mass 7C kin. fit
  //RESOLIST[1]: 3pi mass 7C kin. fit
  //RESOLIST[2]: track mass 7C kin. fit
  //RESOLIST[3]: pi+ pi- invariant mass
  //RESOLIST[4]: betapi0

  ALLCHAIN_CUT.Branch("Br_PULLIST", &PULLIST, "Br_PULLIST[100]/D");

  //ALLCHAIN_CUT.Branch("Br_MOM4TRKPLUS", &MOM4TRKPLUS, "MOM4TRKPLUS[4]/D");

  ALLCHAIN_CUT.Branch("Br_evt_indx", &evt_indx, "Br_evt_indx/I");
  ALLCHAIN_CUT.Branch("Br_E1", &pho_E1, "Br_E1/D");
  ALLCHAIN_CUT.Branch("Br_px1", &pho_px1, "Br_px1/D");
  ALLCHAIN_CUT.Branch("Br_py1", &pho_py1, "Br_py1/D");
  ALLCHAIN_CUT.Branch("Br_pz1", &pho_pz1, "Br_pz1/D");

  ALLCHAIN_CUT.Branch("Br_E2", &pho_E2, "Br_E2/D");
  ALLCHAIN_CUT.Branch("Br_px2", &pho_px2, "Br_px2/D");
  ALLCHAIN_CUT.Branch("Br_py2", &pho_py2, "Br_py2/D");
  ALLCHAIN_CUT.Branch("Br_pz2", &pho_pz2, "Br_pz2/D");

  ALLCHAIN_CUT.Branch("Br_E3", &pho_E3, "Br_E3/D");
  ALLCHAIN_CUT.Branch("Br_px3", &pho_px3, "Br_px3/D");
  ALLCHAIN_CUT.Branch("Br_py3", &pho_py3, "Br_py3/D");
  ALLCHAIN_CUT.Branch("Br_pz3", &pho_pz3, "Br_pz3/D");

  ALLCHAIN_CUT.Branch("Br_ppl_E", &ppl_E, "Br_ppl_E/D");
  ALLCHAIN_CUT.Branch("Br_ppl_px", &ppl_px, "Br_ppl_px/D");
  ALLCHAIN_CUT.Branch("Br_ppl_py", &ppl_py, "Br_ppl_py/D");
  ALLCHAIN_CUT.Branch("Br_ppl_pz", &ppl_pz, "Br_ppl_pz/D");

  ALLCHAIN_CUT.Branch("Br_pmi_E", &pmi_E, "Br_pmi_E/D");
  ALLCHAIN_CUT.Branch("Br_pmi_px", &pmi_px, "Br_pmi_px/D");
  ALLCHAIN_CUT.Branch("Br_pmi_py", &pmi_py, "Br_pmi_py/D");
  ALLCHAIN_CUT.Branch("Br_pmi_pz", &pmi_pz, "Br_pmi_pz/D");

  // mc truth
  ALLCHAIN_CUT.Branch("Br_E1_true", &pho_E1_true, "Br_E1_true/D");
  ALLCHAIN_CUT.Branch("Br_px1_true", &pho_px1_true, "Br_px1_true/D");
  ALLCHAIN_CUT.Branch("Br_py1_true", &pho_py1_true, "Br_py1_true/D");
  ALLCHAIN_CUT.Branch("Br_pz1_true", &pho_pz1_true, "Br_pz1_true/D");

  ALLCHAIN_CUT.Branch("Br_E2_true", &pho_E2_true, "Br_E2_true/D");
  ALLCHAIN_CUT.Branch("Br_px2_true", &pho_px2_true, "Br_px2_true/D");
  ALLCHAIN_CUT.Branch("Br_py2_true", &pho_py2_true, "Br_py2_true/D");
  ALLCHAIN_CUT.Branch("Br_pz2_true", &pho_pz2_true, "Br_pz2_true/D");

  ALLCHAIN_CUT.Branch("Br_E3_true", &pho_E3_true, "Br_E3_true/D");
  ALLCHAIN_CUT.Branch("Br_px3_true", &pho_px3_true, "Br_px3_true/D");
  ALLCHAIN_CUT.Branch("Br_py3_true", &pho_py3_true, "Br_py3_true/D");
  ALLCHAIN_CUT.Branch("Br_pz3_true", &pho_pz3_true, "Br_pz3_true/D");

  //
  TTree ALLCHAIN_TEST ("ALLCHAIN_TEST", "recreate"); ALLCHAIN_TEST.SetAutoSave(0);
  ALLCHAIN_TEST.Branch("Br_bkg_indx", &bkg_indx, "Br_bkg_indx/I");
  ALLCHAIN_TEST.Branch("Br_IM3pi_7C", &IM3pi_7C, "Br_IM3pi_7C/D");
  ALLCHAIN_TEST.Branch("Br_IMisrpho_miss", &IMisrpho_miss, "Br_IMisrpho_miss/D");
  ALLCHAIN_TEST.Branch("Br_angle_ppl_3piboost", &angle_ppl_3piboost, "Br_angle_ppl_3piboost/D");
  
  ///
  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    //if (jentry != 228610) continue;
    //if (jentry > 1) continue; 
    
    //cout << "jentry = " << jentry << endl;

    evnt_sum ++;
    
    if (bit_select == 1) {// pre-filtred
      evnt_pre_filtred ++;
    }
    else {
      evnt_pre_pass ++;
    }

    /// beam
    Beam.SetPxPyPzE(bpx, bpy, bpz, bene);
    //TVector3 Boost3vector = -Beam.BoostVector();
    
    //cout << "beam (bpx, bpy, bpz, bene) = " << "(" << bpx << ", " << bpy << ", " << bpz << ", " << bene << ") \n";

    /// define isr 3pi signals
    double nb_pho_radiv = 0, nb_pi = 0, nb_pi0pho = 0;
    int pi0gam1_ntmc = 0, pi0gam2_ntmc = 0; 
    TVector3 MC_vect;
    TLorentzVector pi0MC_TLvect, piplusMC_TLvect, piminusMC_TLvect, threepi_TLvect, isrpho_TLvect, finalstate_TLvect, pho_radiv1_TLvect, pho_radiv2_TLvect, pi0gam1_TLvect, pi0gam2_TLvect;
   
      
    /// paticle infomation
    //cout << "\n\n" << endl;
    for (int i = 0; i < ntmc; i ++) {// loop over emc
      MC_vect.SetXYZ(pxmc[i], pymc[i], pzmc[i]);

      //if (phid == 5) {
      //evnt_phi5 ++;
      //cout << i << ": kineid = " << kineid << ", phid = " << phid << ", pidmc = " << pidmc[i] << ", virmom = " << virmom[i] << ", mother = " << mother[indv[i] - 1] << ", pxmc = " << pxmc[i] << ", pymc = " << pymc[i] << ", pzmc = " << pzmc[i] << "\n";
      //}
      
      // radiative photons (isr str3)
      if (pidmc[i] == 1 && virmom[i] == 1 && mother[indv[i] - 1] == 50 && nb_pho_radiv == 0) {// first radiative photon
	nb_pho_radiv ++;
	pho_radiv1_TLvect = GetLorentzVector(MC_vect, 0.);
      }
      else if (pidmc[i] == 1 && virmom[i] == 1 && mother[indv[i] - 1] == 50 && nb_pho_radiv == 1) {// second radiative photon
	nb_pho_radiv ++;
	pho_radiv2_TLvect = GetLorentzVector(MC_vect, 0.);
      }
      
      // radiative photons (isr phok5)
      else if (pidmc[i] == 1 && virmom[i] == 0 && mother[indv[i] - 1] == 50 && nb_pho_radiv == 0) {// first radiative photon
	nb_pho_radiv ++;
	pho_radiv1_TLvect = GetLorentzVector(MC_vect, 0.);
      }
      else if (pidmc[i] == 1 && virmom[i] == 0 && mother[indv[i] - 1] == 50 && nb_pho_radiv == 1) {// second radiative photon
	nb_pho_radiv ++;
	pho_radiv2_TLvect = GetLorentzVector(MC_vect, 0.);
      }
      
      // pions (isr str3)
      else if (pidmc[i] == 7 && virmom[i] == 54 && mother[indv[i] - 1] == 50) {// pi0
	nb_pi ++;
	pi0MC_TLvect = GetLorentzVector(MC_vect, massneupion);
	//cout << "vect.X = " << MC_vect.X() << endl;
	//cout << "TLvect.X = " << pi0MC_TLvect.X() << endl;
      }
      else if (pidmc[i] == 8 && virmom[i] == 54 && mother[indv[i] - 1] == 50) {// piplus
	nb_pi ++;
	piplusMC_TLvect = GetLorentzVector(MC_vect, masschpion);
      }
      else if (pidmc[i] == 9 && virmom[i] == 54 && mother[indv[i] - 1] == 50) {// piminus
	nb_pi ++;
	//cout << "vect.X = " << MC_vect.X() << endl;
	piminusMC_TLvect = GetLorentzVector(MC_vect, masschpion);
      }

      // pions (isr phok5)
      else if (pidmc[i] == 7 && virmom[i] == 0 && mother[indv[i] - 1] == 50) {// pi0
	nb_pi ++;
	pi0MC_TLvect = GetLorentzVector(MC_vect, massneupion);
	//cout << "vect.X = " << MC_vect.X() << endl;
	//cout << "TLvect.X = " << pi0MC_TLvect.X() << endl;
      }
      else if (pidmc[i] == 8 && virmom[i] == 0 && mother[indv[i] - 1] == 50) {// piplus
	nb_pi ++;
	piplusMC_TLvect = GetLorentzVector(MC_vect, masschpion);
      }
      else if (pidmc[i] == 9 && virmom[i] == 0 && mother[indv[i] - 1] == 50) {// piminus
	nb_pi ++;
	//cout << "vect.X = " << MC_vect.X() << endl;
	piminusMC_TLvect = GetLorentzVector(MC_vect, masschpion);
      }

      // pions (eta->3pi)
      else if (pidmc[i] == 7 && virmom[i] == 17 && mother[indv[i] - 1] == 50) {// pi0
	nb_pi ++;
	pi0MC_TLvect = GetLorentzVector(MC_vect, massneupion);
	//cout << "vect.X = " << MC_vect.X() << endl;
	//cout << "TLvect.X = " << pi0MC_TLvect.X() << endl;
      }
      else if (pidmc[i] == 8 && virmom[i] == 17 && mother[indv[i] - 1] == 50) {// piplus
	nb_pi ++;
	piplusMC_TLvect = GetLorentzVector(MC_vect, masschpion);
      }
      else if (pidmc[i] == 9 && virmom[i] == 17 && mother[indv[i] - 1] == 50) {// piminus
	nb_pi ++;
	//cout << "vect.X = " << MC_vect.X() << endl;
	piminusMC_TLvect = GetLorentzVector(MC_vect, masschpion);
      }
      
      // pi0 photons
      else if (pidmc[i] == 1 && virmom[i] == 0 && mother[indv[i] - 1] == 7 && nb_pi0pho == 0) {
	nb_pi0pho ++;
	pi0gam1_TLvect = GetLorentzVector(MC_vect, 0.);
	pi0gam1_ntmc = i;
	//cout << pxmc[i] << endl;
	
      }
      else if (pidmc[i] == 1 && virmom[i] == 0 && mother[indv[i] - 1] == 7 && nb_pi0pho == 1) {
	nb_pi0pho ++;
	pi0gam2_TLvect = GetLorentzVector(MC_vect, 0.);
	pi0gam2_ntmc = i;

	//cout << pxmc[i] << endl;
	
      }

    }// end loop emc

    /// MC info
    finalstate_TLvect = piminusMC_TLvect + piplusMC_TLvect + pi0MC_TLvect + pho_radiv1_TLvect + pho_radiv2_TLvect;
    threepi_TLvect = piminusMC_TLvect + piplusMC_TLvect + pi0gam1_TLvect + pi0gam2_TLvect;
    //threepi_TLvect = piminusMC_TLvect + piplusMC_TLvect + pi0MC_TLvect;
    isrpho_TLvect = Beam - (pi0MC_TLvect + piplusMC_TLvect + piminusMC_TLvect);

    IM_3pi = threepi_TLvect.M();//
    //IM_3pi = (Beam - (pho_radiv1_TLvect + pho_radiv2_TLvect)).M(); //threepi_TLvect.M();
    IM_pi0 = (pi0gam1_TLvect + pi0gam2_TLvect).M();//(Beam - (pho_radiv1_TLvect + pho_radiv2_TLvect)).M(); //threepi_TLvect.M();
    ppIM_true = (piplusMC_TLvect + piminusMC_TLvect).M();
    
    E_radiv1 = pho_radiv1_TLvect.E();
    E_radiv2 = pho_radiv2_TLvect.E();
    E_pho_isr = isrpho_TLvect.E();
    Esum = finalstate_TLvect.E();
    EPI0GAM[0] = pi0gam1_TLvect.E();
    EPI0GAM[1] = pi0gam2_TLvect.E();
    EPI0NTMC[0] = pi0gam1_ntmc; 
    EPI0NTMC[1] = pi0gam2_ntmc;
    E_pi0gam1 = pi0gam1_TLvect.E();
    E_pi0gam2 = pi0gam2_TLvect.E();
    betapi0_true = (pi0MC_TLvect.Vect()).Mag() / pi0MC_TLvect.E();
    deltaE_true = DeltaE(piplusMC_TLvect, piminusMC_TLvect);
    trkmass_true = Trkmass(piplusMC_TLvect, piminusMC_TLvect); 
    Angle_pho_isr = TMath::Cos(isrpho_TLvect.Theta());
    angle_pi0gam12_true = pi0gam1_TLvect.Angle(pi0gam2_TLvect.Vect())*TMath::RadToDeg();
    //cout << "deltaE_true = " << deltaE_true << endl;
    //cout << "IM_3pi = " << IM_3pi << ", E_pho_isr = " << E_pho_isr << ", Angle_pho_isr = " << Angle_pho_isr << "\n";
    //cout << "nb_pho_radiv = " << nb_pho_radiv << ", nb_pi = " << nb_pi << ", nb_pi0pho = " << nb_pi0pho << "\n";
    //cout << "pi0gam1 E = " << E_pi0gam1 << ", pi0gam2 E = " << E_pi0gam2 << endl;
      
    /// define ISR signal

    if (nb_pho_radiv > 0 && nb_pi == 3 && nb_pi0pho == 2) {
      sig_type = 1;
      evnt_sig ++;
      
      //cout << phid << endl;
      //cout << nb_pho_radiv << endl;
      //cout << pho_radiv1_TLvect.E() << "," << pho_radiv2_TLvect.E() << ", checked " << E_radiv1 << ", " << E_radiv2 << "\n";
      
      //cout << "jentry = " << jentry << ", finalstate E = " << finalstate_TLvect.E() << endl;
      //if (phid == 5) cout << "nb_pho_radiv = " << nb_pho_radiv << ", nb_pi = " << nb_pi << ", nb_pi0pho = " << nb_pi0pho << endl;
      //TrSample.Fill();
    }
    else {
      sig_type = 0;
      evnt_bkg ++;
    }
    
    // check etagam->3pigamma
    //if (phid == 5) {
    //}

    //cout << sig_type << endl;
    ALLCHAIN_GEN.Fill();

    /// cluster info
    TVector3 Pos_clust(0., 0., 0.);
    int promptnb = 0;
    double Emax_clust_tmp = 0., rt_clust_tmp = 0.;
    double Esum_clust_tmp = 0.;
    
    for (int i = 0; i < nclu; i ++) {// looping over clusters

      //Pos_clust.SetXYZ(xcl[i] - bx, ycl[i] - by, zcl[i] - bz); // cluster position vector
      Pos_clust.SetXYZ(xcl[i], ycl[i], zcl[i]); // cluster position vector
      
      Angle_clust = Pos_clust.Theta() * TMath::RadToDeg(); // cluster polar angle
      Tof_clust = tcl[i] - Pos_clust.Mag() / speedc; // cluster tof
      Sigma2_T_clust = TMath::Power(tres * TMath::Sqrt(1000. / enecl[i]), 2) + TMath::Power(tcst, 2);
      Sigma_T_clust = TMath::Sqrt(Sigma2_T_clust);

      Esum_clust_tmp += enecl[i];
      // maximum cluster energy
      if (enecl[i] > Emax_clust_tmp) {
	Emax_clust_tmp = enecl[i];
      }
      //cout << "enecl[" << i << "] = " << enecl[i] << ", Emax_clust = " << Emax_clust << "\n";

      //cout << promptnb << ", " << Pnum1[i]-1 << ", pxmc=" << pxmc[Pnum1[i]-1] << ", pymc=" << pymc[Pnum1[i]-1] << ", pzmc=" << pzmc[Pnum1[i]-1] << ", E= " << TMath::Sqrt(pxmc[Pnum1[i]-1]*pxmc[Pnum1[i]-1]+pymc[Pnum1[i]-1]*pymc[Pnum1[i]-1]+pzmc[Pnum1[i]-1]*pzmc[Pnum1[i]-1]) << "\n";
      
	
      // select prompt clusters
      if (charged_new[i] == 0 && intime[i] && enecl[i] > egammamin && Angle_clust > minangle && Angle_clust < (180. - minangle) && TMath::Abs(Tof_clust) < nb_sigma_T_clust * Sigma_T_clust && mmclu[i] != 5) {

	hTof_clust -> Fill(Tof_clust);

	//cout << egammamin << endl;
	//cout << i << ": charged_new = " << charged_new[i] << ", intime = " << intime[i] << ", cluster angle = " << Angle_clust << ", cluster tof = " << Tof_clust << ", mmclu = " << mmclu[i] << ", nb_sigma_T_clust = " << nb_sigma_T_clust << ", sigma T cluster = " << Sigma_T_clust << "\n";

	beta_clu = Pos_clust.Mag() / tcl[i] / speedc;
	tclusdiff = Tof_clust;

	// fill prompt cluster info: E, X, Y, Z and T
	E_list[promptnb] = enecl[i]; 
	X_list[promptnb] = xcl[i]; 
	Y_list[promptnb] = ycl[i];
	Z_list[promptnb] = zcl[i];
	T_list[promptnb] = tcl[i];

	pho_indx[promptnb] = Pnum1[i] - 1;
	clu_indx[promptnb] = i;
	pid_indx[promptnb] = Pid1[i];

	if (kineid != -999) {
	  E_true_list[promptnb] = TMath::Sqrt(pxmc[Pnum1[i]-1] * pxmc[Pnum1[i]-1] + pymc[Pnum1[i]-1] * pymc[Pnum1[i]-1] + pzmc[Pnum1[i]-1] * pzmc[Pnum1[i]-1]);
	}
	
	// sigma square
	rt_clust_tmp = TMath::Sqrt(TMath::Power(xcl[i], 2) + TMath::Power(ycl[i], 2));
	
	if (rt_clust_tmp > 200.) {// hits on barrel
	  Sigma2_X_list[promptnb] = xyres * xyres;
	  Sigma2_Y_list[promptnb] = xyres * xyres;
	  Sigma2_Z_list[promptnb] = TMath::Power(zres * TMath::Sqrt(1000. / enecl[i]), 2);
	}
	else {// hits on end-cup
	  Sigma2_X_list[promptnb] = xyres * xyres;
	  Sigma2_Y_list[promptnb] = TMath::Power(zres * TMath::Sqrt(1000. / enecl[i]), 2);
	  Sigma2_Z_list[promptnb] = xyres * xyres;
	}
	Sigma2_E_list[promptnb] = TMath::Power(eres * TMath::Sqrt(1000. * enecl[i]), 2);
	Sigma2_T_list[promptnb] = TMath::Power(tres * TMath::Sqrt(1000. / enecl[i]), 2) + TMath::Power(tcst, 2);
	promptnb ++;
      }// end prompt clusters section
      
      
    }// end looping over clusters
    Emax_clust = Emax_clust_tmp;
    Esum_clust = Esum_clust_tmp;
    //cout << "promptnb = " << promptnb << "\n";
      
    /// track info
    TVectorD trkvect(3);
    trkvect = Getpiontrnb();
    trknb_prompt = trkvect(0);
    trkindx1 = trkvect(1);
    trkindx2 = trkvect(2);
    Bool_t ifbroken = IfBroken(trkindx1, trkindx2);
    
    // CUT5: cut vertex in confined region
    //if (Zv > Zvmax || Rhov > Rhovmax) continue;
      
    //cout << "jentry: " << jentry << ", Zv_max = " << Zv << ", Rhov_max = " << Rhov << "\n" << endl;
  
    /// select classified event
    
    /*
    Bool_t bTagged = false;
    Bool_t chrad = false;
    for(Int_t NrEC=0; NrEC<necls; NrEC++){
      if(eclstream[NrEC]==pstrnb){
	//only interested in the tag if the stream is 4
	bTagged=true;
	if((ecltagnum[NrEC] &  1 )==1) chrad = true;
      }
    }
    
    
    if(!bTagged) continue;
    */

    // CUT0, bit cut
    if(bit_select == 1) continue; 
    evnt_pre ++;

    //cout << "!!!!!!!!!!!!!!!!!!!!!!!!!" << promptnb << endl;
    hprompt_distr -> Fill(promptnb);
    
    // Event classification

    if (necls == 0) {// UFO
      evnt_unstr ++;
    }
    else {// streamed
      evnt_cls ++;
      for (int NrEC=0; NrEC<necls; NrEC++) {     

	//if (eclstream[NrEC] == 2) {
	hstr_distr -> Fill(eclstream[NrEC]);
	//}
	//cout<<eclstream[NrEC]<<endl;
	//if (eclstream[NrEC]==4) {

	//}
      }
    }

    // CUT1, trigger
    if (IfTriggered()) continue; // SYST. CHECK _EMC: pass EMC trigger,  _DC(): pass DC trigger, _EMCandDC(): pass EMC AND DC triggers     
    trigger_indx = 1;
    evnt_trig ++;
    
    // CUT2, FILFO
    if (IfFilfoed()) continue; // SYST. CHECK, bit 20, standard FILFO cut
    filfo_indx = 1;
    evnt_filfo ++;

    // CUT3, ksl stream
    if (!IfStreamed(pstrnb)) continue; 
    evtcls_indx = 1;
    evnt_str ++;
    //cout << pstrnb << endl;

    double Zv_tmp = 0., Rhov_tmp = 0.;
    int iv_ip = 0, iv_indx = -1;
   
    for (int k = 0; k < nv; k++) {// begin loop
      Zv_tmp = TMath::Abs(zv[k] - bz);
      Rhov_tmp = TMath::Sqrt((xv[k] - bx) * (xv[k] - bx) + (yv[k] - by) * (yv[k] - by));
      //cout << "Zv_tmp = " << Zv_tmp << ", Rhov_tmp = " << Rhov_tmp << endl;
      if (Zv_tmp <= Zvmax && Rhov_tmp <= Rhovmax) { // fiducial volume
	iv_ip ++;
	iv_indx = k;
      }
    }

    // fill histos
    h_iv_ip -> Fill(iv_ip);
    h_iv_indx -> Fill(iv_indx);

    //nv .vs. iv_ip
    h_nv_ip -> Fill(nv, iv_ip);

    //cout << "iv_ip = " << iv_ip << ", iv_indx = " << iv_indx << endl;
    
    //if (iv_ip != 1) continue; // select events with only 1 vertex in the fiducial volume
    //fiduial_indx = 1;
    
    //h_iv_ip_1 -> Fill(iv_ip);
    h_iv_indx_1 -> Fill(iv_indx);

    // CUT4, two tracks and 3 prompt photons
    if (trknb_prompt != 2 || ifbroken ) continue; evnt_trk ++; 
    if (promptnb != 3 ) continue; evnt_class ++;

    //nv .vs. iv_ip after the track and prompt photon selection
    h_nv_ip_1 -> Fill(nv, iv_ip);

    //x/y/zpca(1) .vs. x/y/zpaca(2) for the two tracks
    h_xpca -> Fill(xpca[trkindx1],xpca[trkindx2]);
    h_ypca -> Fill(ypca[trkindx1],ypca[trkindx2]);
    h_zpca -> Fill(zpca[trkindx1],zpca[trkindx2]);
    
    //ALLCHAIN_STR2.Fill();

    /// 7C kinematical fit
    const int Row = 15, Col = 15, row = 7;
    const int nfloop_7C = 5;
	  
    TVectorD inputvect_7C(Row), sigma2vect_7C(Row), etakinfit_min_7C(Row);
    TVectorD sigma2vectorkinfit_min_7C(Row), pullkinfit(Row);

    for (int nr1 = 0; nr1 < promptnb - 2; nr1 ++) {// select 3 photon loop
      for (int nr2 = nr1 + 1; nr2 < promptnb - 1; nr2 ++) {
	for (int nr3 = nr2 + 1; nr3 < promptnb; nr3 ++) {
	  
	  int step = 0;

	  double lagravaluekinfittildeloop = 1e14, lagravaluekinfitloop_temp = 1e14, lagravaluekinfitloop = 1e14;
    
	  TVectorD inputvect(Row), sigma2vector(Row), sigma2vectorkinfitloop(Row);
	  TVectorD etatildekinfitloop(Row), etakinfitloop_temp(Row), etakinfitloop(Row), diffvectorkinfitloop(Row);
	  TVectorD gkinfitloop(row), lambdavectorloop(row), pullsvectorloop(Row);
  
	  TMatrixD Getakinfitloop(row, Col), Getakinfittransloop(Col, row), Vmatrixloop(Row, Col), Skinfitloop(row, row), InSkinfitloop(row, row), Vkinfitloop(Row, Col);

	  inputvect = FillInputvector_7C(Row, nr1, nr2, nr3, trkindx1, trkindx2);
	  //cout << "!!!! inspect inputvector before kin. fit. " << endl;
	  //cout << "nr1 = " << nr1 << ", nr2 = " << nr2 << ", nr3 = " << nr3 << endl;
	  //inputvect.Print();
	  sigma2vector = FillSigma2vector_7C(Row, nr1, nr2, nr3, trkindx1, trkindx2); //sigma2vector.Print();

	  // initialization
	  etakinfitloop_temp = inputvect, etakinfitloop = inputvect;
	  Vmatrixloop = CovMatrix(sigma2vector, Row, Col);
  
	  while (step < nfloop_7C) {
	    etatildekinfitloop = etakinfitloop_temp; //etatildekinfitloop.Print();
	    etakinfitloop_temp = etakinfitloop;

	    diffvectorkinfitloop = inputvect - etatildekinfitloop; //diffvectorkinfitloop.Print();
	    gkinfitloop = Gfunc_7C(etatildekinfitloop, row, Beam, trkindx1, trkindx2); //gkinfitloop.Print();
	    
	    lagravaluekinfittildeloop=lagravaluekinfitloop_temp;
	    lagravaluekinfitloop_temp=lagravaluekinfitloop;
	    
	    Getakinfitloop = Getafunc_7C(etatildekinfitloop, row, Col); //Getakinfitloop.Print();
	    Getakinfittransloop = Mtrans(Getakinfitloop);
	    Skinfitloop = Getakinfitloop * (Vmatrixloop * Getakinfittransloop); //Skinfitloop.Print();
	    InSkinfitloop = MInvert(Skinfitloop);
	    lambdavectorloop = Lambdavector(row, InSkinfitloop, Getakinfitloop, gkinfitloop, diffvectorkinfitloop); //lambdavectorloop.Print();

	    lagravaluekinfitloop = lambdavectorloop * (Getakinfitloop * diffvectorkinfitloop + gkinfitloop); //cout << lagravaluekinfitloop << endl;
	    
	    lagravaluekinfitloop_temp = lagravaluekinfitloop;
	    lagravaluekinfitloop = lagravaluekinfittildeloop;
	    
	    etakinfitloop = etakinfitfunc(Vmatrixloop, Getakinfittransloop, lambdavectorloop, inputvect, Row);
	    etakinfitloop_temp = etakinfitloop; //etakinfitloop_temp.Print();
	    etakinfitloop = etatildekinfitloop;
	    
	    // statistics
	    Vkinfitloop = Vmatrixloop - Vmatrixloop * (Getakinfittransloop * InSkinfitloop * Getakinfitloop) * Vmatrixloop; //Vkinfitloop.Print();
	    sigma2vectorkinfitloop = Fillsigma2vectorkinfit(Row, Vkinfitloop, Col, Row); //sigma2vectorkinfitloop.Print();
	    pullsvectorloop = Fillpullsvector(Row, sigma2vector, sigma2vectorkinfitloop, inputvect, etakinfitloop_temp); //pullsvectorloop.Print();
	    
	    step ++;
	  }

	  lagvalue_min_7C = lagravaluekinfitloop_temp;
	  inputvect_7C = inputvect;
	  sigma2vect_7C = sigma2vector;
	  etakinfit_min_7C = etakinfitloop_temp; 
	  pvalue=TMath::Prob(lagvalue_min_7C, 7);
	  sigma2vectorkinfit_min_7C = sigma2vectorkinfitloop;
	  pullkinfit = pullsvectorloop;
	  
	  
	  //cout << "lagvalue_min_7C = " <<  lagvalue_min_7C << "\n";
	  //etakinfit_min_7C.Print();
	}

      }

    }

    	  
    //cout << pullkinfit(0) << ", " << pullkinfit(5) << ", " << pullkinfit(
    /// check kin.fit 7c output
    //cout << "sigma2vect_7C info:" << "\n"; sigma2vect_7C.Print();

    /// fill reconstructed final state particles: pi+, pi-, pi0, pi0 photons, isr photon candidate

    // fill pi+ and pi-
    //TLorentzVector TLVector_pmi = Gettrack4vectorkinfit(curpca[trkindx1], cotpca[trkindx1], phipca[trkindx1]);
    //TLorentzVector TLVector_ppl = Gettrack4vectorkinfit(curpca[trkindx2], cotpca[trkindx2], phipca[trkindx2]);
    TLorentzVector TLVector_pmi = Gettrack4vectorkinfit(trkindx1);
    TLorentzVector TLVector_ppl = Gettrack4vectorkinfit(trkindx2);

    //cout << TLVector_pmi.M() << endl;

    /*
    cout << "\nMC TRUE" << "\n"
	 << "trk + (px, py, pz) = (" << piplusMC_TLvect.Px() << ", " << piplusMC_TLvect.Py() << ", " << piplusMC_TLvect.Pz() << ")" << "\n"
	 << "trk - (px, py, pz) = (" << piminusMC_TLvect.Px() << ", " << piminusMC_TLvect.Py() << ", " << piminusMC_TLvect.Pz() << ")" << "\n"
	 << "RECONSTRUCTED" << "\n"
	 << "trk + (px, py, pz) = (" << TLVector_ppl.Px() << ", " << TLVector_ppl.Py() << ", " << TLVector_ppl.Pz() << ")" << "\n"
	 << "trk - (px, py, pz) = (" << TLVector_pmi.Px() << ", " << TLVector_pmi.Py() << ", " << TLVector_pmi.Pz() << ")" << "\n";
    */
    
    // select pi0 photons and fill all photon final states
    TLorentzVector pionphoton1_tmp, pionphoton2_tmp;
    //cout << "TLVector_pmi info:" << "\n"; TLVector_pmi.Print();
    
    TVectorD inputvect_permut(Row), inputvectfit_permut(Row);
    TVectorD sigma2vect_permut(Row);
    TVectorD inputvect_ordered(Row), inputvect_fitted_ordered(Row);
    
    double mgg_tmp = 0., mgg_tmp_checked = 0., mgg_diff_tmp = 0.;
    double error_mggsq = 0., chi2mgg_tmp = 0., chi2mgg_min = 1e14;
    double m3pi_tmp = 0., m3pi_diff_tmp = 0.;
    double error_m3pisq = 0., chi2m3pi_tmp = 0., chi2m3pi_min = 1e14;
    int nr1 = 0, nr2 = 0, nr3 = 0;
    int isrgam_indx = 0, pi0gam1_indx = 0, pi0gam2_indx = 0;
    const int error_type_indx = 0; // 1 for kin.fitted error

    //cout << "!!!! inspect inputvector after kin. fit. " << endl;
    //etakinfit_min_7C.Print();
    //cout << "!!!! inspect sigma2 inputvector after kin. fit. " << endl;
    //sigma2vectorkinfit_min_7C.Print();
	  
    for (int i = 0; i < 3; i ++) {// loop over number of photon candidates

      // determine permutation indices
      nr1 = i % 3; nr2 = (i + 1) % 3; nr3 = (i + 2) % 3;
      
      // input vectors
      inputvect_permut = Fillpermutvector(Row, inputvect_7C, nr1, nr2, nr3); //inputvect_permut.Print();
      sigma2vect_permut = Fillpermutvector(Row, sigma2vect_7C, nr1, nr2, nr3); //sigma2vect_permut.Print();

      //inputvect_permut = Fillpermutvector(Row, etakinfit_min_7C, nr1, nr2, nr3); //inputvect_permut.Print();
      //sigma2vect_permut = Fillpermutvector(Row, sigma2vectorkinfit_min_7C, nr1, nr2, nr3); //sigma2vect_permut.Print();

      // filling pi0 photon candidates in permutation
      pionphoton1_tmp = Getphoton4vector(inputvect_permut(5), inputvect_permut(6), inputvect_permut(7), inputvect_permut(8));
      pionphoton2_tmp = Getphoton4vector(inputvect_permut(10),inputvect_permut(11),inputvect_permut(12),inputvect_permut(13));
      //pionphoton1_tmp.Print();
      //pionphoton2_tmp.Print();

      // calculate mgg and chi2mgg
      mgg_tmp = TMath::Sqrt((pionphoton1_tmp + pionphoton2_tmp).M2());
      mgg_tmp_checked = TMath::Sqrt(2. * pionphoton1_tmp.E() * pionphoton2_tmp.E() * (1. - TMath::Cos(pionphoton1_tmp.Angle(pionphoton2_tmp.Vect()))));
      mgg_diff_tmp = mgg_tmp - massneupion;
      error_mggsq = TMath::Power(mgg_tmp / 2, 2) * (sigma2vect_permut(5) / TMath::Power(inputvect_permut(5), 2) + sigma2vect_permut(10) / TMath::Power(inputvect_permut(10), 2));
      //error_mggsq = 2.855 * 2.855;
      chi2mgg_tmp = TMath::Power(mgg_diff_tmp, 2) / error_mggsq;
      
      m3pi_tmp = inputvect_permut(0);
      m3pi_diff_tmp = m3pi_tmp - 209.;
      //error_m3pisq = TMath::Power(m3pi_tmp / 2, 2) * (sigma2vect_permut(0)); //2.63 * 2.63;
      chi2m3pi_tmp = chi2mgg_tmp + TMath::Power(m3pi_diff_tmp, 2.) / sigma2vect_permut(0);
      
      //cout << "(nr1, nr2, nr3) = (" << nr1 << ", " << nr2 << ", " << nr3 << "), mgg squared = " << mgg_tmp << "+/-" << error_mggsq << ", checked = " << mgg_tmp_checked << ", mgg_diff = " << mgg_diff_tmp << ", chi2mgg_tmp = " << chi2mgg_tmp << ", error type = " << error_type_indx << endl;
      //cout << "(nr1, nr2, nr3) = (" << nr1 << ", " << nr2 << ", " << nr3 << "), m3pi squared = " << m3pi_tmp << "+/-" << error_m3pisq << ", m3pi_diff = " << m3pi_diff_tmp << ", chi2m3pi_tmp = " << chi2m3pi_tmp << endl;

      //cout << "(nr1, nr2, nr3) = (" << nr1 << ", " << nr2 << ", " << nr3 << "), chi2mgg_tmp = " << chi2mgg_tmp << endl;
      // fill kin.fit values
      inputvectfit_permut = Fillpermutvector(Row, etakinfit_min_7C, nr1, nr2, nr3);

      //if (chi2m3pi_tmp < chi2m3pi_min) {
      if (chi2mgg_tmp < chi2mgg_min) {
	chi2mgg_min = chi2mgg_tmp;
	inputvect_ordered = inputvect_permut;
	inputvect_fitted_ordered = inputvectfit_permut;
	isrgam_indx = nr1;
	pi0gam1_indx = nr2;
	pi0gam2_indx = nr3;
	//cout << "min (nr1, nr2, nr3) = (" << nr1 << ", " << nr2 << ", " << nr3 << "), chi2mgg = " << chi2mgg_tmp << "\n";
      }
	
    }// end loop

    chi2mgg_min_7C = chi2mgg_min;
    //hist1d -> Fill(chi2mgg_min_7C);

    //cout << "massomega = " << massomega << endl;
    //if (chi2mgg_min_7C > 10.) continue;

    //cout << nr1 << ", " << nr2 << ", " << nr3 << endl;

    //cout << "!!!! chi2mgg_min = " << chi2mgg_min << endl;
    
    // without kin.fit
    TLorentzVector TLVector_isrpho = Getphoton4vector(inputvect_ordered(0),inputvect_ordered(1), inputvect_ordered(2), inputvect_ordered(3));
    TLorentzVector TLVector_pi0pho1 = Getphoton4vector(inputvect_ordered(5), inputvect_ordered(6), inputvect_ordered(7), inputvect_ordered(8));
    TLorentzVector TLVector_pi0pho2 = Getphoton4vector(inputvect_ordered(10), inputvect_ordered(11), inputvect_ordered(12), inputvect_ordered(13));
    TLorentzVector TLvector_isrpho_miss = Beam - (TLVector_pi0pho1 + TLVector_pi0pho2 + TLVector_ppl + TLVector_pmi);

    IMisrpho_miss = TLvector_isrpho_miss.M2();
    
    //cout << (TLVector_pi0pho1 + TLVector_pi0pho2).M() << endl;
    
    // masses
    MASSLIST[0] = (TLVector_pi0pho1 + TLVector_pi0pho2).M(); // pi0 mass no fit
    MASSLIST[1] = (TLVector_pi0pho1 + TLVector_pi0pho2 + TLVector_ppl + TLVector_pmi).M(); // 3pi mass no fit
    
    /*
    cout << "mass without kin.fit" << "\n"
	 << "pi0 mass = " << MASSLIST[0] << "\n"
	 << "3pi mass = " << MASSLIST[1] << "\n";
    */
    
    // fill 7C kin.fitted values
    //inputvect_fitted_ordered.Print();
    TLorentzVector TLVector_isrpho_kinfit7C = Getphoton4vector(inputvect_fitted_ordered(0),inputvect_fitted_ordered(1), inputvect_fitted_ordered(2), inputvect_fitted_ordered(3));
    TLorentzVector TLVector_pi0pho1_kinfit7C = Getphoton4vector(inputvect_fitted_ordered(5), inputvect_fitted_ordered(6), inputvect_fitted_ordered(7), inputvect_fitted_ordered(8));
    TLorentzVector TLVector_pi0pho2_kinfit7C = Getphoton4vector(inputvect_fitted_ordered(10), inputvect_fitted_ordered(11), inputvect_fitted_ordered(12), inputvect_fitted_ordered(13));
    TLorentzVector TLVector_pi0gg12_kinfit7C = TLVector_pi0pho1_kinfit7C + TLVector_pi0pho2_kinfit7C;
    TLorentzVector TLVector_3pi_kinfit7C = TLVector_pi0pho1_kinfit7C + TLVector_pi0pho2_kinfit7C + TLVector_ppl + TLVector_pmi;
    TLorentzVector TLVector_2pi_kinfit7C = TLVector_ppl + TLVector_pmi;

    IM3pi_7C = (TLVector_pi0pho1_kinfit7C + TLVector_pi0pho2_kinfit7C + TLVector_ppl + TLVector_pmi).M();
    IM_3pi_nofit = (TLVector_pi0pho1 + TLVector_pi0pho2 + TLVector_ppl + TLVector_pmi).M();

    

    //if (TLVector_pi0pho1_kinfit7C.M2() < 0. || TLVector_pi0pho2_kinfit7C.M2() < 0. || TLVector_isrpho_kinfit7C.M2() < 0.) continue;
    //cout << "TLVector_pi0pho1_kinfit7C.M2() = " << TLVector_pi0pho1_kinfit7C.M2() << endl;
    
    //TLorentzVector TLVector_ppl_boost = TLVector_ppl;
    //TLVector_ppl_boost.Boost(Boost3vector);

    //TLorentzVector TLVector_pmi_boost = TLVector_pmi;
    //TLVector_pmi_boost.Boost(Boost3vector);

    //double pplmomsqr_boost = (TLVector_ppl_boost.Vect()).Mag2();
    //double pmimomsqr_boost = (TLVector_pmi_boost.Vect()).Mag2();

    // boost 3pi system
    TVector3 Boost3vector_3pi = -TLVector_3pi_kinfit7C.BoostVector();
    TVector3 Boost3vector_2pi = -TLVector_2pi_kinfit7C.BoostVector();
    
    TLorentzVector TLVector_ppl_boost = TLVector_ppl;
    TLorentzVector TLVector_pmi_boost = TLVector_pmi;
    TLVector_ppl_boost.Boost(Boost3vector_2pi);
    TLVector_pmi_boost.Boost(Boost3vector_2pi);
    
    angle_ppl_3piboost = TLVector_ppl_boost.Angle(TLVector_pi0gg12_kinfit7C.Vect())*TMath::RadToDeg(); // angle between pi+ in the rest 3pi system with repect to isr photon cadidate
    //angle_ppl_3piboost = TLVector_ppl_boost.Angle(TLVector_pmi_boost.Vect())*TMath::RadToDeg();
    //cout << angle_ppl_3piboost << endl;
    
    betapi0 = (TLVector_pi0gg12_kinfit7C.Vect()).Mag() / TLVector_pi0gg12_kinfit7C.E();
    //cout << "betapi0 = " << betapi0 << endl;

    // get IM3pi12, IM3pi13 and IM3pi23, with kin.fit photons:1, 2 and 3 in energy order E1 > E2 > E3
    // define 4-mom for pho1, pho2 and pho3
    TVector3 pho1_vect, pho2_vect, pho3_vect; 
    TLorentzVector TLVector_pho1, TLVector_pho2, TLVector_pho3;
    TLorentzVector TLVector_3pi12, TLVector_pho13, TLVector_pho23;

    // list of 4_mom of all three kin.fit photons
    TLorentzVector TLVector_pho_list[3] = {TLVector_isrpho_kinfit7C, TLVector_pi0pho1_kinfit7C, TLVector_pi0pho2_kinfit7C};
    
    double pho1E = 0., pho3E = 1E6, pho2E = 0.;

    for (int i = 0; i < 3; i ++) {// find min and max

      if (TLVector_pho_list[i].E() > pho1E) {

	pho1E = TLVector_pho_list[i].E();
	
	TLVector_pho1 = TLVector_pho_list[i];
	
      }

      if (TLVector_pho_list[i].E() < pho3E) {

	pho3E = TLVector_pho_list[i].E();
	
	TLVector_pho3 = TLVector_pho_list[i];
	
      }

      //cout << "phoE_list[" << i << "] = " << TLVector_pho_list[i].E() << ", pho1E = " << TLVector_pho1.E()  << ", pho3E = "<< TLVector_pho3.E() << "\n\n";

      
    }

    for (int i = 0; i < 3; i ++) {//find middle

      if (TLVector_pho_list[i].E() != pho1E && TLVector_pho_list[i].E() != pho3E ) {

	TLVector_pho2 = TLVector_pho_list[i];

	//cout << "pho2E = " << TLVector_pho2.E() << endl;
	
      }
      
      
    }

    IM3pi_pi12 = (TLVector_pho1 + TLVector_pho2 + TLVector_ppl + TLVector_pmi).M();
    IM3pi_pi13 = (TLVector_pho1 + TLVector_pho3 + TLVector_ppl + TLVector_pmi).M();
    IM3pi_pi23 = (TLVector_pho2 + TLVector_pho3 + TLVector_ppl + TLVector_pmi).M();

    /*
    cout << "pho1E = " << TLVector_pho1.E() << ", pho2E = " << TLVector_pho2.E() << ", pho3E = "<< TLVector_pho3.E() << "\n"
	 << "IM3pi_pi12 = " << IM3pi_pi12 << ", IM3pi_pi13 = " << IM3pi_pi13 << ", IM3pi_pi23 = " << IM3pi_pi23 << ", IM3pi_7C = " << IM3pi_7C << "\n\n";
    */
    
    //
    TVector3 pi0gam1_vect, pi0gam2_vect, pi0gam3_vect; // selected pi0 photon, 3pi mc true four vectors
    TLorentzVector TLVector_pi0pho1_true, TLVector_pi0pho2_true, TLVector_isrpho3_true;
    TLorentzVector TLVector_3pi_true;

    pi0gam1_vect = TVector3(0., 0., 0.);
    pi0gam2_vect = TVector3(0., 0., 0.);
    pi0gam3_vect = TVector3(0., 0., 0.);
    TLVector_pi0pho1_true = TLorentzVector(0,0,0,0);
    TLVector_pi0pho2_true = TLorentzVector(0,0,0,0);
    TLVector_isrpho3_true = TLorentzVector(0,0,0,0);
    TLVector_3pi_true = TLorentzVector(0,0,0,0);
    
    if (kineid != -999) {
      pi0gam1_vect.SetXYZ(pxmc[pho_indx[pi0gam1_indx]], pymc[pho_indx[pi0gam1_indx]], pzmc[pho_indx[pi0gam1_indx]]);
      pi0gam2_vect.SetXYZ(pxmc[pho_indx[pi0gam2_indx]], pymc[pho_indx[pi0gam2_indx]], pzmc[pho_indx[pi0gam2_indx]]);
      pi0gam3_vect.SetXYZ(pxmc[pho_indx[isrgam_indx]], pymc[pho_indx[isrgam_indx]], pzmc[pho_indx[isrgam_indx]]);

    }

    TLVector_pi0pho1_true = GetLorentzVector(pi0gam1_vect, 0.);
    TLVector_pi0pho2_true = GetLorentzVector(pi0gam2_vect, 0.);
    TLVector_isrpho3_true = GetLorentzVector(pi0gam3_vect, 0.);
    
    TLVector_3pi_true = piminusMC_TLvect + piplusMC_TLvect + TLVector_pi0pho1_true + TLVector_pi0pho2_true;

    IM3pi_true = TLVector_3pi_true.M();
    
    
    // energy
    ENERGYLIST[0] = TLVector_isrpho_kinfit7C.E();
    ENERGYLIST[1] = TLVector_pi0pho1_kinfit7C.E(); 
    ENERGYLIST[2] = DeltaE(TLVector_ppl, TLVector_pmi);//(TLVector_ppl_boost.Vect() + TLVector_pmi_boost.Vect()).Mag() - (Beam.M() - TMath::Sqrt(masschpion * masschpion + pplmomsqr_boost) - TMath::Sqrt(masschpion * masschpion + pmimomsqr_boost));
    ENERGYLIST[3] = TLVector_pi0pho2_kinfit7C.E(); 
    ENERGYLIST[4] = TLVector_isrpho.E();
    //cout << "!!!!!!!!!!!!!!!!! " << ENERGYLIST[3] << endl;

    // photons 4 vectors
    // pi0 photon1
    MOM4PHO1[0] = TLVector_pi0pho1_kinfit7C.E();
    MOM4PHO1[1] = TLVector_pi0pho1_kinfit7C.X();
    MOM4PHO1[2] = TLVector_pi0pho1_kinfit7C.Y();
    MOM4PHO1[3] = TLVector_pi0pho1_kinfit7C.Z();

    // pi0 photon2
    MOM4PHO2[0] = TLVector_pi0pho2_kinfit7C.E();
    MOM4PHO2[1] = TLVector_pi0pho2_kinfit7C.X();
    MOM4PHO2[2] = TLVector_pi0pho2_kinfit7C.Y();
    MOM4PHO2[3] = TLVector_pi0pho2_kinfit7C.Z();

    // pi0 photon3
    MOM4PHO3[0] = TLVector_isrpho_kinfit7C.E();
    MOM4PHO3[1] = TLVector_isrpho_kinfit7C.X();
    MOM4PHO3[2] = TLVector_isrpho_kinfit7C.Y();
    MOM4PHO3[3] = TLVector_isrpho_kinfit7C.Z();
    
    // track 4 vectors
    MOM4TRKPLUS[0] = TLVector_ppl.E();
    MOM4TRKPLUS[1] = TLVector_ppl.X();
    MOM4TRKPLUS[2] = TLVector_ppl.Y();
    MOM4TRKPLUS[3] = TLVector_ppl.Z();

    MOM4TRKMINS[0] = TLVector_pmi.E();
    MOM4TRKMINS[1] = TLVector_pmi.X();
    MOM4TRKMINS[2] = TLVector_pmi.Y();
    MOM4TRKMINS[3] = TLVector_pmi.Z();

    /*
    cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
	 << "ppl (E, px, py, pz) = (" << MOM4TRKPLUS[0] << ", " << MOM4TRKPLUS[1] << ", " << MOM4TRKPLUS[2] << ", " << MOM4TRKPLUS[3] << ")\n";
    */

    
    // masses
    trkmass = Trkmass(TLVector_ppl, TLVector_pmi);
    ppIM = (TLVector_ppl + TLVector_pmi).M();
    //betapi0 = (TLVector_pi0gg12_kinfit7C.Vect()).Mag() / TLVector_pi0gg12_kinfit7C.E();
    ppIM_beta = ((TLVector_ppl + TLVector_pmi).Vect()).Mag() / (TLVector_ppl + TLVector_pmi).E();
    //(TLVector_pi0pho1_kinfit7C + TLVector_pi0pho2_kinfit7C + TLVector_ppl + TLVector_pmi).M();
    IM3pi_beta = ((TLVector_pi0pho1_kinfit7C + TLVector_pi0pho2_kinfit7C + TLVector_ppl + TLVector_pmi).Vect()).Mag() / (TLVector_pi0pho1_kinfit7C + TLVector_pi0pho2_kinfit7C + TLVector_ppl + TLVector_pmi).E();
    
    IM_pi0_7C = TLVector_pi0gg12_kinfit7C.M();
    IM_pi0_nofit = (TLVector_pi0pho1 + TLVector_pi0pho2).M();
    MASSLIST[2] = IM_pi0_7C; //TLVector_pi0gg12_kinfit7C.M(); // pi0 mass
    MASSLIST[3] = IM3pi_7C; //(TLVector_pi0pho1_kinfit7C + TLVector_pi0pho2_kinfit7C + TLVector_ppl + TLVector_pmi).M(); // 3pi mass 
    MASSLIST[4] = trkmass; //Trkmass(TLVector_ppl, TLVector_pmi); // track mass
    MASSLIST[5] = ppIM; //(TLVector_ppl + TLVector_pmi).M(); // pi+ pi- invaraint mass
    MASSLIST[6] = (TLVector_pi0pho1_kinfit7C + TLVector_isrpho_kinfit7C).M(); // pi0 mass, gamma 13
    MASSLIST[7] = (TLVector_pi0pho2_kinfit7C + TLVector_isrpho_kinfit7C).M(); // pi0 mass, gamma 23
    MASSLIST[8] = (TLVector_pi0pho1_kinfit7C + TLVector_isrpho_kinfit7C + TLVector_ppl + TLVector_pmi).M(); // 3pi mass, gamma13 
    MASSLIST[9] = (TLVector_pi0pho2_kinfit7C + TLVector_isrpho_kinfit7C + TLVector_ppl + TLVector_pmi).M(); // 3pi mass, gamma23 
    MASSLIST[10] = (TLVector_ppl + TLVector_pmi).M2() * 1e-6; //invariant mass square of pi+ pi- in GeV^{2}
    MASSLIST[11] = (TLVector_ppl + TLVector_pi0gg12_kinfit7C).M2() * 1e-6; //invariant mass square of pi+ pi0 in GeV^{2}
    MASSLIST[12] = IMisrpho_miss; // invariant mass of recoil mass against final hadron system 
    MASSLIST[13] = ppIM_beta; // ppIM beta 
    MASSLIST[14] = IM3pi_beta; // IM3pi_beta
    //cout << "isrpho missing mass = " << MASSLIST[12] << endl;
    
    // angle
    ANGLELIST[0] = TLVector_pi0pho1_kinfit7C.Angle(TLVector_pi0pho2_kinfit7C.Vect())*TMath::RadToDeg(); // gamma12
    ANGLELIST[1] = TLVector_pi0pho2_kinfit7C.Angle(TLVector_isrpho_kinfit7C.Vect())*TMath::RadToDeg(); // gamma23
    ANGLELIST[2] = TMath::Cos(TLVector_isrpho_kinfit7C.Theta()); // cos isr
    ANGLELIST[3] = TLVector_pi0pho1_kinfit7C.Angle(TLVector_isrpho_kinfit7C.Vect())*TMath::RadToDeg(); // gamma13
    ANGLELIST[4] = TLVector_pmi.Angle(TLVector_ppl.Vect())*TMath::RadToDeg(); // pi+ pi-
    ANGLELIST[5] = angle_ppl_3piboost;
    /*
    cout << "pho1: px = " << TLVector_pi0pho1_kinfit7C.X() << ", py = " << TLVector_pi0pho1_kinfit7C.Y() << ", pz = " << TLVector_pi0pho1_kinfit7C.Z() << ", E = " << TLVector_pi0pho1_kinfit7C.E() << "\n"
	 << "pho3: px = " << TLVector_isrpho_kinfit7C.X() << ", py = " << TLVector_isrpho_kinfit7C.Y() << ", pz = " << TLVector_isrpho_kinfit7C.Z() << ", E = " << TLVector_isrpho_kinfit7C.E() << "\n"
	 << "angle = " << ANGLELIST[3] << "\n\n";
    */
    
    //cout << ANGLELIST[3] << endl;
    //cout << "ANGLELIST[0] = " << ANGLELIST[0] << endl;
    //cout << "ANGLELIST[2] = " << ANGLELIST[2] << endl;
    
    /*
    cout << "mass with 7C kin.fit" << "\n"
	 << "pi0 mass = " << MASSLIST[2] << "\n"
	 << "3pi mass = " << MASSLIST[3] << "\n"
	 << "trk mass = " << MASSLIST[4] << "\n";
    */
    
    // resolutions
    RESOLIST[0] = TLVector_pi0gg12_kinfit7C.M() - massneupion; // Mpi
    RESOLIST[1] = IM_3pi; // M3pi

    //cout << "RESOLIST[0] = " << RESOLIST[0] << "\n";

    // pull
    PULLIST[0] = pullkinfit(0);
    PULLIST[1] = pullkinfit(1);
    PULLIST[2] = pullkinfit(2);
    PULLIST[3] = pullkinfit(3);
    PULLIST[4] = pullkinfit(4);
    PULLIST[5] = pullkinfit(5);
    PULLIST[6] = pullkinfit(6);
    PULLIST[7] = pullkinfit(7);
    PULLIST[8] = pullkinfit(8);
    PULLIST[9] = pullkinfit(9);
    PULLIST[10] = pullkinfit(10);
    PULLIST[11] = pullkinfit(11);
    PULLIST[12] = pullkinfit(12);
    PULLIST[13] = pullkinfit(13);
    PULLIST[14] = pullkinfit(14);
    
    //cout << Row << pullkinfit(0) << endl;
    test_value = ENERGYLIST[2];//DeltaE(TLVector_ppl_boost, TLVector_pmi_boost);

    //cout << "trkmass = " << ANGLELIST[0] << ", true = " << angle_pi0gam12_true << endl;

    ///
    /*
    TLorentzVector TVect_clust_test1; 
    if (IM_3pi < 820. && IM_3pi > 760. && IM3pi_7C < 750. && IM3pi_7C > 500.  ) {
    //if (E_radiv2 == 0.) {
      hsmearmatr -> Fill(IM_3pi, IM3pi_7C);
      //htestI -> Fill(Angle_pho_isr);
      htestI -> Fill(ppIM_true);
      
      cout << "typeI" << "\n"
	   << "Esum = " << Esum << "\n"
	   << "ppIM = " << ppIM_true << "\n"
	   << "Angle_pho_isr = " << Angle_pho_isr << "\n"
	   << "E_radiv1 = " << E_radiv1 << ", E_radiv2 = " << E_radiv2 << "\n";

      for (int i = 0; i < ntmc; i ++) {
	cout << "ntmc=" << i << ", pidmc=" << pidmc[i] << ", pxmc=" << pxmc[i] << ", pymc=" << pymc[i] << ", pzmc=" << pzmc[i] << "\n";
      }
      cout << "RECON\n";
      //Getphoton4vector(double E, double x, double y, double z)
      for (int i = 0; i < nclu; i ++) {
	TVect_clust_test1 = Getphoton4vector(enecl[i], xcl[i] - bx, ycl[i] - by, zcl[i] - bz); //SetPxPyPzE(xcl[i]-bx, ycl[i]-by, zcl[i]-bz, enecl[i]);
	//cout << xcl[i] - bx << ", " << ycl[i] - by << ", " <<  zcl[i] - bz << ", " << TVect_clust_test.M() << endl;
	cout << "nclu=" << i << ", px=" << TVect_clust_test1.X() << ", py=" << TVect_clust_test1.Y() << ", pz=" << TVect_clust_test1.Z() << ", E=" << TVect_clust_test1.E() << "\n";
      }
      
    }
    else {

      //htestII -> Fill(Angle_pho_isr);
      htestII -> Fill(ppIM_true);

      cout << "typeII" << "\n"
	   << "Esum = " << Esum << "\n"
	   << "ppIM = " << ppIM_true << "\n"
	   << "Angle_pho_isr = " << Angle_pho_isr << "\n"
	   << "E_radiv1 = " << E_radiv1 << ", E_radiv2 = " << E_radiv2 << "\n";

      for (int i = 0; i < ntmc; i ++) {
	cout << "ntmc=" << i << ", pidmc=" << pidmc[i] << ", pxmc=" << pxmc[i] << ", pymc=" << pymc[i] << ", pzmc=" << pzmc[i] << "\n";
      }

      cout << "RECON\n";
      //Getphoton4vector(double E, double x, double y, double z)
      for (int i = 0; i < nclu; i ++) {
	TVect_clust_test1 = Getphoton4vector(enecl[i], xcl[i] - bx, ycl[i] - by, zcl[i] - bz); //SetPxPyPzE(xcl[i]-bx, ycl[i]-by, zcl[i]-bz, enecl[i]);
	//cout << xcl[i] - bx << ", " << ycl[i] - by << ", " <<  zcl[i] - bz << ", " << TVect_clust_test.M() << endl;
	cout << "nclu=" << i << ", px=" << TVect_clust_test1.X() << ", py=" << TVect_clust_test1.Y() << ", pz=" << TVect_clust_test1.Z() << ", E=" << TVect_clust_test1.E() << "\n";
      }
      
    }
    cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n";
    */
    
    /// fill trees
    if (lagvalue_min_7C > 100.) continue; // cut 6
    //if (IM_3pi > 820. || IM_3pi < 750.) continue;
    //if (IM3pi_7C > 730. || IM3pi_7C < 440.) continue;
    evnt_final ++;
    //hsmearmatr_trk -> Fill(ppIM_true, ppIM);
    //hsmearmatr_clust -> Fill(IM_pi0, IM_pi0_7C);
    //hsmearmatr -> Fill(IM_3pi, MASSLIST[1]);
    
    //cout << "!!!!!!!!!!!!!1" << IM_3pi << ", " << IM3pi_7C << endl;
    //cout << "nclu = " << nclu << ", nclumc = " << nclumc << "\n"
    //	 << "mcflag = " << mcflag << "\n"
    //	 << "nt = " << nt << ", ntmc = " << ntmc << ", ntfmc = " << ntfmc << "\n\n";

    //IM3pi_7C = MASSLIST[3];

    //cout << "SELECTED PHOTONS\n";
    double clu_Emin = 1e13;
    int clu_indx_Emin = 0;
    int bkg_counter = 1;
    for (int i = 0; i < 3; i ++) {
      if (pho_indx[i] == -1) {
	bkg_counter = bkg_counter + 1;
      }
      
      if (E_list[i] < clu_Emin) {
	clu_Emin = E_list[i];
	clu_indx_Emin = pho_indx[i];
      }
    }
    
    //for (int i = 0; i < nclumc; i ++) {
    //if (Pnum1[i] - 1 == -1) bkg_indx ++; 
    //}

    //int clust_mult = 1; // cluster multiplicity
    if (pho_indx[0] == pho_indx[1] || pho_indx[0] == pho_indx[2] || pho_indx[1] == pho_indx[2]) {
      bkg_counter = bkg_counter + 1;
      //clust_mult ++;
    }

    //bkg_indx = clust_mult;
    bkg_indx = bkg_counter;

    //if (bkg_indx == 0) {
      //hIM3pi_bkg_indx0 -> Fill(IM3pi_true);
    //}
    //else {
      //hIM3pi_bkg_indx_rest -> Fill(IM3pi_true);
    //}

    //if (bkg_indx == 0) continue;
    
    //if (phid == 5) {
    //cout << "\nMC TRUE " << jentry << "\n";
    //cout << IM3pi_7C << endl;
    //if (IM3pi_7C > 750.) {
    double E_ntmc;
    TString pid_type = "";
    for (int i = 0; i < ntmc; i ++) {
      if (pidmc[i] == 1) {
	pid_type = "gamma";
      }
      else if (pidmc[i] == 2) {
	pid_type = "e+";
      }
      else if (pidmc[i] == 3) {
	pid_type = "e-";
      }
      else if (pidmc[i] == 5) {
	pid_type = "mu+";
      }
      else if (pidmc[i] == 6) {
	pid_type = "mu-";
      }
      else if (pidmc[i] == 7) {
	pid_type = "pi0";
      }
      else if (pidmc[i] == 8) {
	pid_type = "pi+";
      }
      else if (pidmc[i] == 9) {
	pid_type = "pi-";
      }
      else {
	pid_type = "unknown";
      }
      
      E_ntmc = TMath::Sqrt(pxmc[i] * pxmc[i] + pymc[i] * pymc[i] + pzmc[i] * pzmc[i]); 
      //cout << "ntmc=" << i << ", pidmc(pid_type)=" << pidmc[i] << "(" << pid_type << "), pxmc=" << pxmc[i] << ", pymc=" << pymc[i] << ", pzmc=" << pzmc[i] << ", E= " << E_ntmc << "\n";
    }

    //cout << "\nTRUE\n";
    int Pnum_sum = 0;
    for (int i = 0; i < nclumc; i ++) {
      //if (i < 3 && Pnum1[i] - 1 == -1) {
      //if (i == 0 && Pnum1[i] - 1 == -1) {
      if (Pnum1[i] - 1 == -1) {
	Pnum_sum ++;
	//cout << Pnum1[i] - 1 << endl;
      }
      if (Pid1[i] == 1) {
	pid_type = "gamma";
      }
      else {
	pid_type = "others";
      }
      
      //cout << "nclumc="<< i << ", Npar=" << Npar[i] << ", Pnum1=" << Pnum1[i] - 1 << ", Pid1(pid_type)=" << Pid1[i] << "(" << pid_type << "), pxmc=" << pxmc[Pnum1[i]-1] << ", pymc=" << pymc[Pnum1[i]-1] << ", pzmc=" << pzmc[Pnum1[i]-1] << ", E= " << TMath::Sqrt(pxmc[Pnum1[i]-1]*pxmc[Pnum1[i]-1]+pymc[Pnum1[i]-1]*pymc[Pnum1[i]-1]+pzmc[Pnum1[i]-1]*pzmc[Pnum1[i]-1]) << "\n";
      //Pnum2="<< Pnum2[i] - 1 << ", Pid2=" << Pid2[i] << ", Pnum3=" << Pnum3[i] - 1 << ", Pid3=" << Pid3[i] << "\n";
    }

    //cout << "\nRECON\n";
    TLorentzVector TVect_clust_test; //Getphoton4vector(double E, double x, double y, double z)
    for (int i = 0; i < nclu; i ++) {
      TVect_clust_test = Getphoton4vector(enecl[i], xcl[i] - bx, ycl[i] - by, zcl[i] - bz); //SetPxPyPzE(xcl[i]-bx, ycl[i]-by, zcl[i]-bz, enecl[i]);
      //cout << xcl[i] - bx << ", " << ycl[i] - by << ", " <<  zcl[i] - bz << ", " << TVect_clust_test.M() << endl;
      
      //cout << "nclu=" << i << ", px=" << TVect_clust_test.X() << ", py=" << TVect_clust_test.Y() << ", pz=" << TVect_clust_test.Z() << ", E=" << TVect_clust_test.E() << "\n";

    }
    
    
    //cout << "Pnum_sum = " << Pnum_sum << endl;
    
    //if (bkg_counter > 0) {
    //hist1d -> Fill(IM3pi_7C);
    //ALLCHAIN_TEST.Fill();
    
    //cout << "E input: (" << etakinfit_min_7C(0) << ", " << etakinfit_min_7C(5) << ", " << etakinfit_min_7C(10) << ")\n";
    //etakinfit_min_7C.Print();
    //inputvect_7C.Print();
    //cout << "indx ordered: (" << isrgam_indx << ", " << pi0gam1_indx << ", " << pi0gam2_indx << ")\n";
    //cout << "E input ordered: (" << inputvect_ordered(0) << ", " << inputvect_ordered(5) << ", " << inputvect_ordered(10) << ")\n";
    //inputvect_ordered.Print();

    //
    //cout << "\nCheck pi0 photons" << endl;
    int recon_indx_tmp = 0;

    /*
    cout << "pi0gam1_indx = " << pi0gam1_indx << ", pi0gam2_indx = " << pi0gam2_indx << "\n"
	 << "E_true_list[pi0gam1_indx] = " << E_true_list[pi0gam1_indx] << ", (px, py, pz) = (" << pxmc[pho_indx[pi0gam1_indx]] << ", " << pymc[pho_indx[pi0gam1_indx]] << ", " << pzmc[pho_indx[pi0gam1_indx]] << "), pho_indx = " << pho_indx[pi0gam1_indx] << "\n"
	 << "E_true_list[pi0gam2_indx] = " << E_true_list[pi0gam2_indx] << ", (px, py, pz) = (" << pxmc[pho_indx[pi0gam2_indx]] << ", " << pymc[pho_indx[pi0gam2_indx]] << ", " << pzmc[pho_indx[pi0gam2_indx]] << "), pho_indx = " << pho_indx[pi0gam2_indx] << "\n";
    */
    
    //cout << "TLVector_pimins_true: " << endl; piminusMC_TLvect.Print();
    //cout << "TLVector_piplus_true: " << endl; piplusMC_TLvect.Print();

    PI0PHORESD[0] = TLVector_pi0pho1_true.E() - TLVector_pi0pho1.E();
    PI0PHORESD[1] = TLVector_pi0pho2_true.E() - TLVector_pi0pho2.E();
    PI0PHORESD[2] = TLVector_pi0pho1_true.E() - TLVector_pi0pho1_kinfit7C.E();
    PI0PHORESD[3] = TLVector_pi0pho2_true.E() - TLVector_pi0pho2_kinfit7C.E();

    /*
    cout<< "\n2 prompt photons seleted as pi0 photons:\n"
	<< "pi0pho1_true: (px, py, pz, E) = (" << TLVector_pi0pho1_true.X() << ", " << TLVector_pi0pho1_true.Y() << ", " << TLVector_pi0pho1_true.Z() << ", " << TLVector_pi0pho1_true.E() << ")" << "\n"
	 << "\trecon: (px, py, pz, E) = (" << TLVector_pi0pho1.X() << ", " << TLVector_pi0pho1.Y() << ", " << TLVector_pi0pho1.Z() << ", " << TLVector_pi0pho1.E() << ")" << "\n";
    cout << "pi0pho2_true: (px, py, pz, E) = (" << TLVector_pi0pho2_true.X() << ", " << TLVector_pi0pho2_true.Y() << ", " << TLVector_pi0pho2_true.Z() << ", " << TLVector_pi0pho2_true.E() << ")" << endl; //TLVector_pi0pho2_true.Print();
    */
    
    //cout << "TLVector_3pi_true: (px, py, pz, E) = (" << TLVector_3pi_true.X() << ", " << TLVector_3pi_true.Y() << ", " << TLVector_3pi_true.Z() << ", " << TLVector_3pi_true.E() << ")" << endl; TLVector_3pi_true.Print();

    //cout << "\npiminusMC_TLvect: " << endl; piminusMC_TLvect.Print();
    //cout << "piplusMC_TLvect: " << endl; piplusMC_TLvect.Print();
    //cout << "\nMC pi0pho1_true: (px, py, pz, E) = (" << pi0gam1_TLvect.X() << ", " << pi0gam1_TLvect.Y() << ", " << pi0gam1_TLvect.Z() << ", " << pi0gam1_TLvect.E() << ")" << endl; //pi0gam1_TLvect.Print();
    //cout << "MC pi0pho2_true: (px, py, pz, E) = (" << pi0gam2_TLvect.X() << ", " << pi0gam2_TLvect.Y() << ", " << pi0gam2_TLvect.Z() << ", " << pi0gam2_TLvect.E() << ")" << endl; //pi0gam2_TLvect.Print();

    //cout << "3pi_true.M = " << TLVector_3pi_true.M() << ", IM_3pi (checked) = " << IM_3pi << "(" << (pi0gam1_TLvect + pi0gam2_TLvect + piminusMC_TLvect + piplusMC_TLvect).M() << ") \n\n";
    
    //cout << "\nthreepi_TLvect: " << endl; threepi_TLvect.Print();

    //cout << "NTMC_SEL pi0_pho1 = " << pho_indx[pi0gam1_indx] << ", pi0_pho2 = " << pho_indx[pi0gam2_indx] << "\n";

    /*
    for (int i = 0; i < 2; i ++) {
      
      for (int j = 0; j < 2; j ++) {

	if (j >= i) {

	  cout << "(i, j) = (" << i << ", " << j << ")\n";
	      
	}
	
      }
      
    }
    */

    /*
    cout << "pi0gam indx: " << pi0gam1_indx << ", " << pi0gam2_indx << "\n"
	 << "EPI0NTMC, pho_indx 1: " << EPI0NTMC[0] << ", " << pho_indx[pi0gam1_indx] << "\n"
	 << "EPI0NTMC, pho_indx 2: " << EPI0NTMC[1] << ", " << pho_indx[pi0gam2_indx] << endl;
    */
   
    Bool_t checked[2] = {kFALSE, kFALSE};

    for (int i = 0; i < 2; i ++) {// find the first match

      //cout << "EPI0GAM[" << i << "] = " << EPI0GAM[i] << ", EPI0NTMC[" << i << "] = "<< EPI0NTMC[i] << endl;

      //if (TMath::Abs(EPI0GAM[i] - E_true_list[pi0gam1_indx]) < 1e-3) {
      if (pho_indx[pi0gam1_indx] == EPI0NTMC[i]) {
	recon_indx_tmp ++;
	checked[i] = kTRUE;
	//cout << E_true_list[pi0gam1_indx] << ", pi0gam1_indx = " << pi0gam1_indx << ", recon1_indx = " << recon_indx_tmp << ", first mathched at EPI0NTMC=" << EPI0NTMC[i] << ", pho_indx[pi0gam1_indx]=" << pho_indx[pi0gam1_indx] << "\n";
      }
     
      /*
      //if (TMath::Abs(EPI0GAM[i] - E_true_list[pi0gam2_indx]) < 1e-3) {
      if (pho_indx[pi0gam2_indx] == EPI0NTMC[i]) {
	recon_indx_tmp ++;
	checked[i] = kTRUE;
	cout << E_true_list[pi0gam2_indx] << ", pi0gam2_indx = " << pi0gam2_indx << ", recon2_indx = " << recon_indx_tmp << ", mathched " << checked[i] << " at " << i << "\n";
      }
      */

    }

    for (int i = 0; i < 2; i ++) {// find the second match

      //cout << checked[i] << endl;
	
      if (pho_indx[pi0gam2_indx] == EPI0NTMC[i] && checked[i]==kFALSE) {
	recon_indx_tmp ++;
	checked[i] = kTRUE;
	//cout << E_true_list[pi0gam2_indx] << ", pi0gam2_indx = " << pi0gam2_indx << ", recon2_indx = " << recon_indx_tmp << ", first mathched at EPI0NTMC=" << EPI0NTMC[i] << ", pho_indx[pi0gam2_indx]=" << pho_indx[pi0gam2_indx] << "\n";
      }
      
    }
    
    recon_indx = recon_indx_tmp;
    //cout << "Matched recon pi0 photons = " << recon_indx << endl;

    // photon 1 4-mom: true pi0 photon
    pho_E1 = TLVector_pi0pho1_kinfit7C.E();
    pho_px1 = TLVector_pi0pho1_kinfit7C.X();
    pho_py1 = TLVector_pi0pho1_kinfit7C.Y();
    pho_pz1 = TLVector_pi0pho1_kinfit7C.Z();

    pho_E1_true = TLVector_pi0pho1_true.E();
    pho_px1_true = TLVector_pi0pho1_true.X();
    pho_py1_true = TLVector_pi0pho1_true.Y();
    pho_pz1_true = TLVector_pi0pho1_true.Z();
    
    // photon 2 4-mom: true pi0 photon
    pho_E2 = TLVector_pi0pho2_kinfit7C.E();
    pho_px2 = TLVector_pi0pho2_kinfit7C.X();
    pho_py2 = TLVector_pi0pho2_kinfit7C.Y();
    pho_pz2 = TLVector_pi0pho2_kinfit7C.Z();

    pho_E2_true = TLVector_pi0pho2_true.E();
    pho_px2_true = TLVector_pi0pho2_true.X();
    pho_py2_true = TLVector_pi0pho2_true.Y();
    pho_pz2_true = TLVector_pi0pho2_true.Z();
    
    // photon 3 4-mom: unpaired photon
    pho_E3 = TLVector_isrpho_kinfit7C.E();
    pho_px3 = TLVector_isrpho_kinfit7C.X();
    pho_py3 = TLVector_isrpho_kinfit7C.Y();
    pho_pz3 = TLVector_isrpho_kinfit7C.Z();

    pho_E3_true = TLVector_isrpho3_true.E();
    pho_px3_true = TLVector_isrpho3_true.X();
    pho_py3_true = TLVector_isrpho3_true.Y();
    pho_pz3_true = TLVector_isrpho3_true.Z();
    
    // positive track
    ppl_E = TLVector_ppl.E();
    ppl_px = TLVector_ppl.X();
    ppl_py = TLVector_ppl.Y();
    ppl_pz = TLVector_ppl.Z();

    // negative track
    pmi_E = TLVector_pmi.E();
    pmi_px = TLVector_pmi.X();
    pmi_py = TLVector_pmi.Y();
    pmi_pz = TLVector_pmi.Z();

    if (recon_indx == 0) {
      evnt_recon_type0 ++;
    }
    else if (recon_indx == 1) {
      evnt_recon_type1 ++;
    }
    else {
      evnt_recon_type2 ++;
          
    }

    if (recon_indx == 2 && bkg_indx == 1) {
      //cout << recon_indx << ", " << bkg_indx << endl;
      hmpi0_iden -> Fill((TLVector_pi0pho2 + TLVector_pi0pho1).M()); 
    }
    //if (recon_indx == 1 || recon_indx == 0) continue; 

    
    //if (recon_indx == 2) {
      //hsmearmatr_typeI -> Fill(IM3pi_true, IM3pi_7C);
      //hrecon_typeI -> Fill(IM3pi_true);
      //hIM3pi -> Fill(IM_3pi);
    //}
    //else {
      //hsmearmatr_typeII -> Fill(IM3pi_true, IM3pi_7C);
    //}

    //hsmearmatr_clust -> Fill(IM3pi_true, IM3pi_7C);
    herror_type -> Fill(recon_indx, bkg_indx);

    
    //if (clu_indx_Emin == -1){

    /*
    cout << "\nSummary\n"
	 << "jentry = " << jentry << "\n"
	 << "indx pid: (" << pid_indx[0] << ", " << pid_indx[1] << ", " << pid_indx[2] << "), ntmc: (" << pho_indx[0] << ", " << pho_indx[1] << ", " << pho_indx[2] << "), cluster: (" << clu_indx[0] << ", " << clu_indx[1] << ", " << clu_indx[2] << ") \n"
	 << "energy true: (" << E_true_list[0] << ", " << E_true_list[1] << ", " << E_true_list[2] << "), rec: (" << E_list[0] << ", " << E_list[1] << ", " << E_list[2] << ") \n"
      //<< "energy min: E_min = " << clu_Emin << ", ntmc_indx = " << clu_indx_Emin << "\n"
      //<< "bkg_indx = " << bkg_indx << "\n"
      // << "pi0gam1 E = " << E_pi0gam1 << ", pi0gam2 E = " << E_pi0gam2 << "\n"
      //<< EPI0GAM[0] << ", " << EPI0GAM[1] << "\n"
      
	 << "selected pi0gam indx = (" << pi0gam1_indx << ", " << pi0gam2_indx << "), energy = (" << E_true_list[pi0gam1_indx] << ", " << E_true_list[pi0gam2_indx] << ") \n"
	 << "bkg_indx = " << bkg_indx << "\n\n";
    */
      //<< "IM3pi_7C = " << IM3pi_7C << "\n\n";
    //}

    
    //if (recon_indx == 2) {
    //cout << jentry << endl;
    //break;
    //}
    
    
    //inputvect_ordered.Print();
    //}
    
    //}
    //}
    //cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n" << endl;
    
    
    /*
      if (phid == 5) {
      //cout << IM3pi_7C << ", Pnum_sum = " << Pnum_sum << endl;
      if (IM3pi_7C > 750.) {
      cout << "combinatorial" << ", " << Pnum_sum << endl;
      for (int i = 0; i < nclumc; i ++) {
      cout << "nclumc="<< i << ", Npar=" << Npar[i] << ", Pnum1=" << Pnum1[i] - 1 << ", Pid1=" << Pid1[i] << ", pxmc=" << pxmc[Pnum1[i]-1] << ", pymc=" << pymc[Pnum1[i]-1] << ", pzmc=" << pzmc[Pnum1[i]-1] << ", p2_sqrt= " << TMath::Sqrt(pxmc[Pnum1[i]-1]*pxmc[Pnum1[i]-1]+pymc[Pnum2[i]-1]*pymc[Pnum2[i]-1]+pzmc[Pnum1[i]-1]*pzmc[Pnum1[i]-1]) << ", Pnum2="<< Pnum2[i] - 1 << ", Pid2=" << Pid2[i] << ", Pnum3=" << Pnum3[i] - 1 << ", Pid3=" << Pid3[i] << "\n";
      }
      }
      else {
      //cout << "eta peak" << ", " << Pnum_sum << endl;
      for (int i = 0; i < nclumc; i ++) {
      cout << "nclumc="<< i << ", Npar=" << Npar[i] << ", Pnum1=" << Pnum1[i] - 1 << ", Pid1=" << Pid1[i] << ", pxmc=" << pxmc[Pnum1[i]-1] << ", pymc=" << pymc[Pnum1[i]-1] << ", pzmc=" << pzmc[Pnum1[i]-1] << ", p2_sqrt= " << TMath::Sqrt(pxmc[Pnum1[i]-1]*pxmc[Pnum1[i]-1]+pymc[Pnum2[i]-1]*pymc[Pnum2[i]-1]+pzmc[Pnum1[i]-1]*pzmc[Pnum1[i]-1]) << ", Pnum2="<< Pnum2[i] - 1 << ", Pid2=" << Pid2[i] << ", Pnum3=" << Pnum3[i] - 1 << ", Pid3=" << Pid3[i] << "\n";
      }
      }
      cout << "\n" << endl;
      
      if (Pnum_sum > 0){
      //cout << Pnum_sum << endl;
      hist1d -> Fill(IM3pi_7C);
      }
      //IM3pi_7C = MASSLIST[3];
      
      }
    */
    
    ALLCHAIN_CUT.Fill();
  }
  
  // save trees
  //TrSample.Write();
  //ALLCHAIN_TEST.Write();
  ALLCHAIN_GEN.Write();
  ALLCHAIN_STR2.Write();
  ALLCHAIN_CUT.Write();

  // summary
  cout << "SUMMARY" << "\n"
       << "evnt_sum = " << evnt_sum << ": evnt_sig = " << evnt_sig << ", evnt_bkg = " << evnt_bkg << "\n"
       << "evnt_sum = " << evnt_sum << ": evnt_pre_filtred = " << evnt_pre_filtred << ", passed = " << evnt_sum - evnt_pre_filtred << "\n"
       << "etagam evnts: evnt_phi5 = " << evnt_phi5 << "\n"
 
       << "Pre-selected = " << evnt_pre << ", streamed = " << evnt_cls << ", unstr = "<< evnt_unstr << "\n"
       << "evnt_trig = " << evnt_trig << ", trigger_indx = " << trigger_indx << "\n"
       << "evnt_filfo = " << evnt_filfo << ", filfo_indx = " << filfo_indx << "\n"
       << "evnt str2 = " << evnt_str << ", evtcls_indx = " << evtcls_indx << "\n"
       << "" << ", fiduial_indx = " << fiduial_indx << "\n"
       << "evnt_trk = " << evnt_trk << "\n"
       << "evnt_class = " << evnt_class << "\n"
       << "evnt_final = " << evnt_final << "\n";

  cout << "test_indx = " << test_indx << endl;

  cout << "find 0 correct pi0 photons = " << evnt_recon_type0 << "\n"
       << "find 1 correct pi0 photons = " << evnt_recon_type1 << "\n"
       << "find 2 correct pi0 photons = " << evnt_recon_type2 << "\n";

  cout << "PARMETERS" << "\n"
       << "egammamin = " << egammamin << "\n"
       << "Zvmax = " << Zvmax << "\n"
       << "Rhovmax = " << Rhovmax << "\n"
       << "nb_sigma_T_clust = " << nb_sigma_T_clust << "\n";
       
    
  //printf("speed of light c");
}

// funcitons
TLorentzVector MyClass::GetLorentzVector(TVector3 vector, double mass) {

  TLorentzVector tvector(0.,0.,0.,0.);
  double E = TMath::Sqrt(TMath::Power(mass, 2) + vector.Mag2());
  tvector.SetPxPyPzE(vector(0), vector(1), vector(2), E);
  //cout << tvector.M() << endl;
  
  return tvector;
}

bool MyClass::IfStreamed(int pstrnb) {
  Bool_t bTagged = kFALSE, chrad = kFALSE;

  for (int NrEC = 0; NrEC < necls; NrEC ++){

    if (eclstream[NrEC] == pstrnb) {

      bTagged = kTRUE;

      if ((ecltagnum[NrEC] &  1 ) == 1) chrad = kTRUE;

    }

  }

  //if (bTagged) return kTRUE;
  if (bTagged) return kTRUE;
  else return kFALSE;
}

Bool_t MyClass::IfTriggered() {
  Bool_t passtrigger=kFALSE;

  if (trgtype == 6 || trgtype == 2 || trgtype == 4) passtrigger = kTRUE;
  //if (trgtype == 6) passtrigger = kTRUE;

  if (passtrigger) return kFALSE;

  else return kTRUE;
}

Bool_t MyClass::IfTriggered_EMC() {
  Bool_t passtrigger=kFALSE;

  if (trgtype == 2 || trgtype == 6) passtrigger = kTRUE;

  if (passtrigger) return kFALSE;

  else return kTRUE;
}

Bool_t MyClass::IfTriggered_DC() {
  Bool_t passtrigger=kFALSE;

  if (trgtype == 4 || trgtype == 6) passtrigger = kTRUE;

  if (passtrigger) return kFALSE;

  else return kTRUE;
}

Bool_t MyClass::IfTriggered_EMCandDC() {
  Bool_t passtrigger=kFALSE;

  if (trgtype == 6) passtrigger = kTRUE;

  if (passtrigger) return kFALSE;

  else return kTRUE;
}

  

Bool_t MyClass::IfFilfoed() {
  Bool_t passfilfo=kFALSE;

  if (((eclfilfo & ( 1 << 20 )) >> 20) == 1) passfilfo=kTRUE;

  if (passfilfo) return kFALSE;

  else return kTRUE;

}

Bool_t MyClass::IfFilfoed_28() {
  Bool_t passfilfo=kFALSE;

  if (((eclfilfo & ( 1 << 28 )) >> 28) == 1) passfilfo=kTRUE;

  if (passfilfo) return kFALSE;

  else return kTRUE;

}


TVectorD MyClass::Getpiontrnb() {
  TVectorD indexvector(3);
  int svar_pos = 0, svar_neg = 0, index1 = 0, index2 = 0;
  Float_t distance1 = 1000000.0, distance2 = 100000.0;
  Float_t dist = 0.0;

  for (int k = 0; k < nt; k++) {// begin loop
    dist = sqrt((bx-xpca[k])*(bx-xpca[k]) + (by-ypca[k])*(by-ypca[k]) +(bz-zpca[k])*(bz-zpca[k]));
    if (curpca[k] < 0 && dist < distance1) {
      distance1 = dist;
      index1 = k;
      svar_neg=1;
    }
    else if (curpca[k] > 0 && dist < distance2) {
      distance2 = dist;
      index2 = k;
      svar_pos=1;
    }
  }// end loop

  indexvector(0) = svar_neg+svar_pos;
  indexvector(1) = index1;
  indexvector(2) = index2; //indexvector.Print();

  return indexvector;
}

bool MyClass::IfBroken(int idx1, int idx2) {
  if (TMath::Abs(cot[idx1] + cot[idx2]) < 0.1 && TMath::Abs((cur[idx1] + cur[idx2]) / cur[idx1]) < 0.2) {
    //cout << TMath::Abs(cot[idx1] + cot[idx2]) << endl;
    return kTRUE;
  }
  else return kFALSE;
}

TVectorD MyClass::FillInputvector_7C(int size, int index1, int index2, int index3, int indtr1, int indtr2) {
  TVectorD vector(size);
  double inputarray[15] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double r1sq=0., r2sq=0., r3sq=0., pp1sq=0., pp2sq=0.;
  double r1=0., r2=0., r3=0., px1=0., px2=0., px3=0., py1=0., py2=0., py3=0., pz1=0., pz2=0., pz3=0.;
  double Epp1=0., Epp2=0., pp1=0., pp2=0., ppx1=0., ppx2=0., ppy1=0., ppy2=0., ppz1=0., ppz2=0., pp1tran=0., pp2tran=0.;
  double x1=0., x2=0., x3=0., y1=0., y2=0., y3=0., z1=0., z2=0., z3=0., E1=0., E2=0., E3=0., t1=0., t2=0., t3=0.;

  // calulate x1, x2, x3, y1, y2, y3, z1, z2, z3
  // calulate E1, E2, E3, t1, t2, t3
  // calulate px1, px2, px3
  x1=X_list[index1], x2=X_list[index2], x3=X_list[index3];
  y1=Y_list[index1], y2=Y_list[index2], y3=Y_list[index3];
  z1=Z_list[index1], z2=Z_list[index2], z3=Z_list[index3];

  E1=E_list[index1], E2=E_list[index2], E3=E_list[index3];
  t1=T_list[index1], t2=T_list[index2], t3=T_list[index3];

  // measurements array
  inputarray[0]=E1, inputarray[1]=x1, inputarray[2]=y1, inputarray[3]=z1, inputarray[4]=t1; //cout<<inputarray[4]<<endl;
  inputarray[5]=E2, inputarray[6]=x2, inputarray[7]=y2, inputarray[8]=z2, inputarray[9]=t2;
  inputarray[10]=E3, inputarray[11]=x3, inputarray[12]=y3, inputarray[13]=z3, inputarray[14]=t3;
  return vector.Use(size,inputarray);
}

TVectorD MyClass::FillSigma2vector_7C(int size, int index1, int index2, int index3, int indtr1, int indtr2) {
  TVectorD vector(size);
  double sigma2array[15]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  // fill sigma2 array
  sigma2array[0]=Sigma2_E_list[index1], sigma2array[1]=Sigma2_X_list[index1], sigma2array[2]=Sigma2_Y_list[index1], sigma2array[3]=Sigma2_Z_list[index1], sigma2array[4]=Sigma2_T_list[index1];
  sigma2array[5]=Sigma2_E_list[index2], sigma2array[6]=Sigma2_X_list[index2], sigma2array[7]=Sigma2_Y_list[index2], sigma2array[8]=Sigma2_Z_list[index2], sigma2array[9]=Sigma2_T_list[index2];
  sigma2array[10]=Sigma2_E_list[index3], sigma2array[11]=Sigma2_X_list[index3], sigma2array[12]=Sigma2_Y_list[index3], sigma2array[13]=Sigma2_Z_list[index3], sigma2array[14]=Sigma2_T_list[index3];

  //sigma2array[15]=sigcurv[indtr1], sigma2array[16]=sigcot[indtr1], sigma2array[17]=sigphi[indtr1];
  //sigma2array[18]=sigcurv[indtr2], sigma2array[19]=sigcot[indtr2], sigma2array[20]=sigphi[indtr2];

  return vector.Use(size,sigma2array);
}

TVectorD MyClass::Gfunc_7C(TVectorD etavector, int rownb, TLorentzVector beam, int indtr1, int indtr2) {
  TVectorD vector(rownb);
  //double Beam_E=benergy;
  //beam.Print();

  double r1sq=0., r2sq=0., r3sq=0., pp1sq=0., pp2sq=0.;
  double r1=0., r2=0., r3=0., px1=0., px2=0., px3=0., py1=0., py2=0., py3=0., pz1=0., pz2=0., pz3=0.;
  double Epp1=0., Epp2=0., pp1=0., pp2=0., ppx1=0., ppx2=0., ppy1=0., ppy2=0., ppz1=0., ppz2=0., pp1tran=0., pp2tran=0.;
  double curv1=0., cotan1=0., phi1=0, curv2=0., cotan2=0., phi2=0;
  double x1=0., x2=0., x3=0., y1=0., y2=0., y3=0., z1=0., z2=0., z3=0., E1=0., E2=0., E3=0., t1=0., t2=0., t3=0.;
  double garray[7] = {0.,0.,0.,0.,0.,0.,0.};
  double costheta = 0., massfactor = 0.;

  E1=etavector(0), x1=etavector(1), y1=etavector(2), z1=etavector(3), t1=etavector(4);
  E2=etavector(5), x2=etavector(6), y2=etavector(7), z2=etavector(8), t2=etavector(9);
  E3=etavector(10), x3=etavector(11), y3=etavector(12), z3=etavector(13), t3=etavector(14);

  // calculate squared r1sq, r2sq, r3sq and calculate r1, r2, r3
  r1sq=TMath::Power(x1,2.)+TMath::Power(y1,2.)+TMath::Power(z1,2.);
  r2sq=TMath::Power(x2,2.)+TMath::Power(y2,2.)+TMath::Power(z2,2.);
  r3sq=TMath::Power(x3,2.)+TMath::Power(y3,2.)+TMath::Power(z3,2.);

  r1=TMath::Sqrt(r1sq), r2=TMath::Sqrt(r2sq), r3=TMath::Sqrt(r3sq);

  //
  costheta = (x1*x2+y1*y2+z1*z2)/(r1*r2);
  //massfactor = (massneupion*massneupion)/(2*(1-costheta));
  //std::cout<<"E1 = "<<E1<<", E2 = "<<etavector(5)<<std::endl;
  //std::cout<<"costheta = "<<costheta<<", massfactor = "<<massfactor<<std::endl;
  //std::cout<<"constraint E2 = "<<massfactor/E1<<std::endl;

  //E2=massfactor/E1;


  px1=E1*x1/r1, px2=E2*x2/r2, px3=E3*x3/r3;
  py1=E1*y1/r1, py2=E2*y2/r2, py3=E3*y3/r3;
  pz1=E1*z1/r1, pz2=E2*z2/r2, pz3=E3*z3/r3;

  px1=E1*x1/r1, px2=E2*x2/r2, px3=E3*x3/r3;
  py1=E1*y1/r1, py2=E2*y2/r2, py3=E3*y3/r3;
  pz1=E1*z1/r1, pz2=E2*z2/r2, pz3=E3*z3/r3;

  // calculate squared ppx1sq, ppx2s, pp1tran, pp2tran, pp1, pp2, Epp1, Epp2
  // ppi1.SetX(ptransvSmeared * TMath::Cos(phipca[index]));
  //ppi1.SetY(ptransvSmeared * TMath::Sin(phipca[index]));
  //ppi1.SetZ(ptransvSmeared * cotpca[index]);

  curv1=curpca[indtr1], cotan1=cotpca[indtr1], phi1=phipca[indtr1];
  curv2=curpca[indtr2], cotan2=cotpca[indtr2], phi2=phipca[indtr2];
  //cout << "curv1 = " << curv1 << ", cotan1 = " << cotan1 << endl;
  //cout << "curv2 = " << curv2 << ", cotan2 = " << cotan2 << endl;
  //cout<<cotan2-cotpca[indtr2]<<endl;
  //cout<<phi2-phipca[indtr2]<<endl;
  pp1tran=1000.*1./TMath::Abs(curv1);
  pp2tran=1000.*1./TMath::Abs(curv2);

  ppx1=pp1tran*TMath::Cos(phi1);
  ppx2=pp2tran*TMath::Cos(phi2);

  ppy1=pp1tran*TMath::Sin(phi1);
  ppy2=pp2tran*TMath::Sin(phi2);

  ppz1=pp1tran*cotan1;
  ppz2=pp2tran*cotan2;

  pp1sq=TMath::Power(ppx1,2.)+TMath::Power(ppy1,2.)+TMath::Power(ppz1,2.);
  pp2sq=TMath::Power(ppx2,2.)+TMath::Power(ppy2,2.)+TMath::Power(ppz2,2.);

  Epp1=TMath::Sqrt(masschpion*masschpion+pp1sq);
  Epp2=TMath::Sqrt(masschpion*masschpion+pp2sq);

  pp1=TMath::Sqrt(pp1sq);
  pp2=TMath::Sqrt(pp2sq);

  TLorentzVector tvector1(0.,0.,0.,0.), tvector2(0.,0.,0.,0.);
  tvector1.SetPxPyPzE(px1,py1,pz1,E1);
  tvector2.SetPxPyPzE(px2,py2,pz2,E2);

  double openangle_loop=tvector1.Angle(tvector2.Vect());
  double cosangle_loop=TMath::Cos(openangle_loop);

  garray[0]=E1+E2+E3+Epp2+Epp1-beam.E();
  //cout<<"!!!!!!!!!!!!!!Beam E="<<Beam.X()<<"!!!!!!!!!!"<<endl;
  garray[1]=px1+px2+px3+ppx1+ppx2-beam.X();
  garray[2]=py1+py2+py3+ppy1+ppy2-beam.Y();
  garray[3]=pz1+pz2+pz3+ppz1+ppz2-beam.Z();
  garray[4]=r1-speedc*t1;
  garray[5]=r2-speedc*t2; //cout<<garray[5]<<endl;
  garray[6]=r3-speedc*t3;

  //cout << "E1 = " << E1 << ", E2 = " << E2 << ", E3 = " << E3 << endl;
  //cout << "Epp2 = " << Epp2 << ", Epp1 = " << Epp1 << ", beamE = " << beam.E() << endl;
  //cout << "g_new" << endl;
  //cout << "g1 = " << garray[0] << ", g2 = " << garray[1] << ", g3 = " << garray[2] << ", g4 = " << garray[3] << ", g5 = " << garray[4] << ", g6 = " << garray[5] << ", g7 = " << garray[6] << endl;
  //cout << "r1 = " << r1 << ", t1 = " << t1 << endl;

  vector.Use(rownb,garray);

  return vector;
}

TMatrixD MyClass::Getafunc_7C(TVectorD etavector, int rownb, int colnb) {
  TMatrixD matrix_temp(rownb, colnb);
  //cout << rownb << endl;
  double dg1deta[15] = {1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.};
  double dg2deta[15] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double dg3deta[15] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double dg4deta[15] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double dg5deta[15] = {0.,0.,0.,0.,-speedc,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  double dg6deta[15] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,-speedc,0.,0.,0.,0.,0.};
  double dg7deta[15] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-speedc};
  double r1sq=0., r2sq=0., r3sq=0., pp1sq=0., pp2sq=0.;
  double r1=0., r2=0., r3=0., px1=0., px2=0., px3=0., py1=0., py2=0., py3=0., pz1=0., pz2=0., pz3=0.;
  double Epp1=0., Epp2=0., pp1=0., pp2=0., ppx1=0., ppx2=0., ppy1=0., ppy2=0., ppz1=0., ppz2=0., pp1tran=0., pp2tran=0.;
  double curv1=0., cotan1=0., phi1=0, curv2=0., cotan2=0., phi2=0;
  double x1=0., x2=0., x3=0., y1=0., y2=0., y3=0., z1=0., z2=0., z3=0., E1=0., E2=0., E3=0., t1=0., t2=0., t3=0.;
  TVectorD dg1(colnb), dg2(colnb), dg3(colnb), dg4(colnb), dg5(colnb), dg6(colnb), dg7(colnb);
  TVectorD dgarray[7]={dg1, dg2, dg3, dg4, dg5, dg6, dg7};
  double hatE1=0., hatE2=0., hatE3=0, hatx1=0., hatx2=0., hatx3=0., haty1=0., haty2=0., haty3=0., hatz1=0., hatz2=0., hatz3=0., costheta = 0., massfactor = 0.;


  // calulate x1, x2, x3, y1, y2, y3, z1, z2, z3
  // calulate E1, E2, E3, t1, t2, t3
  // calulate px1, px2, px3
  // calculate squared r1sq, r2sq, r3sq and calculate r1, r2, r3
  x1=etavector(1), x2=etavector(6), x3=etavector(11);
  y1=etavector(2), y2=etavector(7), y3=etavector(12);
  z1=etavector(3), z2=etavector(8), z3=etavector(13);

  t1=etavector(4), t2=etavector(9), t3=etavector(14);

  r1sq=TMath::Power(x1,2.)+TMath::Power(y1,2.)+TMath::Power(z1,2.);
  r2sq=TMath::Power(x2,2.)+TMath::Power(y2,2.)+TMath::Power(z2,2.);
  r3sq=TMath::Power(x3,2.)+TMath::Power(y3,2.)+TMath::Power(z3,2.);

  r1=TMath::Sqrt(r1sq), r2=TMath::Sqrt(r2sq), r3=TMath::Sqrt(r3sq);

  E1=etavector(0);
  E2=etavector(5);
  E3=etavector(10);
  //
  costheta = (x1*x2+y1*y2+z1*z2)/(r1*r2);
  massfactor = (massneupion*massneupion)/(2*(1-costheta));
  //std::cout<<"E1 = "<<E1<<", E2 = "<<etavector(5)<<std::endl;
  //std::cout<<"costheta = "<<costheta<<", massfactor = "<<massfactor<<std::endl;
  //std::cout<<"constraint E2 = "<<massfactor/E1<<std::endl;

  //E2=massfactor/E1;


  px1=E1*(x1)/r1, px2=E2*(x2)/r2, px3=E3*(x3)/r3;
  py1=E1*(y1)/r1, py2=E2*(y2)/r2, py3=E3*(y3)/r3;
  pz1=E1*(z1)/r1, pz2=E2*(z2)/r2, pz3=E3*(z3)/r3;


  // calculate more variables
  hatE1=E1/r1, hatE2=E2/r2, hatE3=E3/r3, hatx1=(x1)/r1, hatx2=(x2)/r2, hatx3=(x3)/r3, haty1=(y1)/r1, haty2=(y2)/r2, haty3=(y3)/r3, hatz1=(z1)/r1, hatz2=(z2)/r2, hatz3=(z3)/r3;


  // fill dgdeta array

  dg2deta[0]=hatx1, dg2deta[1]=hatE1*(TMath::Power(haty1,2)+TMath::Power(hatz1,2)), dg2deta[2]=-hatE1*hatx1*haty1, dg2deta[3]=-hatE1*hatz1*hatx1;
  dg2deta[5]=hatx2, dg2deta[6]=hatE2*(TMath::Power(haty2,2)+TMath::Power(hatz2,2)), dg2deta[7]=-hatE2*hatx2*haty2, dg2deta[8]=-hatE2*hatz2*hatx2;
  dg2deta[10]=hatx3, dg2deta[11]=hatE3*(TMath::Power(haty3,2)+TMath::Power(hatz3,2)), dg2deta[12]=-hatE3*hatx3*haty3, dg2deta[13]=-hatE3*hatz3*hatx3;


  dg3deta[0]=haty1, dg3deta[1]=dg2deta[2], dg3deta[2]=hatE1*(TMath::Power(hatx1,2)+TMath::Power(hatz1,2)), dg3deta[3]=-hatE1*hatz1*haty1;
  dg3deta[5]=haty2, dg3deta[6]=dg2deta[7], dg3deta[7]=hatE2*(TMath::Power(hatx2,2)+TMath::Power(hatz2,2)), dg3deta[8]=-hatE2*hatz2*haty2;
  dg3deta[10]=haty3, dg3deta[11]=dg2deta[12], dg3deta[12]=hatE3*(TMath::Power(hatx3,2)+TMath::Power(hatz3,2)), dg3deta[13]=-hatE3*hatz3*haty3;


  dg4deta[0]=hatz1, dg4deta[1]=dg2deta[3], dg4deta[2]=dg3deta[3], dg4deta[3]=hatE1*(TMath::Power(hatx1,2)+TMath::Power(haty1,2));
  dg4deta[5]=hatz2, dg4deta[6]=dg2deta[8], dg4deta[7]=dg3deta[8], dg4deta[8]=hatE2*(TMath::Power(hatx2,2)+TMath::Power(haty2,2));
  dg4deta[10]=hatz3, dg4deta[11]=dg2deta[13], dg4deta[12]=dg3deta[13], dg4deta[13]=hatE3*(TMath::Power(hatx3,2)+TMath::Power(haty3,2));


  dg5deta[1]=hatx1, dg5deta[2]=haty1, dg5deta[3]=hatz1;

  dg6deta[6]=hatx2, dg6deta[7]=haty2, dg6deta[8]=hatz2;

  dg7deta[11]=hatx3, dg7deta[12]=haty3, dg7deta[13]=hatz3;


  // fill dgarray
  dgarray[0]=dg1.Use(15,dg1deta);
  dgarray[1]=dg2.Use(15,dg2deta);
  dgarray[2]=dg3.Use(15,dg3deta);
  dgarray[3]=dg4.Use(15,dg4deta); //cout << "error" << endl;
  dgarray[4]=dg5.Use(15,dg5deta);
  dgarray[5]=dg6.Use(15,dg6deta);
  dgarray[6]=dg7.Use(15,dg7deta);

  // fill G matrix
  TMatrixDRow(matrix_temp,0)=dg1, TMatrixDRow(matrix_temp,1)=dg2, TMatrixDRow(matrix_temp,2)=dg3; TMatrixDRow(matrix_temp,3)=dg4, TMatrixDRow(matrix_temp,4)=dg5, TMatrixDRow(matrix_temp,5)=dg6, TMatrixDRow(matrix_temp,6)=dg7;


  return matrix_temp;
}

TMatrixD MyClass::CovMatrix(TVectorD vector, int row, int col) {
  TMatrixD covmatrix_temp(row,col);

  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      if (i==j) covmatrix_temp(i,j)=vector(i);
    }
  }
  return covmatrix_temp;
}

TMatrixD MyClass::Mtrans(TMatrixD ma) {
  TMatrixD ma_temp = ma;
  return ma_temp.T();
}

TMatrixD MyClass::MInvert(TMatrixD ma) {
  TMatrixD ma_temp = ma;
  return ma_temp.Invert();
}

TVectorD MyClass::Lambdavector(int size, TMatrixD inveSmatrix, TMatrixD dgmatrix, TVectorD constr, TVectorD diffv) {
  TVectorD vector(size), rvector(size); //rvector.Print();
  rvector=constr+dgmatrix*diffv; //rvector.Print();
  vector=inveSmatrix*rvector;

  return vector;

}

TVectorD MyClass::etakinfitfunc(TMatrixD Vmatrix, TMatrixD dgmatrixTrans, TVectorD lambdavector, TVectorD etatilde, int rownb) {
  TVectorD vector(rownb);
  vector=etatilde-(Vmatrix*dgmatrixTrans)*lambdavector;
  return vector;
}

TVectorD MyClass::Fillsigma2vectorkinfit(int size, TMatrixD matrix, int matrixcol, int matrixrow) {
  TVectorD vector(size);

  for (int i = 0; i < matrixrow; i++) {
    for (int j = 0; j< matrixcol; j++) {
      if (i==j) vector[i]=matrix[i][j];
    }
  }

  return vector;
}

TVectorD MyClass::Fillpullsvector(int size, TVectorD sigma2vector_old, TVectorD sigma2vector_new, TVectorD inputvector_old, TVectorD inputvector_new) {
  double Array[21] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  TVectorD vector(size); //vector.Print();

  for (int i=0; i<size; i++) {
    //if((sigma2vector_old-sigma2vector_new)(i)<=0) printf("<=0!!!=%lf\n",(sigma2vector_old-sigma2vector_new)(i));
    Array[i]=((inputvector_old-inputvector_new)(i))/TMath::Sqrt((sigma2vector_old-sigma2vector_new)(i));
    //if (TMath::Abs(vector(i)) < 0.001) cout<<vector(i)<<endl;
    //vector(i)=(inputvector_old-inputvector_new)(i);
  }
  vector.Use(size,Array);
  return vector;
}

TVectorD MyClass::GetTVectorD(int size) {
  TVectorD vector(size);

  return vector;

}
  
//TLorentzVector MyClass::Gettrack4vectorkinfit(double curv, double cotan, double phi) {
TLorentzVector MyClass::Gettrack4vectorkinfit(int index) {

  TLorentzVector tvector_smeared(0.,0.,0.,0.);
  TLorentzVector tvector1(0.,0.,0.,0.);
  double E=0., psq=0., ptran=0., px=0., py=0., pz=0.;
  Double_t smearing_factor = 0.1501;
  Double_t scale_factor = 0.000367;
  double curv = curpca[index];
  double phi = phipca[index];
  double cotan = cotpca[index];
  
  ptran = 1000. * 1. / TMath::Abs(curv);
  Double_t ptransvSmeared = ptran*(1.+scale_factor);
  //ptransvSmeared = generator->Gaus(ptransvSmeared,(TMath::Sqrt(sigcurv[index])/1000.)*smearing_factor*ptransvSmeared*ptransvSmeared);
  ptransvSmeared = generator->Gaus(ptransvSmeared,(TMath::Sqrt(sigcurv[index]))*smearing_factor*ptransvSmeared);
  
  /*
  cout << "scale_factor = " << scale_factor << ", smearing_factor = " << smearing_factor << "\n"
       << "ptran = " << ptran << ", ptran*(1.+scale_factor) = " << ptran*(1.+scale_factor) << ", ptransvSmeared = " << ptransvSmeared << "\n"
       << "TMath::Sqrt(sigcurv[index])) = " << TMath::Sqrt(sigcurv[index]) << ", smearing = " << smearing_factor << "\n";
  */
  
  //px=ptran*TMath::Cos(phi); py=ptran*TMath::Sin(phi), pz=ptran*cotan;
  //psq=TMath::Power(ptran,2)+TMath::Power(pz,2);
  //E=TMath::Sqrt(masschpion*masschpion+psq);
  //cout<<TMath::Abs(curv)<<endl;
  //std::cout<<"E = "<<E<<endl;

  //tvector.SetPxPyPzE(px,py,pz,E);

  tvector1.SetX(ptran * TMath::Cos(phi));
  tvector1.SetY(ptran * TMath::Sin(phi));
  tvector1.SetZ(ptran * cotan);
  tvector1.SetE(sqrt(masschpion*masschpion + tvector1.Vect().Mag2()));

  tvector_smeared.SetX(ptransvSmeared * TMath::Cos(phi));
  tvector_smeared.SetY(ptransvSmeared * TMath::Sin(phi));
  tvector_smeared.SetZ(ptransvSmeared * cotan);
  tvector_smeared.SetE(sqrt(masschpion*masschpion + tvector_smeared.Vect().Mag2()));
  
  //std::cout << "check ptran = " << (tvector1.Vect()).Mag() << endl;
  //std::cout<<"===================="<<endl;
  //cout << "ptran = " << ptran << endl;

  /*
    cout << "no smearing (px, py, pz, E) = (" << tvector1.X() << "," << tvector1.Y() << ", " << tvector1.Z() << ", " << tvector1.E() << ")" << "\n"
    << "smeared (px, py, pz, E) = (" << tvector_smeared.X() << "," << tvector_smeared.Y() << ", " << tvector_smeared.Z() << ", " << tvector_smeared.E() << ")" << "\n\n";
  */    


  //return tvector1;
  return tvector_smeared;
}

TVectorD MyClass::Fillpermutvector(int size, TVectorD input, int index1, int index2, int index3) {
  TVectorD vector(size);
  TVectorD testvector(size);

  for (int i=0;i<15;i++) {
    testvector(i)=i;
  }
  //input.Print();
  //index1=1, index2=0, index3=2;
  //printf("%d,%d,%d \n",5*index1,5*index2,5*index3);
  //cout<<testvector(index1)<<endl;
  vector(0)=input(5*index1), vector(1)=input(5*index1+1), vector(2)=input(5*index1+2), vector(3)=input(5*index1+3), vector(4)=input(5*index1+4);
  vector(5)=input(5*index2), vector(6)=input(5*index2+1), vector(7)=input(5*index2+2), vector(8)=input(5*index2+3), vector(9)=input(5*index2+4);
  vector(10)=input(5*index3), vector(11)=input(5*index3+1), vector(12)=input(5*index3+2), vector(13)=input(5*index3+3), vector(14)=input(5*index3+4);


  //testvector.Print();
  //input.Print();
  //vector.Print();
  //printf("=========================\n");

  return vector;
}

TLorentzVector MyClass::Getphoton4vector(double E, double x, double y, double z) {
  //given a cluster index returns the 4-mom of a photon
  TVector3 gamma(x-bx,y-by,z-bz);
  Double_t scale1;
  scale1=E/gamma.Mag();
  TLorentzVector gamma4mom(scale1*gamma, E);
  //cout << gamma4mom.M() << endl;
  return gamma4mom;

}

double MyClass::DeltaE(TLorentzVector bestppl, TLorentzVector bestpmi) {
  double deltaE = 0.;
  
  TVector3 Boost3vector = -Beam.BoostVector();
    
  TLorentzVector bestppl_boost = bestppl;
  bestppl_boost.Boost(Boost3vector);

  TLorentzVector bestpmi_boost = bestpmi;
  bestpmi_boost.Boost(Boost3vector);

  double pplmomsqr_boost = (bestppl_boost.Vect()).Mag2();
  double pmimomsqr_boost = (bestpmi_boost.Vect()).Mag2();

  deltaE = (bestppl_boost.Vect() + bestpmi_boost.Vect()).Mag() - (Beam.M() - TMath::Sqrt(masschpion * masschpion + pplmomsqr_boost) - TMath::Sqrt(masschpion * masschpion + pmimomsqr_boost));

  //cout << "pplmomsqr_boost = " << pplmomsqr_boost << "\n"; 
  
  return deltaE;
}

double MyClass::Trkmass(TLorentzVector bestppl, TLorentzVector bestpmi) {
  TVector3 ppl3mom, pmi3mom, p3mom_miss;
  double ppl3mom_mag = 0., pmi3mom_mag = 0., p3mom_mag_miss = 0.;
  double sqrtS = 0., mtrk2 = 0., mtrk = 0.;
  double sqr2_denom = 0., sqr2_nom = 0.;

  sqrtS = Beam.E();
  //cout << "Beam.E = " << Beam.E() << ", Beam.rho = " << Beam.Rho() << endl;
  //Beam.Print();
  //cout << "Beam trans = " << Beam.Px() << Beam. endl;

  ppl3mom = bestppl.Vect();
  ppl3mom_mag = ppl3mom.Mag();

  pmi3mom = bestpmi.Vect();
  pmi3mom_mag = pmi3mom.Mag();

  p3mom_miss = Beam.Vect() - ppl3mom - pmi3mom;
  p3mom_mag_miss = p3mom_miss.Mag();

  //cout << "pmi3mom_mag = " << pmi3mom_mag << ", bestpmi.rho = " << bestpmi.Rho() << endl;

  sqr2_denom = TMath::Power((ppl3mom_mag * ppl3mom_mag - pmi3mom_mag * pmi3mom_mag + (sqrtS - p3mom_mag_miss) * (sqrtS - p3mom_mag_miss)),2);
  sqr2_nom = (p3mom_mag_miss - sqrtS) * (p3mom_mag_miss - sqrtS);

  mtrk2 = sqr2_denom / (4 * sqr2_nom) - ppl3mom_mag * ppl3mom_mag;

  //if(mtrk2 < 0.){
  //mtrk2 = 0.;
    //cout << "mtrk2 = " << mtrk2 << endl;
  //}
  
  mtrk = TMath::Sqrt(mtrk2);

  //cout << "sqrt(s) = " << Beam.E() << endl;
  //bestppl.Print();
  //bestpmi.Print();
  //cout << "p+ = " << ppl3mom_mag << ", p- = " << pmi3mom_mag  << endl;
  //cout << "p3mom_mag_miss = " << p3mom_mag_miss << endl;
  //cout << "trkmass = " << sqr2_denom << ", sqr2_nom = " << sqr2_nom << endl;
  //cout << "TRKMASS2 = " << mtrk2 << ", TRKMASS = " << mtrk << endl;

  //cout << "p3mom_mag_miss = " << p3mom_miss.Mag() << endl;
  //cout << "p3mom_mag_miss check = " << TMath::Sqrt(p3mom_miss.Mag2()) << endl;

  return mtrk;
}
