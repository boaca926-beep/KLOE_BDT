#include "../header/sm_para.h"
#include "../header/path.h"
#include "../header/sfw1d.h"
#include "../header/sfw2d.txt"
#include "../header/method.h"

int sfw1d(){
  
  cout << "Tuning signal MC IM3pi ..." << endl;

  /// Get histos
  //TList *H1D_list = new TList();
  TFile *f_input = new TFile(outputHist + "hist.root");
  TFile *f_cut = new TFile(outputCut + "tree_pre.root");  
  TFile *f_output = new TFile(outputSfw1D + "sfw1d.root", "recreate");
  
  TList *HSFW1D = (TList *) f_input -> Get("HSFW1D");
  //gethist(f_input, HSFW2D, 12850);
  checkList(HSFW1D);
  
  // smeared signal IM3pi distr.
  TRandom *generator = new TRandom();
  rnd=new TRandom3();

  getObj(f_cut);

  f_cut -> GetObject("TISR3PI_SIG", TISR3PI_SIG);
 
  double m3pi_true = 0., m3pi_corr = 0.;

  TH1D * hsig = new TH1D("hsig", "", fit_bin, fit_min, fit_max);
  
  for (Int_t irow = 0; irow < TISR3PI_SIG -> GetEntries(); irow++) {// loop chain

    TISR3PI_SIG -> GetEntry(irow);

    m3pi_true = TISR3PI_SIG -> GetLeaf("Br_IM3pi_true") -> GetValue(0);

    m3pi_corr = DetectorEvent(TMath::Abs(m3pi_true));

    //cout << m3pi_corr << endl;
    
    hsig -> Fill(m3pi_corr);
    
  }

  // get data and MC histos
  TH1D *hIM3pi_eeg = (TH1D *) HSFW1D -> FindObject("h1d_IM3pi_sfw_TEEG"); //hIM3pi_eeg -> Draw();
  TH1D *hIM3pi_data = (TH1D *) HSFW1D -> FindObject("h1d_IM3pi_sfw_TDATA");
  //TH1D *hIM3pi_isr3pi = (TH1D *) HSFW1D -> FindObject(2); // reconstructed
  TH1D * hIM3pi_omegapi = (TH1D *) HSFW1D -> FindObject("h1d_IM3pi_sfw_TOMEGAPI");
  TH1D * hIM3pi_kpm = (TH1D *) HSFW1D -> FindObject("h1d_IM3pi_sfw_TKPM");
  TH1D * hIM3pi_ksl = (TH1D *) HSFW1D -> FindObject("h1d_IM3pi_sfw_TKSL");
  TH1D * hIM3pi_3pigam = (TH1D *) HSFW1D -> FindObject("h1d_IM3pi_sfw_T3PIGAM");
  TH1D * hIM3pi_rhopi = (TH1D *) HSFW1D -> FindObject("h1d_IM3pi_sfw_TRHOPI");
  TH1D * hIM3pi_etagam = (TH1D *) HSFW1D -> FindObject("h1d_IM3pi_sfw_TETAGAM");
  TH1D * hIM3pi_bkgrest = (TH1D *) HSFW1D -> FindObject("h1d_IM3pi_sfw_TBKGREST");
  
  TH1D * hIM3pi_mcrest;

  hIM3pi_mcrest = (TH1D*) hIM3pi_bkgrest -> Clone();
  hIM3pi_mcrest -> Add(hIM3pi_kpm, 1.);
  hIM3pi_mcrest -> Add(hIM3pi_rhopi, 1.);
  hIM3pi_mcrest -> SetName("hIM3pi_mcrest");

  
  // scale IM3pi histos

  cout << "SFW2D Fit Results: Scaling Factors" << "\n"
       << "1: eeg     = " << eeg_sfw << "\n"
       << "2: isr3pi  = " << isr3pi_sfw << "\n"
       << "3: omegapi = " << omegapi_sfw << "\n"
       << "4: etagam  = " << etagam_sfw  << "\n"
       << "5: ksl     = " << ksl_sfw     << "\n"
       << "6: mcrest  = " << mcrest_sfw  << "\n"; 

  TList *H1D_list_scaled = new TList();

  // data
  H1D_list_scaled -> Add(hIM3pi_data);
  
  // eeg
  TH1D * hIM3pi_eeg_sc = (TH1D*) hIM3pi_eeg -> Clone();
  hIM3pi_eeg_sc -> Scale(eeg_sfw);
  hIM3pi_eeg_sc -> SetName("hIM3pi_eeg_sc");
  H1D_list_scaled -> Add(hIM3pi_eeg_sc);

  // isr3pi
  TH1D * hsig_sc = (TH1D*) hsig -> Clone();
  hsig_sc -> Scale(isr3pi_sfw);
  hsig_sc -> SetName("hsig_sc");
  H1D_list_scaled -> Add(hsig);

  // omegapi
  TH1D * hIM3pi_omegapi_sc = (TH1D*) hIM3pi_omegapi -> Clone();
  hIM3pi_omegapi_sc -> Scale(omegapi_sfw);
  hIM3pi_omegapi_sc -> SetName("hIM3pi_omegapi_sc");
  H1D_list_scaled -> Add(hIM3pi_omegapi_sc);

  // etagam
  TH1D * hIM3pi_etagam_sc = (TH1D*) hIM3pi_etagam -> Clone();
  hIM3pi_etagam_sc -> Scale(etagam_sfw);
  hIM3pi_etagam_sc -> SetName("hIM3pi_etagam_sc");
  H1D_list_scaled -> Add(H1D_list_scaled);
			 
  // ksl
  TH1D * hIM3pi_ksl_sc = (TH1D*) hIM3pi_ksl -> Clone();
  hIM3pi_ksl_sc -> Scale(ksl_sfw);
  hIM3pi_ksl_sc -> SetName("hIM3pi_ksl_sc");
  H1D_list_scaled -> Add(hIM3pi_ksl_sc);
			 
  // mcrest
  TH1D * hIM3pi_mcrest_sc = (TH1D*) hIM3pi_mcrest -> Clone();
  hIM3pi_mcrest_sc -> Scale(mcrest_sfw);
  hIM3pi_mcrest_sc -> SetName("hIM3pi_mcrest_sc");
  H1D_list_scaled -> Add(hIM3pi_mcrest_sc);

  // bkgsum
  TH1D* hIM3pi_bkgsum_sc = (TH1D*) hIM3pi_eeg_sc -> Clone();
  hIM3pi_bkgsum_sc -> Add(hIM3pi_omegapi_sc, 1.);
  hIM3pi_bkgsum_sc -> Add(hIM3pi_ksl_sc, 1.);
  hIM3pi_bkgsum_sc -> Add(hIM3pi_etagam_sc, 1.);
  hIM3pi_bkgsum_sc -> Add(hIM3pi_mcrest_sc, 1.);
  hIM3pi_bkgsum_sc -> SetName("hIM3pi_bkgsum_sc");

  /// fill in the tree
  double nb_data = 0.;
  double nb_eeg_sc = 0.;
  double nb_ksl_sc = 0.;
  double nb_omegapi_sc = 0.;
  double nb_isr3pi = 0.;
  double nb_etagam_sc = 0.;
  double nb_mcrest_sc = 0.;
  double nb_bkgsum_sc = 0.;

  TSFW -> SetAutoSave(0);
  
  TSFW -> Branch("Br_nb_data", &nb_data, "Br_nb_data/D");
  TSFW -> Branch("Br_nb_eeg_sc", &nb_eeg_sc, "Br_nb_eeg_sc/D");
  TSFW -> Branch("Br_nb_ksl_sc", &nb_ksl_sc, "Br_nb_ksl_sc/D");
  TSFW -> Branch("Br_nb_omegapi_sc", &nb_omegapi_sc, "Br_nb_omegapi_sc/D");
  TSFW -> Branch("Br_nb_isr3pi", &nb_isr3pi, "Br_nb_isr3pi/D");
  TSFW -> Branch("Br_nb_etagam_sc", &nb_etagam_sc, "Br_nb_etagam_sc/D");
  TSFW -> Branch("Br_nb_mcrest_sc", &nb_mcrest_sc, "Br_nb_mcrest_sc/D");
  TSFW -> Branch("Br_nb_bkgsum_sc", &nb_bkgsum_sc, "Br_nb_bkgsum_sc/D");
  
  for (int i = 1; i <= hIM3pi_data -> GetNbinsX(); i ++ ) {

    // data
    nb_data = hIM3pi_data -> GetBinContent(i);

    // eeg
    nb_eeg_sc = hIM3pi_eeg_sc -> GetBinContent(i);

    // ksl
    nb_ksl_sc = hIM3pi_ksl_sc -> GetBinContent(i);

    // omegapi
    nb_omegapi_sc = hIM3pi_omegapi_sc -> GetBinContent(i);

    // isr3pi
    nb_isr3pi = hsig -> GetBinContent(i);

    // etagam
    nb_etagam_sc = hIM3pi_etagam_sc -> GetBinContent(i);

    // mcrest
    nb_mcrest_sc = hIM3pi_mcrest_sc -> GetBinContent(i);

    // bkgsum
    nb_bkgsum_sc = hIM3pi_bkgsum_sc -> GetBinContent(i);

    if (i == 19) {
      cout << i << "\n"
	   << "nb_data = " << nb_data << "\n"
	   << "nb_eeg_sc = " << nb_eeg_sc << "\n"
	   << "nb_ksl_sc = " << nb_ksl_sc << "\n"
	   << "nb_omegapi_sc = " << nb_omegapi_sc << "\n"
	   << "nb_isr3pi = " << nb_isr3pi << "\n"
	   << "nb_etagam_sc = " << nb_etagam_sc << "\n"
	   << "nb_mcrest_sc = " << nb_mcrest_sc << "\n"
	   << "nb_bkgsum_sc = " << nb_bkgsum_sc << "\n";
    }
    
    TSFW -> Fill();
    
  }

  /// fit
  const int para_nb = 1;

  TMinuit* gMinuit = new TMinuit(para_nb); // maximum number of arameters in ()
  gMinuit -> SetFCN(fcn_sfw);

  Double_t arglist[10];
  Int_t ierflg = 0;

  gMinuit -> mnparm(0, "sig_sfw", isr3pi_sfw, 0.01, 0., 10., ierflg);

  // for max likelihood = 0.5, for chisq = 1.0
  gMinuit -> SetErrorDef(1.0);

  // 1 standard
  // 2 try to improve minimum (slower)
  arglist[0] = 1;
  gMinuit -> mnexcm("SET STR", arglist, 1, ierflg);

  // ready for minimization step
  arglist[0] = 500;
  gMinuit -> mnexcm("MIGRAD", arglist, 1, ierflg); // fit 

  // fit results

  int nvpar = 0, nparx = 0, icstat = 0;
  double amin = 0., edm = 0., errdef = 0;
  double sig_sfw = 0., sig_sfw_err = 0.;

  gMinuit -> GetParameter(0, sig_sfw, sig_sfw_err);
  gMinuit -> mnstat(amin, edm, errdef, nvpar, nparx, icstat);

  cout << "frac = " << frac << endl;
  cout << "sig_sfw = " << sig_sfw << "+/-" << sig_sfw_err << endl;
  cout << "chisq = " << amin << ", fit_indx = " << fit_indx << ", chisq rel. = " << amin / (fit_indx - 1)  << "\n";

  ofstream myfile;
  TString myfile_nm = "../header/sfw1d.txt";
  myfile.open(myfile_nm);
  myfile << "const double sig_sfw = " << sig_sfw << ";\n";
  myfile.close();
  
  
  // save

  TTree* TRESULT = new TTree("TRESULT", "recreate");
  TRESULT -> SetAutoSave(0);

  TRESULT -> Branch("Br_sig_sfw", &sig_sfw, "Br_sig_sfw/D");
  TRESULT -> Fill();
  
  hsig -> Write();
  H1D_list_scaled -> Write("HIM3pi", 1);
  TSFW -> Write();
  TRESULT -> Write();
  
  f_output -> Close();
  
  return 0;

}
