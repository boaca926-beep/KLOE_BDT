#include "../header_bdt/sm_para.h"
#include "../header_bdt/path.h"
//#include "../header/sfw2d.txt"
//#include "../header/sfw1d.txt" 
#include "../header_bdt/method.h"
#include "../header_bdt/omega_fit.h"
#include "../header_bdt/cut_para.h"

int omega_fit(){

  gErrorIgnoreLevel = kError;
  
  cout << "Extract omega parameters ..." << endl;

  // scaling factors
  //getObj(f_sfw2d);
  //getObj(f_sfw1d);

  // sfw1d
  TTree *TSFW1D = (TTree*)f_sfw1d -> Get("TRESULT");

  for (Int_t irow = 0; irow < TSFW1D -> GetEntries(); irow++) {// loop chain

    TSFW1D -> GetEntry(irow);

    sig_sfw = TSFW1D -> GetLeaf("Br_sig_sfw") -> GetValue(0);
    
  }

  // sfw2d
  TTree *TSFW2D = (TTree*)f_sfw2d -> Get("TRESULT");

  for (Int_t irow = 0; irow < TSFW2D -> GetEntries(); irow++) {// loop chain

    TSFW2D -> GetEntry(irow);
    
    eeg_sfw = TSFW2D -> GetLeaf("Br_eeg_sfw") -> GetValue(0);
    isr3pi_sfw = TSFW2D -> GetLeaf("Br_isr3pi_sfw") -> GetValue(0);
    omegapi_sfw = TSFW2D -> GetLeaf("Br_omegapi_sfw") -> GetValue(0);
    etagam_sfw = TSFW2D -> GetLeaf("Br_etagam_sfw") -> GetValue(0);
    ksl_sfw = TSFW2D -> GetLeaf("Br_ksl_sfw") -> GetValue(0);
    mcrest_sfw = TSFW2D -> GetLeaf("Br_mcrest_sfw") -> GetValue(0);

  }

  cout << "Saling factors \n"
       << "eeg_sfw = " << eeg_sfw << "\n"
       << "isr3pi_sfw = " << isr3pi_sfw << "\n"
       << "omegapi_sfw = " << omegapi_sfw << "\n" 
       << "etagam_sfw = " << etagam_sfw << "\n" 
       << "ksl_sfw = " << ksl_sfw << "\n"
       << "mcrest_sfw = " << mcrest_sfw << "\n"
       << "sig_sfw = " << sig_sfw << endl;

  
  // data amd bkg histos
  //getObj(f_cut);
  //getObj(f_hist);

  //checkList(HIM3pi_fit);
  //checkList(HSIG);

  // data
  cout << "gsf = " << gsf << ", exp_type = " << exp_type << endl;
  hdata = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_" + exp_type) -> Clone();
  //hdata -> Draw();
  hdata -> SetName("hdata");
  binsize = hdata -> GetNbinsX();
  hmin = hdata -> GetXaxis() -> GetXmin();
  hmax = hdata -> GetXaxis() -> GetXmax();

  // bkg MC
  eeg_hist = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TEEG") -> Clone();
  //eeg_hist -> Draw();
  omegapi_hist = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TOMEGAPI") -> Clone();
  ksl_hist = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TKSL") -> Clone();
  etagam_hist = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TETAGAM") -> Clone();
  hkpm = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TKPM") -> Clone();
  hrhopi = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TRHOPI") -> Clone();
  
  mcrest_hist = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TBKGREST") -> Clone(); // bkgsum
  mcrest_hist -> Add(hkpm, 1.);
  mcrest_hist -> Add(hrhopi, 1.);

  // efficiency
  get_efficy();
  
  // smearing matrix and smearing signal IM3pi MC true distr.
  TRandom *generator = new TRandom();
  rnd=new TRandom3();

  
  //cout << "Indx = " << type_indx << ", name = " << hdata -> GetName() << endl;
  cout << "nb_data = " << hdata -> Integral(1, binsize) << endl;

  hcorrmatrix = getCorrMatrix(binsize, hmin, hmax);
  hcorrmatrix -> SetName("hcorrmatrix");
  hcorrmatrix -> Sumw2();
  
  hsmearmatr = getSmearMatrix(hcorrmatrix, binsize, hmin, hmax);
  hsmearmatr -> SetName("hsmearmatr");
  hsmearmatr -> Sumw2();

  hsmearmatr_err = getSmearMatrixErr(hcorrmatrix, binsize, hmin, hmax);
  hsmearmatr_err -> SetName("hsmearmatr_err");
  hsmearmatr_err -> Sumw2();

  hsmearmatr_efficy = getSmearMatrixEfficy(hcorrmatrix, binsize, hmin, hmax);
  hsmearmatr_efficy -> SetName("hsmearmatr_efficy");
  hsmearmatr_efficy -> Sumw2();

  hsmearmatr_efficy_err = getSmearMatrixEfficy_Err(hcorrmatrix, binsize, hmin, hmax);
  hsmearmatr_efficy_err -> SetName("hsmearmatr_efficy_err");
  hsmearmatr_efficy_err -> Sumw2();
  
  getSmearSig(); // smeared signal

  // MC normalization
  MCNorm(); 
  scaleGSF(); 

  
  // calculate MC yields from sfw2d
  
  // fill the output tree
  int binsize = hdata -> GetNbinsX();
  double xmin = hdata -> GetXaxis() -> GetXmin();
  double xmax = hdata -> GetXaxis() -> GetXmax();

  cout << "binszie = " << binsize << ", x range(" << xmin << ", " << xmax << ") MeV/c^2\n";
  
  double m3pi = 0., m3pi_lower = 0.,  m3pi_upper = 0.;
  double evnt_data = 0., evnt_data_err = 0.;
  double evnt_eeg = 0., evnt_eeg_err = 0.;
  double evnt_ksl = 0., evnt_ksl_err = 0.;
  double evnt_omegapi = 0., evnt_omegapi_err = 0.;
  double evnt_etagam = 0., evnt_etagam_err = 0.;
  double evnt_mcrest = 0., evnt_mcrest_err = 0.;
  double evnt_bkgsum = 0.;
  double nb_isr3pi_obs = 0., nb_isr3pi_obs_err = 0.;
  double efficy = 0., efficy_err = 0.;
  double W0_full = 0.;
  double isrlumi = 0., isrlumi_apprx = 0.;

  TRESULT -> SetAutoSave(0);

  TRESULT -> Branch("Br_m3pi", &m3pi, "Br_m3pi/D");
  TRESULT -> Branch("Br_efficy", &efficy, "Br_efficy/D");
  TRESULT -> Branch("Br_efficy_err", &efficy_err, "Br_efficy_err/D");
  TRESULT -> Branch("Br_isrlumi", &isrlumi, "Br_isrlumi/D");
  TRESULT -> Branch("Br_isrlumi_apprx", &isrlumi_apprx, "Br_isrlumi_apprx/D");
  TRESULT -> Branch("Br_nb_isr3pi_obs", &nb_isr3pi_obs, "Br_nb_isr3pi_obs/D");
  TRESULT -> Branch("Br_nb_isr3pi_obs_err", &nb_isr3pi_obs_err, "Br_nb_isr3pi_obs_err/D");
  TRESULT -> Branch("Br_W0_full", &W0_full, "Br_W0_full/D");
  
  //hefficy -> Draw();

  for (int i = 1; i <= binsize; i ++ ) {

    m3pi = hdata -> GetBinCenter(i);
    m3pi_lower = hdata -> GetBinLowEdge(i);
    m3pi_upper = hdata -> GetBinLowEdge(i + 1);

    // data
    evnt_data = hdata -> GetBinContent(i);
    evnt_data_err = hdata -> GetBinError(i);

    // eeg
    evnt_eeg = heeg -> GetBinContent(i);
    //evnt_eeg_err = GetBkgErr(nb_data_sum, eeg_hist -> GetBinContent(i), nb_eeg_sum, feeg, feeg_err);
    evnt_eeg_err = heeg -> GetBinError(i);
    
    // ksl
    evnt_ksl = hksl -> GetBinContent(i);
    //evnt_ksl_err = GetBkgErr(nb_data_sum, ksl_hist -> GetBinContent(i), nb_ksl_sum, fksl, fksl_err);
    evnt_ksl_err = hksl -> GetBinError(i);

    // omegapi
    evnt_omegapi = homegapi -> GetBinContent(i);
    //evnt_omegapi_err = GetBkgErr(nb_data_sum, omegapi_hist -> GetBinContent(i), nb_omegapi_sum, fomegapi, fomegapi_err);
    evnt_omegapi_err = homegapi -> GetBinError(i);

    // etagam
    evnt_etagam = hetagam -> GetBinContent(i);
    //evnt_etagam_err = GetBkgErr(nb_data_sum, etagam_hist -> GetBinContent(i), nb_etagam_sum, fetagam, fetagam_err);
    evnt_etagam_err = hetagam -> GetBinError(i);

    // mcrest
    evnt_mcrest = hmcrest -> GetBinContent(i);
    //evnt_mcrest_err = GetBkgErr(nb_data_sum, mcrest_hist -> GetBinContent(i), nb_mcrest_sum, fmcrest, fmcrest_err);
    evnt_mcrest_err = hmcrest -> GetBinError(i);

    // bkgsum
    evnt_bkgsum = evnt_eeg + evnt_ksl + evnt_omegapi + evnt_etagam + evnt_mcrest;

    // isr3pi_obs
    nb_isr3pi_obs = evnt_data - evnt_bkgsum; //h1d_IM3pi_T3PIOBS -> GetBinContent(i);
    nb_isr3pi_obs_err = TMath::Sqrt(evnt_data + evnt_eeg_err * evnt_eeg_err + evnt_ksl_err * evnt_ksl_err + evnt_omegapi_err * evnt_omegapi_err + evnt_etagam_err * evnt_etagam_err + evnt_mcrest_err * evnt_mcrest_err);

    // efficiency
    efficy = hefficy -> GetBinContent(i);
    efficy_err = hefficy -> GetBinError(i);

    // isr lumi
    W0_full = Get_W0_full(&m3pi);

    isrlumi = GetISRLumi_exact(m3pi_lower, m3pi_upper);
    isrlumi_apprx = GetISRLumi_apprx(m3pi, m3pi_lower, m3pi_upper, W0_full);

    /*
    if (i == 119) {

      // summary
      
      cout << "bin indx = " << i << ", m3pi = " << m3pi << ", m3pi_lower = " << m3pi_lower << ", m3pi_upper = " << m3pi_upper << "\n"
	   << "evnt_eeg_err = " << evnt_eeg_err << ", nb_data_sum = " << nb_data_sum << ", event_eeg = " << eeg_hist -> GetBinContent(i) << ", nb_eeg_sum = " << nb_eeg_sum << ", feeg = " << feeg << "+/-" << feeg_err << "\n"
	   << "evnt_ksl_err = " << evnt_ksl_err << ", nb_data_sum = " << nb_data_sum << ", event_ksl = " << ksl_hist -> GetBinContent(i) << ", nb_ksl_sum = " << nb_ksl_sum << ", fksl = " << fksl << "+/-" << fksl_err << "\n"
	   << "evnt_omegapi_err = " << evnt_omegapi_err << ", nb_data_sum = " << nb_data_sum << ", event_omegapi = " << omegapi_hist -> GetBinContent(i) << ", nb_omegapi_sum = " << nb_omegapi_sum << ", fomegapi = " << fomegapi << "+/-" << fomegapi_err << "\n"
	   << "evnt_etagam_err = " << evnt_etagam_err << ", nb_data_sum = " << nb_data_sum << ", event_etagam = " << etagam_hist -> GetBinContent(i) << ", nb_etagam_sum = " << nb_etagam_sum << ", fetagam = " << fetagam << "+/-" << fetagam_err << "\n"
	   << "evnt_mcrest_err = " << evnt_mcrest_err << ", nb_data_sum = " << nb_data_sum << ", event_mcrest = " << mcrest_hist -> GetBinContent(i) << ", nb_mcrest_sum = " << nb_mcrest_sum << ", fmcrest = " << fmcrest << "+/-" << fmcrest_err << "\n\n"

      	   << "\tevnt_data = " << evnt_data << " +/- " << evnt_data_err << "\n"
	   << "\tevnt_eeg = " << evnt_eeg << " +/- " << evnt_eeg_err << "\n"
	   << "\tevnt_ksl = " << evnt_ksl << " +/- " << evnt_ksl_err << "\n"
	   << "\tevnt_omgapi = " << evnt_omegapi << " +/- " << evnt_omegapi_err << "\n"
	   << "\tevnt_etagam = " << evnt_etagam << " +/- " << evnt_etagam_err << "\n"
	   << "\tevnt_mcrest = " << evnt_mcrest << " +/- " << evnt_mcrest_err << "\n"
	   << "\tevnt_bkgsum = " << evnt_bkgsum << "\n"
	   << "\tnb_isr3pi_obs = " << nb_isr3pi_obs << "+/-" << nb_isr3pi_obs_err << "\n"
	   << "\tefficy = " << efficy << "+/-" << efficy_err << "\n"
	   << "\tLumi_tot " << Lumi_tot << ", gsf = " << gsf << ", Lumi_int = " << Lumi_int << "\n"
	   << "\tisrlumi_apprx = " << isrlumi_apprx << "\n";

      
    }    
    */
    
    TRESULT -> Fill();

  }

  ///
  cout << "\nPerform cross section fit ..." << endl;

  const int para_nb_crx3pi = 3; //number of parameters

  fun_bw = new TF1("func_bw", Func_bw, 0., 2000., 3);
  
  TMinuit * gMinuit_crx3pi = new TMinuit(para_nb_crx3pi); // maximum number of parameters in ()
  gMinuit_crx3pi -> SetFCN(fcn_crx3pi); 

  Double_t arglist_crx3pi[10];
  Int_t ierflg_crx3pi = 0;

  // Set print level
  //gMinuit_crx3pi -> SetPrintLevel(-1);

  double sigma_max = 1606.53; //nb
  double BB = GetBB_new(sigma_max, mass_omega); //sigma_max * 1e-9 * (mass_omega * mass_omega) / 12. / pi;
  cout << GetCrx3piMax_new(6.7e-5, mass_omega) << " [nb]\n";

  cout << "mass_omega = " << mass_omega << ", Gam_omega = " << Gam_omega << ", BB = " << BB << endl;
  TF1 *f_bw_test = new TF1("f_bw_test", Func_bw, 0., 2000., 3);
  f_bw_test -> SetParameters(mass_omega, Gam_omega, BB);
  cout << f_bw_test -> Eval(mass_omega) << endl;
    
  gMinuit_crx3pi -> mnparm(0, "M_V",   mass_omega, 0.01, 0., 1000., ierflg_crx3pi);
  gMinuit_crx3pi -> mnparm(1, "Gam_V", Gam_omega,  0.01, 0., 100.,  ierflg_crx3pi);
  gMinuit_crx3pi -> mnparm(2, "BB",    BB,         0.01, 0., 1.,    ierflg_crx3pi);

  // No warnings
  //gMinuit_crx3pi -> mnexcm("SET NOW", arglist_crx3pi, 1, ierflg_crx3pi);

  // for max likelihood = 0.5, for chisq = 1.0
  gMinuit_crx3pi -> SetErrorDef(0.5);

  // minimization strategy
  //arglist_crx3pi[0] = 2;
  //gMinuit_crx3pi -> mnexcm("SET STR",arglist,1,ierflg);
  
  // ready for minimization step
  arglist_crx3pi[0] = 500;
  // perform the fit !!
  gMinuit_crx3pi -> mnexcm("MIGRAD", arglist_crx3pi, 1, ierflg_crx3pi);

  // get fit results
  double Mass_omega_fit = 0., Mass_omega_err_fit = 0.; 
  double Gam_omega_fit = 0., Gam_omega_err_fit = 0.; 
  double BB_fit = 0., BB_err_fit = 0.; 
  double crx3pi_max_fit = 0.;
  
  double amin_crx3pi = 0., edm_crx3pi = 0., errdef_crx3pi = 0.;
  int nvpar_crx3pi = 0, nparx_crx3pi = 0, icstat_crx3pi = 0;

  gMinuit_crx3pi -> GetParameter(0, Mass_omega_fit, Mass_omega_err_fit);
  gMinuit_crx3pi -> GetParameter(1, Gam_omega_fit, Gam_omega_err_fit);
  gMinuit_crx3pi -> GetParameter(2, BB_fit, BB_err_fit);
  gMinuit_crx3pi -> mnstat(amin_crx3pi, edm_crx3pi, errdef_crx3pi, nvpar_crx3pi, nparx_crx3pi, icstat_crx3pi);
  
  crx3pi_max_fit = GetCrx3piMax_new(BB_fit, Mass_omega_fit);
  int ndf_crx3pi = residul_size_crx3pi - para_nb_crx3pi;

  ofstream myfile;
  TString myfile_nm = "omega_fit.txt";
  myfile.open(myfile_nm);
  
  myfile << "\nFit Results Summary\n"
	 << "gsf = " << gsf << ", Lumi_int = " << Lumi_int << "\n" 
	 << "xmax = " << xmax << ", xmin = " << xmin << ", binsize = " << binsize << "\n"
	 << "chi2_sum_crx3pi (FCN) = " << chi2_sum_crx3pi << "(" << amin_crx3pi << "), ndf = " << ndf_crx3pi << ", chi2_sum / ndf = " << chi2_sum_crx3pi / ndf_crx3pi << "\n"
	 << "Mass_omega_fit [MeV/c^2] = " << Mass_omega_fit << " +/- " << Mass_omega_err_fit << "\n"
	 << "Gam_omega_fit [MeV/c^2] = " << Gam_omega_fit << " +/- " << Gam_omega_err_fit << "\n"
	 << "BB_fit[10^5] = " << BB_fit * 1e5 << " +/- " << BB_err_fit * 1e5 << "\n"
	 << "crx3pi_max [nb] = " << crx3pi_max_fit << "\n"
	 << "amin_crx3pi = " << amin_crx3pi << "\n\n";
  

  myfile.close();
  
  /// save fit result tree
  int bin_indx = 0;
  double m3pi_fit = 0., m3pi_fit_err = 0.;
  double n3pi_obs_fit = 0., n3pi_obs_fit_err = 0.;
  double n3pi_fit = 0., n3pi_fit_err = 0.;
  double n3pi_diff = 0., n3pi_diff_err = 0.;

  double OMEGA_PARA[3] = {Mass_omega_fit, Gam_omega_fit, BB_fit};
  double OMEGA_PARA_ERR[3] = {Mass_omega_err_fit, Gam_omega_err_fit, BB_err_fit};
  
  cout << Lumi_int << endl;
    
  TCRX3PI -> SetAutoSave(0);

  TCRX3PI -> Branch("Br_Lumi_int", &Lumi_int, "Br_Lumi_int/D");
  TCRX3PI -> Branch("Br_chi2_sum_crx3pi", &chi2_sum_crx3pi, "Br_chi2_sum_crx3pi/D");
  TCRX3PI -> Branch("Br_ndf", &ndf_crx3pi, "Br_ndf/I");
  TCRX3PI -> Branch("Br_gsf", &gsf, "Br_gsf/D");
  TCRX3PI -> Branch("Br_cut_value", &cut_value, "Br_cut_value/D");

  
  //TCRX3PI -> Branch("Br_Mass_omega_fit", &Mass_omega_fit, "Br_Mass_omega_fit/D");
  //TCRX3PI -> Branch("Br_Mass_omega_err_fit", &Mass_omega_err_fit, "Br_Mass_omega_err_fit/D");
  
  //TCRX3PI -> Branch("Br_Gam_omega_fit", &Gam_omega_fit, "Br_Gam_omega_fit/D");
  //TCRX3PI -> Branch("Br_Gam_omega_err_fit", &Gam_omega_err_fit, "Br_Gam_omega_err_fit/D");
  
  //TCRX3PI -> Branch("Br_BB_fit", &BB_fit, "Br_BB_fit/D");
  //TCRX3PI -> Branch("Br_BB_err_fit", &BB_err_fit, "Br_BB_err_fit/D");
  
  TCRX3PI -> Branch("Br_OMEGA_PARA", &OMEGA_PARA, "Br_OMEGA_PARA[3]/D");
  TCRX3PI -> Branch("Br_OMEGA_PARA_ERR", &OMEGA_PARA_ERR, "Br_OMEGA_PARA_ERR[3]/D");
  
  TCRX3PI -> Branch("Br_bin_indx", &bin_indx, "Br_bin_indx/I");
  TCRX3PI -> Branch("Br_M3PI", &M3PI, "Br_M3PI[1000]/D");

  TCRX3PI -> Branch("Br_EFFICY", &EFFICY, "Br_EFFICY[1000]/D");
  TCRX3PI -> Branch("Br_EFFICY_ERR", &EFFICY_ERR, "Br_EFFICY_ERR[1000]/D");

  TCRX3PI -> Branch("Br_ISRLUMI", ISRLUMI, "Br_ISRLUMI[1000]/D");
  TCRX3PI -> Branch("Br_CRX3PI_BW", CRX3PI_BW, "Br_CRX3PI_BW[1000]/D");
    
  TCRX3PI -> Branch("Br_N3PI_OBS", &N3PI_OBS, "Br_N3PI_OBS[1000]/D");
  TCRX3PI -> Branch("Br_N3PI_OBS_ERR", &N3PI_OBS_ERR, "Br_N3PI_OBS_ERR[1000]/D");

  TCRX3PI -> Branch("Br_N3PI_SMEAR", &N3PI_SMEAR, "Br_N3PI_SMEAR[1000]/D");
  TCRX3PI -> Branch("Br_N3PI_SMEAR_ERR", &N3PI_SMEAR_ERR, "Br_N3PI_SMEAR_ERR[1000]/D");

  TCRX3PI -> Branch("Br_N3PI_DIFF", &N3PI_DIFF, "Br_N3PI_DIFF[1000]/D");
  TCRX3PI -> Branch("Br_N3PI_DIFF_ERR", &N3PI_DIFF_ERR, "Br_N3PI_DIFF_ERR[1000]/D");
  
  TCRX3PI -> Fill();
  
  for (int i = 0; i < residul_size_crx3pi; i ++) {

    bin_indx = BIN_INDX[i];
    
    m3pi_fit = M3PI[bin_indx];
      
    n3pi_obs_fit = N3PI_OBS[bin_indx];
    n3pi_obs_fit_err = N3PI_OBS_ERR[bin_indx];

    n3pi_fit = N3PI_SMEAR[bin_indx];
    n3pi_fit_err = N3PI_SMEAR_ERR[bin_indx];

    n3pi_diff = N3PI_DIFF[bin_indx];
    n3pi_diff_err = N3PI_DIFF_ERR[bin_indx];

    cout << i << ", " << bin_indx << ", " << m3pi_fit << ", " << n3pi_obs_fit << "+/-" << n3pi_obs_fit_err << ", " << n3pi_fit << "+/-" << n3pi_fit_err << endl;
    TCRX3PI -> Fill();
    
  }

  /// save histos
  hcorrmatrix -> Write();
  hsmearmatr -> Write();
  hsmearmatr_err -> Write();
  hsmearmatr_efficy -> Write();
  hefficy -> Write();

  hsig -> Write("hsig"); // SIG
  hdata -> Write("hdata"); // DATA
  heeg -> Write("heeg"); // scale EEG
  homegapi -> Write("homegapi"); // scale OMEGAPI
  hetagam -> Write("hetagam"); // scale ETAGAM
  hksl -> Write("hksl"); // scale KSL 
  hmcrest -> Write("hmcrest"); // scale MCREST

  hsig_true -> Write();
  hsig_gen -> Write();
  
  TCRX3PI -> Write();
  TRESULT -> Write();
  
  f_out -> Close();
  
  return 0;
  
}
