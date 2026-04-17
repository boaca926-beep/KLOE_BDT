#include "../header_bdt/sm_para.h"
#include "../header_bdt/path.h"
#include "../header_bdt/sfw2d.h"
#include "../header_bdt/method.h"
#include "../header_bdt/cut_para.h"

int sfw2d(){
  /// get histos
  TFile *f_input = new TFile(outputHist + "hist.root");  
  TFile *f_output = new TFile(outputSfw2D + "sfw2d.root", "recreate");

  cout << cut_value << endl;
  
  //getObj(f_input);
  
  TList *HSFW2D = (TList *) f_input -> Get("HSFW2D");
  checkList(HSFW2D);
  
  TH2D *h2d_sfw_TEEG = (TH2D *) HSFW2D -> FindObject("h2d_sfw_TEEG");
  TH2D *h2d_sfw_TDATA = (TH2D *) HSFW2D -> FindObject("h2d_sfw_TDATA");
  TH2D *h2d_sfw_TISR3PI_SIG = (TH2D *) HSFW2D -> FindObject("h2d_sfw_TISR3PI_SIG");
  TH2D *h2d_sfw_TOMEGAPI = (TH2D *) HSFW2D -> FindObject("h2d_sfw_TOMEGAPI");
  TH2D *h2d_sfw_TKPM = (TH2D *) HSFW2D -> FindObject("h2d_sfw_TKPM");
  TH2D *h2d_sfw_TKSL = (TH2D *) HSFW2D -> FindObject("h2d_sfw_TKSL");
  TH2D *h2d_sfw_T3PIGAM = (TH2D *) HSFW2D -> FindObject("h2d_sfw_T3PIGAM");
  TH2D *h2d_sfw_TRHOPI = (TH2D *) HSFW2D -> FindObject("h2d_sfw_TRHOPI");
  TH2D *h2d_sfw_TETAGAM = (TH2D *) HSFW2D -> FindObject("h2d_sfw_TETAGAM");
  TH2D *h2d_sfw_TBKGREST = (TH2D *) HSFW2D -> FindObject("h2d_sfw_TBKGREST");

  TH2D * h2d_sfw_MCREST;
  TH2D * h2d_sfw_MCSUM;
  
  // h2d_sfw_MCREST
  h2d_sfw_MCREST = (TH2D*) h2d_sfw_TBKGREST -> Clone();
  h2d_sfw_MCREST -> Add(h2d_sfw_TKPM, 1.);
  h2d_sfw_MCREST -> Add(h2d_sfw_TRHOPI, 1.);
  h2d_sfw_MCREST -> SetName("h2d_sfw_MCREST");

  // h2d_sfw_MCSUM
  h2d_sfw_MCSUM = (TH2D*) h2d_sfw_TEEG -> Clone();
  h2d_sfw_MCSUM -> Add(h2d_sfw_TOMEGAPI, 1.);
  h2d_sfw_MCSUM -> Add(h2d_sfw_TKSL, 1.);
  h2d_sfw_MCSUM -> Add(h2d_sfw_TETAGAM, 1.);
  h2d_sfw_MCSUM -> Add(h2d_sfw_TISR3PI_SIG, 1.);
  h2d_sfw_MCSUM -> Add(h2d_sfw_MCREST, 1.);
  h2d_sfw_MCSUM -> SetName("h2d_sfw_MCSUM");

  // define and initialize variable variables
  double nb_data = 0.;
  double nb_eeg = 0.;
  double nb_omegapi = 0.;
  double nb_ksl = 0.;
  double nb_etagam = 0.;
  double nb_isr3pi = 0.;
  double nb_mcrest = 0.;
  double nb_mc = 0., nb_mcsum = 0.;

  // fractions
  double feeg = 0.,      feeg_err = 0.;
  double fisr3pi = 0.,   fisr3pi_err = 0.;
  double fomegapi = 0.,  fomegapi_err = 0.;
  double fetagam = 0.,   fetagam_err = 0.;
  double fksl = 0.,      fksl_err = 0.;
  double fmcrest = 0.,   fmcrest_err = 0.;

  // scaling factors
  double isr3pi_sfw = 0,  isr3pi_sfw_err = 0; 
  double omegapi_sfw = 0, omegapi_sfw_err = 0; 
  double etagam_sfw = 0,  etagam_sfw_err = 0; 
  double ksl_sfw = 0,     ksl_sfw_err = 0; 
  double mcrest_sfw  = 0, mcrest_sfw_err = 0; 

  // define tree
  TSFW2D -> SetAutoSave(0);

  TSFW2D -> Branch("Br_nb_data", &nb_data, "Br_nb_data/D");
  TSFW2D -> Branch("Br_nb_eeg", &nb_eeg, "Br_nb_eeg/D");
  TSFW2D -> Branch("Br_nb_ksl", &nb_ksl, "Br_nb_ksl/D");
  TSFW2D -> Branch("Br_nb_omegapi", &nb_omegapi, "Br_nb_omegapi/D");
  TSFW2D -> Branch("Br_nb_etagam", &nb_etagam, "Br_nb_etagam/D");
  TSFW2D -> Branch("Br_nb_isr3pi", &nb_isr3pi, "Br_nb_isr3pi/D");
  TSFW2D -> Branch("Br_nb_mcrest", &nb_mcrest, "Br_nb_mcrest/D");

  for (int i = 1; i <= h2d_sfw_TDATA -> ProjectionX() -> GetNbinsX(); i ++ ) {

    for (int j = 1; j <= h2d_sfw_TDATA -> ProjectionY() -> GetNbinsX(); j ++ ) {

      // data
      nb_data = h2d_sfw_TDATA -> GetBinContent(i, j);
      nb_data_sum += nb_data;

      //cout << i << ", " << j << ", " << h2d_sfw_TDATA -> GetBinContent(i, j) << endl;

      // eeg
      nb_eeg = h2d_sfw_TEEG -> GetBinContent(i, j);
      nb_eeg_sum += nb_eeg;

      // omegapi
      nb_omegapi = h2d_sfw_TOMEGAPI -> GetBinContent(i, j);
      nb_omegapi_sum += nb_omegapi;

      // ksl
      nb_ksl = h2d_sfw_TKSL -> GetBinContent(i, j);
      nb_ksl_sum += nb_ksl;

      // etagam
      nb_etagam = h2d_sfw_TETAGAM -> GetBinContent(i, j);
      nb_etagam_sum += nb_etagam;

      // isr3pi
      nb_isr3pi = h2d_sfw_TISR3PI_SIG -> GetBinContent(i, j);
      nb_isr3pi_sum += nb_isr3pi;

      // mcrest
      nb_mcrest = h2d_sfw_MCREST -> GetBinContent(i, j);
      nb_mcrest_sum += nb_mcrest;

      // mcsum
      nb_mc = h2d_sfw_MCSUM -> GetBinContent(i, j);
      nb_mcsum += nb_mc;

      
      TSFW2D -> Fill();
  
    }
    
  }

  ofstream myfile;
  TString myfile_nm = "../header/sfw2d.txt";
  myfile.open(myfile_nm);
  myfile << "const double nb_data_sum = " << nb_data_sum << ";\n"
	 << "const double nb_eeg_sum = " << nb_eeg_sum << ";\n"
	 << "const double nb_omegapi_sum = " << nb_omegapi_sum << ";\n"
	 << "const double nb_ksl_sum = " << nb_ksl_sum << ";\n"
	 << "const double nb_etagam_sum = " << nb_etagam_sum << ";\n"
	 << "const double nb_isr3pi_sum = " << nb_isr3pi_sum << ";\n"
	 << "const double nb_mcrest_sum = " << nb_mcrest_sum << ";\n"
	 << "const double nb_mcsum = " << nb_mcsum << ";\n\n";

  /// Fitting

  double feeg_init     = nb_eeg_sum     / nb_mcsum;
  double fisr3pi_init  = nb_isr3pi_sum  / nb_mcsum;
  double fomegapi_init = nb_omegapi_sum / nb_mcsum;
  double fksl_init     = nb_ksl_sum     / nb_mcsum;
  double fmcrest_init  = nb_mcrest_sum  / nb_mcsum;
  double fetagam_init  = nb_etagam_sum  / nb_mcsum;

  const int para_nb_sfw2d = 6;
  //const TString para_nb_str = Form("%d", para_nb);

  TMinuit *gMinuit = new TMinuit(para_nb_sfw2d); // maximum number of parameters in ()
  gMinuit -> SetFCN(fcn_sfw2d);
  
  Double_t arglist[10];
  Int_t ierflg = 0;

  // Set print level
  //gMinuit -> SetPrintLevel(-1);
  
  gMinuit -> mnparm(0, "feeg_ML",     feeg_init,     0.01, 0., 10.,ierflg);
  gMinuit -> mnparm(1, "fisr3pi_ML",  fisr3pi_init,  0.01, 0., 10.,ierflg);
  gMinuit -> mnparm(2, "fomegapi_ML", fomegapi_init, 0.01, 0., 10.,ierflg);
  gMinuit -> mnparm(3, "fetagam_ML",  fetagam_init,  0.01, 0., 10.,ierflg);
  gMinuit -> mnparm(4, "fksl_ML",     fksl_init,     0.01, 0., 10.,ierflg);
  gMinuit -> mnparm(5, "fmcrest_ML",  fmcrest_init,  0.01, 0., 10.,ierflg);

  gMinuit -> SetErrorDef(0.5);

  // ready for minimization step
  arglist[0] = 500;
  gMinuit -> mnexcm("MIGRAD", arglist, 1, ierflg); // fit sfw2d

  // get sfw2d fit results

  gMinuit -> GetParameter(0, feeg, feeg_err);
  gMinuit -> GetParameter(1, fisr3pi, fisr3pi_err);
  gMinuit -> GetParameter(2, fomegapi, fomegapi_err);
  gMinuit -> GetParameter(3, fetagam, fetagam_err);
  gMinuit -> GetParameter(4, fksl, fksl_err);
  gMinuit -> GetParameter(5, fmcrest, fmcrest_err);

  double eeg_sfw = getscale(nb_data_sum, feeg, nb_eeg_sum);
  double eeg_sfw_err = GetScalError(nb_data_sum, nb_eeg_sum, feeg, feeg_err);

  isr3pi_sfw = getscale(nb_data_sum, fisr3pi, nb_isr3pi_sum);
  isr3pi_sfw_err = GetScalError(nb_data_sum, nb_isr3pi_sum, fisr3pi, fisr3pi_err);

  omegapi_sfw = getscale(nb_data_sum, fomegapi, nb_omegapi_sum);
  omegapi_sfw_err = GetScalError(nb_data_sum, nb_omegapi_sum, fomegapi, fomegapi_err);

  etagam_sfw = getscale(nb_data_sum, fetagam, nb_etagam_sum);
  etagam_sfw_err = GetScalError(nb_data_sum, nb_etagam_sum, fetagam, fetagam_err);

  ksl_sfw = getscale(nb_data_sum, fksl, nb_ksl_sum);
  ksl_sfw_err = GetScalError(nb_data_sum, nb_ksl_sum, fksl, fksl_err);

  mcrest_sfw  = getscale(nb_data_sum, fmcrest, nb_mcrest_sum);
  mcrest_sfw_err = GetScalError(nb_data_sum, nb_mcrest_sum, fmcrest, fmcrest_err);

  // write in the output file
  myfile << "const double feeg = " << feeg << ";\n"
	 << "const double fisr3pi = " << fisr3pi << ";\n"
	 << "const double fomegapi = " << fomegapi << ";\n"
	 << "const double fetagam = " << fetagam << ";\n"
	 << "const double fksl = " << fksl << ";\n"   
 	 << "const double fmcrest = " << fmcrest << ";\n\n"

	 << "const double feeg_err = " << feeg_err << ";\n"
	 << "const double fisr3pi_err = " << fisr3pi_err << ";\n"
	 << "const double fomegapi_err = " << fomegapi_err << ";\n"
	 << "const double fetagam_err = " << fetagam_err << ";\n"
	 << "const double fksl_err = " << fksl_err << ";\n"
	 << "const double fmcrest_err = " << fmcrest_err << ";\n\n"
  
	 << "const double eeg_sfw = " << eeg_sfw << ";\n"
	 << "const double isr3pi_sfw = " << isr3pi_sfw << ";\n"
	 << "const double omegapi_sfw = " << omegapi_sfw << ";\n"
	 << "const double etagam_sfw = " << etagam_sfw << ";\n"
	 << "const double ksl_sfw = " << ksl_sfw << ";\n"    
 	 << "const double mcrest_sfw = " << mcrest_sfw << ";\n";
  myfile.close();

  //
  cout << "\nSFW2D Fit Results: Fractions (initial)[%]" << "\n"
    
       << "1: eeg     = " << feeg * 100. << "(" << feeg_init * 100. << ") +/- " << feeg_err * 100.     << "\n"
       << "2: isr3pi  = " << fisr3pi * 100. << "(" << fisr3pi_init * 100. << ") +/- " << fisr3pi_err * 100. << "\n"
       << "3: omegapi = " << fomegapi * 100. << "(" << fomegapi_init * 100. << ") +/- " << fomegapi_err * 100. << "\n"
       << "4: etagam  = " << fetagam * 100. << "(" << fetagam_init * 100. << ")+/- " << fetagam_err * 100. << "\n"
       << "5: ksl     = " << fksl * 100. << "(" << fksl_init * 100. << ") +/- " << fksl_err * 100. << "\n"
       << "6: mcrest  = " << fmcrest * 100. << "(" << fmcrest_init * 100. << ") +/- " << fmcrest_err * 100. << "\n"
       << "       sum = " << (feeg + fisr3pi + fomegapi + fetagam + fksl + fmcrest) * 100. << "(" << (feeg_init + fisr3pi_init + fomegapi_init + fetagam_init + fksl_init + fmcrest_init) * 100. << ")\n\n";

  //double para_nb_sfw2d = 6.;
  int ndf_sfw2d = residul_size_sfw2d - para_nb_sfw2d;
  double p_value_chi2 = TMath::Prob(chi2_sfw2d_sum, ndf_sfw2d);
    
  cout << "SFW2D Fit Results: Scaling Factors" << "\n"
       << "1: eeg     = " << eeg_sfw     << " +/- " << eeg_sfw_err << "\n"
       << "2: isr3pi  = " << isr3pi_sfw << " +/- " << isr3pi_sfw_err << "\n"
       << "3: omegapi = " << omegapi_sfw << " +/- " << omegapi_sfw_err << "\n"
       << "4: etagam  = " << etagam_sfw  << " +/- " << etagam_sfw_err << "\n"
       << "5: ksl     = " << ksl_sfw     << " +/- " << ksl_sfw_err << "\n"
       << "6: mcrest  = " << mcrest_sfw  << " +/- " << mcrest_sfw_err << "\n\n";

  cout << "SFW2D Fit Results: Chi2, ndf ..." << "\n"
       << "residul_size = " << residul_size_sfw2d << ", para_nb = " << para_nb_sfw2d << ", ndf = " << ndf_sfw2d << "\n"
       << "chi2_sum = " << chi2_sfw2d_sum << ", chi2_sum / ndf = " << chi2_sfw2d_sum / ndf_sfw2d << "\n"
       << "p-value = " << p_value_chi2 << "\n";

  TTree* TRESULT = new TTree("TRESULT", "recreate");
  TRESULT -> SetAutoSave(0);

  const int sf_nb = 6;
  double SF[sf_nb] = {eeg_sfw,
		      isr3pi_sfw,
		      omegapi_sfw,
		      etagam_sfw,
		      ksl_sfw,
		      mcrest_sfw};
  
  double SF_ERR[sf_nb] = {eeg_sfw_err,
			  isr3pi_sfw_err,
			  omegapi_sfw_err,
			  etagam_sfw_err,
			  ksl_sfw_err,
			  mcrest_sfw_err};
  
  // scaling factors
  TRESULT -> Branch("Br_cut_value", &cut_value, "Br_cut_value/D");
  
  TRESULT -> Branch("Br_SF", &SF, "Br_SF[6]/D");
  TRESULT -> Branch("Br_SF_ERR", &SF_ERR, "Br_SF_ERR[6]/D");
  
  TRESULT -> Branch("Br_eeg_sfw", &eeg_sfw, "Br_eeg_sfw/D");
  TRESULT -> Branch("Br_eeg_sfw_err", &eeg_sfw_err, "Br_eeg_sfw_err/D");
  
  TRESULT -> Branch("Br_isr3pi_sfw", &isr3pi_sfw, "Br_isr3pi_sfw/D");
  TRESULT -> Branch("Br_isr3pi_sfw_err", &isr3pi_sfw_err, "Br_isr3pi_sfw_err/D");

  TRESULT -> Branch("Br_omegapi_sfw", &omegapi_sfw, "Br_omegapi_sfw/D");
  TRESULT -> Branch("Br_omegapi_sfw_err", &omegapi_sfw_err, "Br_omegapi_sfw_err/D");

  TRESULT -> Branch("Br_etagam_sfw", &etagam_sfw, "Br_etagam_sfw/D");
  TRESULT -> Branch("Br_etagam_sfw_err", &etagam_sfw_err, "Br_etagam_sfw_err/D");

  TRESULT -> Branch("Br_ksl_sfw", &ksl_sfw, "Br_ksl_sfw/D");
  TRESULT -> Branch("Br_ksl_sfw_err", &ksl_sfw_err, "Br_ksl_sfw_err/D");

  TRESULT -> Branch("Br_mcrest_sfw", &mcrest_sfw, "Br_mcrest_sfw/D");
  TRESULT -> Branch("Br_mcrest_sfw_err", &mcrest_sfw_err, "Br_mcrest_sfw_err/D");

  // fractions
  TRESULT -> Branch("Br_feeg", &feeg, "Br_feeg/D");
  TRESULT -> Branch("Br_feeg_err", &feeg_err, "Br_feeg_err/D");

  TRESULT -> Branch("Br_fisr3pi", &fisr3pi, "Br_fisr3pi/D");
  TRESULT -> Branch("Br_fisr3pi_err", &fisr3pi_err, "Br_fisr3pi_err/D");

  TRESULT -> Branch("Br_fomegapi", &fomegapi, "Br_fomegapi/D");
  TRESULT -> Branch("Br_fomegapi_err", &fomegapi_err, "Br_fomegapi_err/D");

  TRESULT -> Branch("Br_fetagam", &fetagam, "Br_fetagam/D");
  TRESULT -> Branch("Br_fetagam_err", &fetagam_err, "Br_fetagam_err/D");

  TRESULT -> Branch("Br_fksl", &fksl, "Br_fksl/D");
  TRESULT -> Branch("Br_fksl_err", &fksl_err, "Br_fksl_err/D");

  TRESULT -> Branch("Br_fmcrest", &fmcrest, "Br_fmcrest/D");
  TRESULT -> Branch("Br_fmcrest_err", &fmcrest_err, "Br_fmcrest_err/D");

  // event yield
  TRESULT -> Branch("Br_nb_data_sum", &nb_data_sum, "Br_nb_data_sum/D");
  TRESULT -> Branch("Br_nb_eeg_sum", &nb_eeg_sum, "Br_nb_eeg_sum/D");
  TRESULT -> Branch("Br_nb_isr3pi_sum", &nb_isr3pi_sum, "Br_nb_isr3pi_sum/D");
  TRESULT -> Branch("Br_nb_omegapi_sum", &nb_omegapi_sum, "Br_nb_omegapi_sum/D");
  TRESULT -> Branch("Br_nb_etagam_sum", &nb_etagam_sum, "Br_nb_etagam_sum/D");
  TRESULT -> Branch("Br_nb_ksl_sum", &nb_ksl_sum, "Br_nb_ksl_sum/D");
  TRESULT -> Branch("Br_nb_mcrest_sum", &nb_mcrest_sum, "Br_nb_mcrest_sum/D");
  
  
  TRESULT -> Fill();

  // save
  
  TRESULT -> Write();
  TSFW2D -> Write();
  
  f_output -> Close();
  
  return 0;
  
}
