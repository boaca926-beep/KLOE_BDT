TH1D *hufo;
TH1D *hdata;
TH1D *heeg, *eeg_hist;
TH1D *hksl, *ksl_hist;
TH1D *homegapi, *omegapi_hist;
TH1D *hetagam, *etagam_hist;
TH1D *hmcrest, *mcrest_hist;
TH1D *hkpm;
TH1D *hrhopi;

TH1D *hsig;
TH1D *hefficy;

TH2D *hcorrmatrix;
TH2D *hsmearmatr;
TH2D *hsmearmatr_err;
TH2D *hsmearmatr_efficy;
TH2D *hsmearmatr_efficy_err;

TF1 *fun_bw;

const int list_size = 1000;
int binsize = 0;
double hmin = 0.;
double hmax = 0.;

double M3PI[list_size]; 
double EFFICY[list_size]; 
double EFFICY_ERR[list_size]; 
double ISRLUMI[list_size]; 
double CRX3PI_BW[list_size]; 

double N3PI_SMEAR[list_size], N3PI_SMEAR_ERR[list_size];
double N3PI_OBS[list_size], N3PI_OBS_ERR[list_size];
double N3PI_PRE[list_size];
double N3PI_DIFF[list_size], N3PI_DIFF_ERR[list_size];
double BIN_INDX[list_size];
double M3PI_FIT[list_size], M3PI_FIT_ERR[list_size]; 

double chi2_sum_crx3pi = 0.;
int residul_size_crx3pi = 0;

TFile *f_hist = new TFile(outputHist + "hist.root"); 
TFile *f_cut = new TFile(outputCut + "tree_pre.root");  
TFile *f_sfw2d = new TFile(outputSfw2D + "sfw2d.root");  
TFile *f_sfw1d = new TFile(outputSfw1D + "sfw1d.root");  

TFile *f_out = new TFile(outputOmega + "omega_fit.root", "recreate");

TTree* TRESULT = new TTree("TRESULT", "recreate");
TTree* TCRX3PI = new TTree("TCRX3PI", "recreate");
TTree* TINPUT = new TTree("TINPUT", "recreate");
  
TList *HIM3pi_fit = (TList *) f_hist -> Get("HIM3pi_fit");
TList *HSIG = (TList *) f_hist -> Get("HSIG");

TH1D *hsig_true = (TH1D *) HSIG -> FindObject("hsig_true");  
TH1D *hsig_gen = (TH1D *) HSIG -> FindObject("hsig_gen");

// scaling factors
double eeg_sfw = 0.; //1.66579;
double isr3pi_sfw = 0.; //0.0491344;
double omegapi_sfw = 0.; //1.24086;
double etagam_sfw = 0.; //1.01706;
double ksl_sfw = 0.; //1.03746;
double mcrest_sfw = 0.; //4.84081;

double sig_sfw = 0.; //0.0461961;

//double gsf = 1;

//const double crx_scale = 0.3894e6; // 1 [GeV]^{-2}=0.3894 mb=0.3894e6nb

const double fit_min = 760;
const double fit_max = 800;

// methods
double get_sigma2_y(double a, double b, double a_err, double b_err) {

  double term1 = TMath::Power(b * a_err, 2);
  double term2 = TMath::Power(a * b_err, 2);
  
  double sigma2 = term1 + term2;
  
  return sigma2;
}

double Get_W0_full(double * m3pi) {

  double W0_full = 0.;
  double mm = 0., x = 0.;
  
  mm = m3pi[0] * 1e-3;
  x = 1.-TMath::Power(mm/sqrtS, 2.);
  
  W0_full = (alphapi/x)*(TMath::Log(TMath::Power(sqrtS/me,2))-1)*(2-2*x+x*x);

  //cout << mm << endl;
  
  return W0_full;
 
}

double W0_full_fun(double *xx, double *par) {

  double W0_full = 0.;
  double x = 0.;
  
  x = xx[0];

  W0_full = (alphapi/x)*(TMath::Log(TMath::Power(sqrtS/me,2))-1)*(2-2*x+x*x);

  //cout << "(par1, par2, par3) = (" << par[0] << ", " << par[1] << ", " << par[2] << ")" << endl;
  //cout << "x = " << x << ", W0_full = " << W0_full << "\n";
  
  return W0_full;
}

double GetISRLumi_exact(double m3pi_lower, double m3pi_upper) {

  double x_upper = 1. - TMath::Power(m3pi_lower * 1e-3 / sqrtS, 2.);
  double x_lower = 1. - TMath::Power(m3pi_upper * 1e-3 / sqrtS, 2.);

  TF1 *func_intW0 = new TF1("fit_intW0", W0_full_fun, x_lower, x_upper, 3);

  // try to set parameters, if there is any
  func_intW0 -> SetParameters(4, 5, 6); // some dummy variables
  
  double integral = func_intW0 -> Integral(x_lower, x_upper);

  /*
  cout << "theta0 = " << theta0 << ", sqrtS = " << sqrtS << " [GeV/c^2] \n"
       << "m3pi range (" << m3pi_lower << ", " << m3pi_upper << ") [MeV/c^2] \n"
       << "x range (" << x_lower << ", " << x_upper << ") \n"
    //<< "evaluate func_inW0 at m3pi = 1018.82 GeV, x = " << 1. - TMath::Power(1018.82 * 1e-3 / sqrtS, 2.) << ", func_intW0 -> Eval = " << func_intW0 -> Eval(1. - TMath::Power(1018.82 * 1e-3 / sqrtS, 2.)) << "\n" // W0 = 186.672
       << "integral = " << integral << ", Lumi_int = " << Lumi_int << " [nb^-1]\n"
       << "isrlumi = Lumi_int * integral = " << isrlumi << " [nb^-1]\n\n";
  */
  
  return integral * Lumi_int;
  
}

//
double GetISRLumi_apprx(double m3pi, double m3pi_lower, double m3pi_upper, double W0_full) {

  double isrlumi = 0.;
  double Delta_m3pi = (m3pi_upper - m3pi_lower) * 1e-3;
  
  m3pi = m3pi * 1e-3;

  		 
  isrlumi = W0_full * (2. * m3pi / TMath::Power(sqrtS, 2.)) * Lumi_int * Delta_m3pi;
  //isrlumi = W0_full * (2. * m3pi / TMath::Power(sqrtS, 2.)) * 38794.1 * Delta_m3pi;
    
  //cout << "m3pi = " << m3pi << " [GeV/c^2], width = " << Delta_m3pi << " [GeV/c^2], W0_full = " << W0_full << ", isrlumi = " << isrlumi << ", Lumi_int = " << Lumi_int << "\n";

  return isrlumi;
  
}

//
double binomial_err(double nb_true, double nb_gen) {
  double error = 0.;
  double ratio = 0.; 

  if (nb_gen != 0.) {
    ratio = nb_true / nb_gen;
    error = TMath::Sqrt(ratio * (1. - ratio) / nb_gen);
  }
   
  //cout << "true = " << nb_true << ", gen = " << nb_gen << ", ratio = " << ratio << ", error = " << error << endl;

  return error;
}

//
double GetCrx3piMax_new(double BB, double Mass_V) {

  double crx3pi = 12 * pi * BB / Mass_V / Mass_V * 1e6 * crx_scale;

  /*
  cout << "\n\tBB = " << BB << "\n"
       << "\tM_V = " << mass_omega << " [MeV]\n"
       << "\tcrx3pi = " << crx3pi << " [nb]\n";
  */
  
  return crx3pi;
  
}

void getSmearSig() {

  getObj(f_cut);

  double m3pi_true = 0., m3pi_corr = 0.;

  hsig = new TH1D("hsig", "", binsize, hmin, hmax);
  hsig -> Sumw2();

  TTree *TISR3PI_SIG = (TTree *) f_cut -> Get("TISR3PI_SIG");

  for (Int_t irow = 0; irow < TISR3PI_SIG -> GetEntries(); irow++) {// loop chain

    TISR3PI_SIG -> GetEntry(irow);

    m3pi_true = TISR3PI_SIG -> GetLeaf("Br_IM3pi_true") -> GetValue(0);

    m3pi_corr = DetectorEvent(TMath::Abs(m3pi_true));

    //cout << m3pi_corr << endl;
    
    hsig -> Fill(m3pi_corr);
    
  }
  
}

//
void scaleGSF() {

  cout << gsf << endl;
  heeg -> Scale(gsf); // scale EEG
  homegapi -> Scale(gsf); // scale OMEGAPI
  hetagam -> Scale(gsf); // scale ETAGAM
  hksl -> Scale(gsf); // scale KSL 
  hmcrest -> Scale(gsf); // scale MCREST
  hsig -> Scale(gsf); // scale signal
  
  Lumi_int = Lumi_int * gsf; // scale the total luminosity
  
}




//
void MCNorm() {

  // sig
  hsig -> Scale(sig_sfw);
  
  // eeg
  heeg = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TEEG");
  heeg -> Scale(eeg_sfw);
  heeg -> SetName("heeg");
  //heeg -> Draw();
  cout << "nb_eeg = " << heeg -> Integral(1, heeg -> GetNbinsX()) << endl;
  
  // omegapi
  homegapi = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TOMEGAPI");
  homegapi -> Scale(omegapi_sfw);
  homegapi -> SetName("homegapi");
  //homegapi -> Draw();
  cout << "nb_omegapi = " << homegapi -> Integral(1, homegapi -> GetNbinsX()) << endl;
  
  // ksl
  hksl = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TKSL");
  hksl -> Scale(ksl_sfw);
  hksl -> SetName("hksl");
  //hksl -> Draw();
  cout << "nb_ksl = " << hksl -> Integral(1, hksl -> GetNbinsX()) << endl;
  
  // etagam
  hetagam = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TETAGAM");
  hetagam -> Scale(etagam_sfw);
  hetagam -> SetName("hetagam");
  //hetagam -> Draw();
  cout << "nb_etagam = " << hetagam -> Integral(1, hetagam -> GetNbinsX()) << endl;
  
  // mcrest
  TH1D *h1d_IM3pi_TKPM = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TKPM");
  TH1D *h1d_IM3pi_TRHOPI = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TRHOPI");
  TH1D *h1d_IM3pi_TBKGREST = (TH1D *) HIM3pi_fit -> FindObject("h1d_IM3pi_TBKGREST");
  
  hmcrest = (TH1D*) h1d_IM3pi_TBKGREST -> Clone();
  hmcrest -> Add(h1d_IM3pi_TKPM, 1.);
  hmcrest -> Add(h1d_IM3pi_TRHOPI, 1.);
  hmcrest -> Scale(mcrest_sfw);
  hmcrest -> SetName("hmcrest");
  //hmcrest -> Draw();
  cout << "nb_mcrest = " << hmcrest -> Integral(1, hmcrest -> GetNbinsX()) << endl;
  
}

//
void get_efficy() {

  int binsize = hsig_true -> GetNbinsX();
  
  double hmin = hsig_true -> GetXaxis() -> GetXmin();
  double hmax = hsig_true -> GetXaxis() -> GetXmax();

  //cout << "binsize = " << binsize << ", hmin = " << hmin << ", hmax = " << hmax << endl;
  
  hefficy = new TH1D("hefficy", "", binsize, hmin, hmax);
  hefficy -> Sumw2();

  double nb_gen = 0., nb_true = 0.;
  double efficy = 0., efficy_err = 0.;
  
  for (int i = 1; i <= binsize; i ++ ) {

    nb_gen = hsig_gen -> GetBinContent(i);
    
    nb_true = hsig_true -> GetBinContent(i);
    
    efficy = nb_true / nb_gen;

    efficy_err = binomial_err(nb_true, nb_gen);
    
    if (nb_true == 0. || nb_gen == 0.) {

      efficy = 0.;
      
      efficy_err = 0.;
      
    }

    hefficy -> SetBinContent(i, efficy);
    hefficy -> SetBinError(i, efficy_err);

    //cout << "nb_true = " << nb_true << ", nb_gen = " << nb_gen << ", efficy = " << hefficy -> GetBinContent(i) << "+/-" << hefficy -> GetBinError(i) << ", efficy_err = " << efficy_err << endl;
      
    
  }
  
}

//
double Func_bw(double *x, double *par) {

  double bw = 0.;
  double m_3pi = x[0];
  double m_V = par[0]; // Vector meson mass in GeV/c^2
  double W_V = par[1]; // Vector meson mass width in GeV/c^2
  double BB = par[2]; // B_ee * B_3pi
  double sigma_max = GetCrx3piMax_new(BB, m_V);// 12. * pi * BB / (m_V * m_V); // Normalization factor, 12pi * B_ee * B_3pi / m^2_V

  // s-M^{2}_V + i * sqrt{s} * W_V
  double denomin_sqr_real = TMath::Power(m_3pi * m_3pi - m_V * m_V, 2.) + TMath::Power(m_V * W_V, 2.);
  TComplex denomin(m_3pi * m_3pi - m_V * m_V, m_V * W_V);
  TComplex denomin_sqr = denomin.Rho2();
  
  
  bw = sigma_max * (m_3pi * m_3pi) * (W_V * W_V) / denomin_sqr;
  //bw = 1. / denomin_sqr;

  /*
  cout << "\nm_3pi = " << m_3pi << " [GeV/c^2]\n"
    //<< "sqrtS = " << sqrtS << " [GeV/c^2]\n"
       << "para: m_V = " << m_V << " [GeV/c^2], W_V = " << W_V << " [GeV/c^2], BB = " << BB << "\ n"
       << "denomin = (" << denomin.Re() << ", " << denomin.Im() << "), Rho2 (check) = " << denomin.Rho2() << "(" << denomin_sqr_real << ") [GeV/c^2]\n"
       << "sigma_max = " << sigma_max << "\n" 
       << "bw = " << bw << "\n";
  */

  //cout << "sigma_max = " << sigma_max << ", m_V = " << m_V << ", BB = " << BB << endl;

  //cout << "m3pi = " << m_3pi << ", bw = " << bw << endl;

  return bw;//nb
  
}

//
TH2D *getCorrMatrix(int bins, double xmin, double xmax) {

  TTree *TISR3PI_SIG = (TTree*)f_cut -> Get("TISR3PI_SIG");

  // create corrlation matrix
  TH2D *h2d_corr = new TH2D("", "", bins, xmin, xmax, bins, xmin, xmax);
  
  //cout << "correlation matrix is created!!! " << TISR3PI_SIG -> GetName() << endl;

  double IM3pi_true = 0.;
  double IM3pi_corred = 0.;
  double PARA[5] = {frac, smallBias, smallSigma, wideBias, wideSigma};
  
  for (Int_t irow = 0; irow < TISR3PI_SIG -> GetEntries(); irow++) {// loop chain

    TISR3PI_SIG -> GetEntry(irow);
    
    IM3pi_true = TISR3PI_SIG -> GetLeaf("Br_IM3pi_true") -> GetValue(0);
    
    //IM3pi_corred = DetectorEvent(TMath::Abs(IM3pi_true));
    IM3pi_corred = DetectorEvent_fcn(TMath::Abs(IM3pi_true), PARA);
    //cout << IM3pi_corred << ", " << DetectorEvent(TMath::Abs(IM3pi_true), para) << endl;
    
    h2d_corr -> Fill(IM3pi_true, IM3pi_corred);
    
  }
  
  return h2d_corr;
  
}

//
double get_ratio(double nb_true, double nb_gen) {

  double ratio = 0.;
  
  if (nb_gen != 0.) {

    ratio = nb_true / nb_gen;
    
  }
  
  return ratio;
}

//
TH2D* getSmearMatrix(TH2D *h2d_corr, int bins, double xmin, double xmax) {

  TH2D *h2d_smear = new TH2D("", "", bins, xmin, xmax, bins, xmin, xmax);
  
  double z_value = 0.;
  double nb_true = 0.;
  double p_value = 0., p_value_err = 0.;

  int entries = 0;
  
  for (int j = 1; j <= h2d_corr -> ProjectionX() -> GetNbinsX(); j ++ ) {

    for (int i = 1; i <= h2d_corr -> ProjectionY() -> GetNbinsX(); i ++ ) {

      z_value = h2d_corr -> GetBinContent(j, i);
      nb_true = h2d_corr -> ProjectionX() -> GetBinContent(j);
      
      p_value = get_ratio(z_value, nb_true);

      if (p_value != 0.) {

	p_value_err = binomial_err(z_value, nb_true) / p_value; //Relative err

      }
      else {

	p_value_err = 0.;

      }

      h2d_smear -> SetBinContent(j, i, p_value);
      //h2d_smear_matrix_err -> SetBinContent(j, i, p_value_err);
      entries ++;

    }
	
  }
  
  return h2d_smear;
  
}

double Q_err(double pn, double epsilon, double sigma_epsilon, double mu_b, double sigma_b) {

  double q_err = 0.;
  double var_q = 0.;
  double term1 = 0., term2 = 0., term3 = 0., term4 = 0.;

  term1 = TMath::Power(pn * sigma_epsilon, 2);
  term2 = pn * (1 - pn) / mu_b;
  term3 = 1 + TMath::Power(sigma_b / mu_b, 2);
  term4 = epsilon * epsilon + sigma_epsilon * sigma_epsilon;

  var_q = term1 + term2 * term3 * term4;

  if (var_q > 0.) {
    q_err = TMath::Sqrt(var_q);
  }
  else {
    q_err = 0.;
  }

  //cout << "pn = " << pn << ", epsilon = " << epsilon << "+/-" << sigma_epsilon << ", b = "<< mu_b << "+/-" << sigma_b << "\n";
  
  return q_err;
  
}

//
TH2D *getSmearMatrixEfficy(TH2D *h2d_corr, int bins, double xmin, double xmax) {

  TH2D *h2d_smear = new TH2D("", "", bins, xmin, xmax, bins, xmin, xmax);
  
  double z_value = 0.;
  double nb_true = 0., nb_true_err = 0.;
  double nb_gen = 0.;
  double p_value = 0., p_value_err = 0.;
  //double p_value_check = 0.;
  double efficy = 0., efficy_err = 0.;

  int entries = 0;
  
  for (int j = 1; j <= h2d_corr -> ProjectionX() -> GetNbinsX(); j ++ ) {

    nb_true = h2d_corr -> ProjectionX() -> GetBinContent(j);
    nb_true_err = h2d_corr -> ProjectionX() -> GetBinError(j);
    nb_gen = hsig_gen -> GetBinContent(j);
    
    efficy = hefficy -> GetBinContent(j);
    efficy_err = hefficy -> GetBinError(j);
    
    for (int i = 1; i <= h2d_corr -> ProjectionY() -> GetNbinsX(); i ++ ) {

      z_value = h2d_corr -> GetBinContent(j, i);
      p_value = get_ratio(z_value, nb_gen);
      //p_value_check = z_value / nb_true * efficy;
      p_value_err = Q_err(p_value, efficy, efficy_err, nb_true, nb_true_err); 
      
      h2d_smear -> SetBinContent(j, i, p_value);
      entries ++;

    }

    //cout << j << ", nb_true (checked) = " << nb_true << "(" << hsig_true -> GetBinContent(j) << "+/-" << hsig_true -> GetBinError(j) << "), nb_gen = " << nb_gen << ", efficy (checked) = " << efficy << " +/- " << efficy_err << " (" << nb_true / nb_gen << ")" << endl;
	
  }
  
  return h2d_smear;
  
}

//
TH2D* getSmearMatrixErr(TH2D *h2d_corr, int bins, double xmin, double xmax) {

  TH2D *h2d = new TH2D("", "", bins, xmin, xmax, bins, xmin, xmax);
  
  double z_value = 0.;
  double nb_true = 0.;
  double p_value = 0., p_value_err = 0.;

  int entries = 0;
  
  for (int j = 1; j <= h2d_corr -> ProjectionX() -> GetNbinsX(); j ++ ) {

    for (int i = 1; i <= h2d_corr -> ProjectionY() -> GetNbinsX(); i ++ ) {

      z_value = h2d_corr -> GetBinContent(j, i);
      nb_true = h2d_corr -> ProjectionX() -> GetBinContent(j);
      
      p_value = get_ratio(z_value, nb_true);

      if (p_value != 0.) {

	//p_value_err = binomial_err(z_value, nb_true) / p_value; //Relative err
	p_value_err = binomial_err(z_value, nb_true); //Relative err

      }
      else {

	p_value_err = 0.;

      }

      //cout << z_value << ", " << nb_true << ", " << p_value << "+/-" << p_value_err << endl;
      
      h2d -> SetBinContent(j, i, p_value_err);
      entries ++;

    }
	
  }
  
  return h2d;
  
}

//
TH2D *getSmearMatrixEfficy_Err(TH2D *h2d_corr, int bins, double xmin, double xmax) {

  TH2D *h2d_smear = new TH2D("", "", bins, xmin, xmax, bins, xmin, xmax);
  
  double z_value = 0.;
  double nb_true = 0., nb_true_err = 0.;
  double nb_gen = 0.;
  double p_value = 0., p_value_err = 0.;
  //double p_value_check = 0.;
  double efficy = 0., efficy_err = 0.;

  int entries = 0;
  
  for (int j = 1; j <= h2d_corr -> ProjectionX() -> GetNbinsX(); j ++ ) {

    nb_true = h2d_corr -> ProjectionX() -> GetBinContent(j);
    nb_true_err = h2d_corr -> ProjectionX() -> GetBinError(j);
    nb_gen = hsig_gen -> GetBinContent(j);
    
    efficy = hefficy -> GetBinContent(j);
    efficy_err = hefficy -> GetBinError(j);
    
    for (int i = 1; i <= h2d_corr -> ProjectionY() -> GetNbinsX(); i ++ ) {

      z_value = h2d_corr -> GetBinContent(j, i);
      p_value = get_ratio(z_value, nb_gen);
      //p_value_check = z_value / nb_true * efficy;

      if (p_value > 0.) {
	p_value_err = Q_err(p_value, efficy, efficy_err, nb_true, nb_true_err); 
      }
      else {
	p_value_err = 0.;
      }
      
      h2d_smear -> SetBinContent(j, i, p_value_err);
      entries ++;

    }

  }
  
  return h2d_smear;
  
}


void fcn_crx3pi(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){

  //const int list_size = 1000;

  double Mass_V = par[0];
  //double Gam_V = par[1] + Delta_FWHM;
  //double BB = par[2] - Delta_BB;

  double Gam_V = par[1];
  double BB = par[2];

  //double M3PI[list_size];
  //double EFFICY[list_size];
  //double ISRLUMI[list_size];
  //double CRX3PI_BW[list_size];

  //double N3PI_TRUE[list_size];
  //double N3PI_GEN[list_size];
  //double N3PI_CORR[list_size];
  //double N3PI_PRE[list_size];
  //double N3PI_SMEAR[list_size];
  //double N3PI_OBS[list_size];
  //double N3PI_OBS_ERR[list_size];
  //TF1 *fun_bw = new TF1("func_bw", Func_bw, 0., 2000., 3);
  fun_bw -> SetParameters(Mass_V, Gam_V, BB);

  int counter = 0;
  
  for (Int_t irow = 0; irow < TRESULT -> GetEntries(); irow++) {// begin for loop

    TRESULT -> GetEntry(irow);

    M3PI[counter] = TRESULT -> GetLeaf("Br_m3pi") -> GetValue(0);
    EFFICY[counter] = TRESULT -> GetLeaf("Br_efficy") -> GetValue(0);
    EFFICY_ERR[counter] = TRESULT -> GetLeaf("Br_efficy_err") -> GetValue(0);
    ISRLUMI[counter] = TRESULT -> GetLeaf("Br_isrlumi_apprx") -> GetValue(0);
    CRX3PI_BW[counter] = fun_bw -> Eval(M3PI[counter]);

    // updated from TUFO
    N3PI_OBS[counter] = TRESULT -> GetLeaf("Br_nb_isr3pi_obs") -> GetValue(0);
    N3PI_OBS_ERR[counter] = TRESULT -> GetLeaf("Br_nb_isr3pi_obs_err") -> GetValue(0);
	       
    // predicted number of 3pi events after the acceptance correction
    // apply globle scaling factor
    N3PI_PRE[counter] = ISRLUMI[counter] * CRX3PI_BW[counter];

    //cout << "ISRLUMI = " << ISRLUMI[counter] << ", CRX3PI_BW = " << CRX3PI_BW[counter] << "\n";
    //cout << "N3PI_OBS = " << N3PI_OBS[counter] << "\n";

    //cout << "bin indx = " << irow + 1 << ", M3PI = " << M3PI[counter] << ", EFFICY [%] = " << EFFICY[counter] * 100. << " +/- " << EFFICY_ERR[counter] * 100. << ", ISRLUMI = " << ISRLUMI[counter] << ", CRX3PI_BW = " << CRX3PI_BW[counter] << ", N3PI_OBS = " << N3PI_OBS[counter] << " +/- " << N3PI_OBS_ERR[counter] << ", N3PI_PRE = " << N3PI_PRE[counter] << "\n";
    
    counter ++;
    
  }

  //cout << "\n"; 

  
  
  // smearing
  const int bin_indx = 119;
    
  for (int j = 1; j <= hsmearmatr -> ProjectionY() -> GetNbinsX(); j ++ ) {

    double X_sum = 0., X_mean = 0.;
    double w = 0., wX_sum = 0., w_sum = 0.;
    double err2_sum = 0.;
    for (int k = 1; k <= binsize; k ++ ) {

      X_sum += hsmearmatr_efficy -> GetBinContent(k, j) * N3PI_PRE[k - 1];
      err2_sum = get_sigma2_y(hsmearmatr_efficy -> GetBinContent(k, j), N3PI_PRE[k - 1], hsmearmatr_efficy_err -> GetBinContent(k, j), 0.);
      //X_sum += hsmearmatr_efficy -> GetBinContent(k, j);
      if (err2_sum !=0){
	w = 1 / err2_sum;
      }
      else {
	w = 0.;
      }
      wX_sum += w * hsmearmatr_efficy -> GetBinContent(k, j) * N3PI_PRE[k - 1];
      w_sum += w;
    }

    X_mean = wX_sum / w_sum;
    
    double prob_tmp = 0., prob_err_tmp = 0.;
    double prob_efficy_tmp = 0., prob_efficy_err_tmp = 0.;
    double nb_sig = 0., nb_sig_err = 0.;
    double nb_sig_true = 0., nb_sig_true_err = 0.;
    double nb_sig_gen = 0., nb_sig_gen_err = 0.;
    double nb_sig_smeared = 0.; // use semarmatr
    double nb_sig_smeared_efficy = 0., nb_sig_smeared_efficy_err = 0.; // use smearmatr_efficy
    
    double nb_rec_tmp = 0.;
    double nb_rec_tmp_err2 = 0., nb_rec_tmp_err2_sum = 0.;
    double efficy = 0.;
    double sigma2 = 0., sigma2_sum = 0.;
    double var_X = 0., varsum_X = 0.;
    double s2 = 0.;
    
    for (int i = 1; i <= hsmearmatr -> ProjectionX() -> GetNbinsX(); i ++ ) {

      prob_tmp = hsmearmatr -> GetBinContent(i, j);
      prob_err_tmp = hsmearmatr_err -> GetBinError(i, j);

      prob_efficy_tmp = hsmearmatr_efficy -> GetBinContent(i, j);
      prob_efficy_err_tmp = hsmearmatr_efficy_err -> GetBinContent(i, j);
      
      efficy = hefficy -> GetBinContent(i);
      
      nb_sig_true = hsig_true -> GetBinContent(i);
      nb_sig_true_err = hsig_true -> GetBinError(i);
	
      nb_sig_gen = hsig_gen -> GetBinContent(i);
      nb_sig_gen_err = hsig_gen -> GetBinError(i);

      nb_sig_smeared += prob_tmp * nb_sig_true;
      nb_sig_smeared_efficy += prob_efficy_tmp * nb_sig_gen;

      nb_rec_tmp += prob_efficy_tmp * N3PI_PRE[i - 1];
      //nb_rec_tmp_err2_sum += nb_gen * prob_tmp * (1 - prob_tmp) * (N3PI_PRE[i - 1] * N3PI_PRE[i - 1]); 
	
      //cout << "(rec_indx, true_indx) = (" << j << ", " << i << "), N3pi_true = " << N3PI_TRUE[i - 1] << ", prob_tmp = " << prob_tmp << ", nb_rec_tmp = " << nb_rec_tmp << ", N3pi_corr = " << N3PI_CORR[i - 1] << "\n";

      //cout << "(rec_indx, true_indx) = (" << j << ", " << i << "), prob_tmp = " << prob_tmp << ", N3PI_PRE = " << N3PI_PRE[i - 1] << ", nb_rec_tmp = " << nb_rec_tmp << "\n";
      
      //cout << "(rec_indx, true_indx) = (" << j << ", " << i << "), prob_tmp = " << prob_tmp << ", N3PI_PRE = " << N3PI_PRE[i - 1] << ", nb_rec_tmp = " << nb_rec_tmp << "\n";

      if (prob_efficy_tmp != 0.) {
	sigma2 = get_sigma2_y(prob_efficy_tmp, nb_sig_gen, prob_efficy_err_tmp, nb_sig_gen_err);
	nb_rec_tmp_err2 = get_sigma2_y(prob_efficy_tmp, N3PI_PRE[i - 1], prob_efficy_err_tmp, 0.);
	w = 1 / nb_rec_tmp_err2;
	var_X = w * (prob_efficy_tmp * N3PI_PRE[i - 1] - X_mean) * (prob_efficy_tmp * N3PI_PRE[i - 1] - X_mean);
      }
      else {
	sigma2 = 0.;
	nb_rec_tmp_err2 = 0.;
	w = 0.;
	var_X = 0.;
      }

      // sum error2
      sigma2_sum += sigma2;
      s2 += var_X;
      w_sum += w;
      nb_rec_tmp_err2_sum += s2 / (w_sum - 1);

      if (j == bin_indx) {
	//cout << "(rec_indx, true_indx) = (" << j << ", " << i << "), nb_sig_true = " << nb_sig_true << ", prob_tmp = " << prob_tmp << "+/-"<< prob_err_tmp << ", nb_sig_gen = " << nb_sig_gen << ", prob_efficy_tmp = " << prob_efficy_tmp << "+/-" << prob_efficy_err_tmp << ", N3PI_PRE = " << N3PI_PRE[i - 1] << ", nb_rec_tmp_err2_sum = " <<  nb_rec_tmp_err2_sum << endl;
      }
      
    }

    N3PI_SMEAR[j -1] = nb_rec_tmp;
    N3PI_SMEAR_ERR[j -1] = TMath::Sqrt(nb_rec_tmp_err2_sum);

    //cout << N3PI_SMEAR_ERR[j -1] << ", " << TMath::Sqrt(nb_rec_tmp) << endl;
    
    
    //cout << "bin indx = " << j << ", N3PI_SMEAR = " << N3PI_SMEAR[j - 1] << " +/- " << N3PI_SMEAR_ERR[j - 1] << "\n";

    if (j == bin_indx) {
      //cout << "rec bin = " << j << ", m3pi = " << hsig_true -> GetBinCenter(j) << ", nb_sig = " << nb_sig << "+/-" << nb_sig_err << ", nb_sig_smeared = " << nb_sig_smeared << ", nb_sig_smeared_efficy = " << nb_sig_smeared_efficy << "+/-" << nb_sig_smeared_efficy_err <<  ", nb_sig_true = " << hsig_true -> GetBinContent(j) << ", nb_sig_gen = " << hsig_gen -> GetBinContent(j) << ", N3PI_PRE = " << N3PI_PRE[j - 1] << ", N3PI_SMEAR = " << N3PI_SMEAR[j -1] << "+/-" << N3PI_SMEAR_ERR[j -1] << "\n\n";
      //cout << "rec bin = " << j << ", m3pi = " << hsig_true -> GetBinCenter(j) << ", nb_sig_smeared = " << nb_sig_smeared <<  ", nb_sig_true = " << hsig_true -> GetBinContent(j) << ", nb_sig_gen = " << hsig_gen -> GetBinContent(j) << ", N3PI_PRE = " << N3PI_PRE[j - 1] << ", N3PI_SMEAR = " << N3PI_SMEAR[j -1] << "+/-" << N3PI_SMEAR_ERR[j -1] << "\n\n";
    }
    
    //cout << "\n"; 

    
  }

  //cout << "\n"; 

  // mass interval
  
      
  //
  double diff_sqr = 0., err_sqr = 0.;
  double chi2_sum_tmp = 0., chi2_tmp = 0., chi2_tmp_err = 0.;
  
  int fit_indx = 0;
  int range_indx = 0, sub_indx = 0;

  for (int i = 0; i < counter; i ++) {
    
    if (M3PI[i] > fit_min && M3PI[i] < fit_max) {
      
      /*
	if (M3PI[i] > M3PI_LOWER[range_indx] && range_indx < mass_steps) {
	
	range_indx ++;
	
	sub_indx = 0;
	
	}
	
	if (M3PI[i] > M3PI_LOWER[range_indx - 1] && M3PI[i] < M3PI_UPPER[range_indx - 1]) {
	
	sub_indx ++;
	
	}
      */
      
      //cout << sub_indx << ", " << M3PI[i] << ", " << M3PI_LOWER[range_indx - 1] << ", " << M3PI_UPPER[range_indx - 1] << endl;
      
      //cout << h1d_IM3pi_TISR3PI_SIG_TRUE -> GetBinContent(i) << endl;
      
      diff_sqr = TMath::Abs((N3PI_SMEAR[i] - N3PI_OBS[i]) * (N3PI_SMEAR[i] - N3PI_OBS[i]));
      
      err_sqr = N3PI_OBS_ERR[i] * N3PI_OBS_ERR[i] + N3PI_SMEAR_ERR[i] * N3PI_SMEAR_ERR[i];
      
      chi2_tmp = diff_sqr / err_sqr;
      //chi2_tmp_err = 0.1 * chi2_tmp;
      
      chi2_sum_tmp += chi2_tmp;
      
      BIN_INDX[fit_indx] = i;
      
      N3PI_DIFF_ERR[i] = TMath::Sqrt(N3PI_OBS_ERR[i] * N3PI_OBS_ERR[i] + N3PI_SMEAR_ERR[i] * N3PI_SMEAR_ERR[i]);
      N3PI_DIFF[i] = N3PI_OBS[i] - N3PI_SMEAR[i];

      //cout << i << ", " << M3PI[i] <<  BIN_INDX[fit_indx] << ", fit_indx = " << endl;

      /*
      cout << "\t/////////////////////////////////\n"
	   << "\tbin " << i << "\n"
	   << "\tMass_V = " << Mass_V << ", Gam_V = " << Gam_V << ", BB = " << BB << "\n" 
	   << "\tm3pi = "   << M3PI[i] << " [MeV/c^2]\n"
	   << "\tefficy [%] = " << EFFICY[i] * 100. << " +/- " << EFFICY_ERR[i] << "\n"
	   << "\tisrlumi = " << ISRLUMI[i] << "[nb^-1]\n"
	   << "\tcrx3pi_bw = " << CRX3PI_BW[i] << "\n"
	//<< "\tnb_isr3pi_gen = " << N3PI_TRUE[i] << ", efficy = " << EFFICY[i] << ", gen * efficy (true) = " << N3PI_GEN[i] * EFFICY[i] << "(" << N3PI_TRUE[i] << ")\n"
	//<< "\tnb_isr3pi_corr = " << N3PI_CORR[i] << "\n"
	   << "\tnb_isr3pi_obs = " << N3PI_OBS[i] << " +/- " << N3PI_OBS_ERR[i] << "\n"
	   << "\tnb_isr3pi_fit = " << N3PI_SMEAR[i] << " +/- " << N3PI_SMEAR_ERR[i] << "\n"
	   << "\tn3pi_diff = obs - fit = " << N3PI_DIFF[i] << " +/- " << N3PI_DIFF_ERR[i] << "\n"
	   << "\tchi2_sum_tmp = " << chi2_sum_tmp << "\n";
      */
      //cout << "bin " << i << "\n";
      
      fit_indx ++;
      
    }

    //cout << chi2_sum_tmp << endl;

      
  }

  //cout << "\n";
  
  chi2_sum_crx3pi = chi2_sum_tmp,
  residul_size_crx3pi = fit_indx;


  f = chi2_sum_tmp;
 
  
}

//
double GetBB_new(double sigma_max, double Mass_V) {

  double BB = sigma_max * 1e-6 / crx_scale * Mass_V * Mass_V / 12. / pi;

  return BB;

}


//
double GetBkgErr(double NData, double nMC, double NMC, double f, double f_err) {

  double r = nMC / NMC;
  double r_err = 0.; 
  double r_ratio = 0.;
  
  if (nMC != 0.) {
    r = nMC / NMC;
    r_err = binomial_err(nMC, NMC);
    r_ratio = r_err / r;
  }
  
  double f_ratio = f_err / f;
  
  double n_scaled = NData * r * f;

  //cout << "n_scaled = " << n_scaled << endl;
  
  double n_scaled_err = n_scaled * TMath::Sqrt(r_ratio * r_ratio + f_ratio * f_ratio);
  
  return n_scaled_err;
  
}


  
