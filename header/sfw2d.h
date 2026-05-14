double nb_data_sum = 0.;
double nb_eeg_sum = 0.;
double nb_ksl_sum = 0.;
double nb_omegapi_sum = 0.;
double nb_etagam_sum = 0.;
double nb_isr3pi_sum = 0.;
double nb_mcrest_sum = 0.;

double chi2_sfw2d_sum = 0., residul_size_sfw2d = 0.;

TTree* TSFW2D = new TTree("TSFW2D", "recreate");

//
double GetScalError(double N_d, double N, double f, double f_error) {
  double error = 0.;
  double scale = N_d * f / N;
  error = scale * TMath::Sqrt(1. / N_d + 1. / N + TMath::Power(f_error / f, 2));

  //cout << "D = " << N_d << ", N = "<< N << ", f = " << f << "+/-" << f_error << ", scale = " << scale << " +/- " << error << endl;

  return error;

}

//
double getloglh(double n_d, double mu) {
  double value = 0.;

  //value = n_d * TMath::Log(mu) - mu + n_d - n_d * TMath::Log(n_d);
  //value = n_d * TMath::Log(mu) - mu - (n_d * TMath::Log(n_d) - n_d + 0.5 * TMath::Log(n_d) + 0.5 * TMath::Log(2 * TMath::Pi()));
  value = n_d * TMath::Log(mu) - mu;
  //cout << "n_d = " << n_d << endl;

  return value;
}

//
double getscale(double Nd, double fra, double N){
  double value = 0.;

  value =  Nd * fra / N;

  return value;
}

//
void fcn_sfw2d(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){


  int counter = 0;

  double nb_data = 0.;
  double nb_eeg = 0.;
  double nb_ksl = 0.;
  double nb_omegapi = 0.;
  double nb_etagam = 0.;
  double nb_isr3pi = 0.;
  double nb_mcrest = 0.;

  double f_eeg = par[0];
  double f_isr3pi = par[1];
  double f_omegapi = par[2];
  double f_etagam = par[3];
  double f_ksl = par[4];
  double f_mcrest = par[5];

  double eeg_sfw = 0., eeg_mu = 0.;
  double isr3pi_sfw = 0., isr3pi_mu = 0.;
  double omegapi_sfw = 0., omegapi_mu = 0.;
  double etagam_sfw = 0., etagam_mu = 0.;
  double ksl_sfw = 0., ksl_mu = 0.;
  double mcrest_sfw = 0., mcrest_mu = 0.;

  double mu_tmp = 0.;

  double chi2_sum_tmp = 0.;
  double llh_sum = 0.;
    
  //cout << "nb_data_sum = " << nb_data_sum << endl;
  
  for (Int_t irow = 0; irow < TSFW2D -> GetEntries(); irow++) {// fcn loop
    
    TSFW2D -> GetEntry(irow);

    // data
    nb_data = TSFW2D -> GetLeaf("Br_nb_data") -> GetValue(0);

    // eeg
    nb_eeg = TSFW2D -> GetLeaf("Br_nb_eeg") -> GetValue(0);
    eeg_sfw = getscale(nb_data_sum, f_eeg, nb_eeg_sum);
    eeg_mu = nb_eeg * eeg_sfw;

    // ksl
    nb_ksl = TSFW2D -> GetLeaf("Br_nb_ksl") -> GetValue(0);
    ksl_sfw = getscale(nb_data_sum, f_ksl, nb_ksl_sum);
    ksl_mu = nb_ksl * ksl_sfw;

    // omegapi
    nb_omegapi = TSFW2D -> GetLeaf("Br_nb_omegapi") -> GetValue(0);
    omegapi_sfw = getscale(nb_data_sum, f_omegapi, nb_omegapi_sum);
    omegapi_mu = nb_omegapi * omegapi_sfw;

    // isr3pi
    nb_isr3pi = TSFW2D -> GetLeaf("Br_nb_isr3pi") -> GetValue(0);
    isr3pi_sfw = getscale(nb_data_sum, f_isr3pi, nb_isr3pi_sum);
    isr3pi_mu = nb_isr3pi * isr3pi_sfw;
    
    // etagam
    nb_etagam = TSFW2D -> GetLeaf("Br_nb_etagam") -> GetValue(0);
    etagam_sfw = getscale(nb_data_sum, f_etagam, nb_etagam_sum);
    etagam_mu = nb_etagam * etagam_sfw;

    // mcrest
    nb_mcrest = TSFW2D -> GetLeaf("Br_nb_mcrest") -> GetValue(0);
    mcrest_sfw = getscale(nb_data_sum, f_mcrest, nb_mcrest_sum);
    mcrest_mu = nb_mcrest * mcrest_sfw;

    // mu
    mu_tmp = eeg_mu + ksl_mu + etagam_mu + omegapi_mu + mcrest_mu + isr3pi_mu;

    if (mu_tmp > 0. && nb_data > 0.) {
    
      //
      /*
      if (irow == 55) { 
	cout << "indx = " << irow << "\n"
	     << "nb_data(sum) = " << nb_data << "(" << nb_data_sum << ")\n"
	     << "nb_eeg(sum) = " << nb_eeg << "(" << nb_eeg_sum << "), f = " << f_eeg << ", sfw = " << eeg_sfw << ", mu = " << eeg_mu << "\n"
	     << "nb_ksl(sum) = " << nb_ksl << "(" << nb_ksl_sum << "), f = " << f_ksl << ", sfw = " << ksl_sfw << ", mu = " << ksl_mu << "\n"
	     << "nb_omegapi(sum) = " << nb_omegapi << "(" << nb_omegapi_sum << "), f = " << f_omegapi << ", sfw = " << omegapi_sfw << ", mu = " << omegapi_mu << "\n"
	     << "nb_isr3pi(sum) = " << nb_isr3pi << "(" << nb_isr3pi_sum << "), f_isr3pi = " << f_isr3pi << ", sfw = " << isr3pi_sfw << ", mu = " << isr3pi_mu << "\n"
	     << "nb_etagam(sum) = " << nb_etagam << "(" << nb_etagam_sum << "), f_etagam = " << f_etagam << ", sfw = " << etagam_sfw << ", mu = " << etagam_mu << "\n"
	     << "nb_mcrest(sum) = " << nb_mcrest << "(" << nb_mcrest_sum << "), f_mcrest = " << f_mcrest << ", sfw = " << mcrest_sfw << ", mu = " << mcrest_mu << "\n\n";
	
      }
      */
      
      chi2_sum_tmp += (nb_data - mu_tmp) * (nb_data - mu_tmp) / (nb_data + mu_tmp);
      //chi2_sum_tmp += (nb_data - mu_tmp) * (nb_data - mu_tmp) / nb_data;

      //cout << counter << endl;
      
      counter ++;
      
      llh_sum -= 2. * getloglh(nb_data, mu_tmp);
      
    }

  }// end fcn loop 

  chi2_sfw2d_sum = chi2_sum_tmp;
  residul_size_sfw2d = counter;
   
  f = llh_sum;
  //f = chi2_sfw2d_sum;
    
}
