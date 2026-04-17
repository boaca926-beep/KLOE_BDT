const double fit_min = 770;
const double fit_max = 800;
const int fit_bin = 60; 

double fit_indx = 0;

TTree *TISR3PI_SIG = nullptr;
TTree *TSFW = new TTree("TSFW", "recreate");

void fcn_sfw(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){

  double nb_data = 0.;
  double nb_isr3pi = 0., isr3pi_mu = 0.;
  double nb_bkgsum_sc = 0.;
  
  double mu_tmp = 0.;
  double chi2_sum_tmp = 0.;

  double counter = 0;

  double sfw_isr3pi = par[0];

  for (Int_t irow = 0; irow < TSFW -> GetEntries(); irow++) {// fcn loop
    
    TSFW -> GetEntry(irow);

    // data
    nb_data = TSFW -> GetLeaf("Br_nb_data") -> GetValue(0);

    // eeg
    //nb_eeg_sc = TSFW -> GetLeaf("Br_nb_eeg_sc") -> GetValue(0);

    // ksl
    //nb_ksl_sc = TSFW -> GetLeaf("Br_nb_ksl_sc") -> GetValue(0);

    // omegapi
    //nb_omegapi_sc = TSFW -> GetLeaf("Br_nb_omegapi_sc") -> GetValue(0);

    // isr3pi
    nb_isr3pi = TSFW -> GetLeaf("Br_nb_isr3pi") -> GetValue(0);
    isr3pi_mu = nb_isr3pi * sfw_isr3pi;

    // etagam
    //nb_etagam_sc = TSFW -> GetLeaf("Br_nb_etagam_sc") -> GetValue(0);

    // mcrest
    //nb_mcrest_sc = TSFW -> GetLeaf("Br_nb_mcrest_sc") -> GetValue(0);

    // bkgsum
    nb_bkgsum_sc = TSFW -> GetLeaf("Br_nb_bkgsum_sc") -> GetValue(0);

    // mu
    mu_tmp = nb_bkgsum_sc + isr3pi_mu;
    //mu_tmp = nb_eeg_sc + nb_ksl_sc + nb_omegapi_sc + isr3pi_mu + nb_etagam_sc + nb_mcrest_sc;

    if (mu_tmp > 0. && nb_data > 0.) {

      chi2_sum_tmp += (nb_data - mu_tmp) * (nb_data - mu_tmp) / (nb_data + mu_tmp);

      counter ++;

    }

    cout << irow + 1 << "\n"
	 << "nb_data = " << nb_data << ", mu = " << mu_tmp << "\n"
      // << "nb_eeg_sc = " << nb_eeg_sc << "\n"
      // << "nb_ksl_sc = " << nb_ksl_sc << "\n"
      // << "nb_omegapi_sc = " << nb_omegapi_sc << "\n"
	 << "nb_isr3pi = " << nb_isr3pi << ", sfw_isr3pi = " << sfw_isr3pi << ", mu = " << isr3pi_mu << "\n"
      // << "nb_etagam_sc = " << nb_etagam_sc << "\n"
      // << "nb_mcrest_sc = " << nb_mcrest_sc << "\n"
	 << "nb_bkgsum_sc = " << nb_bkgsum_sc << "\n";
    
  }

  fit_indx = counter;
  f = chi2_sum_tmp;

}
