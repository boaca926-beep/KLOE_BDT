#include "../header_bdt/sm_para.h"
#include "../header_bdt/path.h"
#include "../header_bdt/method.h"
#include "../header_bdt/cut_para.h"

#include <TMath.h>
#include <TTree.h>
#include <TMinuit.h>  // Add this

// ----------------------------------------------------------------------
// Global sums (declarations only – defined in sfw2d.C)
// ----------------------------------------------------------------------
double nb_data_sum = 0.;
double nb_eeg_sum = 0.;
double nb_ksl_sum = 0.;
double nb_omegapi_sum = 0.;
double nb_nonReson_sum = 0.;
double nb_isr3pi_sum = 0.;
double nb_mcrest_sum = 0.;

double chi2_sfw2d_sum = 0.;
double residul_size_sfw2d = 0.;

TTree* TSFW2D = new TTree("TSFW2D", "recreate");

// ----------------------------------------------------------------------
// Helper functions (inline to avoid multiple definitions)
// ----------------------------------------------------------------------
inline double GetScalError(double N_d, double N, double f, double f_error) {
  if (N == 0.0) return 0.0;
  double scale = N_d * f / N;
  double error = scale * TMath::Sqrt(1.0/N_d + 1.0/N + TMath::Power(f_error/f, 2));
  return error;
}

inline double getloglh(double n_d, double mu) {
  if (mu <= 0.0) return -1e9;   // invalid, minimiser will avoid
  return n_d * TMath::Log(mu) - mu;
}

inline double getscale(double Nd, double fra, double N) {
  if (N == 0.0) return 0.0;
  return Nd * fra / N;
}

// ----------------------------------------------------------------------
// Fit function for TMinuit
// ----------------------------------------------------------------------
void fcn_sfw2d(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
  // Static branch addresses – set only once
  static bool first = true;
  static Double_t nb_data, nb_eeg, nb_ksl, nb_omegapi, nb_nonReson, nb_isr3pi, nb_mcrest;

  if (first) {
    TSFW2D->SetBranchAddress("Br_nb_data", &nb_data);
    TSFW2D->SetBranchAddress("Br_nb_eeg", &nb_eeg);
    TSFW2D->SetBranchAddress("Br_nb_ksl", &nb_ksl);
    TSFW2D->SetBranchAddress("Br_nb_omegapi", &nb_omegapi);
    TSFW2D->SetBranchAddress("Br_nb_nonReson", &nb_nonReson);
    TSFW2D->SetBranchAddress("Br_nb_isr3pi", &nb_isr3pi); // omega peak
    //add non-reson background
    TSFW2D->SetBranchAddress("Br_nb_mcrest", &nb_mcrest); // include etagam
    first = false;
  }

  double f_eeg     = par[0];
  double f_isr3pi  = par[1]; // omegapeak
  double f_omegapi = par[2];
  double f_nonReson  = par[3]; 
  double f_ksl     = par[4];
  double f_mcrest  = par[5];

  double eeg_sfw = 0., eeg_mu = 0.;
  double isr3pi_sfw = 0., isr3pi_mu = 0.;
  double omegapi_sfw = 0., omegapi_mu = 0.;
  double nonReson_sfw = 0., nonReson_mu = 0.; 
  double ksl_sfw = 0., ksl_mu = 0.;
  double mcrest_sfw = 0., mcrest_mu = 0.;

  double mu_tmp = 0.;
  double chi2_sum_tmp = 0.;
  double llh_sum = 0.;
  Int_t counter = 0;

  for (Int_t irow = 0; irow < TSFW2D->GetEntries(); irow++) {
    TSFW2D->GetEntry(irow);   // fills all branch variables

    // Scale each component using the current fractions (par)
    eeg_sfw     = getscale(nb_data_sum, f_eeg,     nb_eeg_sum);
    eeg_mu      = nb_eeg * eeg_sfw;

    ksl_sfw     = getscale(nb_data_sum, f_ksl,     nb_ksl_sum);
    ksl_mu      = nb_ksl * ksl_sfw;

    omegapi_sfw = getscale(nb_data_sum, f_omegapi, nb_omegapi_sum);
    omegapi_mu  = nb_omegapi * omegapi_sfw;

    isr3pi_sfw  = getscale(nb_data_sum, f_isr3pi,  nb_isr3pi_sum);
    isr3pi_mu   = nb_isr3pi * isr3pi_sfw;

    nonReson_sfw  = getscale(nb_data_sum, f_nonReson,  nb_nonReson_sum);
    nonReson_mu   = nb_nonReson * nonReson_sfw;
    
    mcrest_sfw  = getscale(nb_data_sum, f_mcrest,  nb_mcrest_sum);
    mcrest_mu   = nb_mcrest * mcrest_sfw;

    mu_tmp = eeg_mu + ksl_mu + nonReson_mu + omegapi_mu + mcrest_mu + isr3pi_mu;
   
    if (mu_tmp > 0.0 && nb_data > 0.0) {
      chi2_sum_tmp += (nb_data - mu_tmp) * (nb_data - mu_tmp) / (nb_data + mu_tmp);
      counter++;
      llh_sum -= 2.0 * getloglh(nb_data, mu_tmp);
    }
  }

  chi2_sfw2d_sum = chi2_sum_tmp;
  residul_size_sfw2d = counter;
  f = llh_sum;   // minimise -2 ln L
}

// ----------------------------------------------------------------------
// Main function
// ----------------------------------------------------------------------
int sfw2d_noeta() {
  // Open input and output files
  TFile *f_input  = TFile::Open(outputHist + "hist.root");
  if (!f_input || f_input->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << outputHist << "hist.root" << std::endl;
    return 1;
  }
  TFile *f_output = TFile::Open(outputSfw2D + "sfw2d.root", "recreate");
  if (!f_output || f_output->IsZombie()) {
    std::cerr << "ERROR: Cannot create " << outputSfw2D << "sfw2d.root" << std::endl;
    return 1;
  }

  std::cout << "cut_value = " << cut_value << std::endl;

  // Retrieve the list of 2D histograms
  TList *HSFW2D = (TList*) f_input->Get("HSFW2D");
  if (!HSFW2D) {
    std::cerr << "ERROR: HSFW2D not found in input file" << std::endl;
    return 1;
  }
  checkList(HSFW2D);

  // Get all required histograms
  TH2D *h2d_sfw_TEEG        = (TH2D*) HSFW2D->FindObject("h2d_sfw_TEEG");
  TH2D *h2d_sfw_TDATA       = (TH2D*) HSFW2D->FindObject("h2d_sfw_TDATA");
  TH2D *h2d_sfw_TISR3PI_SIG_peak = (TH2D*) HSFW2D->FindObject("h2d_sfw_TISR3PI_SIG_peak");
  TH2D *h2d_sfw_TISR3PI_SIG_non_reson = (TH2D*) HSFW2D->FindObject("h2d_sfw_TISR3PI_SIG_non_reson");
  TH2D *h2d_sfw_TOMEGAPI    = (TH2D*) HSFW2D->FindObject("h2d_sfw_TOMEGAPI");
  TH2D *h2d_sfw_TKPM        = (TH2D*) HSFW2D->FindObject("h2d_sfw_TKPM");
  TH2D *h2d_sfw_TKSL        = (TH2D*) HSFW2D->FindObject("h2d_sfw_TKSL");
  TH2D *h2d_sfw_T3PIGAM     = (TH2D*) HSFW2D->FindObject("h2d_sfw_T3PIGAM");
  TH2D *h2d_sfw_TRHOPI      = (TH2D*) HSFW2D->FindObject("h2d_sfw_TRHOPI");
  TH2D *h2d_sfw_TETAGAM     = (TH2D*) HSFW2D->FindObject("h2d_sfw_TETAGAM");
  TH2D *h2d_sfw_TBKGREST    = (TH2D*) HSFW2D->FindObject("h2d_sfw_TBKGREST");

  // Verify that critical histograms exist
  if (!h2d_sfw_TDATA || !h2d_sfw_TEEG || !h2d_sfw_TISR3PI_SIG_peak || !h2d_sfw_TISR3PI_SIG_non_reson || !h2d_sfw_TOMEGAPI || !h2d_sfw_TKPM || !h2d_sfw_TKSL || !h2d_sfw_T3PIGAM || !h2d_sfw_TRHOPI || !h2d_sfw_TETAGAM || !h2d_sfw_TBKGREST) {
    std::cerr << "ERROR: Missing required histograms (TDATA or TEEG)" << std::endl;
    return 1;
  }

  // Build combined histograms
  TH2D *h2d_sfw_MCREST = (TH2D*) h2d_sfw_TBKGREST->Clone();
  h2d_sfw_MCREST->Add(h2d_sfw_TKPM, 1.0);
  h2d_sfw_MCREST->Add(h2d_sfw_TRHOPI, 1.0);
  h2d_sfw_MCREST->Add(h2d_sfw_TETAGAM, 1.0);
  h2d_sfw_MCREST->SetName("h2d_sfw_MCREST");

  TH2D *h2d_sfw_MCSUM = (TH2D*) h2d_sfw_TEEG->Clone();
  h2d_sfw_MCSUM->Add(h2d_sfw_TOMEGAPI, 1.0);
  h2d_sfw_MCSUM->Add(h2d_sfw_TKSL, 1.0);
  h2d_sfw_MCSUM->Add(h2d_sfw_TISR3PI_SIG_peak, 1.0);
  h2d_sfw_MCSUM->Add(h2d_sfw_TISR3PI_SIG_non_reson, 1.0);
  h2d_sfw_MCSUM->Add(h2d_sfw_MCREST, 1.0);
  h2d_sfw_MCSUM->SetName("h2d_sfw_MCSUM");

  // Variables for sums and later fitting
  double nb_data     = 0., nb_eeg     = 0., nb_omegapi = 0., nb_ksl = 0.;
  double nb_isr3pi   = 0., nb_nonReson  = 0., nb_mcrest  = 0.;

  nb_data_sum = nb_eeg_sum = nb_omegapi_sum = nb_ksl_sum = 0.;
  nb_nonReson_sum = nb_isr3pi_sum = nb_mcrest_sum = 0.;

  // Fill TSFW2D tree (one entry per 2D bin)
  TSFW2D->SetAutoSave(0);
  TSFW2D->Branch("Br_nb_data",    &nb_data,    "Br_nb_data/D");
  TSFW2D->Branch("Br_nb_eeg",     &nb_eeg,     "Br_nb_eeg/D");
  TSFW2D->Branch("Br_nb_ksl",     &nb_ksl,     "Br_nb_ksl/D");
  TSFW2D->Branch("Br_nb_omegapi", &nb_omegapi, "Br_nb_omegapi/D");
  TSFW2D->Branch("Br_nb_nonReson",  &nb_nonReson,  "Br_nb_nonReson/D");
  TSFW2D->Branch("Br_nb_isr3pi",  &nb_isr3pi,  "Br_nb_isr3pi/D");
  TSFW2D->Branch("Br_nb_mcrest",  &nb_mcrest,  "Br_nb_mcrest/D");

  Int_t nx = h2d_sfw_TDATA->GetXaxis()->GetNbins();
  Int_t ny = h2d_sfw_TDATA->GetYaxis()->GetNbins();

  for (Int_t i = 1; i <= nx; i++) {
    for (Int_t j = 1; j <= ny; j++) {
      nb_data    = h2d_sfw_TDATA->GetBinContent(i, j);
      nb_eeg     = h2d_sfw_TEEG->GetBinContent(i, j);
      nb_omegapi = h2d_sfw_TOMEGAPI->GetBinContent(i, j);
      nb_ksl     = h2d_sfw_TKSL->GetBinContent(i, j);
      nb_nonReson = h2d_sfw_TISR3PI_SIG_non_reson->GetBinContent(i, j);
      nb_isr3pi  = h2d_sfw_TISR3PI_SIG_peak->GetBinContent(i, j);
      nb_mcrest  = h2d_sfw_MCREST->GetBinContent(i, j);

      // Update sums
      nb_data_sum    += nb_data;
      nb_eeg_sum     += nb_eeg;
      nb_omegapi_sum += nb_omegapi;
      nb_ksl_sum     += nb_ksl;
      nb_nonReson_sum  += nb_nonReson;
      nb_isr3pi_sum  += nb_isr3pi;
      nb_mcrest_sum  += nb_mcrest;

      TSFW2D->Fill();
    }
  }

  double nb_mcsum = nb_eeg_sum + nb_omegapi_sum + nb_ksl_sum + nb_nonReson_sum + nb_isr3pi_sum + nb_mcrest_sum;

  // Write sums to a text file (optional, kept from original)
  std::ofstream myfile;
  TString myfile_nm = "../header_bdt/sfw2d.txt";
  myfile.open(myfile_nm.Data());
  myfile << "const double nb_data_sum = " << nb_data_sum << ";\n"
         << "const double nb_eeg_sum = " << nb_eeg_sum << ";\n"
         << "const double nb_omegapi_sum = " << nb_omegapi_sum << ";\n"
         << "const double nb_ksl_sum = " << nb_ksl_sum << ";\n"
         << "const double nb_nonReson_sum = " << nb_nonReson_sum << ";\n"
         << "const double nb_isr3pi_sum = " << nb_isr3pi_sum << ";\n"
         << "const double nb_mcrest_sum = " << nb_mcrest_sum << ";\n"
         << "const double nb_mcsum = " << nb_mcsum << ";\n\n";
  myfile.close();

  // ----- Fitting using TMinuit -----
  const Int_t npar = 6;
  TMinuit *gMinuit = new TMinuit(npar);
  gMinuit->SetFCN(fcn_sfw2d);

  // Initial fractions (from sums)
  double feeg_init     = nb_eeg_sum     / nb_mcsum;
  double fisr3pi_init  = nb_isr3pi_sum  / nb_mcsum;
  double fomegapi_init = nb_omegapi_sum / nb_mcsum;
  double fnonReson_init  = nb_nonReson_sum  / nb_mcsum;
  double fksl_init     = nb_ksl_sum     / nb_mcsum;
  double fmcrest_init  = nb_mcrest_sum  / nb_mcsum;

  Int_t ierflg = 0;
  gMinuit->mnparm(0, "feeg_ML",     feeg_init,     0.01, 0.0, 1.0, ierflg);
  gMinuit->mnparm(1, "fisr3pi_ML",  fisr3pi_init,  0.01, 0.0, 1.0, ierflg);
  gMinuit->mnparm(2, "fomegapi_ML", fomegapi_init, 0.01, 0.0, 1.0, ierflg);
  gMinuit->mnparm(3, "fnonReson_ML",  fnonReson_init,  0.01, 0.0, 1.0, ierflg);
  gMinuit->mnparm(4, "fksl_ML",     fksl_init,     0.01, 0.0, 1.0, ierflg);
  gMinuit->mnparm(5, "fmcrest_ML",  fmcrest_init,  0.01, 0.0, 1.0, ierflg);

  gMinuit->SetErrorDef(1.0);
  Double_t arglist[1] = {500};
  gMinuit->mnexcm("MIGRAD", arglist, 1, ierflg);

  // Retrieve fit results
  double feeg, fisr3pi, fomegapi, fnonReson, fksl, fmcrest;
  double feeg_err, fisr3pi_err, fomegapi_err, fnonReson_err, fksl_err, fmcrest_err;
  gMinuit->GetParameter(0, feeg,     feeg_err);
  gMinuit->GetParameter(1, fisr3pi,  fisr3pi_err);
  gMinuit->GetParameter(2, fomegapi, fomegapi_err);
  gMinuit->GetParameter(3, fnonReson, fnonReson_err);
  gMinuit->GetParameter(4, fksl,     fksl_err);
  gMinuit->GetParameter(5, fmcrest,  fmcrest_err);

  // Compute scaling factors
  double eeg_sfw     = getscale(nb_data_sum, feeg,     nb_eeg_sum);
  double eeg_sfw_err = GetScalError(nb_data_sum, nb_eeg_sum, feeg, feeg_err);
  double isr3pi_sfw     = getscale(nb_data_sum, fisr3pi,  nb_isr3pi_sum);
  double isr3pi_sfw_err = GetScalError(nb_data_sum, nb_isr3pi_sum, fisr3pi, fisr3pi_err);
  double omegapi_sfw     = getscale(nb_data_sum, fomegapi, nb_omegapi_sum);
  double omegapi_sfw_err = GetScalError(nb_data_sum, nb_omegapi_sum, fomegapi, fomegapi_err);
  double nonReson_sfw     = getscale(nb_data_sum, fnonReson,  nb_nonReson_sum);
  double nonReson_sfw_err = GetScalError(nb_data_sum, nb_nonReson_sum, fnonReson, fnonReson_err);
  double ksl_sfw     = getscale(nb_data_sum, fksl,     nb_ksl_sum);
  double ksl_sfw_err = GetScalError(nb_data_sum, nb_ksl_sum, fksl, fksl_err);
  double mcrest_sfw     = getscale(nb_data_sum, fmcrest,  nb_mcrest_sum);
  double mcrest_sfw_err = GetScalError(nb_data_sum, nb_mcrest_sum, fmcrest, fmcrest_err);

  // Append fit results to the text file
  myfile.open(myfile_nm.Data(), std::ios::app);
  myfile << "const double feeg = " << feeg << ";\n"
         << "const double fisr3pi = " << fisr3pi << ";\n"
         << "const double fomegapi = " << fomegapi << ";\n"
         << "const double fnonReson = " << fnonReson << ";\n"
         << "const double fksl = " << fksl << ";\n"
         << "const double fmcrest = " << fmcrest << ";\n\n"
         << "const double feeg_err = " << feeg_err << ";\n"
         << "const double fisr3pi_err = " << fisr3pi_err << ";\n"
         << "const double fomegapi_err = " << fomegapi_err << ";\n"
         << "const double fnonReson_err = " << fnonReson_err << ";\n"
         << "const double fksl_err = " << fksl_err << ";\n"
         << "const double fmcrest_err = " << fmcrest_err << ";\n\n"
         << "const double eeg_sfw = " << eeg_sfw << ";\n"
         << "const double isr3pi_sfw = " << isr3pi_sfw << ";\n"
         << "const double omegapi_sfw = " << omegapi_sfw << ";\n"
         << "const double nonReson_sfw = " << nonReson_sfw << ";\n"
         << "const double ksl_sfw = " << ksl_sfw << ";\n"
         << "const double mcrest_sfw = " << mcrest_sfw << ";\n";
  myfile.close();

  // Print results to screen
  std::cout << "\nSFW2D Fit Results: Fractions (initial)[%]\n"
            << "1: eeg     = " << feeg*100 << "(" << feeg_init*100 << ") +/- " << feeg_err*100 << "\n"
            << "2: isr3pi  = " << fisr3pi*100 << "(" << fisr3pi_init*100 << ") +/- " << fisr3pi_err*100 << "\n"
            << "3: omegapi = " << fomegapi*100 << "(" << fomegapi_init*100 << ") +/- " << fomegapi_err*100 << "\n"
            << "4: nonReson  = " << fnonReson*100 << "(" << fnonReson_init*100 << ") +/- " << fnonReson_err*100 << "\n"
            << "5: ksl     = " << fksl*100 << "(" << fksl_init*100 << ") +/- " << fksl_err*100 << "\n"
            << "6: mcrest  = " << fmcrest*100 << "(" << fmcrest_init*100 << ") +/- " << fmcrest_err*100 << "\n"
            << "       sum = " << (feeg+fisr3pi+fomegapi+fnonReson+fksl+fmcrest)*100
            << "(" << (feeg_init+fisr3pi_init+fomegapi_init+fnonReson_init+fksl_init+fmcrest_init)*100 << ")\n\n";

  std::cout << "SFW2D Fit Results: Scaling Factors\n"
            << "1: eeg     = " << eeg_sfw     << " +/- " << eeg_sfw_err << "\n"
            << "2: isr3pi  = " << isr3pi_sfw  << " +/- " << isr3pi_sfw_err << "\n"
            << "3: omegapi = " << omegapi_sfw << " +/- " << omegapi_sfw_err << "\n"
            << "4: nonReson  = " << nonReson_sfw  << " +/- " << nonReson_sfw_err << "\n"
            << "5: ksl     = " << ksl_sfw     << " +/- " << ksl_sfw_err << "\n"
            << "6: mcrest  = " << mcrest_sfw  << " +/- " << mcrest_sfw_err << "\n\n";

  Int_t ndf = residul_size_sfw2d - npar;
  double p_value = TMath::Prob(chi2_sfw2d_sum, ndf);
  std::cout << "SFW2D Fit Results: Chi2, ndf ...\n"
            << "residul_size = " << residul_size_sfw2d << ", npar = " << npar << ", ndf = " << ndf << "\n"
            << "chi2_sum = " << chi2_sfw2d_sum << ", chi2/ndf = " << chi2_sfw2d_sum/ndf << "\n"
            << "p-value = " << p_value << "\n";

  // Write final results tree
  TTree *TRESULT = new TTree("TRESULT", "recreate");
  TRESULT->Branch("Br_cut_value", &cut_value, "Br_cut_value/D");

  double SF[6]     = {eeg_sfw,     isr3pi_sfw,     omegapi_sfw,     nonReson_sfw,     ksl_sfw,     mcrest_sfw};
  double SF_ERR[6] = {eeg_sfw_err, isr3pi_sfw_err, omegapi_sfw_err, nonReson_sfw_err, ksl_sfw_err, mcrest_sfw_err};
  TRESULT->Branch("Br_SF",       SF,      "Br_SF[6]/D");
  TRESULT->Branch("Br_SF_ERR",   SF_ERR,  "Br_SF_ERR[6]/D");

  TRESULT->Branch("Br_eeg_sfw",     &eeg_sfw,     "Br_eeg_sfw/D");
  TRESULT->Branch("Br_eeg_sfw_err", &eeg_sfw_err, "Br_eeg_sfw_err/D");
  TRESULT->Branch("Br_isr3pi_sfw",     &isr3pi_sfw,     "Br_isr3pi_sfw/D");
  TRESULT->Branch("Br_isr3pi_sfw_err", &isr3pi_sfw_err, "Br_isr3pi_sfw_err/D");
  TRESULT->Branch("Br_omegapi_sfw",     &omegapi_sfw,     "Br_omegapi_sfw/D");
  TRESULT->Branch("Br_omegapi_sfw_err", &omegapi_sfw_err, "Br_omegapi_sfw_err/D");
  TRESULT->Branch("Br_nonReson_sfw",     &nonReson_sfw,     "Br_nonReson_sfw/D");
  TRESULT->Branch("Br_nonReson_sfw_err", &nonReson_sfw_err, "Br_nonReson_sfw_err/D");
  TRESULT->Branch("Br_ksl_sfw",     &ksl_sfw,     "Br_ksl_sfw/D");
  TRESULT->Branch("Br_ksl_sfw_err", &ksl_sfw_err, "Br_ksl_sfw_err/D");
  TRESULT->Branch("Br_mcrest_sfw",     &mcrest_sfw,     "Br_mcrest_sfw/D");
  TRESULT->Branch("Br_mcrest_sfw_err", &mcrest_sfw_err, "Br_mcrest_sfw_err/D");

  TRESULT->Branch("Br_feeg",     &feeg,     "Br_feeg/D");
  TRESULT->Branch("Br_feeg_err", &feeg_err, "Br_feeg_err/D");
  TRESULT->Branch("Br_fisr3pi",     &fisr3pi,     "Br_fisr3pi/D");
  TRESULT->Branch("Br_fisr3pi_err", &fisr3pi_err, "Br_fisr3pi_err/D");
  TRESULT->Branch("Br_fomegapi",     &fomegapi,     "Br_fomegapi/D");
  TRESULT->Branch("Br_fomegapi_err", &fomegapi_err, "Br_fomegapi_err/D");
  TRESULT->Branch("Br_fnonReson",     &fnonReson,     "Br_fnonReson/D");
  TRESULT->Branch("Br_fnonReson_err", &fnonReson_err, "Br_fnonReson_err/D");
  TRESULT->Branch("Br_fksl",     &fksl,     "Br_fksl/D");
  TRESULT->Branch("Br_fksl_err", &fksl_err, "Br_fksl_err/D");
  TRESULT->Branch("Br_fmcrest",     &fmcrest,     "Br_fmcrest/D");
  TRESULT->Branch("Br_fmcrest_err", &fmcrest_err, "Br_fmcrest_err/D");

  TRESULT->Branch("Br_nb_data_sum",    &nb_data_sum,    "Br_nb_data_sum/D");
  TRESULT->Branch("Br_nb_eeg_sum",     &nb_eeg_sum,     "Br_nb_eeg_sum/D");
  TRESULT->Branch("Br_nb_isr3pi_sum",  &nb_isr3pi_sum,  "Br_nb_isr3pi_sum/D");
  TRESULT->Branch("Br_nb_omegapi_sum", &nb_omegapi_sum, "Br_nb_omegapi_sum/D");
  TRESULT->Branch("Br_nb_nonReson_sum",  &nb_nonReson_sum,  "Br_nb_nonReson_sum/D");
  TRESULT->Branch("Br_nb_ksl_sum",     &nb_ksl_sum,     "Br_nb_ksl_sum/D");
  TRESULT->Branch("Br_nb_mcrest_sum",  &nb_mcrest_sum,  "Br_nb_mcrest_sum/D");

  TRESULT->Fill();

  // Write everything to output file
  f_output->cd();
  TRESULT->Write();
  TSFW2D->Write();
  f_output->Close();
  f_input->Close();

  if (gMinuit) delete gMinuit;

  return 0;
}
