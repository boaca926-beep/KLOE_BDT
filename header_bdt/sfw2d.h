#ifndef SFW2D_H
#define SFW2D_H

#include <TMath.h>
#include <TTree.h>

// Global sums (filled by sfw2d.C before fitting)
extern double nb_data_sum;
extern double nb_eeg_sum;
extern double nb_ksl_sum;
extern double nb_omegapi_sum;
extern double nb_etagam_sum;
extern double nb_isr3pi_sum;
extern double nb_mcrest_sum;

// Global fit results (updated by fcn_sfw2d)
extern double chi2_sfw2d_sum;
extern double residul_size_sfw2d;

// The tree that stores per‑bin yields
extern TTree* TSFW2D;

// ----------------------------------------------------------------------
// Helper functions
// ----------------------------------------------------------------------

inline double GetScalError(double N_d, double N, double f, double f_error) {
  // Scaling factor = N_d * f / N, with error propagation
  if (N == 0.0) return 0.0;                // avoid division by zero
  double scale = N_d * f / N;
  double error = scale * TMath::Sqrt(1.0/N_d + 1.0/N + TMath::Power(f_error/f, 2));
  return error;
}

inline double getloglh(double n_d, double mu) {
  // Poisson log‑likelihood (constant terms omitted for minimisation)
  if (mu <= 0.0) return -1e9;             // invalid argument
  return n_d * TMath::Log(mu) - mu;
}

inline double getscale(double Nd, double fra, double N) {
  if (N == 0.0) return 0.0;
  return Nd * fra / N;
}

// ----------------------------------------------------------------------
// Fit function for TMinuit
// ----------------------------------------------------------------------
void fcn_sfw2d(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

#endif
