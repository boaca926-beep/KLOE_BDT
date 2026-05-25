#include "../header_bdt/sm_para.h"
#include "../header_bdt/path.h"
#include "../header_bdt/method.h"
#include "../header_bdt/omega_fit.h"
#include "../header_bdt/cut_para.h"

int omega_fit(){

  gErrorIgnoreLevel = kError;
  
  cout << "Extract omega parameters ..." << endl;

  // efficiency
  get_efficy();
  //hsig_true->Draw();
  hefficy->Draw();

  // --- Check efficiency in omega region (760-820 MeV) ---
  double omega_low = 400.; //760.0;
  double omega_high = 900.; //820.0;
  
  int bin_low = hefficy->FindBin(omega_low);
  int bin_high = hefficy->FindBin(omega_high);
  
  double eff_sum = 0.0;
  double err2_sum = 0.0;
  int n_bins = 0;
  
  for (int bin = bin_low; bin <= bin_high; ++bin) {
    double eff = hefficy->GetBinContent(bin);
    double err = hefficy->GetBinError(bin);
    eff_sum += eff;
    err2_sum += err * err;
    n_bins++;
  }
  
  double eff_avg = (n_bins > 0) ? eff_sum / n_bins : 0.0;
  double eff_avg_err = (n_bins > 0) ? TMath::Sqrt(err2_sum) / n_bins : 0.0;
  
  cout << "Efficiency in ω region (" << omega_low << "-" << omega_high << " MeV/c²):" << endl;
  cout << "  Average = " << eff_avg << " ± " << eff_avg_err << endl;
  
  // Optional: print bin-by-bin values
  for (int bin = bin_low; bin <= bin_high; ++bin) {
    double mass = hefficy->GetXaxis()->GetBinCenter(bin);
    double eff = hefficy->GetBinContent(bin);
    double err = hefficy->GetBinError(bin);
    cout << "  M = " << mass << " MeV/c² → efficiency = " << eff << " ± " << err << endl;
  }

  return 0;
  
}
