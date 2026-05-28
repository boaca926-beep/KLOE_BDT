#include "../header/sm_para.h"
#include "../header/plot.h"
#include "../header/efficy.h"
#include "../header/plot_efficy.h"
#include "../header/method.h"
#include "../header/graph.h"

double ratioErr(double a, double sigma_a, double b, double sigma_b) {// ratio = a / b

  if (b == 0) {

    cout << "Division by zero!" << endl;
    return 0;

  }

  double c = a / b;
  double rela_err = TMath::Sqrt(TMath::Power(sigma_a / a, 2) + TMath::Power(sigma_b / b, 2));
  double sigma_c = c * rela_err;

  //cout << TMath::Power(sigma_a / a, 2) + TMath::Power(sigma_b / b, 2) << ", a = " << a << ", b = " << b << endl;
   
  return sigma_c;
  
}

double get_ratio(double a, double b) {// ratio = a / b

  if (b == 0) {

    cout << "Division by zero!" << endl;
    return 0;

  }

  double c = a / b;

  return c;
  
}

int plot_efficy() {

  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFitFormat("6.4g");

  cout << "Plot efficiency ..." << endl;

  cout << "input path: " << input_folder << endl;

  TFile *f_input = new TFile(input_folder + "/efficy.root");
  //getObj(f_input);

  TGraphErrors* gf_efficy_sig = (TGraphErrors*)f_input -> Get("gf_efficy_TISR3PI_SIG");
  gf_efficy_sig -> SetLineColor(kBlue);
  gf_efficy_sig -> SetMarkerSize(.8);
  gf_efficy_sig -> SetMarkerStyle(20);
  
  TGraphErrors* gf_efficy_ufo = (TGraphErrors*)f_input -> Get("gf_efficy_TUFO");
  gf_efficy_ufo -> SetLineColor(kBlack);
  gf_efficy_ufo -> SetMarkerSize(.8);
  gf_efficy_ufo -> SetMarkerStyle(22);
  
  // calcualte ratio
  int nPoints = gf_efficy_sig -> GetN();
  double *x_efficy_sig = gf_efficy_sig -> GetX();
  double *x_efficy_sig_err = gf_efficy_sig -> GetEX();

  double *y_efficy_sig = gf_efficy_sig -> GetY();
  double *y_efficy_sig_err = gf_efficy_sig -> GetEY();
  
  double *y_efficy_ufo = gf_efficy_ufo -> GetY();
  double *y_efficy_ufo_err = gf_efficy_ufo -> GetEY();
  
  double RATIO[nPoints], RATIO_ERR[nPoints];

  cout << "nPoints = " << nPoints << endl;
  
  for (int i = 0; i < nPoints; i ++) {

    if (y_efficy_sig[i] == 0. || y_efficy_ufo[i] == 0.) {
        RATIO[i] = 0.;
	RATIO_ERR[i] = 0.;
    }
    else {
      RATIO[i] = get_ratio(y_efficy_ufo[i], y_efficy_sig[i]); //y_efficy_ufo[i] / y_efficy_sig[i];
      RATIO_ERR[i] = ratioErr(y_efficy_ufo[i], y_efficy_ufo_err[i], y_efficy_sig[i], y_efficy_sig_err[i]);
  
    }
    
    //cout << "point " << i << ", mass =" << x_efficy_sig[i] << ", efficy_sig = " << y_efficy_sig[i] << "+/-" << y_efficy_sig_err[i] << ", efficy_ufo = " << y_efficy_ufo[i] << "+/-" << y_efficy_ufo_err[i] << ", efficy ratio (ufo/sig) = " << RATIO[i] << "+/-" << RATIO_ERR[i] << endl;
    
  }

  TGraphErrors *gf_ratio = get_graph_syst(x_efficy_sig, RATIO, x_efficy_sig_err, RATIO_ERR, nPoints);
  gf_ratio -> SetName("gf_ratio");
  //gf_ratio -> Draw("AP");
  
  /*
  double a = 10., sigma_a = 1.;
  double b = 2., sigma_b = .5;
  double c = a / b;
  
  double sigma_c = ratioErr(a, sigma_a, b, sigma_b);
  cout << "c = " << c << " +/- " << sigma_c << endl;
  */
  
  // plot
  TCanvas *cv_efficy = plotting_efficy("cv_efficy", "Efficiency Comparsion", gf_efficy_sig, gf_efficy_ufo, gf_ratio);

  //gf_efficy_sig -> Draw("AP");
  //gf_efficy_ufo -> Draw("P");

  // save
  cv_efficy -> SaveAs(input_folder + "/cv_efficy.pdf");
  
  cout << input_folder << endl;
  
  return 0;
  
}
