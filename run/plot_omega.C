#include "../header/sm_para.h"
#include "../header/plot.h"
#include "../header/graph.h"
#include "../header/plot_omega.h"
#include "../header/method.h"

//int plot_omega(TString data_type = "DATA"){// fit_DATA.root
int plot_omega(){
  
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetFitFormat("6.4g");

  cout << input_file << endl;
  f_input = new TFile(input_file);
  getObj(f_input);

  f_input -> GetObject("hdata", hdata);
  f_input -> GetObject("heeg", heeg);
  f_input -> GetObject("homegapi", homegapi);
  f_input -> GetObject("hetagam", hetagam);
  f_input -> GetObject("hksl", hksl);
  f_input -> GetObject("hsig", hsig);
  f_input -> GetObject("hmcrest", hmcrest);
  f_input -> GetObject("hsmearmatr", hsmearmatr);
  f_input -> GetObject("hcorrmatrix", hcorrmatrix); 
  f_input -> GetObject("TCRX3PI", TCRX3PI);

  //TRESULT = (TTree*)TCRX3PI->Clone("TRESULT");
  
  
  //TRESULT = (TTree*)TCRX3PI->Clone("TRESULT");
  //if (!TRESULT) {
  //std::cerr << "Error: Could not find the tree 'TCRX3PI'." << std::endl;
  //f_input->Close();
  //}

  //TRESULT->Print();
  
  
  hbkgsum = (TH1D*) hmcrest -> Clone();
  hbkgsum -> Add(hksl, 1.);
  hbkgsum -> Add(hetagam, 1.);
  hbkgsum -> Add(homegapi, 1.);
  hbkgsum -> Add(heeg, 1.);
  
  hmcsum = (TH1D*) hbkgsum -> Clone();
  hmcsum -> Add(hsig, 1.);
  Delta_m3pi = getbinwidth(hdata);
  //cout << Delta_m3pi << endl;

  // Plot correlation and smearing matrices
  //plot_matr();

  // IM3pi distri.
  plot_IM3pi();

  // fit result
  plot_fit();

  TRESULT -> SetAutoSave(0);

  TRESULT -> Branch("Br_OMEGA_PARA", &OMEGA_PARA, "Br_OMEGA_PARA[3]/D");
  TRESULT -> Branch("Br_OMEGA_PARA_ERR", &OMEGA_PARA_ERR, "Br_OMEGA_PARA_ERR[3]/D");
  TRESULT -> Branch("Br_chi2_sum", &chi2_sum, "Br_chi2_sum/D");
  TRESULT -> Branch("Br_Lumi_int_fit", &Lumi_int_fit, "Br_Lumi_int_fit/D");
  TRESULT -> Branch("Br_ndf", &ndf, "Br_ndf/I");
  
  TRESULT -> Fill();
    
  cout << OMEGA_PARA[0] << "+/-" << OMEGA_PARA_ERR[0]<< ", " << OMEGA_PARA[1] << "+/-" << OMEGA_PARA_ERR[1] << ", " << OMEGA_PARA[2] << "+/-" << OMEGA_PARA_ERR[2] << ", " << chi2_sum << ", " << Lumi_int_fit << ", " << ndf << endl;
       
  // save
  //f_input -> GetObject("TCRX3PI", TRESULT);
  
  TFile *f_out = new TFile(outputPlot + "plot_omega_" + file_type + ".root", "recreate");
  cout << f_out -> GetName() << endl;
  
  gf_n3pi_obs -> Write();
  gf_n3pi_fit -> Write();
  gf_n3pi_diff -> Write();
  hdata -> Write();
  TRESULT -> Write();
  
  // pdf
  cv_n3pi_fit -> SaveAs(outputPlot + "crx3pi_" + file_type + ".pdf");
  cv_compr -> SaveAs(outputPlot + "n3pi_" + file_type + ".pdf");

  hmcsum -> Write();
  f_out -> Close();
    
  /*
  ofstream myfile;
  TString myfile_nm = "../header/omega_" + file_type + ".txt";
  myfile.open(myfile_nm);
  myfile << "const double Mass_omega_fit_" + file_type + "=" << Mass_omega_fit << ";\n"
	 << "const double Mass_omega_err_fit_" + file_type + "=" << Mass_omega_err_fit << ";\n"
    	 << "const double Gam_omega_fit_" + file_type + "=" << Gam_omega_fit << ";\n"
    	 << "const double Gam_omega_err_fit_" + file_type + "=" << Gam_omega_err_fit << ";\n"
    	 << "const double BB_fit_" + file_type + "=" << BB_fit << ";\n"
	 << "const double BB_err_fit_" + file_type + "=" << BB_err_fit << ";\n"
    	 << "const double Lumi_int_fit_" + file_type + "=" << Lumi_int_fit << ";\n"
	 << "const double chi2_sum_" + file_type + "=" << chi2_sum << ";\n"
	 << "const int ndf_" + file_type + "=" << ndf << ";\n";
    	 
  myfile.close();
  */
  
  return 0;

}

