// Full RooFit implementation for 2D template fit 
// using extended likelihood with yield parameters.

#include "../header_bdt/path.h"
#include <fstream>
#include <TLegend.h>

using namespace RooFit;

// ----------------------------------------------------------------------
// Helper: add a tiny constant to all bins of a TH2 to avoid zero PDF values
// ----------------------------------------------------------------------
void addTinyConstant(TH2* h, double eps = 1e-12) {
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        for (int j = 1; j <= h->GetNbinsY(); ++j) {
            double old = h->GetBinContent(i, j);
            h->SetBinContent(i, j, old + eps);
            double err = h->GetBinError(i, j);
            h->SetBinError(i, j, std::sqrt(err*err + eps));
        }
    }
}

// ----------------------------------------------------------------------
int sfw2d_rootfit() {
  // ------------------------------------------------------------------
  // 1. Open input file and retrieve the 2D histogram list
  // ------------------------------------------------------------------
  
  TString inputFileName = outputHist + "hist.root";
  cout << inputFileName << endl;
  TFile *fin = TFile::Open(inputFileName);
  if (!fin || fin->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << inputFileName << std::endl;
    return 1;
  }

  TList *HSFW2D = (TList*) fin->Get("HSFW2D");
  if (!HSFW2D) {
    std::cerr << "ERROR: HSFW2D not found in file" << std::endl;
    return 1;
  }

  TH2D *h_data      = (TH2D*) HSFW2D->FindObject("h2d_sfw_TDATA");
  TH2D *h_eeg       = (TH2D*) HSFW2D->FindObject("h2d_sfw_TEEG");
  TH2D *h_isr3pi_peak = (TH2D*) HSFW2D->FindObject("h2d_sfw_TISR3PI_SIG_peak");
  TH2D *h_isr3pi_non_reson = (TH2D*) HSFW2D->FindObject("h2d_sfw_TISR3PI_SIG_non_reson");
  TH2D *h_omegapi   = (TH2D*) HSFW2D->FindObject("h2d_sfw_TOMEGAPI");
  TH2D *h_etagam    = (TH2D*) HSFW2D->FindObject("h2d_sfw_TETAGAM");
  TH2D *h_ksl       = (TH2D*) HSFW2D->FindObject("h2d_sfw_TKSL");
  TH2D *h_kpm       = (TH2D*) HSFW2D->FindObject("h2d_sfw_TKPM");
  TH2D *h_rhopi     = (TH2D*) HSFW2D->FindObject("h2d_sfw_TRHOPI");
  TH2D *h_bkgrest   = (TH2D*) HSFW2D->FindObject("h2d_sfw_TBKGREST");

  if (!h_data || !h_eeg || !h_isr3pi_peak || !h_isr3pi_non_reson || !h_omegapi || !h_etagam || !h_ksl || !h_kpm || !h_rhopi || !h_bkgrest) {
    std::cerr << "ERROR: Missing one or more required histograms" << std::endl;
    return 1;
  }

  // Summary
  double nb_data = h_data->Integral();
  double nb_eeg = h_eeg->Integral();
  double nb_isr3pi_peak = h_isr3pi_peak->Integral();
  double nb_isr3pi_non_reson = h_isr3pi_non_reson->Integral();
  double nb_omegapi = h_omegapi->Integral();
  double nb_etagam = h_etagam->Integral();
  double nb_ksl = h_ksl->Integral();
  double nb_kpm = h_kpm->Integral();
  double nb_rhopi = h_rhopi->Integral();
  double nb_bkgrest = h_bkgrest->Integral();

  // Build combined "mcrest" histogram (same as in sfw2d.C)
  TH2D *h_mcrest = (TH2D*) h_bkgrest->Clone("h_mcrest");
  h_mcrest->Add(h_kpm, 1.0);
  h_mcrest->Add(h_rhopi, 1.0);
  h_mcrest->Add(h_etagam, 1.0);

  double nb_mcrest = h_mcrest->Integral();
  
  cout << "nb_data = " << nb_data << "\n"
       << "MC: " << "\n"
       << "1. eeg = " << nb_eeg << "\n"
       << "2. isr3pi_peak = " << nb_isr3pi_peak << "\n"
       << "3. isr3pi_non_reson = " << nb_isr3pi_non_reson << "\n"
       << "4. omegapi = " << nb_omegapi << "\n"
       << "5. ksl = " << nb_ksl << "\n"
       << "6. mcrest = " << nb_mcrest << ", checked = " << nb_kpm + nb_rhopi + nb_bkgrest + nb_etagam << "\n"
       << "\tnb_kpm = " << nb_kpm << "\n"
       << "\tnb_rhopi = " << nb_rhopi << "\n"
       << "\tnb_bkgrest = " << nb_bkgrest << "\n"
       << "\tnb_etagam = " << nb_etagam << "\n";
       

  // ------------------------------------------------------------------
  // 2. Add a tiny constant to all MC histograms (to avoid zero PDF values)
  //    Do NOT modify the data histogram!
  // ------------------------------------------------------------------
  const double eps = 1e-12;
  addTinyConstant(h_eeg, eps);
  addTinyConstant(h_isr3pi_peak, eps);
  addTinyConstant(h_isr3pi_non_reson, eps);
  addTinyConstant(h_omegapi, eps);
  addTinyConstant(h_ksl, eps);
  addTinyConstant(h_mcrest, eps);
  // h_data remains unchanged

  // ------------------------------------------------------------------
  // 3. Define the observables (x = your 1D fit variable, ppIM;
  //    y = energy of unpaired photon). Adjust axis ranges to match your histograms.
  // ------------------------------------------------------------------
  RooRealVar x("x", "Fit variable",
	       h_data->GetXaxis()->GetXmin(),
	       h_data->GetXaxis()->GetXmax());
  RooRealVar y("y", "E_{#gamma_{3}} [MeV]",
	       h_data->GetYaxis()->GetXmin(),
	       h_data->GetYaxis()->GetXmax());

  // ------------------------------------------------------------------
  // 4. Convert TH2 histograms to RooDataHist (binned data sets)
  // ------------------------------------------------------------------
  RooDataHist dataHist("dataHist", "Data", RooArgSet(x, y), h_data);
  RooDataHist eegHist("eegHist", "EEG MC", RooArgSet(x, y), h_eeg);
  RooDataHist isr3piHist("isr3piHist", "ISR3pi MC", RooArgSet(x, y), h_isr3pi_peak);
  RooDataHist nonResonHist("nonResonHist", "Non-reson MC", RooArgSet(x, y), h_isr3pi_non_reson);
  RooDataHist omegapiHist("omegapiHist", "OmegaPi MC", RooArgSet(x, y), h_omegapi);
  RooDataHist kslHist("kslHist", "KSL MC", RooArgSet(x, y), h_ksl);
  RooDataHist mcrestHist("mcrestHist", "MC Rest", RooArgSet(x, y), h_mcrest);

  // ------------------------------------------------------------------
  // 5. Create PDFs from templates (RooHistPdf)
  // ------------------------------------------------------------------
  RooHistPdf eegPdf("eegPdf", "EEG PDF", RooArgSet(x, y), eegHist);
  RooHistPdf isr3piPdf("isr3piPdf", "ISR3pi PDF", RooArgSet(x, y), isr3piHist);
  RooHistPdf omegapiPdf("omegapiPdf", "OmegaPi PDF", RooArgSet(x, y), omegapiHist);
  RooHistPdf nonResonPdf("nonResonPdf", "Non-reson PDF", RooArgSet(x, y), nonResonHist);
  RooHistPdf kslPdf("kslPdf", "KSL PDF", RooArgSet(x, y), kslHist);
  RooHistPdf mcrestPdf("mcrestPdf", "MC Rest PDF", RooArgSet(x, y), mcrestHist);
  
  // ------------------------------------------------------------------
  // 6. Get total MC yields for initial values of the fit parameters
  // ------------------------------------------------------------------
  double nb_eeg_sum     = h_eeg->Integral();
  double nb_isr3pi_sum  = h_isr3pi_peak->Integral();
  double nb_nonReson_sum  = h_isr3pi_non_reson->Integral();
  double nb_omegapi_sum = h_omegapi->Integral();
  double nb_ksl_sum     = h_ksl->Integral();
  double nb_mcrest_sum  = h_mcrest->Integral();

  // ------------------------------------------------------------------
  // 7. Create yield parameters (extended likelihood will use these as coefficients)
  // ------------------------------------------------------------------

  // Calculate scale factor to bring MC total to data size
  double mc_total = nb_isr3pi_sum + nb_nonReson_sum + nb_eeg_sum + nb_omegapi_sum + nb_ksl_sum + nb_mcrest_sum;
  double scale_to_data = nb_data / mc_total;  // ~0.1
  
  cout << "\n=== PARAMETER INITIALIZATION ===" << endl;
  cout << "Data events: " << nb_data << endl;
  cout << "MC total events: " << mc_total << endl;
  cout << "Scale factor to data: " << scale_to_data << endl;
  cout << "================================\n" << endl;

  // Signal components (ISR3pi) - expected to dominate the fit
  // Start with 60% of data for peak, 10% for non-resonant based on physics expectations
  RooRealVar N_isr3pi("N_isr3pi", "ISR3pi peak yield", 
                      nb_data * 0.60,      // initial: ~10224
                      0.0,                 // lower bound
                      nb_data * 1.2);      // upper bound: ~20448
  
  RooRealVar N_nonReson("N_nonReson", "ISR3pi non-resonant yield", 
                        nb_data * 0.10,    // initial: ~1704
                        0.0,
                        nb_data * 0.8);    // upper bound: ~13632
  
  // Background components - scale MC expectations to data size
  RooRealVar N_eeg("N_eeg", "EEG yield", 
                   nb_eeg_sum * scale_to_data,   // initial: ~32.1
                   0.0, 
                   nb_eeg_sum * scale_to_data * 20.0);
  
  RooRealVar N_omegapi("N_omegapi", "OmegaPi yield", 
                       nb_omegapi_sum * scale_to_data,   // initial: ~191
                       0.0, 
                       nb_omegapi_sum * scale_to_data * 20.0);
  
  RooRealVar N_ksl("N_ksl", "KSL yield", 
                   nb_ksl_sum * scale_to_data,   // initial: ~320
                   0.0, 
                   nb_ksl_sum * scale_to_data * 20.0);
  
  RooRealVar N_mcrest("N_mcrest", "MC Rest yield", 
                      nb_mcrest_sum * scale_to_data,   // initial: ~8.6
                      0.0, 
                      nb_mcrest_sum * scale_to_data * 30.0);

  //N_mcrest.setConstant(kTRUE);
  cout << "NOTE: N_mcrest fixed to " << N_mcrest.getVal() << " (low statistics component)" << endl;
  
  // Print initial values for verification
  cout << "Initial parameter values (scaled to data):" << endl;
  cout << "  N_isr3pi     = " << N_isr3pi.getVal() << " (MC input: " << nb_isr3pi_sum << ")" << endl;
  cout << "  N_nonReson   = " << N_nonReson.getVal() << " (MC input: " << nb_nonReson_sum << ")" << endl;
  cout << "  N_eeg        = " << N_eeg.getVal() << " (MC input: " << nb_eeg_sum << ")" << endl;
  cout << "  N_omegapi    = " << N_omegapi.getVal() << " (MC input: " << nb_omegapi_sum << ")" << endl;
  cout << "  N_ksl        = " << N_ksl.getVal() << " (MC input: " << nb_ksl_sum << ")" << endl;
  cout << "  N_mcrest     = " << N_mcrest.getVal() << " (MC input: " << nb_mcrest_sum << ")" << endl;
  
  double initial_total = N_isr3pi.getVal() + N_nonReson.getVal() + N_eeg.getVal() + 
    N_omegapi.getVal() + N_ksl.getVal() + N_mcrest.getVal();
  cout << "  Initial total = " << initial_total << " (target: " << nb_data << ")" << endl;
  cout << "================================\n" << endl;
  

  // ------------------------------------------------------------------
  // 8. Build the total PDF as a sum of PDFs weighted by yields
  //    This automatically makes the fit extended (total yield = sum of yields)
  // ------------------------------------------------------------------
  RooAddPdf totalPdf("totalPdf", "Total model",
		     RooArgList(eegPdf, isr3piPdf, omegapiPdf, nonResonPdf, kslPdf, mcrestPdf),
		     RooArgList(N_eeg, N_isr3pi, N_omegapi, N_nonReson, N_ksl, N_mcrest));

  // ------------------------------------------------------------------
  // 9. Perform extended binned maximum likelihood fit
  // ------------------------------------------------------------------
  RooFitResult *fitRes = totalPdf.fitTo(dataHist,
					SumW2Error(kTRUE),   // include MC statistical errors
					Save(kTRUE),
					PrintLevel(1));
  fitRes->Print("v");

  // Print correlations between floating parameters only
  TMatrixDSym corr = fitRes->correlationMatrix();
  cout << "\n=== CORRELATION MATRIX (floating parameters) ===" << endl;
  cout << "Correlation(N_isr3pi, N_nonReson) = " << corr(0, 1) << endl;
  cout << "Correlation(N_isr3pi, N_eeg) = " << corr(0, 2) << endl;
  cout << "Correlation(N_isr3pi, N_omegapi) = " << corr(0, 3) << endl;
  cout << "Correlation(N_isr3pi, N_ksl) = " << corr(0, 4) << endl;

  // ------------------------------------------------------------------
  // 10. Print fitted yields and derived scaling factors
  // ------------------------------------------------------------------
  std::cout << "\n=== Fitted yields (absolute normalisation) ===" << std::endl;
  std::cout << "EEG        : " << N_eeg.getVal() << " +/- " << N_eeg.getError() << std::endl;
  std::cout << "ISR3pi     : " << N_isr3pi.getVal() << " +/- " << N_isr3pi.getError() << std::endl;
  std::cout << "Non-reson  : " << N_nonReson.getVal() << " +/- " << N_nonReson.getError() << std::endl;
  std::cout << "OmegaPi    : " << N_omegapi.getVal() << " +/- " << N_omegapi.getError() << std::endl;
  std::cout << "KSL        : " << N_ksl.getVal() << " +/- " << N_ksl.getError() << std::endl;
  std::cout << "MC Rest    : " << N_mcrest.getVal() << " +/- " << N_mcrest.getError() << std::endl;

  // Calculate total fitted yield
  double total_fitted = N_eeg.getVal() + N_isr3pi.getVal() + N_nonReson.getVal() + 
                        N_omegapi.getVal() + N_ksl.getVal() + N_mcrest.getVal();
  std::cout << "\nTotal fitted: " << total_fitted << " (data: " << nb_data << ")" << std::endl;

  // Compute scaling factors (Data / scaled MC expectation)
  // Scaled MC = MC integral * scale_to_data
  auto computeSF = [scale_to_data](double N_fit, double N_mc, double err_fit, double err_mc=0.0) {
    double N_mc_scaled = N_mc * scale_to_data;
    double sf = N_fit / N_mc_scaled;
    double rel_err = std::hypot(err_fit/N_fit, err_mc/N_mc);
    return std::make_pair(sf, sf * rel_err);
  };
  
  auto [sf_eeg, err_eeg]         = computeSF(N_eeg.getVal(), nb_eeg_sum, N_eeg.getError());
  auto [sf_isr3pi, err_isr3pi]   = computeSF(N_isr3pi.getVal(), nb_isr3pi_sum, N_isr3pi.getError());
  auto [sf_nonReson, err_nonReson] = computeSF(N_nonReson.getVal(), nb_nonReson_sum, N_nonReson.getError());
  auto [sf_omegapi, err_omegapi] = computeSF(N_omegapi.getVal(), nb_omegapi_sum, N_omegapi.getError());
  auto [sf_ksl, err_ksl]         = computeSF(N_ksl.getVal(), nb_ksl_sum, N_ksl.getError());
  auto [sf_mcrest, err_mcrest]   = computeSF(N_mcrest.getVal(), nb_mcrest_sum, N_mcrest.getError());

  std::cout << "\n=== Scaling factors (Data / (MC * scale_to_data)) ===" << std::endl;
  std::cout << "EEG        : " << sf_eeg << " +/- " << err_eeg << std::endl;
  std::cout << "ISR3pi     : " << sf_isr3pi << " +/- " << err_isr3pi << std::endl;
  std::cout << "Non-reson  : " << sf_nonReson << " +/- " << err_nonReson << std::endl;
  std::cout << "OmegaPi    : " << sf_omegapi << " +/- " << err_omegapi << std::endl;
  std::cout << "KSL        : " << sf_ksl << " +/- " << err_ksl << std::endl;
  std::cout << "MC Rest    : " << sf_mcrest << " +/- " << err_mcrest << std::endl;

  // Print scaled MC expectations for reference
  std::cout << "\n=== Scaled MC expectations (MC * scale_to_data) ===" << std::endl;
  std::cout << "EEG        : " << nb_eeg_sum * scale_to_data << std::endl;
  std::cout << "ISR3pi     : " << nb_isr3pi_sum * scale_to_data << std::endl;
  std::cout << "Non-reson  : " << nb_nonReson_sum * scale_to_data << std::endl;
  std::cout << "OmegaPi    : " << nb_omegapi_sum * scale_to_data << std::endl;
  std::cout << "KSL        : " << nb_ksl_sum * scale_to_data << std::endl;
  std::cout << "MC Rest    : " << nb_mcrest_sum * scale_to_data << std::endl;

  // Append fit results to the text file
  std::ofstream myfile;
  TString myfile_nm = "../header_bdt/sfw2d_bdt.txt";
  myfile.open(myfile_nm.Data(), std::ios::out | std::ios::trunc);
  myfile << "const double eeg_sfw = " << sf_eeg << ";\n"
	 << "const double isr3pi_sfw = " << sf_isr3pi << ";\n"
    	 << "const double nonReson_sfw = " << sf_nonReson << ";\n"
	 << "const double omegapi_sfw = " << sf_omegapi << ";\n"
	 << "const double ksl_sfw = " << sf_ksl << ";\n"
	 << "const double mcrest_sfw = " << sf_mcrest << ";\n";
  myfile.close();

  cout << "\nScaling factors saved to: " << myfile_nm << endl;

  // ------------------------------------------------------------------
  // 11. Project the fit onto the y axis (integrate over x)
  // ------------------------------------------------------------------
  TCanvas *c1 = new TCanvas("c1", "Eisr projection", 800, 600);
  RooPlot *EisrFrame = y.frame(Title("Projection onto E_{#gamma_{3}}"));
  dataHist.plotOn(EisrFrame, MarkerSize(0.8), MarkerColor(kBlack));
  totalPdf.plotOn(EisrFrame, LineColor(kBlue), LineWidth(2));
  
  // Add individual components
  totalPdf.plotOn(EisrFrame, Components(isr3piPdf), LineStyle(kDashed), LineColor(kRed));
  totalPdf.plotOn(EisrFrame, Components(nonResonPdf), LineStyle(kDashed), LineColor(kOrange));
  totalPdf.plotOn(EisrFrame, Components(RooArgSet(eegPdf, omegapiPdf, kslPdf, mcrestPdf)), 
                  LineStyle(kDashed), LineColor(kGreen));
  
  // Add legend
  TLegend* leg = new TLegend(0.65, 0.65, 0.89, 0.89);
  leg->AddEntry(EisrFrame->getObject(0), "Data", "LP");
  leg->AddEntry((TObject*)0, "Total Model", "L");
  leg->AddEntry((TObject*)0, "ISR3pi Peak", "L");
  leg->AddEntry((TObject*)0, "ISR3pi Non-reson", "L");
  leg->AddEntry((TObject*)0, "Backgrounds", "L");
  leg->Draw();
  
  EisrFrame->Draw();
  c1->SaveAs("sfw2d_projection.pdf");
  
  //delete leg;
  //delete c1;
  
  // ------------------------------------------------------------------
  // 12. Save results to ROOT file
  // ------------------------------------------------------------------
  TFile *fout = new TFile("sfw2d_fit_result.root", "RECREATE");
  fitRes->Write();
  totalPdf.Write();
  dataHist.Write();
  EisrFrame->Write();
  fout->Close();
  
  cout << "\n=== FIT COMPLETED SUCCESSFULLY ===" << endl;
  cout << "Output files:" << endl;
  cout << "  - sfw2d_projection.pdf" << endl;
  cout << "  - sfw2d_fit_result.root" << endl;
  cout << "  - ../header_bdt/sfw2d_bdt.txt (appended)" << endl;
  
  fin->Close();
  return 0;
}
