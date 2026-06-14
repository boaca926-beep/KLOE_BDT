// Full RooFit implementation for 2D template fit 
// using extended likelihood with yield parameters.
// Features: Gaussian constraints, MINOS errors, correlation analysis

#include "../header_bdt/path.h"
#include <fstream>
#include <TLegend.h>
#include <ctime>

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
int sfw2d_rootfit_improved() {
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

  // Build combined "mcrest" histogram
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
       << "6. mcrest = " << nb_mcrest << "\n";
       

  // ------------------------------------------------------------------
  // 2. Add tiny constant to MC histograms
  // ------------------------------------------------------------------
  const double eps = 1e-12;
  addTinyConstant(h_eeg, eps);
  addTinyConstant(h_isr3pi_peak, eps);
  addTinyConstant(h_isr3pi_non_reson, eps);
  addTinyConstant(h_omegapi, eps);
  addTinyConstant(h_ksl, eps);
  addTinyConstant(h_mcrest, eps);

  // ------------------------------------------------------------------
  // 3. Define observables
  // ------------------------------------------------------------------
  RooRealVar x("x", "Fit variable",
	       h_data->GetXaxis()->GetXmin(),
	       h_data->GetXaxis()->GetXmax());
  RooRealVar y("y", "E_{#gamma_{3}} [MeV]",
	       h_data->GetYaxis()->GetXmin(),
	       h_data->GetYaxis()->GetXmax());

  // ------------------------------------------------------------------
  // 4. Convert to RooDataHist
  // ------------------------------------------------------------------
  RooDataHist dataHist("dataHist", "Data", RooArgSet(x, y), h_data);
  RooDataHist eegHist("eegHist", "EEG MC", RooArgSet(x, y), h_eeg);
  RooDataHist isr3piHist("isr3piHist", "ISR3pi MC", RooArgSet(x, y), h_isr3pi_peak);
  RooDataHist nonResonHist("nonResonHist", "Non-reson MC", RooArgSet(x, y), h_isr3pi_non_reson);
  RooDataHist omegapiHist("omegapiHist", "OmegaPi MC", RooArgSet(x, y), h_omegapi);
  RooDataHist kslHist("kslHist", "KSL MC", RooArgSet(x, y), h_ksl);
  RooDataHist mcrestHist("mcrestHist", "MC Rest", RooArgSet(x, y), h_mcrest);

  // ------------------------------------------------------------------
  // 5. Create PDFs
  // ------------------------------------------------------------------
  RooHistPdf eegPdf("eegPdf", "EEG PDF", RooArgSet(x, y), eegHist);
  RooHistPdf isr3piPdf("isr3piPdf", "ISR3pi PDF", RooArgSet(x, y), isr3piHist);
  RooHistPdf omegapiPdf("omegapiPdf", "OmegaPi PDF", RooArgSet(x, y), omegapiHist);
  RooHistPdf nonResonPdf("nonResonPdf", "Non-reson PDF", RooArgSet(x, y), nonResonHist);
  RooHistPdf kslPdf("kslPdf", "KSL PDF", RooArgSet(x, y), kslHist);
  RooHistPdf mcrestPdf("mcrestPdf", "MC Rest PDF", RooArgSet(x, y), mcrestHist);
  
  // ------------------------------------------------------------------
  // 6. Get MC yields
  // ------------------------------------------------------------------
  double nb_eeg_sum     = h_eeg->Integral();
  double nb_isr3pi_sum  = h_isr3pi_peak->Integral();
  double nb_nonReson_sum  = h_isr3pi_non_reson->Integral();
  double nb_omegapi_sum = h_omegapi->Integral();
  double nb_ksl_sum     = h_ksl->Integral();
  double nb_mcrest_sum  = h_mcrest->Integral();

  // ------------------------------------------------------------------
  // 7. Calculate scale factor
  // ------------------------------------------------------------------
  double mc_total = nb_isr3pi_sum + nb_nonReson_sum + nb_eeg_sum + nb_omegapi_sum + nb_ksl_sum + nb_mcrest_sum;
  double scale_to_data = nb_data / mc_total;
  
  cout << "\n=== PARAMETER INITIALIZATION ===" << endl;
  cout << "Data events: " << nb_data << endl;
  cout << "MC total events: " << mc_total << endl;
  cout << "Scale factor to data: " << scale_to_data << endl;
  cout << "================================\n" << endl;

  // ------------------------------------------------------------------
  // 8. Create yield parameters with Gaussian constraints
  // ------------------------------------------------------------------
  
  // Signal component - wide prior (let it float freely)
  RooRealVar N_isr3pi("N_isr3pi", "ISR3pi peak yield", 
                      nb_data * 0.60, 0.0, nb_data * 1.2);
  
  // Non-resonant - tighter constraint (factor 2 uncertainty)
  RooRealVar N_nonReson("N_nonReson", "ISR3pi non-resonant yield", 
                        nb_data * 0.10, 0.0, nb_data * 0.5);
  
  // Background components with Gaussian constraints (factor 2-3 uncertainty)
  double eeg_scaled = nb_eeg_sum * scale_to_data;
  RooRealVar N_eeg("N_eeg", "EEG yield", eeg_scaled, 0.0, eeg_scaled * 10.0);
  RooGaussian eeg_constraint("eeg_constraint", "EEG constraint", N_eeg, 
                             RooConst(eeg_scaled), RooConst(eeg_scaled * 2.0));
  
  double omegapi_scaled = nb_omegapi_sum * scale_to_data;
  RooRealVar N_omegapi("N_omegapi", "OmegaPi yield", omegapi_scaled, 0.0, omegapi_scaled * 10.0);
  RooGaussian omegapi_constraint("omegapi_constraint", "OmegaPi constraint", N_omegapi,
                                  RooConst(omegapi_scaled), RooConst(omegapi_scaled * 2.0));
  
  double ksl_scaled = nb_ksl_sum * scale_to_data;
  RooRealVar N_ksl("N_ksl", "KSL yield", ksl_scaled, 0.0, ksl_scaled * 10.0);
  RooGaussian ksl_constraint("ksl_constraint", "KSL constraint", N_ksl,
                              RooConst(ksl_scaled), RooConst(ksl_scaled * 2.0));
  
  // MC Rest - low statistics, fix to scaled value
  RooRealVar N_mcrest("N_mcrest", "MC Rest yield", 
                      nb_mcrest_sum * scale_to_data, 0.0, nb_mcrest_sum * scale_to_data * 5.0);
  //N_mcrest.setConstant(kTRUE);
  
  // Print initial values
  cout << "Initial parameter values:" << endl;
  cout << "  N_isr3pi     = " << N_isr3pi.getVal() << " (no constraint)" << endl;
  cout << "  N_nonReson   = " << N_nonReson.getVal() << " (tight constraint)" << endl;
  cout << "  N_eeg        = " << N_eeg.getVal() << " (Gaussian constraint: " << eeg_scaled << " +/- " << eeg_scaled*2 << ")" << endl;
  cout << "  N_omegapi    = " << N_omegapi.getVal() << " (Gaussian constraint: " << omegapi_scaled << " +/- " << omegapi_scaled*2 << ")" << endl;
  cout << "  N_ksl        = " << N_ksl.getVal() << " (Gaussian constraint: " << ksl_scaled << " +/- " << ksl_scaled*2 << ")" << endl;
  cout << "  N_mcrest     = " << N_mcrest.getVal() << " (fixed)" << endl;
  
  double initial_total = N_isr3pi.getVal() + N_nonReson.getVal() + N_eeg.getVal() + 
                         N_omegapi.getVal() + N_ksl.getVal() + N_mcrest.getVal();
  cout << "  Initial total = " << initial_total << " (target: " << nb_data << ")" << endl;
  cout << "================================\n" << endl;

  // ------------------------------------------------------------------
  // 9. Build total PDF
  // ------------------------------------------------------------------
  RooAddPdf totalPdf("totalPdf", "Total model",
		     RooArgList(eegPdf, isr3piPdf, omegapiPdf, nonResonPdf, kslPdf, mcrestPdf),
		     RooArgList(N_eeg, N_isr3pi, N_omegapi, N_nonReson, N_ksl, N_mcrest));

  // ------------------------------------------------------------------
  // 10. Perform fit with Gaussian constraints
  // ------------------------------------------------------------------
  RooArgSet constraints(eeg_constraint, omegapi_constraint, ksl_constraint);
  
  cout << "\n=== Performing fit with Gaussian constraints ===" << endl;
  
  RooFitResult *fitRes = totalPdf.fitTo(dataHist,
                                       SumW2Error(kTRUE),
                                       ExternalConstraints(constraints),
                                       Save(kTRUE),
                                       PrintLevel(1),
                                       Minimizer("Minuit2", "migrad"),
                                       Strategy(2),
                                       NumCPU(4, kTRUE),
                                       Optimize(2));
  
  // Run HESSE for better errors
  totalPdf.fitTo(dataHist,
                SumW2Error(kTRUE),
                ExternalConstraints(constraints),
                Save(kTRUE),
                PrintLevel(-1),
                Minimizer("Minuit2", "hesse"));
  
  // ------------------------------------------------------------------
  // 11. Run MINOS for asymmetric errors on signal yield
  // ------------------------------------------------------------------
  cout << "\n=== Running MINOS for asymmetric errors on signal ===" << endl;
  
  RooFitResult *minosRes = totalPdf.fitTo(dataHist,
                                         SumW2Error(kTRUE),
                                         ExternalConstraints(constraints),
                                         Save(kTRUE),
                                         PrintLevel(1),
                                         Minimizer("Minuit2", "minos"),
                                         NumCPU(4, kTRUE));
  
  // ------------------------------------------------------------------
  // 12. Print correlation matrix
  // ------------------------------------------------------------------
  TMatrixDSym corr = fitRes->correlationMatrix();
  cout << "\n=== CORRELATION MATRIX ===" << endl;
  cout << "Correlation(N_isr3pi, N_nonReson) = " << corr(0, 1) << endl;
  cout << "Correlation(N_isr3pi, N_eeg) = " << corr(0, 2) << endl;
  cout << "Correlation(N_isr3pi, N_omegapi) = " << corr(0, 3) << endl;
  cout << "Correlation(N_isr3pi, N_ksl) = " << corr(0, 4) << endl;
  
  // Check for high correlations (>0.7)
  if (std::abs(corr(0, 1)) > 0.7) {
    cout << "WARNING: High correlation between N_isr3pi and N_nonReson!" << endl;
  }
  if (std::abs(corr(2, 3)) > 0.7) {
    cout << "WARNING: High correlation between N_eeg and N_omegapi!" << endl;
  }

  // ------------------------------------------------------------------
  // 13. Print fitted yields
  // ------------------------------------------------------------------
  std::cout << "\n=== Fitted yields ===" << std::endl;
  std::cout << "EEG        : " << N_eeg.getVal() << " +/- " << N_eeg.getError() << std::endl;
  std::cout << "ISR3pi     : " << N_isr3pi.getVal() << " +/- " << N_isr3pi.getError() << std::endl;
  std::cout << "Non-reson  : " << N_nonReson.getVal() << " +/- " << N_nonReson.getError() << std::endl;
  std::cout << "OmegaPi    : " << N_omegapi.getVal() << " +/- " << N_omegapi.getError() << std::endl;
  std::cout << "KSL        : " << N_ksl.getVal() << " +/- " << N_ksl.getError() << std::endl;
  std::cout << "MC Rest    : " << N_mcrest.getVal() << " +/- " << N_mcrest.getError() << std::endl;

  // MINOS asymmetric errors for signal
  double signal_err_lo = N_isr3pi.getErrorLo();
  double signal_err_hi = N_isr3pi.getErrorHi();
  std::cout << "\nISR3pi MINOS asymmetric errors: +" << signal_err_hi << " / -" << signal_err_lo << std::endl;

  // Calculate total fitted yield
  double total_fitted = N_eeg.getVal() + N_isr3pi.getVal() + N_nonReson.getVal() + 
                        N_omegapi.getVal() + N_ksl.getVal() + N_mcrest.getVal();
  std::cout << "\nTotal fitted: " << total_fitted << " (data: " << nb_data << ")" << std::endl;
  std::cout << "Difference: " << (total_fitted - nb_data) << " events" << std::endl;

  // ------------------------------------------------------------------
  // 14. Compute scaling factors
  // ------------------------------------------------------------------
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

  std::cout << "\n=== Scaling factors ===" << std::endl;
  std::cout << "EEG        : " << sf_eeg << " +/- " << err_eeg << std::endl;
  std::cout << "ISR3pi     : " << sf_isr3pi << " +/- " << err_isr3pi << std::endl;
  std::cout << "Non-reson  : " << sf_nonReson << " +/- " << err_nonReson << std::endl;
  std::cout << "OmegaPi    : " << sf_omegapi << " +/- " << err_omegapi << std::endl;
  std::cout << "KSL        : " << sf_ksl << " +/- " << err_ksl << std::endl;
  std::cout << "MC Rest    : " << sf_mcrest << " +/- " << err_mcrest << std::endl;

  // ------------------------------------------------------------------
  // 15. Save results to file
  // ------------------------------------------------------------------
  std::ofstream myfile;
  TString myfile_nm = "../header_bdt/sfw2d_bdt.txt";
  myfile.open(myfile_nm.Data(), std::ios::out | std::ios::trunc);
  myfile << "// Fit results with Gaussian constraints and MINOS errors\n";
  myfile << "// Date: " << TDatime().AsString() << "\n\n";
  myfile << "const double eeg_sfw = " << sf_eeg << ";\n"
	 << "const double isr3pi_sfw = " << sf_isr3pi << ";\n"
    	 << "const double nonReson_sfw = " << sf_nonReson << ";\n"
	 << "const double omegapi_sfw = " << sf_omegapi << ";\n"
	 << "const double ksl_sfw = " << sf_ksl << ";\n"
	 << "const double mcrest_sfw = " << sf_mcrest << ";\n\n";
  myfile << "// MINOS asymmetric errors for ISR3pi: +" << signal_err_hi << " / -" << signal_err_lo << "\n";
  myfile.close();
  cout << "\nResults saved to: " << myfile_nm << endl;

  // ------------------------------------------------------------------
  // 16. Create projection plot
  // ------------------------------------------------------------------
  TCanvas *c1 = new TCanvas("c1", "Projection", 800, 600);
  RooPlot *frame = y.frame(Title("Projection onto E_{#gamma_{3}}"));
  dataHist.plotOn(frame, MarkerSize(0.8), MarkerColor(kBlack));
  totalPdf.plotOn(frame, LineColor(kBlue), LineWidth(2));
  totalPdf.plotOn(frame, Components(isr3piPdf), LineStyle(kDashed), LineColor(kRed));
  totalPdf.plotOn(frame, Components(nonResonPdf), LineStyle(kDashed), LineColor(kOrange));
  totalPdf.plotOn(frame, Components(RooArgSet(eegPdf, omegapiPdf, kslPdf, mcrestPdf)), 
                  LineStyle(kDashed), LineColor(kGreen));
  
  TLegend* leg = new TLegend(0.65, 0.65, 0.89, 0.89);
  leg->AddEntry(frame->getObject(0), "Data", "LP");
  leg->AddEntry((TObject*)0, "Total Model", "L");
  leg->AddEntry((TObject*)0, "ISR3pi Peak", "L");
  leg->AddEntry((TObject*)0, "ISR3pi Non-reson", "L");
  leg->AddEntry((TObject*)0, "Backgrounds", "L");
  leg->Draw();
  
  frame->Draw();
  c1->SaveAs("sfw2d_projection.pdf");
  
  // ------------------------------------------------------------------
  // 17. Save to ROOT file
  // ------------------------------------------------------------------
  TFile *fout = new TFile("sfw2d_fit_result.root", "RECREATE");
  fitRes->Write("fitresult");
  if (minosRes) minosRes->Write("minosResult");
  totalPdf.Write();
  dataHist.Write();
  frame->Write();
  fout->Close();
  
  cout << "\n=== FIT COMPLETED SUCCESSFULLY ===" << endl;
  cout << "Output files:" << endl;
  cout << "  - sfw2d_projection.pdf" << endl;
  cout << "  - sfw2d_fit_result.root" << endl;
  cout << "  - ../header_bdt/sfw2d_bdt.txt" << endl;
  
  delete leg;
  delete c1;
  delete fout;
  fin->Close();
  
  return 0;
}
