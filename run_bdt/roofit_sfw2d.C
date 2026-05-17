// fit_3pi_mass_fixed.C
// Full RooFit implementation for 2D template fit (fit variable + 3π mass)
// using extended likelihood with yield parameters.

#include <TFile.h>
#include <TH2D.h>
#include <TList.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <RooFit.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <iostream>
#include <cmath>
#include "../header_bdt/path.h"

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
int roofit_sfw2d() {
    // ------------------------------------------------------------------
    // 1. Open input file and retrieve the 2D histogram list
    // ------------------------------------------------------------------
    // Adjust the file name to your actual path (same as in sfw2d.C)
    TString inputFileName = outputHist + "hist.root";   // modify as needed
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

    // Get all needed 2D histograms (same names as in sfw2d.C)
    TH2D *h_data      = (TH2D*) HSFW2D->FindObject("h2d_sfw_TDATA");
    TH2D *h_eeg       = (TH2D*) HSFW2D->FindObject("h2d_sfw_TEEG");
    TH2D *h_isr3pi    = (TH2D*) HSFW2D->FindObject("h2d_sfw_TISR3PI_SIG");
    TH2D *h_omegapi   = (TH2D*) HSFW2D->FindObject("h2d_sfw_TOMEGAPI");
    TH2D *h_etagam    = (TH2D*) HSFW2D->FindObject("h2d_sfw_TETAGAM");
    TH2D *h_ksl       = (TH2D*) HSFW2D->FindObject("h2d_sfw_TKSL");
    TH2D *h_kpm       = (TH2D*) HSFW2D->FindObject("h2d_sfw_TKPM");
    TH2D *h_rhopi     = (TH2D*) HSFW2D->FindObject("h2d_sfw_TRHOPI");
    TH2D *h_bkgrest   = (TH2D*) HSFW2D->FindObject("h2d_sfw_TBKGREST");

    if (!h_data || !h_eeg || !h_isr3pi || !h_omegapi || !h_etagam || !h_ksl) {
        std::cerr << "ERROR: Missing one or more required histograms" << std::endl;
        return 1;
    }

    // Build combined "mcrest" histogram (same as in sfw2d.C)
    TH2D *h_mcrest = (TH2D*) h_bkgrest->Clone("h_mcrest");
    h_mcrest->Add(h_kpm, 1.0);
    h_mcrest->Add(h_rhopi, 1.0);

    // ------------------------------------------------------------------
    // 2. Add a tiny constant to all MC histograms (to avoid zero PDF values)
    //    Do NOT modify the data histogram!
    // ------------------------------------------------------------------
    const double eps = 1e-12;
    addTinyConstant(h_eeg, eps);
    addTinyConstant(h_isr3pi, eps);
    addTinyConstant(h_omegapi, eps);
    addTinyConstant(h_etagam, eps);
    addTinyConstant(h_ksl, eps);
    addTinyConstant(h_mcrest, eps);
    // h_data remains unchanged

    // ------------------------------------------------------------------
    // 3. Define the observables (x = your 1D fit variable, e.g. chi2 or deltaE;
    //    y = 3π mass). Adjust axis ranges to match your histograms.
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
    RooDataHist isr3piHist("isr3piHist", "ISR3pi MC", RooArgSet(x, y), h_isr3pi);
    RooDataHist omegapiHist("omegapiHist", "OmegaPi MC", RooArgSet(x, y), h_omegapi);
    RooDataHist etagamHist("etagamHist", "EtaGamma MC", RooArgSet(x, y), h_etagam);
    RooDataHist kslHist("kslHist", "KSL MC", RooArgSet(x, y), h_ksl);
    RooDataHist mcrestHist("mcrestHist", "MC Rest", RooArgSet(x, y), h_mcrest);

    // ------------------------------------------------------------------
    // 5. Create PDFs from templates (RooHistPdf)
    // ------------------------------------------------------------------
    RooHistPdf eegPdf("eegPdf", "EEG PDF", RooArgSet(x, y), eegHist);
    RooHistPdf isr3piPdf("isr3piPdf", "ISR3pi PDF", RooArgSet(x, y), isr3piHist);
    RooHistPdf omegapiPdf("omegapiPdf", "OmegaPi PDF", RooArgSet(x, y), omegapiHist);
    RooHistPdf etagamPdf("etagamPdf", "EtaGamma PDF", RooArgSet(x, y), etagamHist);
    RooHistPdf kslPdf("kslPdf", "KSL PDF", RooArgSet(x, y), kslHist);
    RooHistPdf mcrestPdf("mcrestPdf", "MC Rest PDF", RooArgSet(x, y), mcrestHist);

    // ------------------------------------------------------------------
    // 6. Get total MC yields for initial values of the fit parameters
    // ------------------------------------------------------------------
    double nb_eeg_sum     = h_eeg->Integral();
    double nb_isr3pi_sum  = h_isr3pi->Integral();
    double nb_omegapi_sum = h_omegapi->Integral();
    double nb_etagam_sum  = h_etagam->Integral();
    double nb_ksl_sum     = h_ksl->Integral();
    double nb_mcrest_sum  = h_mcrest->Integral();

    // ------------------------------------------------------------------
    // 7. Create yield parameters (extended likelihood will use these as coefficients)
    //    Initial values = MC integral; allow the fit to vary up to 10× the MC integral.
    // ------------------------------------------------------------------
    RooRealVar N_eeg("N_eeg", "EEG yield", nb_eeg_sum, 0.0, nb_eeg_sum*10);
    RooRealVar N_isr3pi("N_isr3pi", "ISR3pi yield", nb_isr3pi_sum, 0.0, nb_isr3pi_sum*10);
    RooRealVar N_omegapi("N_omegapi", "OmegaPi yield", nb_omegapi_sum, 0.0, nb_omegapi_sum*10);
    RooRealVar N_etagam("N_etagam", "EtaGamma yield", nb_etagam_sum, 0.0, nb_etagam_sum*10);
    RooRealVar N_ksl("N_ksl", "KSL yield", nb_ksl_sum, 0.0, nb_ksl_sum*10);
    RooRealVar N_mcrest("N_mcrest", "MC Rest yield", nb_mcrest_sum, 0.0, nb_mcrest_sum*10);

    // ------------------------------------------------------------------
    // 8. Build the total PDF as a sum of PDFs weighted by yields
    //    This automatically makes the fit extended (total yield = sum of yields)
    // ------------------------------------------------------------------
    RooAddPdf totalPdf("totalPdf", "Total model",
                       RooArgList(eegPdf, isr3piPdf, omegapiPdf, etagamPdf, kslPdf, mcrestPdf),
                       RooArgList(N_eeg, N_isr3pi, N_omegapi, N_etagam, N_ksl, N_mcrest));

    // ------------------------------------------------------------------
    // 9. Perform extended binned maximum likelihood fit
    // ------------------------------------------------------------------
    RooFitResult *fitRes = totalPdf.fitTo(dataHist,
                                          SumW2Error(kTRUE),   // include MC statistical errors
                                          Save(kTRUE),
                                          PrintLevel(1));
    fitRes->Print("v");

    // ------------------------------------------------------------------
    // 10. Print fitted yields and derived scaling factors
    // ------------------------------------------------------------------
    std::cout << "\n=== Fitted yields (absolute normalisation) ===" << std::endl;
    std::cout << "EEG        : " << N_eeg.getVal() << " +/- " << N_eeg.getError() << std::endl;
    std::cout << "ISR3pi     : " << N_isr3pi.getVal() << " +/- " << N_isr3pi.getError() << std::endl;
    std::cout << "OmegaPi    : " << N_omegapi.getVal() << " +/- " << N_omegapi.getError() << std::endl;
    std::cout << "EtaGamma   : " << N_etagam.getVal() << " +/- " << N_etagam.getError() << std::endl;
    std::cout << "KSL        : " << N_ksl.getVal() << " +/- " << N_ksl.getError() << std::endl;
    std::cout << "MC Rest    : " << N_mcrest.getVal() << " +/- " << N_mcrest.getError() << std::endl;

    // Compute scaling factors (Data/MC normalisation) and their errors
    auto computeSF = [](double N_fit, double N_mc, double err_fit, double err_mc=0.0) {
        double sf = N_fit / N_mc;
        double rel_err = std::hypot(err_fit/N_fit, err_mc/N_mc);
        return std::make_pair(sf, sf * rel_err);
    };
    auto [sf_eeg, err_eeg]     = computeSF(N_eeg.getVal(), nb_eeg_sum, N_eeg.getError());
    auto [sf_isr3pi, err_isr3pi] = computeSF(N_isr3pi.getVal(), nb_isr3pi_sum, N_isr3pi.getError());
    auto [sf_omegapi, err_omegapi] = computeSF(N_omegapi.getVal(), nb_omegapi_sum, N_omegapi.getError());
    auto [sf_etagam, err_etagam] = computeSF(N_etagam.getVal(), nb_etagam_sum, N_etagam.getError());
    auto [sf_ksl, err_ksl]     = computeSF(N_ksl.getVal(), nb_ksl_sum, N_ksl.getError());
    auto [sf_mcrest, err_mcrest] = computeSF(N_mcrest.getVal(), nb_mcrest_sum, N_mcrest.getError());

    std::cout << "\n=== Scaling factors (Data / MC integral) ===" << std::endl;
    std::cout << "EEG        : " << sf_eeg << " +/- " << err_eeg << std::endl;
    std::cout << "ISR3pi     : " << sf_isr3pi << " +/- " << err_isr3pi << std::endl;
    std::cout << "OmegaPi    : " << sf_omegapi << " +/- " << err_omegapi << std::endl;
    std::cout << "EtaGamma   : " << sf_etagam << " +/- " << err_etagam << std::endl;
    std::cout << "KSL        : " << sf_ksl << " +/- " << err_ksl << std::endl;
    std::cout << "MC Rest    : " << sf_mcrest << " +/- " << err_mcrest << std::endl;

    // ------------------------------------------------------------------
    // 11. Project the fit onto the Eisr axis (integrate over x)
    // ------------------------------------------------------------------
    TCanvas *c1 = new TCanvas("c1", "Eisr projection", 800, 600);
    RooPlot *EisrFrame = y.frame(Title("Projection onto M_{3#pi}"));
    dataHist.plotOn(EisrFrame);
    totalPdf.plotOn(EisrFrame, ProjWData(dataHist));   // project using data as integration reference
    // Optionally add individual components (use different line styles)
    totalPdf.plotOn(EisrFrame, Components(eegPdf), LineStyle(kDashed), LineColor(kRed));
    totalPdf.plotOn(EisrFrame, Components(isr3piPdf), LineStyle(kDotted), LineColor(kBlue));
    totalPdf.plotOn(EisrFrame, Components(omegapiPdf), LineStyle(kDashDotted), LineColor(kGreen));
    // ... add others if desired
    EisrFrame->Draw();
    c1->SaveAs("fit_Eisr_projection.pdf");

    // ------------------------------------------------------------------
    // 12. Save the workspace and results
    // ------------------------------------------------------------------
    TFile *fout = new TFile("fit_2d_result.root", "RECREATE");
    fitRes->Write();
    totalPdf.Write();
    dataHist.Write();
    fout->Close();

    // Cleanup
    delete c1;
    fin->Close();
    return 0;
}
