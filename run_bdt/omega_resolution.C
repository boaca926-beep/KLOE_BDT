// ============================================================================
// omega_resolution.C
//
// Fit the ω peak in M3π with a Gaussian + constant for both MC and
// background‑subtracted data. Extracts the Gaussian sigma (resolution).
// Accepts fits even with large χ² (due to natural Breit‑Wigner width).
//
// Usage:
//   .x omega_resolution.C("tuning_false")
//   .x omega_resolution.C("tuning_true")
// ============================================================================

#include "../header_method/method.h"
#include "../header_plot/plot.h"
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>

struct FitResult {
    TString name;
    double mean, mean_err;
    double sigma, sigma_err;
    double chi2_ndf;
    int entries;
};

// ----------------------------------------------------------------------------
// Fit Gaussian + constant (or pure Gaussian) – accepts large χ²
// ----------------------------------------------------------------------------
bool fitGaussian(TH1D *h, FitResult &res, const TString &label = "") {
    if (!h || h->GetEntries() < 50) return false;

    // Estimate initial parameters
    double peak = h->GetMaximum();
    double mean0 = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
    double rms = h->GetRMS();
    double sigma_guess = (rms > 0) ? rms * 0.6 : 3.0;
    if (sigma_guess < 0.5) sigma_guess = 1.0;
    // Clamp mean to reasonable range
    if (mean0 < 775) mean0 = 782;
    if (mean0 > 790) mean0 = 782;

    // Strategies: core fits only
    const int nStrategies = 4;
    struct Strategy {
        double xmin, xmax;
        TString func;
        int npar;
        double sigma_min, sigma_max;
    } strategies[nStrategies] = {
      {775.0, 790.0, "gaus(0) + pol0(3)", 4, 0.1, 20.0},
      {775.0, 790.0, "gaus(0)", 3, 0.1, 20.0},
      {775.0, 790.0, "gaus(0) + pol0(3)", 4, 0.1, 20.0},
      {775.0, 790.0, "gaus(0)", 3, 0.1, 20.0}
    };

    for (int s = 0; s < nStrategies; ++s) {
        double xmin = strategies[s].xmin;
        double xmax = strategies[s].xmax;
        TString funcStr = strategies[s].func;
        int npar = strategies[s].npar;
        double sigma_min = strategies[s].sigma_min;
        double sigma_max = strategies[s].sigma_max;

        TF1 *fit = new TF1("omega_fit", funcStr, xmin, xmax);
        fit->SetParameter(0, peak);
        fit->SetParameter(1, mean0);
        fit->SetParameter(2, sigma_guess);
        for (int i = 3; i < npar; ++i) fit->SetParameter(i, 0.0);

        fit->SetParLimits(0, 0, peak * 5);
        fit->SetParLimits(1, mean0 - 2.0, mean0 + 2.0);
        fit->SetParLimits(2, sigma_min, sigma_max);
        if (funcStr.Contains("pol0")) {
            fit->SetParLimits(3, -peak/2, peak/2);
        }

        Int_t status = h->Fit(fit, "RQS");
        // Ignore chi2, just check convergence and physical sigma
        double sigma_val = fit->GetParameter(2);
        if ((status == 0 || status == 1) && sigma_val > sigma_min + 0.01 && sigma_val < sigma_max - 0.01) {
            res.mean = fit->GetParameter(1);
            res.mean_err = fit->GetParError(1);
            res.sigma = sigma_val;
            res.sigma_err = fit->GetParError(2);
            double chi2 = fit->GetChisquare();
            int ndf = fit->GetNDF();
            res.chi2_ndf = (ndf > 0) ? chi2 / ndf : 0;
            res.entries = h->GetEntries();
            std::cout << label << " Fit succeeded (strategy " << s
                      << ", range " << xmin << "-" << xmax << ", " << funcStr
                      << "): sigma = " << res.sigma << " +/- " << res.sigma_err
                      << ", chi2/ndf = " << res.chi2_ndf << std::endl;
            delete fit;
            return true;
        }
        delete fit;
    }

    std::cerr << label << " All fit strategies failed." << std::endl;
    return false;
}

// ----------------------------------------------------------------------------
// Main macro
// ----------------------------------------------------------------------------
int omega_resolution(const TString tuning_type = "raw_false",
                     const TString var_nm = "m3pi_bdt",
                     const TString var_symb = "M_{3#pi} [MeV/c^{2}]") {

    const TString tree_file_nm = "../" + tuning_type + "_" + var_nm + "/hist.root";
    const TString out_dir = "../omega_resolution_" + tuning_type;

    gErrorIgnoreLevel = kError;
    TGaxis::SetMaxDigits(4);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0.8);
    TH1::SetDefaultSumw2();

    gSystem->mkdir(out_dir, kTRUE);

    TFile* tree_file = new TFile(tree_file_nm);
    if (!tree_file || tree_file->IsZombie()) {
        std::cerr << "ERROR: Cannot open " << tree_file_nm << std::endl;
        return 1;
    }

    checkFile(tree_file); // optional

    // ---- Retrieve histograms ----
    TH1D *hist_eeg = (TH1D*)tree_file->Get("hist_eeg_sc");
    TH1D *hist_signal = (TH1D*)tree_file->Get("hist_isr3pi_sc");
    TH1D *hist_omegapi = (TH1D*)tree_file->Get("hist_omegapi_sc");
    TH1D *hist_nonreson = (TH1D*)tree_file->Get("hist_nonreson_sc");
    TH1D *hist_ksl = (TH1D*)tree_file->Get("hist_ksl_sc");
    TH1D *hist_mcrest = (TH1D*)tree_file->Get("hist_mcrest_sc");
    TH1D *hist_data = (TH1D*)tree_file->Get("hist_data");

    if (!hist_data || !hist_signal) {
        std::cerr << "ERROR: Missing histograms." << std::endl;
        tree_file->Close();
        return 1;
    }

    // ---- Build background sum and subtract ----
    TH1D *hist_bkg_sum = (TH1D*)hist_data->Clone("hist_bkg_sum");
    hist_bkg_sum->Reset();
    if (hist_eeg)     hist_bkg_sum->Add(hist_eeg, 1.0);
    if (hist_omegapi) hist_bkg_sum->Add(hist_omegapi, 1.0);
    if (hist_nonreson)hist_bkg_sum->Add(hist_nonreson, 1.0);
    if (hist_ksl)     hist_bkg_sum->Add(hist_ksl, 1.0);
    if (hist_mcrest)  hist_bkg_sum->Add(hist_mcrest, 1.0);

    TH1D *hist_data_sub = (TH1D*)hist_data->Clone("hist_data_sub");
    hist_data_sub->Add(hist_bkg_sum, -1.0);

    // ------------------------------------------------------------------
    // Fit MC signal and background-subtracted data with Gaussian
    // ------------------------------------------------------------------
    FitResult resMC, resData;
    bool okMC = fitGaussian(hist_signal, resMC, "MC");

    // Try background-subtracted first
    bool okData = fitGaussian(hist_data_sub, resData, "Data (bkg sub)");
    TH1D *h_data_used = (okData) ? hist_data_sub : nullptr;

    // Fallback: if background-subtracted fails, fit raw data
    if (!okData) {
        std::cout << "Background-subtracted fit failed. Trying raw data..." << std::endl;
        okData = fitGaussian(hist_data, resData, "Data (raw)");
        if (okData) h_data_used = hist_data;
    }

    if (!okMC || !okData) {
        std::cerr << "ERROR: Fit results invalid; cannot compute ratio." << std::endl;
        tree_file->Close();
        return 1;
    }

    // ---- Compute ratio ----
    double R = resData.sigma / resMC.sigma;
    double R_err = R * TMath::Sqrt(
        TMath::Power(resData.sigma_err / resData.sigma, 2) +
        TMath::Power(resMC.sigma_err / resMC.sigma, 2)
    );

    std::cout << "\n========================================" << std::endl;
    std::cout << "Omega Gaussian resolution results:" << std::endl;
    std::cout << "MC:   sigma = " << resMC.sigma << " +/- " << resMC.sigma_err << " MeV" << std::endl;
    std::cout << "Data: sigma = " << resData.sigma << " +/- " << resData.sigma_err << " MeV" << std::endl;
    std::cout << "Ratio R = " << R << " +/- " << R_err << std::endl;
    std::cout << "========================================" << std::endl;

    // ---- Write parameters to header file ----
    std::ofstream myfile;
    TString myfile_nm = "../trackmass_scan/omega_params.h";
    gSystem->mkdir("../trackmass_scan", kTRUE);
    myfile.open(myfile_nm.Data());
    myfile << "// Omega Gaussian resolution parameters\n";
    myfile << "const double SIGMA_MC_OMEGA   = " << resMC.sigma << ";\n";
    myfile << "const double SIGMA_MC_OMEGA_ERR = " << resMC.sigma_err << ";\n";
    myfile << "const double SIGMA_DATA_OMEGA = " << resData.sigma << ";\n";
    myfile << "const double SIGMA_DATA_OMEGA_ERR = " << resData.sigma_err << ";\n";
    myfile << "const double R_OMEGA = " << R << ";\n";
    myfile << "const double R_OMEGA_ERR = " << R_err << ";\n";
    myfile.close();
    std::cout << "Parameters written to " << myfile_nm << std::endl;

    // ---- Draw and save individual plots ----
    auto drawHist = [&](TH1D *h, const FitResult &res, const TString &title, const TString &fname) {
        TCanvas *c = new TCanvas("c", title, 900, 700);
        gPad->SetBottomMargin(0.15);
        gPad->SetLeftMargin(0.15);
        h->SetLineColor(kBlack);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.6);
        h->GetXaxis()->SetTitle(var_symb);
        h->GetYaxis()->SetTitle("Events");
        h->GetYaxis()->SetRangeUser(0.01, h->GetMaximum() * 1.6);
        h->Draw("E0");
        TF1 *fit = (TF1*)h->GetFunction("omega_fit");
        if (fit) {
            fit->SetLineColor(kRed);
            fit->SetLineWidth(2);
            fit->Draw("same");
        }
        TPaveText *pt = new TPaveText(0.35, 0.65, 0.85, 0.85, "NDC");
        pt->SetFillColor(0);
        pt->SetBorderSize(0);
        pt->SetTextAlign(12);
        pt->SetTextSize(0.04);
	cout << res.sigma << endl;
        pt->AddText(Form("#sigma = %.3f #pm %.3f MeV/c^{2}", res.sigma, res.sigma_err));
        pt->AddText(Form("#chi^{2}/NDF = %.2f", res.chi2_ndf));
        pt->Draw();
        c->SaveAs(fname);
        delete c;
    };

    drawHist(hist_signal, resMC, "MC #omega peak", out_dir + "/omega_fit_MC.pdf");
    drawHist(h_data_used, resData, "Data #omega peak", out_dir + "/omega_fit_Data.pdf");

    // ---- Comparison plot: MC vs Data overlaid (raw yields, no normalization) ----
    TCanvas *c_comp = new TCanvas("c_comp", "MC vs Data #omega peak", 900, 900);
    gPad->SetBottomMargin(0.15);
    gPad->SetLeftMargin(0.15);

    // Clone histograms (without scaling) to avoid modifying originals
    TH1D *h_mc_draw = (TH1D*)hist_signal->Clone("h_mc_draw");
    TH1D *h_data_draw = (TH1D*)h_data_used->Clone("h_data_draw");

    // Set styles for MC (black markers)
    h_mc_draw->SetMarkerColor(kBlue);
    h_mc_draw->SetLineColor(kBlue);
    h_mc_draw->SetMarkerStyle(20);
    h_mc_draw->SetMarkerSize(0.8);

    // Set styles for Data (blue markers)
    h_data_draw->SetMarkerColor(kBlack);
    h_data_draw->SetLineColor(kBlack);
    h_data_draw->SetMarkerStyle(21);
    h_data_draw->SetMarkerSize(0.8);

    // Set Y range to cover the maximum of both histograms
    double y_max = TMath::Max(h_mc_draw->GetMaximum(), h_data_draw->GetMaximum());
    h_mc_draw->SetMarkerStyle(20);
    h_mc_draw->SetMarkerSize(0.6);
    h_mc_draw->GetYaxis()->SetTitle("Events");
    h_mc_draw->GetYaxis()->SetRangeUser(0.01, y_max * 1.2);
    h_mc_draw->GetYaxis()->CenterTitle();
    h_mc_draw->GetYaxis()->SetTitleSize(0.05);
    h_mc_draw->GetYaxis()->SetTitleOffset(1.4);
    h_mc_draw->GetYaxis()->SetLabelSize(0.04);
    h_mc_draw->GetXaxis()->SetTitle(var_symb);
    h_mc_draw->GetXaxis()->SetTitleSize(0.05);
    h_mc_draw->GetXaxis()->SetTitleOffset(1.2);
    h_mc_draw->GetXaxis()->SetLabelSize(0.04);
    h_mc_draw->GetXaxis()->CenterTitle();
    
    // Draw MC first, then Data on top
    h_mc_draw->Draw("hist");
    h_data_draw->Draw("E0 same");

    // Retrieve fit functions (raw amplitudes, no scaling)
    TF1 *fit_mc = (TF1*)hist_signal->GetFunction("omega_fit")->Clone("fit_mc_clone");
    if (fit_mc) {
        fit_mc->SetLineColor(kRed);
        fit_mc->SetLineWidth(2);
        fit_mc->Draw("same");
    }

    TF1 *fit_data = (TF1*)h_data_used->GetFunction("omega_fit")->Clone("fit_data_clone");
    if (fit_data) {
        fit_data->SetLineColor(kGreen + 2);
        fit_data->SetLineWidth(2);
        fit_data->Draw("same");
    }

    // Legend with sigma values
    TLegend *leg = new TLegend(0.40, 0.65, 0.7, 0.85);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(h_mc_draw, "MC", "lep");
    leg->AddEntry(h_data_draw, "Data (bkg sub)", "lep");
    leg->AddEntry(fit_mc, Form("MC #sigma = %.3f MeV", resMC.sigma), "l");
    leg->AddEntry(fit_data, Form("Data #sigma = %.3f MeV", resData.sigma), "l");
    leg->Draw();

    c_comp->SaveAs(out_dir + "/omega_comparison.pdf");
    delete c_comp;

    tree_file->Close();
    return 0;
}
