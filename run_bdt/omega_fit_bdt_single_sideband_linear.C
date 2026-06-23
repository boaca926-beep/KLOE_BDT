// omega_fit_bdt.C – template fit with proper memory management
// All masses in MeV
// Non-resonant template + Chebyshev 2nd order background with single sideband

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TFitResult.h>
#include <iostream>
#include <cmath>
#include <vector>

#include "../header_bdt/sfw2d.txt"
#include "../header_bdt/omega_fit_bdt.h"

// Global pointers for templates (needed by fit function)
TH1D *gSigTemplate = nullptr;
TH1D *gBkgTemplate = nullptr;

// Fit function: α * signal_template + β * nonReson_template + Chebyshev polynomial (2nd order)
Double_t template_chebyshev2_with_nonreson(Double_t *x, Double_t *par) {
    int bin_sig = gSigTemplate->FindBin(x[0]);
    int bin_bkg = gBkgTemplate->FindBin(x[0]);
    Double_t sig = gSigTemplate->GetBinContent(bin_sig);
    Double_t nonReson = gBkgTemplate->GetBinContent(bin_bkg);
    // Chebyshev: T0=1, T1=x, T2=2x²-1
    Double_t X = x[0];
    Double_t cheb = par[2] + par[3]*X + par[4]*(2*X*X - 1);
    return par[0] * sig + par[1] * nonReson + cheb;
}

void omega_fit_bdt_single_sideband_linear() {
    // ------------------------------------------------------------------
    // 1. Open tree file
    // ------------------------------------------------------------------
    TFile *ftree = TFile::Open(treeFile);
    if (!ftree || ftree->IsZombie()) {
        std::cerr << "ERROR: cannot open " << treeFile << std::endl;
        return;
    }

    TTree *tdata = (TTree*) ftree->Get("TDATA");
    if (!tdata) { std::cerr << "ERROR: TDATA not found." << std::endl; return; }

    double low, high;
    int nbins = NBINS;

    low  = MASS_MIN;          // e.g., 600.0 MeV
    high = MASS_MAX;          // e.g., 1000.0 MeV

    std::cout << "Mass range [" << low << ", " << high << "] MeV/c²\n";

    // ------------------------------------------------------------------
    // 2. Data histogram (detached from file)
    // ------------------------------------------------------------------
    TH1D *h_data = new TH1D("h_data", "", nbins, low, high);
    h_data->Sumw2();
    h_data->SetDirectory(0);
    double mtest;
    tdata->SetBranchAddress("Br_m3pi_bdt", &mtest);
    for (Long64_t i = 0; i < tdata->GetEntries(); ++i) {
        tdata->GetEntry(i);
        h_data->Fill(mtest);
    }
    std::cout << "Data integral: " << h_data->Integral() << std::endl;

    // ------------------------------------------------------------------
    // 3. Load scaled MC components (detach each histogram)
    // ------------------------------------------------------------------
    auto makeScaledHist = [&](const char* tname, double scale, int color = 1, int style = 1) -> TH1D* {
        TTree *t = (TTree*) ftree->Get(tname);
        if (!t) return nullptr;
        if (!t->GetBranch("Br_m3pi_bdt")) return nullptr;
        TH1D *h = new TH1D(Form("h_%s", tname), "", nbins, low, high);
        h->Sumw2();
        h->SetDirectory(0);
        double val;
        t->SetBranchAddress("Br_m3pi_bdt", &val);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            h->Fill(val);
        }
        h->Scale(scale);
        h->SetLineColor(color);
        h->SetLineStyle(style);
        h->SetLineWidth(2);
        return h;
    };

    // Same color codes and line styles as correct_and_plot.C
    TH1D *h_eeg       = makeScaledHist("TEEG", eeg_sfw * 2., 6, 7);
    TH1D *h_isr3pi    = makeScaledHist("TISR3PI_SIG_PEAK", isr3pi_sfw, 4, 2);
    TH1D *h_nonReson  = makeScaledHist("TISR3PI_SIG_NON_RESON", nonReson_sfw, 2, 3);
    TH1D *h_omegapi   = makeScaledHist("TOMEGAPI", omegapi_sfw, 7, 5);
    TH1D *h_ksl       = makeScaledHist("TKSL", ksl_sfw, 28, 4);

    // MC Rest (includes TETAGAM)
    TH1D *h_mcrest = new TH1D("h_mcrest", "", nbins, low, high);
    h_mcrest->Sumw2();
    h_mcrest->SetDirectory(0);
    for (const char* name : {"TKPM", "TRHOPI", "TBKGREST", "TETAGAM"}) {
        TTree *t = (TTree*) ftree->Get(name);
        if (!t) continue;
        if (!t->GetBranch("Br_m3pi_bdt")) continue;
        double val;
        t->SetBranchAddress("Br_m3pi_bdt", &val);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            h_mcrest->Fill(val);
        }
    }
    h_mcrest->Scale(mcrest_sfw);
    h_mcrest->SetLineColor(37);
    h_mcrest->SetLineStyle(6);
    h_mcrest->SetLineWidth(2);

    if (!h_isr3pi) { std::cerr << "ERROR: No ISR3pi histogram." << std::endl; return; }
    if (!h_nonReson) { std::cerr << "ERROR: No nonReson histogram." << std::endl; return; }

    // Summary
    std::cout << "EEG integral: " << h_eeg->Integral() << ", sfw = " << eeg_sfw * 2. << "\n"
              << "ISR3PI integral: " << h_isr3pi->Integral() << ", sfw = " << isr3pi_sfw << "\n"
              << "OMEGAPI integral: " << h_omegapi->Integral() << ", sfw = " << omegapi_sfw << "\n"
              << "NonReson integral: " << h_nonReson->Integral() << ", sfw = " << nonReson_sfw << "\n"
              << "KSL integral: " << h_ksl->Integral() << ", sfw = " << ksl_sfw << "\n"
              << "MCREST integral: " << h_mcrest->Integral() << ", sfw = " << mcrest_sfw << "\n"
              << std::endl;

    // ------------------------------------------------------------------
    // 4. Subtract backgrounds -> h_data_isr (detached)
    // ------------------------------------------------------------------
    TH1D *h_data_isr = (TH1D*) h_data->Clone("h_data_isr");
    h_data_isr->SetDirectory(0);

    // Propagate errors properly
    for (int bin = 1; bin <= h_data_isr->GetNbinsX(); ++bin) {
        double data_val = h_data->GetBinContent(bin);
        double data_err = h_data->GetBinError(bin);
        double bkg_total = 0;
        double bkg_err_sq = 0;

        if (h_eeg) {
            bkg_total += h_eeg->GetBinContent(bin);
            bkg_err_sq += pow(h_eeg->GetBinError(bin), 2);
        }
        if (h_omegapi) {
            bkg_total += h_omegapi->GetBinContent(bin);
            bkg_err_sq += pow(h_omegapi->GetBinError(bin), 2);
        }
        if (h_ksl) {
            bkg_total += h_ksl->GetBinContent(bin);
            bkg_err_sq += pow(h_ksl->GetBinError(bin), 2);
        }
        if (h_mcrest) {
            bkg_total += h_mcrest->GetBinContent(bin);
            bkg_err_sq += pow(h_mcrest->GetBinError(bin), 2);
        }
        // nonReson is NOT subtracted (it's used in the fit)

        double new_val = data_val - bkg_total;
        double new_err = sqrt(data_err * data_err + bkg_err_sq);

        h_data_isr->SetBinContent(bin, new_val);
        h_data_isr->SetBinError(bin, new_err);

        if (h_data_isr->GetBinContent(bin) < 0) {
            h_data_isr->SetBinContent(bin, 0);
        }
    }

    // ------------------------------------------------------------------
    // 5. Normalize templates to unit area
    // ------------------------------------------------------------------
    double sig_int = h_isr3pi->Integral();
    double bkg_int = h_nonReson->Integral();
    if (sig_int > 0) h_isr3pi->Scale(1.0 / sig_int);
    if (bkg_int > 0) h_nonReson->Scale(1.0 / bkg_int);

    // Set global pointers for the fit function
    gSigTemplate = h_isr3pi;
    gBkgTemplate = h_nonReson;

    // ------------------------------------------------------------------
    // 6. Single sideband Chebyshev fit (2nd order) + non-resonant template
    // ------------------------------------------------------------------
    double peak_low = 740.0;
    double peak_high = 820.0;
    
    // SINGLE SIDEBAND ABOVE 800 MeV
    double sb_low = 850.0;    // MeV
    double sb_high = 950.0;   // MeV

    std::cout << "\nSingle sideband Chebyshev fit (2nd order) + non-resonant template:\n"
              << "  Sideband: " << sb_low << "-" << sb_high << " MeV\n"
              << "  Fit range: " << peak_low << "-" << peak_high << " MeV\n";

    // Clone and zero peak region
    TH1D *h_side = (TH1D*) h_data_isr->Clone("h_side");
    for (int bin = 1; bin <= h_side->GetNbinsX(); ++bin) {
        double x = h_side->GetBinCenter(bin);
        if (x >= peak_low && x <= peak_high) h_side->SetBinContent(bin, 0);
    }

    // 2nd‑order Chebyshev polynomial: c0 + c1*x + c2*(2x²-1)
    TF1 *bkg_cheb = new TF1("bkg_cheb", 
        "[0] + [1]*x + [2]*(2*x*x - 1)", 
        sb_low, sb_high);

    // Fit over the single sideband (Chebyshev only)
    h_side->Fit(bkg_cheb, "QN", "", sb_low, sb_high);

    double cheb_c0 = bkg_cheb->GetParameter(0);
    double cheb_c1 = bkg_cheb->GetParameter(1);
    double cheb_c2 = bkg_cheb->GetParameter(2);
    std::cout << "  Chebyshev (initial): c0 = " << cheb_c0 << ", c1 = " << cheb_c1 
              << ", c2 = " << cheb_c2 << "\n";

    // Total fit function (5 parameters: alpha, beta, c0, c1, c2)
    TF1 *total_func = new TF1("total_func", template_chebyshev2_with_nonreson, low, high, 5);
    total_func->SetParameters(1000, 1000, cheb_c0, cheb_c1, cheb_c2);
    total_func->SetParNames("alpha", "beta", "c0", "c1", "c2");

    // Parameter limits with fallback
    auto setLim = [](TF1 *f, int ipar, double val) {
        double lim = 0.3 * std::abs(val);
        if (lim < 1e-6) lim = 1.0;
        f->SetParLimits(ipar, val - lim, val + lim);
    };
    setLim(total_func, 2, cheb_c0);
    setLim(total_func, 3, cheb_c1);
    setLim(total_func, 4, cheb_c2);
    // Limit alpha and beta to reasonable ranges
    total_func->SetParLimits(0, 0, 200000);  // alpha
    total_func->SetParLimits(1, 0, 100000);  // beta

    std::cout << "\nFitting from " << peak_low << " to " << peak_high << " MeV\n";

    TFitResultPtr r = h_data_isr->Fit(total_func, "RQS", "", peak_low, peak_high);

    double alpha = total_func->GetParameter(0);
    double beta = total_func->GetParameter(1);
    double poly_c0 = total_func->GetParameter(2);
    double poly_c1 = total_func->GetParameter(3);
    double poly_c2 = total_func->GetParameter(4);

    double chi2 = r->Chi2();
    int ndf = r->Ndf();
    double chi2_ndf = chi2 / ndf;

    std::cout << "\nFit results: α = " << alpha << ", β = " << beta
              << ", Chebyshev: c0 = " << poly_c0 << ", c1 = " << poly_c1 
              << ", c2 = " << poly_c2 << "\n";
    std::cout << "Fit quality: χ² = " << chi2 << ", ndf = " << ndf 
              << ", χ²/ndf = " << chi2_ndf << std::endl;

    if (chi2_ndf > 5.0) {
        std::cerr << "\nWARNING: Poor fit quality (χ²/ndf = " << chi2_ndf << ")" << std::endl;
    }

    // ------------------------------------------------------------------
    // 7. Create signal & background histograms (scaled)
    // ------------------------------------------------------------------
    TH1D *h_signal = (TH1D*) h_isr3pi->Clone("h_signal");
    h_signal->SetDirectory(0);
    h_signal->Scale(alpha);
    h_signal->SetLineColor(kBlue);
    h_signal->SetLineWidth(2);

    TH1D *h_nonReson_scaled = (TH1D*) h_nonReson->Clone("h_nonReson_scaled");
    h_nonReson_scaled->SetDirectory(0);
    h_nonReson_scaled->Scale(beta);
    h_nonReson_scaled->SetLineColor(kOrange + 1);
    h_nonReson_scaled->SetLineWidth(2);
    h_nonReson_scaled->SetLineStyle(5);

    TH1D *h_background_cheb = (TH1D*) h_data_isr->Clone("h_background_cheb");
    h_background_cheb->SetDirectory(0);
    h_background_cheb->Reset();
    for (int bin = 1; bin <= h_background_cheb->GetNbinsX(); ++bin) {
        double x = h_background_cheb->GetBinCenter(bin);
        double val = poly_c0 + poly_c1*x + poly_c2*(2*x*x - 1);
        if (val < 0) val = 0;
        h_background_cheb->SetBinContent(bin, val);
        h_background_cheb->SetBinError(bin, TMath::Sqrt(val));
    }
    h_background_cheb->SetLineColor(kRed);
    h_background_cheb->SetLineStyle(3);
    h_background_cheb->SetLineWidth(2);

    // Total background (nonReson + Chebyshev)
    TH1D *h_background_total = (TH1D*) h_nonReson_scaled->Clone("h_background_total");
    h_background_total->Add(h_background_cheb);
    h_background_total->SetLineColor(kMagenta);
    h_background_total->SetLineStyle(4);
    h_background_total->SetLineWidth(2);

    // ------------------------------------------------------------------
    // 8. Compute correction weights (using h_signal as desired shape)
    // ------------------------------------------------------------------
    TH1D *h_weight = (TH1D*) h_data_isr->Clone("h_weight");
    h_weight->SetDirectory(0);
    h_weight->Reset();
    for (int bin = 1; bin <= h_weight->GetNbinsX(); ++bin) {
        double x = h_weight->GetBinCenter(bin);
        double mc_val = h_isr3pi->GetBinContent(bin);
        double desired = h_signal->GetBinContent(bin);
        if (x >= peak_low && x <= peak_high && mc_val > 0 && desired > 0) {
            double ratio = desired / mc_val;
            if (ratio > 2.5) ratio = 2.5;
            if (ratio < 0.4) ratio = 0.4;
            h_weight->SetBinContent(bin, ratio);
        } else {
            h_weight->SetBinContent(bin, 1.0);
        }
    }
    // Smooth weights
    TH1D *h_weight_smooth = (TH1D*) h_weight->Clone("h_weight_smooth");
    h_weight_smooth->SetDirectory(0);
    for (int bin = 2; bin <= h_weight_smooth->GetNbinsX() - 1; ++bin) {
        double w_avg = (h_weight->GetBinContent(bin - 1) +
                        h_weight->GetBinContent(bin) +
                        h_weight->GetBinContent(bin + 1)) / 3.0;
        h_weight_smooth->SetBinContent(bin, w_avg);
    }

    // ------------------------------------------------------------------
    // 9. Apply correction to ISR3pi MC (use original scaled MC)
    // ------------------------------------------------------------------
    TH1D *h_isr3pi_orig = makeScaledHist("TISR3PI_SIG_PEAK", isr3pi_sfw);
    if (!h_isr3pi_orig) {
        std::cerr << "ERROR: Cannot reload ISR3pi for correction" << std::endl;
        return;
    }

    TH1D *h_isr3pi_corrected = (TH1D*) h_isr3pi_orig->Clone("h_isr3pi_corrected");
    h_isr3pi_corrected->SetDirectory(0);
    for (int bin = 1; bin <= h_isr3pi_corrected->GetNbinsX(); ++bin) {
        double w = h_weight_smooth->GetBinContent(bin);
        double old = h_isr3pi_corrected->GetBinContent(bin);
        h_isr3pi_corrected->SetBinContent(bin, old * w);
        double err = h_isr3pi_corrected->GetBinError(bin);
        h_isr3pi_corrected->SetBinError(bin, err * w);
    }
    // Renormalise
    double orig_int = h_isr3pi_orig->Integral(peak_low, peak_high);
    double new_int = h_isr3pi_corrected->Integral(peak_low, peak_high);
    if (new_int > 0 && orig_int > 0) {
        double renorm = orig_int / new_int;
        h_isr3pi_corrected->Scale(renorm);
        std::cout << "Renormalisation factor: " << renorm << std::endl;
    }
    h_isr3pi_corrected->SetLineColor(kGreen);
    h_isr3pi_corrected->SetLineWidth(2);

    // ------------------------------------------------------------------
    // 10. Build total MC sum for plotting
    // ------------------------------------------------------------------
    std::vector<TH1D*> comps;
    comps.push_back(h_eeg);
    comps.push_back(h_omegapi);
    comps.push_back(h_ksl);
    comps.push_back(h_mcrest);
    comps.push_back(h_signal);
    comps.push_back(h_nonReson_scaled);
    comps.push_back(h_background_cheb);

    TH1D *h_mc_total = (TH1D*) h_mcrest->Clone("h_mc_total");
    h_mc_total->Reset();
    h_mc_total->Sumw2();
    for (auto h : comps) if (h) h_mc_total->Add(h);
    h_mc_total->SetLineColor(kRed);
    h_mc_total->SetLineWidth(2);

    // ------------------------------------------------------------------
    // 11. Pull distribution (full range)
    // ------------------------------------------------------------------
    TH1D *h_pull = new TH1D("h_pull", "", nbins, low, high);
    for (int bin = 1; bin <= nbins; ++bin) {
        double d = h_data->GetBinContent(bin);
        double m = h_mc_total->GetBinContent(bin);
        double err = std::sqrt(d + m);
        if (err > 0) h_pull->SetBinContent(bin, (d - m) / err);
        else h_pull->SetBinContent(bin, 0);
    }
    h_pull->SetMarkerStyle(20);
    h_pull->SetMarkerSize(0.6);
    h_pull->SetLineWidth(0);

    // ------------------------------------------------------------------
    // 12. Main plotting
    // ------------------------------------------------------------------
    TCanvas *c1 = new TCanvas("c1", "3π mass projection (nonReson + Chebyshev 2nd order)", 1200, 700);
    c1->SetBottomMargin(0.13);
    c1->SetLeftMargin(0.12);

    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.12);
    pad1->Draw();
    pad1->cd();

    double max_val = h_data->GetMaximum();
    max_val = std::max(max_val, h_mc_total->GetMaximum());
    for (auto h : comps) if (h) max_val = std::max(max_val, h->GetMaximum());
    h_data->GetYaxis()->SetRangeUser(0, max_val * 1.2);
    double bin_width = h_data->GetBinWidth(1);

    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.6);
    h_data->Draw("E0");
    h_mc_total->Draw("hist same");
    for (auto h : comps) if (h) h->Draw("hist same");

    h_data->GetYaxis()->SetTitle(Form("Events / [%.1f MeV/c^{2}]", bin_width));
    h_data->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
    h_data->GetYaxis()->CenterTitle();
    h_data->GetXaxis()->SetTitleSize(0.06);
    h_data->GetYaxis()->SetTitleSize(0.07);
    h_data->GetYaxis()->SetTitleOffset(0.7);
    h_data->GetYaxis()->SetLabelSize(0.06);
    h_data->GetYaxis()->SetNdivisions(505);
    h_data->GetXaxis()->SetLabelOffset(0.2);

    // Legend
    TLegend *leg = new TLegend(0.50, 0.15, 0.9, 0.9);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(h_data, "Data", "lep");
    leg->AddEntry(h_mc_total, "Total MC", "l");
    leg->AddEntry(h_eeg, "EEG", "l");
    leg->AddEntry(h_omegapi, "#omega#pi^{0}", "l");
    leg->AddEntry(h_ksl, "K_{S}K_{L}", "l");
    leg->AddEntry(h_signal, "Corrected #omega peak", "l");
    leg->AddEntry(h_nonReson_scaled, "Non-resonant template", "l");
    leg->AddEntry(h_background_cheb, "Chebyshev background (2nd order)", "l");
    leg->AddEntry(h_mcrest, "Others", "l");
    leg->Draw();

    // Add fit info to plot
    TPaveText *pt = new TPaveText(0.15, 0.78, 0.45, 0.92, "NDC");
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->SetTextSize(0.04);
    pt->AddText(Form("#chi^{2}/ndf = %.2f", chi2_ndf));
    pt->AddText(Form("#alpha = %.0f, #beta = %.0f", alpha, beta));
    pt->Draw();

    // Pull pad
    c1->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.12);
    pad2->Draw();
    pad2->cd();
    gPad->SetGrid();

    h_pull->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
    h_pull->GetXaxis()->SetTitleSize(0.15);
    h_pull->GetXaxis()->SetTitleOffset(1.0);
    h_pull->GetXaxis()->SetLabelSize(0.15);
    h_pull->GetYaxis()->SetTitle("Pull");
    h_pull->GetYaxis()->SetTitleSize(0.2);
    h_pull->GetYaxis()->SetTitleOffset(0.2);
    h_pull->GetYaxis()->SetLabelSize(0.15);
    h_pull->GetYaxis()->SetRangeUser(-5, 5);
    h_pull->GetXaxis()->SetNdivisions(505);
    h_pull->GetYaxis()->SetNdivisions(505);
    h_pull->GetXaxis()->CenterTitle();
    h_pull->GetYaxis()->CenterTitle();
    h_pull->Draw("P");

    TLine *line = new TLine(low, 0, high, 0);
    line->SetLineStyle(2);
    line->Draw();

    c1->SaveAs(output_path + "omega_nonreson_chebyshev2_fit.pdf");

    // ------------------------------------------------------------------
    // 13. Background-subtracted ω signal
    // ------------------------------------------------------------------
    TH1D *h_signal_data = (TH1D*) h_data->Clone("h_signal_data");
    h_signal_data->Add(h_eeg, -1.0);
    h_signal_data->Add(h_omegapi, -1.0);
    h_signal_data->Add(h_ksl, -1.0);
    h_signal_data->Add(h_mcrest, -1.0);
    h_signal_data->Add(h_nonReson_scaled, -1.0);
    h_signal_data->Add(h_background_cheb, -1.0);
    for (int bin = 1; bin <= h_signal_data->GetNbinsX(); ++bin)
        if (h_signal_data->GetBinContent(bin) < 0) h_signal_data->SetBinContent(bin, 0);

    TCanvas *c_bkg = new TCanvas("c_bkg", "Background-subtracted ω signal", 1200, 700);
    c_bkg->SetBottomMargin(0.13);
    c_bkg->SetLeftMargin(0.12);

    h_signal_data->SetMarkerStyle(20);
    h_signal_data->SetMarkerSize(0.6);
    h_signal_data->GetYaxis()->SetTitle(Form("Events / [%.1f MeV/c^{2}]", bin_width));
    h_signal_data->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
    h_signal_data->GetYaxis()->CenterTitle();
    h_signal_data->GetXaxis()->CenterTitle();
    h_signal_data->Draw("E0");
    h_signal->SetLineColor(kBlue);
    h_signal->Draw("hist same");

    TLegend *leg2 = new TLegend(0.55, 0.7, 0.9, 0.9);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.04);
    leg2->AddEntry(h_signal_data, "Data - backgrounds", "lep");
    leg2->AddEntry(h_signal, "Corrected #omega peak", "l");
    leg2->Draw();
    c_bkg->SaveAs(output_path + "omega_background_subtracted_nonreson_chebyshev2.pdf");

    // ------------------------------------------------------------------
    // 14. Correction weights canvas
    // ------------------------------------------------------------------
    TCanvas *c_weight = new TCanvas("c_weight", "Correction weights", 800, 600);
    h_weight_smooth->SetTitle("Correction weight for ISR3pi");
    h_weight_smooth->GetXaxis()->SetTitle("M_{3π} [MeV/c^{2}]");
    h_weight_smooth->GetYaxis()->SetTitle("Weight");
    h_weight_smooth->SetLineColor(kBlue);
    h_weight_smooth->SetLineWidth(2);
    h_weight_smooth->Draw();
    c_weight->SaveAs(output_path + "omega_correction_weights_nonreson_chebyshev2.pdf");

    // ------------------------------------------------------------------
    // 15. Save results
    // ------------------------------------------------------------------
    TFile *fout = new TFile(output_path + "omega_fit_nonreson_chebyshev2.root", "RECREATE");
    h_data_isr->Write();
    h_isr3pi->Write();
    h_nonReson->Write();
    h_signal->Write();
    h_nonReson_scaled->Write();
    h_background_cheb->Write();
    h_background_total->Write();
    h_weight_smooth->Write();
    h_isr3pi_corrected->Write();
    h_mc_total->Write();
    h_pull->Write();
    total_func->Write();
    fout->Close();

    // ------------------------------------------------------------------
    // 16. Clean up
    // ------------------------------------------------------------------
    gSigTemplate = nullptr;
    gBkgTemplate = nullptr;
    total_func->SetParent(0);

    if (total_func) delete total_func;
    if (bkg_cheb) delete bkg_cheb;
    if (h_side) delete h_side;
    if (c1) delete c1;
    if (c_bkg) delete c_bkg;
    if (c_weight) delete c_weight;
    if (h_pull) delete h_pull;
    if (h_signal_data) delete h_signal_data;
    if (h_mc_total) delete h_mc_total;
    if (pt) delete pt;
    if (leg) delete leg;
    if (leg2) delete leg2;
    if (line) delete line;
    if (pad1) delete pad1;
    if (pad2) delete pad2;

    if (ftree) {
        ftree->Close();
        delete ftree;
    }

    if (fout) {
        fout->Close();
        delete fout;
    }

    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Saved " << output_path << "omega_fit_nonreson_chebyshev2.root" << std::endl;
    std::cout << "Saved " << output_path << "omega_nonreson_chebyshev2_fit.pdf" << std::endl;
    std::cout << "Saved " << output_path << "omega_background_subtracted_nonreson_chebyshev2.pdf" << std::endl;
    std::cout << "Saved " << output_path << "omega_correction_weights_nonreson_chebyshev2.pdf" << std::endl;
    std::cout << "Fit: α = " << alpha << ", β = " << beta << std::endl;
    std::cout << "χ²/ndf = " << chi2_ndf << std::endl;
}
