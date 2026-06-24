// omega_fit_bdt_poly2_peak.C – template fit with 2nd order polynomial in signal region
// All masses in MeV
// Signal template + 2nd order polynomial (p0 + p1*x + p2*x²) fitted in peak region

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

// Global pointer for signal template (detached from file)
TH1D *gSigTemplate = nullptr;

// Fit function: α * signal_template + 2nd order polynomial (p0 + p1*x + p2*x²)
Double_t template_poly2_peak(Double_t *x, Double_t *par) {
    if (!gSigTemplate) {
        return 0.0;
    }
    
    int bin_sig = gSigTemplate->FindBin(x[0]);
    Double_t sig = gSigTemplate->GetBinContent(bin_sig);
    Double_t X = x[0];
    
    // 2nd order polynomial: p0 + p1*x + p2*x²
    Double_t poly = par[1] + par[2] * X + par[3] * X * X;
    return par[0] * sig + poly;
}

void omega_fit_bdt_poly2() {
    gErrorIgnoreLevel = kError;
    TGaxis::SetMaxDigits(4);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0.8);
    TH1::SetDefaultSumw2();

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

    // MC Rest
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
        // nonReson is NOT subtracted (not used in this fit)
        if (h_ksl) {
            bkg_total += h_ksl->GetBinContent(bin);
            bkg_err_sq += pow(h_ksl->GetBinError(bin), 2);
        }
        if (h_mcrest) {
            bkg_total += h_mcrest->GetBinContent(bin);
            bkg_err_sq += pow(h_mcrest->GetBinError(bin), 2);
        }
        
        double new_val = data_val - bkg_total;
        double new_err = sqrt(data_err*data_err + bkg_err_sq);
        
        h_data_isr->SetBinContent(bin, new_val);
        h_data_isr->SetBinError(bin, new_err);
        
        if (h_data_isr->GetBinContent(bin) < 0) {
            h_data_isr->SetBinContent(bin, 0);
        }
    }

    // ------------------------------------------------------------------
    // 5. Normalize signal template to unit area
    // ------------------------------------------------------------------
    double sig_int = h_isr3pi->Integral();
    if (sig_int > 0) h_isr3pi->Scale(1.0 / sig_int);
    
    // Set global pointer for the fit function
    gSigTemplate = h_isr3pi;

    // ------------------------------------------------------------------
    // 6. Determine 2nd order polynomial background from signal region
    // ------------------------------------------------------------------
    double peak_low = 740.0;
    double peak_high = 820.0;
    
    // Fit polynomial in the signal region (background underneath the peak)
    double poly_low = 740.0;
    double poly_high = 820.0;

    std::cout << "\nUsing 2nd order polynomial in signal region:\n"
              << "  Range: " << poly_low << "-" << poly_high << " MeV\n";

    // 2nd order polynomial: p0 + p1*x + p2*x²
    TF1 *bkg_poly2 = new TF1("bkg_poly2", "[0] + [1]*x + [2]*x*x", poly_low, poly_high);
    
    // Fit the polynomial in the signal region (background underneath the peak)
    h_data_isr->Fit(bkg_poly2, "RQN", "", poly_low, poly_high);

    double p0 = bkg_poly2->GetParameter(0);
    double p1 = bkg_poly2->GetParameter(1);
    double p2 = bkg_poly2->GetParameter(2);
    std::cout << "Polynomial fit: p0 = " << p0 << ", p1 = " << p1 << ", p2 = " << p2 << "\n";
    std::cout << "Polynomial χ²/ndf = " << bkg_poly2->GetChisquare() / bkg_poly2->GetNDF() << "\n";

    // ------------------------------------------------------------------
    // 7. Template fit with 2nd order polynomial background
    // ------------------------------------------------------------------
    double fit_low = 740.0;
    double fit_high = 820.0;
    
    std::cout << "\nFitting from " << fit_low << " to " << fit_high << " MeV\n";

    // Total fit function (4 parameters: alpha, p0, p1, p2)
    TF1 *total_func = new TF1("total_func", template_poly2_peak, low, high, 4);
    total_func->SetParameters(1000, p0, p1, p2);
    total_func->SetParNames("alpha", "p0", "p1", "p2");
    
    // Parameter limits
    total_func->SetParLimits(0, 0, 200000);  // alpha
    
    // Allow polynomial coefficients to vary within reasonable range
    double lim_p0 = 0.5 * std::abs(p0);
    if (lim_p0 < 1e-6) lim_p0 = 100.0;
    double lim_p1 = 0.5 * std::abs(p1);
    if (lim_p1 < 1e-6) lim_p1 = 1.0;
    double lim_p2 = 0.5 * std::abs(p2);
    if (lim_p2 < 1e-6) lim_p2 = 0.01;
    
    total_func->SetParLimits(1, p0 - lim_p0, p0 + lim_p0);
    total_func->SetParLimits(2, p1 - lim_p1, p1 + lim_p1);
    total_func->SetParLimits(3, p2 - lim_p2, p2 + lim_p2);
    
    // Single fit with quality assessment
    TFitResultPtr r = h_data_isr->Fit(total_func, "RQS", "", fit_low, fit_high);
    
    // Check if fit converged
    if (!r->IsValid()) {
        std::cerr << "WARNING: Fit did not converge!" << std::endl;
    } else {
        std::cout << "Fit converged successfully." << std::endl;
    }
    
    double alpha = total_func->GetParameter(0);
    double poly0 = total_func->GetParameter(1);
    double poly1 = total_func->GetParameter(2);
    double poly2 = total_func->GetParameter(3);
    double chi2 = r->Chi2();
    int ndf = r->Ndf();
    double chi2_ndf = chi2 / ndf;
    
    std::cout << "Fit results: α = " << alpha 
              << ", p0 = " << poly0 << ", p1 = " << poly1 << ", p2 = " << poly2 << std::endl;
    std::cout << "Fit quality: χ² = " << chi2 << ", ndf = " << ndf << ", χ²/ndf = " << chi2_ndf << std::endl;

    // ------------------------------------------------------------------
    // 8. Create signal & background histograms (scaled)
    // ------------------------------------------------------------------
    TH1D *h_signal = (TH1D*) h_isr3pi->Clone("h_signal");
    h_signal->SetDirectory(0);
    h_signal->Scale(alpha);
    h_signal->SetLineColor(kBlue);
    h_signal->SetLineWidth(2);

    TH1D *h_background_poly2 = (TH1D*) h_data_isr->Clone("h_background_poly2");
    h_background_poly2->SetDirectory(0);
    h_background_poly2->Reset();
    for (int bin = 1; bin <= h_background_poly2->GetNbinsX(); ++bin) {
        double x = h_background_poly2->GetBinCenter(bin);
        double val = poly0 + poly1 * x + poly2 * x * x;
        if (val < 0) val = 0;
        h_background_poly2->SetBinContent(bin, val);
        h_background_poly2->SetBinError(bin, TMath::Sqrt(val));
    }
    h_background_poly2->SetLineColor(kRed);
    h_background_poly2->SetLineStyle(3);
    h_background_poly2->SetLineWidth(2);

    // Other backgrounds (from MC)
    TH1D *h_background_mc = (TH1D*) h_mcrest->Clone("h_background_mc");
    h_background_mc->Add(h_eeg);
    h_background_mc->Add(h_omegapi);
    h_background_mc->Add(h_ksl);
    h_background_mc->SetLineColor(kOrange);
    h_background_mc->SetLineStyle(4);
    h_background_mc->SetLineWidth(2);

    // Total background (MC backgrounds + polynomial)
    TH1D *h_background_total = (TH1D*) h_background_mc->Clone("h_background_total");
    h_background_total->Add(h_background_poly2);
    h_background_total->SetLineColor(kMagenta);
    h_background_total->SetLineStyle(4);
    h_background_total->SetLineWidth(2);
    
    // ------------------------------------------------------------------
    // 9. Compute correction weights (using h_signal as desired shape)
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
    for (int bin = 2; bin <= h_weight_smooth->GetNbinsX()-1; ++bin) {
        double w_avg = (h_weight->GetBinContent(bin-1) +
                        h_weight->GetBinContent(bin) +
                        h_weight->GetBinContent(bin+1)) / 3.0;
        h_weight_smooth->SetBinContent(bin, w_avg);
    }

    // ------------------------------------------------------------------
    // 10. Apply correction to ISR3pi MC (use original scaled MC)
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
    double orig_int = h_isr3pi_orig->Integral(fit_low, fit_high);
    double new_int = h_isr3pi_corrected->Integral(fit_low, fit_high);
    if (new_int > 0 && orig_int > 0) {
        double renorm = orig_int / new_int;
        h_isr3pi_corrected->Scale(renorm);
        std::cout << "Renormalisation factor: " << renorm << std::endl;
    }
    h_isr3pi_corrected->SetLineColor(kGreen);
    h_isr3pi_corrected->SetLineWidth(2);

    // ------------------------------------------------------------------
    // 11. Build total MC sum for plotting
    // ------------------------------------------------------------------
    std::vector<TH1D*> comps;
    comps.push_back(h_eeg);
    comps.push_back(h_omegapi);
    comps.push_back(h_ksl);
    comps.push_back(h_mcrest);
    comps.push_back(h_signal);
    comps.push_back(h_background_poly2);

    TH1D *h_mc_total = (TH1D*) h_mcrest->Clone("h_mc_total");
    h_mc_total->Reset();
    h_mc_total->Sumw2();
    for (auto h : comps) if (h) h_mc_total->Add(h);
    h_mc_total->SetLineColor(kRed);
    h_mc_total->SetLineStyle(1);
    h_mc_total->SetLineWidth(2);

    // ------------------------------------------------------------------
    // 12. Pull distribution (full range)
    // ------------------------------------------------------------------
    TH1D *h_pull = new TH1D("h_pull", "", nbins, low, high);
    for (int bin = 1; bin <= nbins; ++bin) {
        double d = h_data->GetBinContent(bin);
        double m = h_mc_total->GetBinContent(bin);
        double err = std::sqrt(d + m);
        if (err > 0) h_pull->SetBinContent(bin, (d - m)/err);
        else h_pull->SetBinContent(bin, 0);
    }
    h_pull->SetMarkerStyle(20);
    h_pull->SetMarkerSize(0.6);
    h_pull->SetLineWidth(0);

    // ------------------------------------------------------------------
    // 13. Main plotting
    // ------------------------------------------------------------------
    TCanvas *c = new TCanvas("c", "3π mass projection (polynomial in peak)", 1200, 700);
    c->SetBottomMargin(0.12);
    c->SetLeftMargin(0.12);

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
    h_data->GetXaxis()->SetLabelOffset(0.1);
    h_data->GetYaxis()->SetTitleSize(0.07);
    h_data->GetYaxis()->SetTitleOffset(0.7);
    h_data->GetYaxis()->SetLabelSize(0.04);
    h_data->GetYaxis()->SetNdivisions(505);

    // Legend
    TLegend *leg = new TLegend(0.15, 0.30, 0.6, 0.9);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(h_data, "Data", "lep");
    leg->AddEntry(h_mc_total, "Total MC", "l");
    leg->AddEntry(h_eeg, "EEG", "l");
    leg->AddEntry(h_omegapi, "#omega#pi^{0}", "l");
    leg->AddEntry(h_ksl, "K_{S}K_{L}", "l");
    leg->AddEntry(h_mcrest, "Others", "l");
    leg->AddEntry(h_signal, "Corrected #omega peak", "l");
    leg->AddEntry(h_background_poly2, "Polynomial background (peak)", "l");
    leg->Draw();

    // Add fit info to plot
    TPaveText *pt = new TPaveText(0.15, 0.78, 0.45, 0.92, "NDC");
    pt->SetFillStyle(0);
    pt->SetBorderSize(0);
    pt->SetTextSize(0.04);
    pt->AddText(Form("#chi^{2}/ndf = %.2f", chi2_ndf));
    pt->AddText(Form("#alpha = %.0f", alpha));
    pt->Draw();

    // Pull pad
    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.3);
    pad2->SetLeftMargin(0.12);
    pad2->Draw();
    pad2->cd();
    gPad->SetGrid();

    h_pull->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
    h_pull->GetXaxis()->SetTitleSize(0.12);
    h_pull->GetXaxis()->SetTitleOffset(1.0);
    h_pull->GetXaxis()->SetLabelSize(0.1);
    h_pull->GetYaxis()->SetTitle("Pull");
    h_pull->GetYaxis()->SetTitleSize(0.2);
    h_pull->GetYaxis()->SetTitleOffset(0.2);
    h_pull->GetYaxis()->SetLabelSize(0.1);
    h_pull->GetYaxis()->SetRangeUser(-5, 5);
    h_pull->GetYaxis()->SetNdivisions(505);
    h_pull->GetXaxis()->CenterTitle();
    h_pull->GetYaxis()->CenterTitle();
    h_pull->Draw("P");

    TLine *line = new TLine(low, 0, high, 0);
    line->SetLineStyle(2);
    line->Draw();

    c->SaveAs(output_path + "omega_poly2_peak_fit.pdf");

    // ------------------------------------------------------------------
    // 14. Background-subtracted ω signal
    // ------------------------------------------------------------------
    TH1D *h_signal_data = (TH1D*) h_data->Clone("h_signal_data");
    h_signal_data->Add(h_eeg, -1.0);
    h_signal_data->Add(h_omegapi, -1.0);
    h_signal_data->Add(h_ksl, -1.0);
    h_signal_data->Add(h_mcrest, -1.0);
    h_signal_data->Add(h_background_poly2, -1.0);
    for (int bin = 1; bin <= h_signal_data->GetNbinsX(); ++bin)
        if (h_signal_data->GetBinContent(bin) < 0) h_signal_data->SetBinContent(bin, 0);

    TCanvas *c2 = new TCanvas("c2", "Background-subtracted ω signal", 1200, 700);
    c2->SetBottomMargin(0.12);
    c2->SetLeftMargin(0.12);
    h_signal_data->SetMarkerStyle(20);
    h_signal_data->SetMarkerSize(0.6);
    h_signal_data->GetYaxis()->SetTitle(Form("Events / [%.1f MeV/c^{2}]", bin_width));
    h_signal_data->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
    h_signal_data->GetYaxis()->CenterTitle();
    h_signal_data->GetXaxis()->CenterTitle();
    h_signal_data->Draw("E0");
    h_signal->SetLineColor(kBlue);
    h_signal->Draw("hist same");

    TLegend *leg2 = new TLegend(0.6, 0.7, 0.9, 0.9);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.04);
    leg2->AddEntry(h_signal_data, "Data - backgrounds", "lep");
    leg2->AddEntry(h_signal, "Corrected #omega peak", "l");
    leg2->Draw();
    c2->SaveAs(output_path + "omega_background_subtracted_poly2_peak.pdf");

    // ------------------------------------------------------------------
    // 15. Correction weights canvas
    // ------------------------------------------------------------------
    TCanvas *c_weight = new TCanvas("c_weight", "Correction weights", 800, 600);
    h_weight_smooth->SetTitle("Correction weight for ISR3pi");
    h_weight_smooth->GetXaxis()->SetTitle("M_{3π} [MeV/c^{2}]");
    h_weight_smooth->GetYaxis()->SetTitle("Weight");
    h_weight_smooth->SetLineColor(kBlue);
    h_weight_smooth->SetLineWidth(2);
    h_weight_smooth->Draw();
    c_weight->SaveAs(output_path + "omega_correction_weights_poly2_peak.pdf");

    // ------------------------------------------------------------------
    // 16. Save results
    // ------------------------------------------------------------------
    TFile *fout = new TFile(output_path + "omega_fit_poly2_peak.root", "RECREATE");
    h_data_isr->Write();
    h_isr3pi->Write();
    h_signal->Write();
    h_background_poly2->Write();
    h_background_mc->Write();
    h_background_total->Write();
    h_weight_smooth->Write();
    h_isr3pi_corrected->Write();
    h_mc_total->Write();
    h_pull->Write();
    total_func->Write();
    fout->Close();

    // ------------------------------------------------------------------
    // 17. Clean up
    // ------------------------------------------------------------------
    delete total_func;
    delete bkg_poly2;

    gSigTemplate = nullptr;
    
    delete h_data_isr;
    delete h_background_mc;
    delete h_background_total;
    delete h_isr3pi_corrected;
    delete h_isr3pi_orig;
    delete h_weight;
    
    delete c;
    delete c2;
    delete c_weight;
    
    if (ftree) {
        ftree->Close();
        delete ftree;
    }
    
    if (fout) {
        fout->Close();
        delete fout;
    }
    
    std::cout << "\n=== Summary ===" << std::endl;
    std::cout << "Saved " << output_path << "omega_fit_poly2_peak.root" << std::endl;
    std::cout << "Saved " << output_path << "omega_poly2_peak_fit.pdf" << std::endl;
    std::cout << "Saved " << output_path << "omega_background_subtracted_poly2_peak.pdf" << std::endl;
    std::cout << "Saved " << output_path << "omega_correction_weights_poly2_peak.pdf" << std::endl;
    std::cout << "Fit: α = " << alpha << std::endl;
    std::cout << "χ²/ndf = " << chi2_ndf << std::endl;
}
