// correct_and_plot_single_sideband_poly.C
// Single lower‑sideband ordinary polynomial fit (2nd order) with full plotting

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLine.h>
#include <TFitResult.h>
#include <iostream>
#include <cmath>
#include <vector>

#include "../header_bdt/sfw2d_bdt.txt"
#include "../header_bdt/correct_omega.h"
#include "../header_bdt/path.h"

TH1D *gSigTemplate = nullptr;

// Fit function: α * signal_template + ordinary polynomial (quadratic)
Double_t template_poly2(Double_t *x, Double_t *par) {
    int bin = gSigTemplate->FindBin(x[0]);
    Double_t sig = gSigTemplate->GetBinContent(bin);
    Double_t X = x[0];
    // Ordinary polynomial: p0 + p1*x + p2*x²
    Double_t poly = par[1] + par[2]*X + par[3]*X*X;
    return par[0] * sig + poly;
}

void correct_and_plot_single_sideband_poly() {
    gErrorIgnoreLevel = kError;
    TGaxis::SetMaxDigits(4);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0.8);
    TH1::SetDefaultSumw2();

    // ------------------------------------------------------------------
    // 1. Open tree file
    // ------------------------------------------------------------------
    TString treeFile = "/home/kloe/Desktop/input_bdt_TDATA_norm/cut/tree_pre_bdt.root";
    TFile *ftree = TFile::Open(treeFile);
    if (!ftree || ftree->IsZombie()) {
        std::cerr << "ERROR: cannot open " << treeFile << std::endl;
        return;
    }

    TTree *tdata = (TTree*) ftree->Get("TDATA");
    if (!tdata) { std::cerr << "ERROR: TDATA not found." << std::endl; return; }

    const double mass_min = 600.0;
    const double mass_max = 1000.0;
    const int nbins = 200;
    std::cout << "Using range [" << mass_min << ", " << mass_max << "] with " << nbins << " bins.\n";

    // ------------------------------------------------------------------
    // 2. Data histogram
    // ------------------------------------------------------------------
    TH1D *h_data = new TH1D("h_data", "", nbins, mass_min, mass_max);
    h_data->Sumw2();
    h_data->SetDirectory(0);
    double m3pi;
    tdata->SetBranchAddress("Br_m3pi_bdt", &m3pi);
    for (Long64_t i = 0; i < tdata->GetEntries(); ++i) {
        tdata->GetEntry(i);
        h_data->Fill(m3pi);
    }
    std::cout << "Data integral: " << h_data->Integral() << std::endl;

    // ------------------------------------------------------------------
    // 3. Load scaled MC background components (with color and style)
    // ------------------------------------------------------------------
    auto makeScaledHist = [&](const char* tname, double scale, int color, int style) -> TH1D* {
        TTree *t = (TTree*) ftree->Get(tname);
        if (!t) return nullptr;
        if (!t->GetBranch("Br_m3pi_bdt")) return nullptr;
        TH1D *h = new TH1D(Form("h_%s", tname), "", nbins, mass_min, mass_max);
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

    TH1D *h_eeg       = makeScaledHist("TEEG",        eeg_sfw,      6, 7);
    TH1D *h_omegapi   = makeScaledHist("TOMEGAPI",    omegapi_sfw,  7, 5);
    TH1D *h_ksl       = makeScaledHist("TKSL",        ksl_sfw,     28, 4);
    TH1D *h_etagam    = makeScaledHist("TETAGAM",     etagam_sfw,   3, 3);
    TH1D *h_isr3pi    = makeScaledHist("TISR3PI_SIG", isr3pi_sfw,   4, 2);

    // MC Rest
    TH1D *h_mcrest = new TH1D("h_mcrest", "", nbins, mass_min, mass_max);
    h_mcrest->Sumw2();
    h_mcrest->SetDirectory(0);
    for (const char* name : {"TKPM", "TRHOPI", "TBKGREST"}) {
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

    // ------------------------------------------------------------------
    // 4. Subtract MC backgrounds -> h_data_isr (for fit)
    // ------------------------------------------------------------------
    TH1D *h_data_isr = (TH1D*) h_data->Clone("h_data_isr");
    h_data_isr->SetDirectory(0);
    if (h_eeg)     h_data_isr->Add(h_eeg, -1.0);
    if (h_omegapi) h_data_isr->Add(h_omegapi, -1.0);
    if (h_ksl)     h_data_isr->Add(h_ksl, -1.0);
    if (h_etagam)  h_data_isr->Add(h_etagam, -1.0);
    if (h_mcrest)  h_data_isr->Add(h_mcrest, -1.0);
    for (int bin = 1; bin <= h_data_isr->GetNbinsX(); ++bin)
        if (h_data_isr->GetBinContent(bin) < 0) h_data_isr->SetBinContent(bin, 0);

    // ------------------------------------------------------------------
    // 5. Create signal template (correctly paired ω)
    // ------------------------------------------------------------------
    TTree *tmc = (TTree*) ftree->Get("TISR3PI_SIG");
    if (!tmc) { std::cerr << "ERROR: TISR3PI_SIG not found." << std::endl; return; }
    int recon_indx_bdt = 0, bkg_indx = 0;
    double m3pi_tmp;
    tmc->SetBranchAddress("Br_recon_indx_bdt", &recon_indx_bdt);
    tmc->SetBranchAddress("Br_bkg_indx", &bkg_indx);
    tmc->SetBranchAddress("Br_m3pi_bdt", &m3pi_tmp);

    TH1D *h_signal_template = new TH1D("h_signal_template", "", nbins, mass_min, mass_max);
    h_signal_template->Sumw2();
    h_signal_template->SetDirectory(0);
    for (Long64_t i = 0; i < tmc->GetEntries(); ++i) {
        tmc->GetEntry(i);
        if (recon_indx_bdt == 2 && bkg_indx == 1)
            h_signal_template->Fill(m3pi_tmp);
    }
    double sig_int = h_signal_template->Integral();
    if (sig_int > 0) h_signal_template->Scale(1.0 / sig_int);
    gSigTemplate = h_signal_template;

    // ------------------------------------------------------------------
    // 6. Single lower‑sideband ordinary polynomial fit (quadratic)
    // ------------------------------------------------------------------
    double peak_low = 740.0;
    double peak_high = 820.0;
    double sb_low  = 670.0;
    double sb_high = 730.0;

    // Clone and zero peak region (so it does not affect sideband fit)
    TH1D *h_side = (TH1D*) h_data_isr->Clone("h_side");
    for (int bin = 1; bin <= h_side->GetNbinsX(); ++bin) {
        double x = h_side->GetBinCenter(bin);
        if (x >= peak_low && x <= peak_high) h_side->SetBinContent(bin, 0);
    }

    // Ordinary quadratic polynomial: p0 + p1*x + p2*x²
    TF1 *bkg_poly = new TF1("bkg_poly", 
        "[0] + [1]*x + [2]*x*x", 
        sb_low, sb_high);

    // Fit only the lower sideband
    h_side->Fit(bkg_poly, "QN", "", sb_low, sb_high);

    double p0 = bkg_poly->GetParameter(0);
    double p1 = bkg_poly->GetParameter(1);
    double p2 = bkg_poly->GetParameter(2);
    std::cout << "Single lower‑sideband ordinary polynomial fit (quadratic):\n"
              << "  Sideband: " << sb_low << "–" << sb_high << " MeV\n"
              << "  Polynomial: p0 = " << p0 << ", p1 = " << p1 
              << ", p2 = " << p2 << "\n";

    // Total fit function (4 parameters: alpha, p0, p1, p2)
    TF1 *total_func = new TF1("total_func", template_poly2, mass_min, mass_max, 4);
    total_func->SetParameters(1000, p0, p1, p2);
    total_func->SetParNames("alpha", "p0", "p1", "p2");

    // Parameter limits with fallback
    auto setLim = [](TF1 *f, int ipar, double val) {
        double lim = 0.3 * std::abs(val);
        if (lim < 1e-6) lim = 1.0;
        f->SetParLimits(ipar, val - lim, val + lim);
    };
    setLim(total_func, 1, p0);
    setLim(total_func, 2, p1);
    setLim(total_func, 3, p2);

    TFitResultPtr r = h_data_isr->Fit(total_func, "RQS", "", peak_low, peak_high);
    double alpha = total_func->GetParameter(0);
    double poly0 = total_func->GetParameter(1);
    double poly1 = total_func->GetParameter(2);
    double poly2 = total_func->GetParameter(3);

    double chi2 = r->Chi2();
    int ndf = r->Ndf();
    double chi2_ndf = chi2 / ndf;
    std::cout << "Fit results: α = " << alpha
              << ", Polynomial: p0 = " << poly0 << ", p1 = " << poly1 
              << ", p2 = " << poly2 << "\n";
    std::cout << "Fit quality: χ² = " << chi2 << ", ndf = " << ndf 
              << ", χ²/ndf = " << chi2_ndf << std::endl;

    // ------------------------------------------------------------------
    // 7. Create corrected signal and background histograms
    // ------------------------------------------------------------------
    TH1D *h_signal = (TH1D*) h_signal_template->Clone("h_signal");
    h_signal->SetDirectory(0);
    h_signal->Scale(alpha);
    h_signal->SetLineColor(kBlue);
    h_signal->SetLineWidth(2);

    TH1D *h_background = (TH1D*) h_data_isr->Clone("h_background");
    h_background->SetDirectory(0);
    h_background->Reset();
    for (int bin = 1; bin <= h_background->GetNbinsX(); ++bin) {
        double x = h_background->GetBinCenter(bin);
        double val = poly0 + poly1*x + poly2*x*x;
        if (val < 0) val = 0;
        h_background->SetBinContent(bin, val);
        h_background->SetBinError(bin, TMath::Sqrt(val));
    }
    h_background->SetLineColor(kRed);
    h_background->SetLineStyle(3);
    h_background->SetLineWidth(2);

    // ------------------------------------------------------------------
    // 8. Compute correction weights and apply to ISR3pi MC
    // ------------------------------------------------------------------
    TH1D *h_weight = (TH1D*) h_data_isr->Clone("h_weight");
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
    TH1D *h_weight_smooth = (TH1D*) h_weight->Clone("h_weight_smooth");
    for (int bin = 2; bin <= h_weight_smooth->GetNbinsX()-1; ++bin) {
        double w_avg = (h_weight->GetBinContent(bin-1) +
                        h_weight->GetBinContent(bin) +
                        h_weight->GetBinContent(bin+1)) / 3.0;
        h_weight_smooth->SetBinContent(bin, w_avg);
    }

    TH1D *h_isr3pi_corrected = (TH1D*) h_isr3pi->Clone("h_isr3pi_corrected");
    for (int bin = 1; bin <= h_isr3pi_corrected->GetNbinsX(); ++bin) {
        double w = h_weight_smooth->GetBinContent(bin);
        double old = h_isr3pi_corrected->GetBinContent(bin);
        h_isr3pi_corrected->SetBinContent(bin, old * w);
        h_isr3pi_corrected->SetBinError(bin, h_isr3pi_corrected->GetBinError(bin) * w);
    }
    double orig_peak = h_isr3pi->Integral(peak_low, peak_high);
    double new_peak = h_isr3pi_corrected->Integral(peak_low, peak_high);
    if (new_peak > 0) h_isr3pi_corrected->Scale(orig_peak / new_peak);
    h_isr3pi_corrected->SetLineColor(kGreen);
    h_isr3pi_corrected->SetLineStyle(1);
    h_isr3pi_corrected->SetLineWidth(2);

    // ------------------------------------------------------------------
    // 9. Build total MC sum for plotting
    // ------------------------------------------------------------------
    std::vector<TH1D*> comps;
    comps.push_back(h_eeg);
    comps.push_back(h_omegapi);
    comps.push_back(h_ksl);
    comps.push_back(h_etagam);
    comps.push_back(h_mcrest);
    comps.push_back(h_signal);
    comps.push_back(h_background);

    TH1D *h_mc_total = (TH1D*) h_mcrest->Clone("h_mc_total");
    h_mc_total->Reset();
    h_mc_total->Sumw2();
    for (auto h : comps) if (h) h_mc_total->Add(h);
    h_mc_total->SetLineColor(kRed);
    h_mc_total->SetLineWidth(2);
    h_mc_total->SetLineStyle(1);

    // ------------------------------------------------------------------
    // 10. Pull distribution
    // ------------------------------------------------------------------
    TH1D *h_pull = new TH1D("h_pull", "", nbins, mass_min, mass_max);
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
    // 11. Plotting
    // ------------------------------------------------------------------
    TCanvas *c = new TCanvas("c", "3π mass projection (ordinary polynomial, single sideband)", 1200, 700);
    c->SetBottomMargin(0.13);
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

    h_data->Draw("E0");
    h_mc_total->Draw("hist same");
    for (auto h : comps) if (h) h->Draw("hist same");
    h_data->GetYaxis()->SetTitle(Form("Events / [%.1f MeV/c^{2}]", bin_width));
    h_data->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
    h_data->GetXaxis()->SetRangeUser(mass_min, mass_max);
    h_data->GetYaxis()->CenterTitle();
    h_data->GetXaxis()->SetTitleSize(0.06);
    h_data->GetYaxis()->SetTitleSize(0.07);
    h_data->GetYaxis()->SetTitleOffset(0.7);
    h_data->GetYaxis()->SetLabelSize(0.06);
    h_data->GetYaxis()->SetNdivisions(505);
    h_data->GetXaxis()->SetLabelOffset(0.2);
    
    // Legend
    TLegend *leg = new TLegend(0.55, 0.25, 0.9, 0.9);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(h_data, "Data", "lep");
    leg->AddEntry(h_mc_total, "Total MC", "l");
    leg->AddEntry(h_eeg, "EEG", "l");
    leg->AddEntry(h_omegapi, "#omega#pi^{0}", "l");
    leg->AddEntry(h_ksl, "K_{S}K_{L}", "l");
    leg->AddEntry(h_etagam, "#eta#gamma", "l");
    leg->AddEntry(h_signal, "Corrected #omega peak", "l");
    leg->AddEntry(h_background, "Non-resonant background (poly)", "l");
    leg->AddEntry(h_mcrest, "Others", "l");
    leg->Draw();

    // Pull pad
    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.35);
    pad2->SetLeftMargin(0.12);
    pad2->Draw();
    pad2->cd();
    gPad->SetGrid();

    h_pull->GetXaxis()->SetTitle("M_{3#pi} [MeV]");
    h_pull->GetXaxis()->SetTitleSize(0.15);
    h_pull->GetXaxis()->SetTitleOffset(1.);
    h_pull->GetXaxis()->SetLabelSize(0.15);
    h_pull->GetYaxis()->SetTitle("Pull");
    h_pull->GetYaxis()->SetTitleSize(0.2);
    h_pull->GetYaxis()->SetTitleOffset(0.2);
    h_pull->GetYaxis()->SetLabelSize(0.15);
    h_pull->GetYaxis()->SetRangeUser(-50, 50);
    h_pull->GetXaxis()->SetNdivisions(505);
    h_pull->GetYaxis()->SetNdivisions(505);
    h_pull->GetXaxis()->CenterTitle();
    h_pull->GetYaxis()->CenterTitle();
    h_pull->Draw("P");

    TLine *line = new TLine(mass_min, 0, mass_max, 0);
    line->SetLineStyle(2);
    line->Draw();

    c->SaveAs(output_path + "combined_fit_and_plot_single_sideband_poly.pdf");
    delete c;

    // ------------------------------------------------------------------
    // 12. Background‑subtracted ω signal (optional)
    // ------------------------------------------------------------------
    TH1D *h_signal_data = (TH1D*) h_data->Clone("h_signal_data");
    h_signal_data->Add(h_eeg, -1.0);
    h_signal_data->Add(h_omegapi, -1.0);
    h_signal_data->Add(h_ksl, -1.0);
    h_signal_data->Add(h_etagam, -1.0);
    h_signal_data->Add(h_mcrest, -1.0);
    h_signal_data->Add(h_background, -1.0);
    for (int bin = 1; bin <= h_signal_data->GetNbinsX(); ++bin)
        if (h_signal_data->GetBinContent(bin) < 0) h_signal_data->SetBinContent(bin, 0);

    TCanvas *c22 = new TCanvas("c22", "Background‑subtracted #omega signal", 1200, 700);
    c22->SetBottomMargin(0.13);
    c22->SetLeftMargin(0.12);

    h_signal_data->SetMarkerStyle(20);
    h_signal_data->SetMarkerSize(0.6);
    h_signal_data->GetXaxis()->SetNdivisions(505);
    h_signal_data->GetYaxis()->SetNdivisions(505);
    h_signal_data->GetYaxis()->SetTitle(Form("Events / [%.1f MeV/c^{2}]", bin_width));
    h_signal_data->GetXaxis()->SetTitle("M_{3#pi} [MeV/c^{2}]");
    h_signal_data->GetXaxis()->SetTitleSize(0.05);
    h_signal_data->GetXaxis()->SetTitleOffset(1.2);
    h_signal_data->GetXaxis()->SetLabelSize(0.06);
    h_signal_data->GetXaxis()->SetLabelOffset(0.01);
    h_signal_data->GetXaxis()->SetTickLength(0.03);
    h_signal_data->GetYaxis()->SetTitleOffset(0.9);
    h_signal_data->GetXaxis()->SetRangeUser(mass_min, mass_max);
    h_signal_data->GetXaxis()->CenterTitle();
    h_signal_data->GetYaxis()->CenterTitle();
    h_signal->SetLineColor(kBlue);
    h_background->SetLineColor(kRed);
    
    h_signal_data->Draw("E0");
    h_signal->Draw("hist same");
    h_background->Draw("hist same");
    
    TLegend *leg2 = new TLegend(0.55, 0.7, 0.9, 0.9);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.04);
    leg2->AddEntry(h_signal_data, "Data - backgrounds", "lep");
    leg2->AddEntry(h_signal, "Corrected #omega peak", "l");
    leg2->AddEntry(h_background, "Non-resonant background (poly)", "l");
    leg2->Draw();

    c22->Update();
    c22->SaveAs(output_path + "background_subtracted_single_sideband_poly.pdf");
    delete c22;
    
    // ------------------------------------------------------------------
    // 13. Save output ROOT file
    // ------------------------------------------------------------------
    TFile *fout = new TFile(output_path + "corrected_and_plotted_single_sideband_poly.root", "RECREATE");
    h_isr3pi_corrected->Write();
    h_signal->Write();
    h_background->Write();
    h_weight_smooth->Write();
    h_data->Write();
    h_mc_total->Write();
    fout->Close();

    // ------------------------------------------------------------------
    // 14. Clean up
    // ------------------------------------------------------------------
    total_func->SetParent(0);
    delete total_func;
    delete bkg_poly;
    delete h_side;
    ftree->Close();
    delete ftree;
    delete fout;
    gSigTemplate = nullptr;

    gSystem->Exit(0);

    std::cout << "\nAll done. Outputs saved to " << output_path << std::endl;
}

int main() {
    correct_and_plot_single_sideband_poly();
    return 0;
}
