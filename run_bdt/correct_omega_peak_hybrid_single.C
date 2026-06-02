// correct_omega_peak_hybrid_lower_linear.C
// Data-driven correction for ω → 3π shape using:
//   - signal template from TISR3PI_SIG (recon_indx_bdt==2 && bkg_indx==1)
//   - 1st‑order (linear) polynomial background, using ONLY the lower sideband (680-730 MeV)
// Outputs: corrected_isr3pi_hybrid_lower_linear.root

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <cmath>

#include "../header_bdt/sfw2d.txt"
#include "../header_bdt/correct_omega.h"
#include "../header_bdt/binning.h"

// Global pointer to signal template (detached from file)
TH1D *gSigTemplate = nullptr;

// Fit function: α * signal_template + linear polynomial
Double_t template_linear(Double_t *x, Double_t *par) {
    int bin = gSigTemplate->FindBin(x[0]);
    Double_t sig = gSigTemplate->GetBinContent(bin);
    Double_t poly = par[1] + par[2] * x[0];   // linear
    return par[0] * sig + poly;
}

void correct_omega_peak_hybrid_single() {
    // ------------------------------------------------------------------
    // 1. Open tree file
    // ------------------------------------------------------------------
    TString treeFile = "/home/kloe/Desktop/input_bdt_TDATA_chain/cut/tree_pre_bdt.root";
    TFile *ftree = TFile::Open(treeFile);
    if (!ftree || ftree->IsZombie()) {
        std::cerr << "ERROR: cannot open " << treeFile << std::endl;
        return;
    }

    TTree *tdata = (TTree*) ftree->Get("TDATA");
    if (!tdata) { std::cerr << "ERROR: TDATA not found." << std::endl; return; }

    // Determine mass unit (MeV or GeV)
    double mtest;
    tdata->SetBranchAddress("Br_m3pi_bdt", &mtest);
    tdata->GetEntry(0);
    bool is_mev = (mtest > 10);
    double low, high;
    int nbins = NBINS;
    
    if (is_mev) {
        low  = MASS_MIN;          // e.g., 600.0 MeV
        high = MASS_MAX;          // e.g., 1000.0 MeV
    } else {
        low  = MASS_MIN / 1000.0; // convert to GeV
        high = MASS_MAX / 1000.0;
    }
    std::cout << "Mass unit: " << (is_mev ? "MeV" : "GeV")
              << " range [" << low << ", " << high << "]\n";

    // ------------------------------------------------------------------
    // 2. Data histogram (detached)
    // ------------------------------------------------------------------
    TH1D *h_data = new TH1D("h_data", "", nbins, low, high);
    h_data->Sumw2();
    h_data->SetDirectory(0);
    tdata->SetBranchAddress("Br_m3pi_bdt", &mtest);
    for (Long64_t i = 0; i < tdata->GetEntries(); ++i) {
        tdata->GetEntry(i);
        h_data->Fill(mtest);
    }
    std::cout << "Data integral: " << h_data->Integral() << std::endl;

    // ------------------------------------------------------------------
    // 3. Load scaled MC components (detach each)
    // ------------------------------------------------------------------
    auto makeScaledHist = [&](const char* tname, double scale) -> TH1D* {
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
        return h;
    };

    TH1D *h_eeg       = makeScaledHist("TEEG", eeg_sfw);
    TH1D *h_omegapi   = makeScaledHist("TOMEGAPI", omegapi_sfw);
    TH1D *h_ksl       = makeScaledHist("TKSL", ksl_sfw);
    TH1D *h_etagam    = makeScaledHist("TETAGAM", etagam_sfw);
    TH1D *h_isr3pi    = makeScaledHist("TISR3PI_SIG", isr3pi_sfw);

    // MC Rest
    TH1D *h_mcrest = new TH1D("h_mcrest", "", nbins, low, high);
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

    if (!h_isr3pi) { std::cerr << "ERROR: No ISR3pi histogram." << std::endl; return; }

    // ------------------------------------------------------------------
    // 4. Subtract backgrounds -> h_data_isr (detached)
    // ------------------------------------------------------------------
    TH1D *h_data_isr = (TH1D*) h_data->Clone("h_data_isr");
    h_data_isr->SetDirectory(0);
    if (h_eeg) h_data_isr->Add(h_eeg, -1.0);
    if (h_omegapi) h_data_isr->Add(h_omegapi, -1.0);
    if (h_ksl) h_data_isr->Add(h_ksl, -1.0);
    if (h_etagam) h_data_isr->Add(h_etagam, -1.0);
    if (h_mcrest) h_data_isr->Add(h_mcrest, -1.0);
    for (int bin = 1; bin <= h_data_isr->GetNbinsX(); ++bin)
        if (h_data_isr->GetBinContent(bin) < 0) h_data_isr->SetBinContent(bin, 0);

    // ------------------------------------------------------------------
    // 5. Create signal template from TISR3PI_SIG classification
    // ------------------------------------------------------------------
    TTree *tmc = (TTree*) ftree->Get("TISR3PI_SIG");
    if (!tmc) { std::cerr << "ERROR: TISR3PI_SIG not found." << std::endl; return; }
    int recon_indx_bdt = 0, bkg_indx = 0;
    double m3pi = 0.;
    tmc->SetBranchAddress("Br_recon_indx_bdt", &recon_indx_bdt);
    tmc->SetBranchAddress("Br_bkg_indx", &bkg_indx);
    tmc->SetBranchAddress("Br_m3pi_bdt", &m3pi);

    TH1D *h_signal_template = new TH1D("h_signal_template", "", nbins, low, high);
    h_signal_template->Sumw2();
    h_signal_template->SetDirectory(0);
    for (Long64_t i = 0; i < tmc->GetEntries(); ++i) {
        tmc->GetEntry(i);
        if (recon_indx_bdt == 2 && bkg_indx == 1)
            h_signal_template->Fill(m3pi);
    }
    // Normalise to unit area
    double sig_int = h_signal_template->Integral();
    if (sig_int > 0) h_signal_template->Scale(1.0 / sig_int);

    gSigTemplate = h_signal_template;

    // ------------------------------------------------------------------
    // 6. Template + linear polynomial fit – using ONLY lower sideband
    // ------------------------------------------------------------------
    double peak_low = is_mev ? 740 : 0.74;
    double peak_high = is_mev ? 820 : 0.82;

    // Lower sideband only (avoid bump above 800 MeV)
    double sb_low = is_mev ? 680 : 0.68;
    double sb_high = is_mev ? 730 : 0.73;

    // Create a temporary histogram for sideband fit (only background)
    TH1D *h_side = (TH1D*) h_data_isr->Clone("h_side");
    // Blank the peak region
    for (int bin = 1; bin <= h_side->GetNbinsX(); ++bin) {
        double x = h_side->GetBinCenter(bin);
        if (x >= peak_low && x <= peak_high) h_side->SetBinContent(bin, 0);
    }

    // Fit linear polynomial to lower sideband only
    TF1 *bkg_linear = new TF1("bkg_linear", "pol1", sb_low, sb_high);
    h_side->Fit(bkg_linear, "QN", "", sb_low, sb_high);
    double p0 = bkg_linear->GetParameter(0);
    double p1 = bkg_linear->GetParameter(1);
    std::cout << "Lower sideband linear fit (680-730 MeV): "
              << p0 << " + " << p1 << "·x\n";

    // Total model: α * signal_template + linear polynomial (3 parameters)
    TF1 *total_func = new TF1("total_func", template_linear, low, high, 3);
    total_func->SetParameters(1000, p0, p1);
    total_func->SetParNames("alpha", "p0", "p1");
    // Allow linear parameters to vary within ±30% of sideband values
    total_func->SetParLimits(1, p0 - 0.3*std::abs(p0), p0 + 0.3*std::abs(p0));
    total_func->SetParLimits(2, p1 - 0.3*std::abs(p1), p1 + 0.3*std::abs(p1));

    h_data_isr->Fit(total_func, "RQ", "", peak_low, peak_high);
    double alpha = total_func->GetParameter(0);
    double poly0 = total_func->GetParameter(1);
    double poly1 = total_func->GetParameter(2);

    std::cout << "Fit results: α = " << alpha << ", linear: "
              << poly0 << " + " << poly1 << "·x\n";

    // ------------------------------------------------------------------
    // 7. Create signal & background histograms
    // ------------------------------------------------------------------
    TH1D *h_signal = (TH1D*) h_signal_template->Clone("h_signal");
    h_signal->SetDirectory(0);
    h_signal->Scale(alpha);

    TH1D *h_background = (TH1D*) h_data_isr->Clone("h_background");
    h_background->SetDirectory(0);
    h_background->Reset();
    for (int bin = 1; bin <= h_background->GetNbinsX(); ++bin) {
        double x = h_background->GetBinCenter(bin);
        double poly_val = poly0 + poly1 * x;
        if (poly_val < 0) poly_val = 0;
        h_background->SetBinContent(bin, poly_val);
        h_background->SetBinError(bin, TMath::Sqrt(poly_val));
    }

    // ------------------------------------------------------------------
    // 8. Compute correction weights (desired shape = h_signal)
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
    // 9. Apply correction to ISR3pi MC
    // ------------------------------------------------------------------
    TH1D *h_isr3pi_corrected = (TH1D*) h_isr3pi->Clone("h_isr3pi_corrected");
    h_isr3pi_corrected->SetDirectory(0);
    for (int bin = 1; bin <= h_isr3pi_corrected->GetNbinsX(); ++bin) {
        double w = h_weight_smooth->GetBinContent(bin);
        double old = h_isr3pi_corrected->GetBinContent(bin);
        h_isr3pi_corrected->SetBinContent(bin, old * w);
        double err = h_isr3pi_corrected->GetBinError(bin);
        h_isr3pi_corrected->SetBinError(bin, err * w);
    }
    // Renormalise to keep integral in peak region unchanged
    double orig_int = h_isr3pi->Integral(peak_low, peak_high);
    double new_int = h_isr3pi_corrected->Integral(peak_low, peak_high);
    if (new_int > 0 && orig_int > 0) {
        double renorm = orig_int / new_int;
        h_isr3pi_corrected->Scale(renorm);
        std::cout << "Renormalisation factor: " << renorm << std::endl;
    }

    // ------------------------------------------------------------------
    // 10. Save results (close input file first)
    // ------------------------------------------------------------------
    ftree->Close();
    delete ftree;

    TFile *fout = new TFile(output_path + "corrected_isr3pi_hybrid_lower_linear.root", "RECREATE");
    h_isr3pi_corrected->Write();
    h_signal->Write();
    h_background->Write();
    h_weight_smooth->Write();
    h_data_isr->Write();
    fout->Close();

    // ------------------------------------------------------------------
    // 11. Visualise
    // ------------------------------------------------------------------
    TCanvas *c = new TCanvas("c", "ω peak correction (lower sideband, linear)", 1000, 600);
    c->Divide(2,1);
    c->cd(1);
    h_data_isr->SetMarkerStyle(20);
    h_data_isr->SetMarkerSize(0.6);
    h_data_isr->Draw("E");
    h_isr3pi->SetLineColor(kBlue);
    h_isr3pi->Draw("SameHist");
    h_isr3pi_corrected->SetLineColor(kRed);
    h_isr3pi_corrected->Draw("SameHist");
    total_func->SetLineColor(kGreen);
    total_func->SetLineStyle(kDashed);
    total_func->Draw("Same");
    h_signal->SetLineColor(kMagenta);
    h_signal->Draw("Same");
    h_background->SetLineColor(kMagenta);
    h_background->SetLineStyle(kDotted);
    h_background->Draw("Same");

    TLegend *leg = new TLegend(0.65,0.6,0.9,0.9);
    leg->AddEntry(h_data_isr, "Data - other backgrounds", "lep");
    leg->AddEntry(h_isr3pi, "Original ISR3pi MC", "l");
    leg->AddEntry(h_isr3pi_corrected, "Corrected ISR3pi MC", "l");
    leg->AddEntry(total_func, "Fit: α·signal_template + linear", "l");
    leg->AddEntry(h_signal, "Fitted signal (α·template)", "l");
    leg->AddEntry(h_background, "Fitted background (linear)", "l");
    leg->Draw();

    c->cd(2);
    h_weight_smooth->SetTitle("Correction weight for ISR3pi");
    h_weight_smooth->GetXaxis()->SetTitle(is_mev ? "M_{3π} [MeV]" : "M_{3π} [GeV]");
    h_weight_smooth->GetYaxis()->SetTitle("Weight");
    h_weight_smooth->Draw();
    c->SaveAs(output_path + "omega_correction_hybrid_lower_linear.pdf");

    gSigTemplate = nullptr;
    delete total_func;
    delete c;

    std::cout << "\nSaved " << output_path + "corrected_isr3pi_hybrid_lower_linear.root and omega_correction_hybrid_lower_linear.pdf\n";

    gSystem->Exit(0);
}
