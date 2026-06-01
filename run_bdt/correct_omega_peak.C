// correct_omega_peak.C
// Safe data-driven correction for ω → 3π shape in ISR3π MC.

/*
  Reads the data and MC trees
Opens the same ROOT file (tree_pre_bdt.root) used in the analysis, containing the Br_m3pi_bdt branch (3π invariant mass). Automatically detects whether the mass is in MeV or GeV.

Creates a data histogram (h_data) from the TDATA tree (no scaling).

Loads all scaled MC components using the latest scaling factors from the 2D fit: ../header_bdt/sfw2d.txt

    EEG, OmegaPi, KSL, EtaGamma, MC Rest, and the original ISR3π.

Subtracts all backgrounds from data to isolate the data‑driven ISR3π shape:
h_data_isr = h_data – (EEG + OmegaPi + KSL + EtaGamma + MC Rest).
Negative bins are set to zero.

Defines the ω peak region (e.g., 740–820 MeV) and sidebands for background estimation.

Fits the ω peak in the subtracted data histogram with a Breit‑Wigner + 2nd‑order polynomial model. The polynomial background is first fitted in the sidebands, then fixed during the peak fit.

Computes correction weights as the ratio (data_fit) / (original ISR3pi MC) bin‑by‑bin within the ω peak region (clamped between 0.4 and 2.5). Outside the peak, weights are set to 1.0.

Smooths the weights using a 3‑bin moving average to avoid statistical fluctuations.

Applies the weights to the original ISR3pi histogram to produce a corrected ISR3pi shape.

Renormalises the corrected histogram so that its integral in the ω peak region matches the original MC integral (preserving total yield).

Saves the corrected histogram (h_isr3pi_corrected) into a file called corrected_isr3pi.root, along with the weight histogram and the subtracted data histogram.

Produces a visualisation (a two‑panel PDF):

    Left panel: subtracted data points vs. original ISR3pi (blue) vs. corrected ISR3pi (red).

    Right panel: the correction weight curve.
*/

// correct_omega_peak.C
// Data-driven correction for ω → 3π shape using Breit‑Wigner signal + polynomial background.
// Outputs: corrected_isr3pi.root (containing h_isr3pi_corrected, h_signal, h_background, h_data_isr, h_weight_smooth)

// correct_omega_peak.C
// Data-driven correction for ω → 3π shape in ISR3π MC.
// Uses Breit‑Wigner signal + 2nd order polynomial background (sideband‑determined).
// Outputs: corrected_isr3pi.root and omega_correction.pdf

// correct_omega_peak.C
// Data-driven correction for ω → 3π shape in ISR3π MC.
// Uses Breit‑Wigner signal + Chebyshev polynomial background (3rd order).
// Outputs: corrected_isr3pi.root and omega_correction.pdf

// correct_omega_peak.C
// Data-driven correction for ω → 3π shape in ISR3π MC.
// Uses Breit‑Wigner signal + Chebyshev polynomial background (3rd order).
// Robust version: fixes background from sidebands, wide limits, fallback.
// Outputs: corrected_isr3pi.root and omega_correction.pdf

// correct_omega_peak.C
// Data-driven correction for ω → 3π shape in ISR3π MC.
// Uses Breit‑Wigner signal + Chebyshev polynomial background (2nd order).
// Robust version: fallback to upper sideband only, wide limits.
// Outputs: corrected_isr3pi.root and omega_correction.pdf

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFitResult.h>
#include <TLine.h>
#include <iostream>
#include <cmath>

#include "../header_bdt/sfw2d.txt"
#include "../header_bdt/correct_omega.h"   // defines output_path as TString
#include "../header_bdt/binning.h"

void correct_omega_peak() {
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
    if (!tdata) {
        std::cerr << "ERROR: TDATA tree not found." << std::endl;
        return;
    }
    if (!tdata->GetBranch("Br_m3pi_bdt")) {
        std::cerr << "ERROR: Branch Br_m3pi_bdt not found in TDATA." << std::endl;
        return;
    }

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
    // 2. Create data histogram (no scaling)
    // ------------------------------------------------------------------
    TH1D *h_data = new TH1D("h_data", "", nbins, low, high);
    h_data->Sumw2();
    tdata->SetBranchAddress("Br_m3pi_bdt", &mtest);
    Long64_t nentries = tdata->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tdata->GetEntry(i);
        h_data->Fill(mtest);
    }
    std::cout << "Data histogram integral: " << h_data->Integral() << std::endl;

    // ------------------------------------------------------------------
    // 3. Load scaled MC components
    // ------------------------------------------------------------------
    auto makeScaledHist = [&](const char* tname, double scale) -> TH1D* {
        TTree *t = (TTree*) ftree->Get(tname);
        if (!t) {
            std::cerr << "Warning: tree " << tname << " not found." << std::endl;
            return nullptr;
        }
        if (!t->GetBranch("Br_m3pi_bdt")) {
            std::cerr << "Warning: branch Br_m3pi_bdt missing in " << tname << std::endl;
            return nullptr;
        }
        TH1D *h = new TH1D(Form("h_%s", tname), "", nbins, low, high);
        h->Sumw2();
        double val;
        t->SetBranchAddress("Br_m3pi_bdt", &val);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            h->Fill(val);
        }
        h->Scale(scale);
        return h;
    };

    TH1D *h_eeg       = makeScaledHist("TEEG",         eeg_sfw);
    TH1D *h_omegapi   = makeScaledHist("TOMEGAPI",     omegapi_sfw);
    TH1D *h_ksl       = makeScaledHist("TKSL",         ksl_sfw);
    TH1D *h_etagam    = makeScaledHist("TETAGAM",      etagam_sfw);
    TH1D *h_isr3pi    = makeScaledHist("TISR3PI_SIG",  isr3pi_sfw);

    // MC Rest: combine TKPM, TRHOPI, TBKGREST
    TH1D *h_mcrest = new TH1D("h_mcrest", "", nbins, low, high);
    h_mcrest->Sumw2();
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

    if (!h_isr3pi) {
        std::cerr << "ERROR: Could not create ISR3pi histogram." << std::endl;
        return;
    }

    // ------------------------------------------------------------------
    // 4. Subtract backgrounds to isolate data‑driven ISR3π shape
    // ------------------------------------------------------------------
    TH1D *h_data_isr = (TH1D*) h_data->Clone("h_data_isr");
    if (h_eeg)     h_data_isr->Add(h_eeg, -1.0);
    if (h_omegapi) h_data_isr->Add(h_omegapi, -1.0);
    if (h_ksl)     h_data_isr->Add(h_ksl, -1.0);
    if (h_etagam)  h_data_isr->Add(h_etagam, -1.0);
    if (h_mcrest)  h_data_isr->Add(h_mcrest, -1.0);
    for (int bin = 1; bin <= h_data_isr->GetNbinsX(); ++bin) {
        if (h_data_isr->GetBinContent(bin) < 0) h_data_isr->SetBinContent(bin, 0);
    }

    // ------------------------------------------------------------------
    // 5. Define ω peak region and sidebands
    // ------------------------------------------------------------------
    double peak_low, peak_high;
    if (is_mev) { peak_low = 740.0; peak_high = 820.0; }
    else        { peak_low = 0.74;  peak_high = 0.82; }

    double sb_low1, sb_high1, sb_low2, sb_high2;
    if (is_mev) {
        sb_low1 = 670.0; sb_high1 = 730.0;
        sb_low2 = 830.0; sb_high2 = 890.0;
    } else {
        sb_low1 = 0.67; sb_high1 = 0.73;
        sb_low2 = 0.83; sb_high2 = 0.89;
    }

    // ------------------------------------------------------------------
    // 6. Fit background with 2nd-order Chebyshev (robust fallback)
    // ------------------------------------------------------------------
    TH1D *h_side = (TH1D*) h_data_isr->Clone("h_side");
    for (int bin = 1; bin <= h_side->GetNbinsX(); ++bin) {
        double x = h_side->GetBinCenter(bin);
        if (x >= peak_low && x <= peak_high) h_side->SetBinContent(bin, 0);
    }

    // 2nd order Chebyshev: c0 + c1*x + c2*(2x^2-1)
    TF1 *cheb_bkg = new TF1("cheb_bkg", 
        "[0] + [1]*x + [2]*(2*x*x - 1)",
        sb_low1, sb_high2);
    
    // Try both sidebands first
    int fitStat = h_side->Fit(cheb_bkg, "QNR", "", sb_low1, sb_high2);
    
    if (fitStat != 0) {
        std::cerr << "Both sidebands fit failed. Trying only upper sideband..." << std::endl;
        cheb_bkg->SetRange(sb_low2, sb_high2);
        fitStat = h_side->Fit(cheb_bkg, "QNR", "", sb_low2, sb_high2);
    }
    
    if (fitStat != 0) {
        std::cerr << "FATAL: Background fit failed even with upper sideband. Exiting." << std::endl;
        return;
    }
    
    double c0 = cheb_bkg->GetParameter(0);
    double c1 = cheb_bkg->GetParameter(1);
    double c2 = cheb_bkg->GetParameter(2);
    std::cout << "Sideband Chebyshev (2nd order) fit: c0 = " << c0 
              << "  c1 = " << c1 << "  c2 = " << c2 << std::endl;

    // Check background estimate under the peak
    double peak_int = h_data_isr->Integral(h_data_isr->FindBin(peak_low), h_data_isr->FindBin(peak_high));
    double bg_int_side = cheb_bkg->Integral(peak_low, peak_high);
    std::cout << "Data integral in peak region: " << peak_int << std::endl;
    std::cout << "Background integral under peak (from sideband fit): " << bg_int_side << std::endl;
    if (bg_int_side >= peak_int) {
        std::cerr << "WARNING: Sideband background exceeds data in peak. You may need to adjust sideband ranges." << std::endl;
    }

    // ------------------------------------------------------------------
    // 7. Total fit: Breit‑Wigner + Chebyshev (fix background for stability)
    // ------------------------------------------------------------------
    // Fix the Chebyshev background to sideband values (most stable)
    TF1 *total = new TF1("total", 
        "[0]*TMath::BreitWigner(x, [1], [2]) + ([3] + [4]*x + [5]*(2*x*x-1))",
        peak_low, peak_high);
    total->SetParNames("Norm", "Mean", "Width", "c0", "c1", "c2");

    // Fix background parameters
    total->SetParameter(3, c0);
    total->SetParameter(4, c1);
    total->SetParameter(5, c2);
    total->FixParameter(3, c0);
    total->FixParameter(4, c1);
    total->FixParameter(5, c2);
    std::cout << "Fixing Chebyshev background to sideband values." << std::endl;

    // Estimate initial signal normalisation
    double est_signal = peak_int - bg_int_side;
    if (est_signal < 100) est_signal = 100;
    double init_mean = is_mev ? 782.0 : 0.782;
    double init_width = is_mev ? 8.5 : 0.0085;
    total->SetParameter(0, est_signal);
    total->SetParameter(1, init_mean);
    total->SetParameter(2, init_width);
    total->SetParLimits(1, init_mean - 15.0, init_mean + 15.0);
    total->SetParLimits(2, 0.5, 30.0);

    // Perform fit
    std::cout << "Performing Breit‑Wigner + Chebyshev(2) fit in region [" 
              << peak_low << ", " << peak_high << "]" << std::endl;
    TFitResultPtr r = h_data_isr->Fit(total, "RQS", "", peak_low, peak_high);
    int fitStatus = r;
    if (fitStatus != 0) {
        std::cerr << "Fit failed (status " << fitStatus 
                  << "). Releasing background parameters..." << std::endl;
        total->ReleaseParameter(3);
        total->ReleaseParameter(4);
        total->ReleaseParameter(5);
        total->SetParLimits(3, c0 - 0.5*fabs(c0), c0 + 0.5*fabs(c0));
        total->SetParLimits(4, c1 - 0.5*fabs(c1), c1 + 0.5*fabs(c1));
        total->SetParLimits(5, c2 - 0.5*fabs(c2), c2 + 0.5*fabs(c2));
        r = h_data_isr->Fit(total, "RQS", "", peak_low, peak_high);
        fitStatus = r;
    }

    bool fit_succeeded = (fitStatus == 0);
    if (!fit_succeeded) {
        std::cerr << "WARNING: Fit still failed. Using sideband Chebyshev as background only." 
                  << std::endl;
    } else {
        std::cout << "Fit succeeded!" << std::endl;
        for (int i=0; i<6; ++i) {
            std::cout << "  " << total->GetParName(i) << " = " << total->GetParameter(i)
                      << " +/- " << total->GetParError(i) << std::endl;
        }
    }

    // ------------------------------------------------------------------
    // 8. Create signal and background histograms
    // ------------------------------------------------------------------
    TH1D *h_signal = (TH1D*) h_data_isr->Clone("h_signal");
    TH1D *h_background = (TH1D*) h_data_isr->Clone("h_background");
    h_signal->Reset();
    h_background->Reset();

    if (fit_succeeded) {
        double norm = total->GetParameter(0);
        double mean = total->GetParameter(1);
        double width = total->GetParameter(2);
        double fc0 = total->GetParameter(3);
        double fc1 = total->GetParameter(4);
        double fc2 = total->GetParameter(5);
        for (int bin = 1; bin <= h_signal->GetNbinsX(); ++bin) {
            double x = h_signal->GetBinCenter(bin);
            double bw = norm * TMath::BreitWigner(x, mean, width);
            double cheb = fc0 + fc1*x + fc2*(2*x*x - 1);
            if (bw < 0) bw = 0;
            if (cheb < 0) cheb = 0;
            h_signal->SetBinContent(bin, bw);
            h_background->SetBinContent(bin, cheb);
            h_signal->SetBinError(bin, TMath::Sqrt(bw));
            h_background->SetBinError(bin, TMath::Sqrt(cheb));
        }
    } else {
        // Fallback: use sideband Chebyshev as background, signal = data - background
        for (int bin = 1; bin <= h_background->GetNbinsX(); ++bin) {
            double x = h_background->GetBinCenter(bin);
            double cheb = c0 + c1*x + c2*(2*x*x - 1);
            if (cheb < 0) cheb = 0;
            h_background->SetBinContent(bin, cheb);
            double data_val = h_data_isr->GetBinContent(bin);
            double sig_val = data_val - cheb;
            if (sig_val < 0) sig_val = 0;
            h_signal->SetBinContent(bin, sig_val);
            h_signal->SetBinError(bin, TMath::Sqrt(sig_val));
        }
    }

    // ------------------------------------------------------------------
    // 9. Compute correction weights (data / MC) within peak region
    // ------------------------------------------------------------------
    TH1D *h_weight = (TH1D*) h_data_isr->Clone("h_weight");
    h_weight->Reset();
    for (int bin = 1; bin <= h_weight->GetNbinsX(); ++bin) {
        double x = h_weight->GetBinCenter(bin);
        double data_val = h_data_isr->GetBinContent(bin);
        double mc_val = h_isr3pi->GetBinContent(bin);
        if (x >= peak_low && x <= peak_high && mc_val > 0) {
            double ratio = data_val / mc_val;
            if (ratio > 2.5) ratio = 2.5;
            if (ratio < 0.4) ratio = 0.4;
            h_weight->SetBinContent(bin, ratio);
        } else {
            h_weight->SetBinContent(bin, 1.0);
        }
    }

    // Smooth weights with 3‑bin moving average
    TH1D *h_weight_smooth = (TH1D*) h_weight->Clone("h_weight_smooth");
    for (int bin = 2; bin <= h_weight_smooth->GetNbinsX()-1; ++bin) {
        double w_avg = (h_weight->GetBinContent(bin-1) +
                        h_weight->GetBinContent(bin) +
                        h_weight->GetBinContent(bin+1)) / 3.0;
        h_weight_smooth->SetBinContent(bin, w_avg);
    }

    // ------------------------------------------------------------------
    // 10. Apply correction to ISR3π MC
    // ------------------------------------------------------------------
    TH1D *h_isr3pi_corrected = (TH1D*) h_isr3pi->Clone("h_isr3pi_corrected");
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
    // 11. Save results to ROOT file
    // ------------------------------------------------------------------
    TFile *fout = new TFile(output_path + "corrected_isr3pi.root", "RECREATE");
    h_isr3pi_corrected->Write();
    h_signal->Write();
    h_background->Write();
    h_weight_smooth->Write();
    h_data_isr->Write();
    fout->Close();

    // ------------------------------------------------------------------
    // 12. Visualisation
    // ------------------------------------------------------------------
    TCanvas *c = new TCanvas("c", "ω peak correction (Breit‑Wigner + Chebyshev 2nd order)", 1200, 600);
    c->Divide(2,1);

    // Left panel: data after background subtraction, MC shapes, and fit
    c->cd(1);
    h_data_isr->SetMarkerStyle(20);
    h_data_isr->SetMarkerSize(0.6);
    h_data_isr->SetStats(0);
    h_data_isr->GetYaxis()->SetTitle(Form("Events / [%.1f %s]", 
        (is_mev ? h_data_isr->GetBinWidth(1) : h_data_isr->GetBinWidth(1)*1000),
        (is_mev ? "MeV/c^{2}" : "GeV/c^{2}")));
    h_data_isr->GetXaxis()->SetTitle(is_mev ? "M_{3#pi} [MeV/c^{2}]" : "M_{3#pi} [GeV/c^{2}]");
    h_data_isr->Draw("E0");
    h_isr3pi->SetLineColor(kBlue);
    h_isr3pi->Draw("hist same");
    h_isr3pi_corrected->SetLineColor(kRed);
    h_isr3pi_corrected->Draw("hist same");
    if (fit_succeeded) {
        total->SetLineColor(kGreen);
        total->SetLineStyle(kDashed);
        total->Draw("same");
    }
    h_signal->SetLineColor(kMagenta);
    h_signal->SetLineStyle(kSolid);
    h_signal->Draw("same");
    h_background->SetLineColor(kMagenta);
    h_background->SetLineStyle(kDotted);
    h_background->Draw("same");

    TLegend *leg = new TLegend(0.55, 0.55, 0.9, 0.9);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(h_data_isr, "Data - other backgrounds", "lep");
    leg->AddEntry(h_isr3pi, "Original ISR3pi MC", "l");
    leg->AddEntry(h_isr3pi_corrected, "Corrected ISR3pi MC", "l");
    if (fit_succeeded) leg->AddEntry(total, "Breit‑Wigner + Chebyshev(2) fit", "l");
    leg->AddEntry(h_signal, "Signal component (Breit‑Wigner)", "l");
    leg->AddEntry(h_background, "Non‑resonant background (Chebyshev)", "l");
    leg->Draw();

    // Right panel: correction weights
    c->cd(2);
    h_weight_smooth->SetTitle("");
    h_weight_smooth->GetXaxis()->SetTitle(is_mev ? "M_{3#pi} [MeV/c^{2}]" : "M_{3#pi} [GeV/c^{2}]");
    h_weight_smooth->GetYaxis()->SetTitle("Correction weight");
    h_weight_smooth->SetLineColor(kBlack);
    h_weight_smooth->SetMarkerStyle(20);
    h_weight_smooth->SetMarkerSize(0.6);
    h_weight_smooth->Draw("PE0");
    TLine *line = new TLine(low, 1.0, high, 1.0);
    line->SetLineStyle(2);
    line->Draw();

    c->SaveAs(output_path + "omega_correction.pdf");
    delete c;

    std::cout << "\nAll done. Saved:\n"
              << "  " << output_path << "corrected_isr3pi.root\n"
              << "  " << output_path << "omega_correction.pdf\n";

    ftree->Close();
    delete ftree;
}

int main() {
    correct_omega_peak();
    return 0;
}
