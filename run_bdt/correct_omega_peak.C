// correct_omega_peak_safe.C
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
// Data-driven correction for ω → 3π shape using Crystal Ball signal + polynomial background.
// Outputs: corrected_isr3pi.root (containing h_isr3pi_corrected, h_signal, h_background, h_data_isr, h_weight_smooth)

// correct_omega_peak.C
// Data-driven correction for ω → 3π shape using Crystal Ball signal + polynomial background.
// Outputs: corrected_isr3pi.root (containing h_isr3pi_corrected, h_signal, h_background, h_data_isr, h_weight_smooth)

// correct_omega_peak.C (minimal adaptation: Breit‑Wigner → Crystal Ball)
// ... (same header comments) ...

// correct_omega_peak.C
// Data-driven correction for ω → 3π shape using Breit‑Wigner signal + polynomial background.
// Outputs: corrected_isr3pi.root (containing h_isr3pi_corrected, h_signal, h_background, h_data_isr, h_weight_smooth)

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <cmath>

#include "../header_bdt/sfw2d.txt"

void correct_omega_peak() {
    // ------------------------------------------------------------------
    // 1. Open tree file and verify branches
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

    // Determine unit (MeV or GeV) from first entry
    double mtest;
    tdata->SetBranchAddress("Br_m3pi_bdt", &mtest);
    tdata->GetEntry(0);
    bool is_mev = (mtest > 10);
    double low, high;
    if (is_mev) { low = 600; high = 1000; }
    else { low = 0.6; high = 1.0; }
    int nbins = 200;
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

    TH1D *h_eeg       = makeScaledHist("TEEG", eeg_sfw);
    TH1D *h_omegapi   = makeScaledHist("TOMEGAPI", omegapi_sfw);
    TH1D *h_ksl       = makeScaledHist("TKSL", ksl_sfw);
    TH1D *h_etagam    = makeScaledHist("TETAGAM", etagam_sfw);
    TH1D *h_isr3pi    = makeScaledHist("TISR3PI_SIG", isr3pi_sfw);

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
    // 4. Subtract backgrounds to isolate ISR3π component
    // ------------------------------------------------------------------
    TH1D *h_data_isr = (TH1D*) h_data->Clone("h_data_isr");
    if (h_eeg) h_data_isr->Add(h_eeg, -1.0);
    if (h_omegapi) h_data_isr->Add(h_omegapi, -1.0);
    if (h_ksl) h_data_isr->Add(h_ksl, -1.0);
    if (h_etagam) h_data_isr->Add(h_etagam, -1.0);
    if (h_mcrest) h_data_isr->Add(h_mcrest, -1.0);
    for (int bin = 1; bin <= h_data_isr->GetNbinsX(); ++bin) {
        if (h_data_isr->GetBinContent(bin) < 0) h_data_isr->SetBinContent(bin, 0);
    }

    // ------------------------------------------------------------------
    // 5. Define ω peak region and sidebands
    // ------------------------------------------------------------------
    double peak_low, peak_high;
    if (is_mev) { peak_low = 740; peak_high = 820; }
    else { peak_low = 0.74; peak_high = 0.82; }

    double sb_low1, sb_high1, sb_low2, sb_high2;
    if (is_mev) {
        sb_low1 = 700; sb_high1 = 740;
        sb_low2 = 820; sb_high2 = 900;
    } else {
        sb_low1 = 0.70; sb_high1 = 0.74;
        sb_low2 = 0.82; sb_high2 = 0.90;
    }

    // ------------------------------------------------------------------
    // 6. Fit with Breit‑Wigner signal + polynomial background
    // ------------------------------------------------------------------
    TH1D *h_fit = (TH1D*) h_data_isr->Clone("h_fit");

    // First, fit a second‑order polynomial to sidebands to estimate background
    TF1 *bkg = new TF1("bkg", "pol2", sb_low1, sb_high2);
    h_fit->Fit(bkg, "QN", "", sb_low1, sb_high1);
    h_fit->Fit(bkg, "QN+", "", sb_low2, sb_high2);
    double b0 = bkg->GetParameter(0);
    double b1 = bkg->GetParameter(1);
    double b2 = bkg->GetParameter(2);

    // Create total model: Breit‑Wigner signal + polynomial
    TF1 *total = new TF1("total", "[0]*TMath::BreitWigner(x, [1], [2]) + pol2(3)", peak_low, peak_high);
    total->SetParNames("Norm", "Mean", "Width", "B0", "B1", "B2");

    double init_mass = is_mev ? 782.0 : 0.782;
    double init_width = is_mev ? 8.5 : 0.0085;
    total->SetParameters(100, init_mass, init_width, b0, b1, b2);

    // Set limits to help convergence
    total->SetParLimits(1, init_mass - 10, init_mass + 10);
    total->SetParLimits(2, 1.0, 20.0);                // width
    // Allow polynomial to vary, but limit near sideband values
    total->SetParLimits(3, b0 - 0.5*std::abs(b0), b0 + 0.5*std::abs(b0));
    total->SetParLimits(4, b1 - 0.5*std::abs(b1), b1 + 0.5*std::abs(b1));
    total->SetParLimits(5, b2 - 0.5*std::abs(b2), b2 + 0.5*std::abs(b2));

    // Perform fit
    int fitStatus = h_fit->Fit(total, "RQ", "", peak_low, peak_high);
    if (fitStatus != 0) {
        std::cerr << "Fit did not converge. Trying with relaxed polynomial limits..." << std::endl;
        total->SetParLimits(3, -1e6, 1e6);
        total->SetParLimits(4, -1e6, 1e6);
        total->SetParLimits(5, -1e6, 1e6);
        fitStatus = h_fit->Fit(total, "RQ", "", peak_low, peak_high);
    }

    std::cout << "Fit status: " << fitStatus << std::endl;
    for (int i=0; i<6; ++i) {
        std::cout << "  " << total->GetParName(i) << " = " << total->GetParameter(i)
                  << " +/- " << total->GetParError(i) << std::endl;
    }

    // ------------------------------------------------------------------
    // 7. Create separate signal (Breit‑Wigner) and background (polynomial) histograms
    // ------------------------------------------------------------------
    TH1D *h_signal = (TH1D*) h_data_isr->Clone("h_signal");
    TH1D *h_background = (TH1D*) h_data_isr->Clone("h_background");
    h_signal->Reset();
    h_background->Reset();

    double bw_norm = total->GetParameter(0);
    double bw_mean = total->GetParameter(1);
    double bw_width = total->GetParameter(2);
    double poly0 = total->GetParameter(3);
    double poly1 = total->GetParameter(4);
    double poly2 = total->GetParameter(5);

    for (int bin = 1; bin <= h_signal->GetNbinsX(); ++bin) {
        double x = h_signal->GetBinCenter(bin);
        double bw_val = bw_norm * TMath::BreitWigner(x, bw_mean, bw_width);
        double poly_val = poly0 + poly1*x + poly2*x*x;
        if (bw_val < 0) bw_val = 0;
        if (poly_val < 0) poly_val = 0;
        h_signal->SetBinContent(bin, bw_val);
        h_background->SetBinContent(bin, poly_val);
        h_signal->SetBinError(bin, TMath::Sqrt(bw_val));
        h_background->SetBinError(bin, TMath::Sqrt(poly_val));
    }

    // ------------------------------------------------------------------
    // 8. Compute correction weights using the data histogram
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

    // Smooth weights
    TH1D *h_weight_smooth = (TH1D*) h_weight->Clone("h_weight_smooth");
    for (int bin = 2; bin <= h_weight_smooth->GetNbinsX()-1; ++bin) {
        double w_avg = (h_weight->GetBinContent(bin-1) +
                        h_weight->GetBinContent(bin) +
                        h_weight->GetBinContent(bin+1)) / 3.0;
        h_weight_smooth->SetBinContent(bin, w_avg);
    }

    // ------------------------------------------------------------------
    // 9. Apply correction to ISR3π MC
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
    // 10. Save results
    // ------------------------------------------------------------------
    TFile *fout = new TFile("corrected_isr3pi.root", "RECREATE");
    h_isr3pi_corrected->Write();
    h_signal->Write();
    h_background->Write();
    h_weight_smooth->Write();
    h_data_isr->Write();
    fout->Close();

    // ------------------------------------------------------------------
    // 11. Visualise
    // ------------------------------------------------------------------
    TCanvas *c = new TCanvas("c", "ω peak correction (Breit‑Wigner)", 1000, 600);
    c->Divide(2,1);
    c->cd(1);
    h_data_isr->SetMarkerStyle(20);
    h_data_isr->SetMarkerSize(0.6);
    h_data_isr->Draw("E");
    h_isr3pi->SetLineColor(kBlue);
    h_isr3pi->Draw("SameHist");
    h_isr3pi_corrected->SetLineColor(kRed);
    h_isr3pi_corrected->Draw("SameHist");
    total->SetLineColor(kGreen);
    total->SetLineStyle(kDashed);
    total->Draw("Same");
    h_signal->SetLineColor(kMagenta);
    h_signal->SetLineStyle(kSolid);
    h_signal->Draw("Same");
    h_background->SetLineColor(kMagenta);
    h_background->SetLineStyle(kDotted);
    h_background->Draw("Same");

    TLegend *leg = new TLegend(0.65,0.6,0.9,0.9);
    leg->AddEntry(h_data_isr, "Data - other backgrounds", "lep");
    leg->AddEntry(h_isr3pi, "Original ISR3pi MC", "l");
    leg->AddEntry(h_isr3pi_corrected, "Corrected ISR3pi MC", "l");
    leg->AddEntry(total, "Breit‑Wigner + polynomial fit", "l");
    leg->AddEntry(h_signal, "ω → 3π (Breit‑Wigner)", "l");
    leg->AddEntry(h_background, "Non‑resonant background", "l");
    leg->Draw();

    c->cd(2);
    h_weight_smooth->SetTitle("Correction weight for ISR3pi");
    h_weight_smooth->GetXaxis()->SetTitle(is_mev ? "M_{3π} [MeV]" : "M_{3π} [GeV]");
    h_weight_smooth->GetYaxis()->SetTitle("Weight");
    h_weight_smooth->Draw();
    c->SaveAs("omega_correction_safe.pdf");

    std::cout << "\nSaved corrected_isr3pi.root and omega_correction_safe.pdf\n";

    ftree->Close();
}
