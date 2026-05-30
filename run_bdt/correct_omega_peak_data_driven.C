// correct_omega_peak_datadriven.C
// Pure data‑driven sideband fit for ω peak correction.
//   - Signal template from correctly paired ω MC (recon_indx_bdt==2 && bkg_indx==1)
//   - Background: 2nd‑order polynomial fitted to data sidebands (680‑730 MeV, 830‑880 MeV)
//   - No MC background subtraction – sidebands taken directly from data.
// Outputs: corrected_isr3pi_datadriven.root (h_isr3pi_corrected, h_signal, h_background, h_data, h_weight_smooth)

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <iostream>
#include <cmath>
#include <vector>

#include "../header_bdt/sfw2d.txt"       // provides eeg_sfw, isr3pi_sfw, etc.
#include "../header_bdt/correct_omega.h" // provides output_path, etc.

// Global pointer to signal template (detached)
TH1D *gSigTemplate = nullptr;

// Fit function: α * signal_template + polynomial
Double_t template_poly(Double_t *x, Double_t *par) {
    int bin = gSigTemplate->FindBin(x[0]);
    Double_t sig = gSigTemplate->GetBinContent(bin);
    Double_t poly = par[1] + par[2] * x[0] + par[3] * x[0] * x[0];
    return par[0] * sig + poly;
}

void correct_omega_peak_datadriven() {
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
    double low = is_mev ? 600 : 0.6;
    double high = is_mev ? 1000 : 1.0;
    int nbins = 200;
    std::cout << "Mass unit: " << (is_mev ? "MeV" : "GeV")
              << " range [" << low << ", " << high << "]\n";

    // ------------------------------------------------------------------
    // 2. Data histogram (raw, no MC subtraction)
    // ------------------------------------------------------------------
    TH1D *h_data = new TH1D("h_data", "", nbins, low, high);
    h_data->Sumw2();
    h_data->SetDirectory(0);
    tdata->SetBranchAddress("Br_m3pi_bdt", &mtest);
    for (Long64_t i = 0; i < tdata->GetEntries(); ++i) {
        tdata->GetEntry(i);
        h_data->Fill(mtest);
    }
    std::cout << "Data integral (full range): " << h_data->Integral() << std::endl;

    // ------------------------------------------------------------------
    // 3. Load original ISR3pi MC (to be corrected)
    // ------------------------------------------------------------------
    TTree *tisr = (TTree*) ftree->Get("TISR3PI_SIG");
    if (!tisr) { std::cerr << "ERROR: TISR3PI_SIG not found." << std::endl; return; }
    TH1D *h_isr3pi_orig = new TH1D("h_isr3pi_orig", "", nbins, low, high);
    h_isr3pi_orig->Sumw2();
    h_isr3pi_orig->SetDirectory(0);
    double val;
    tisr->SetBranchAddress("Br_m3pi_bdt", &val);
    for (Long64_t i = 0; i < tisr->GetEntries(); ++i) {
        tisr->GetEntry(i);
        h_isr3pi_orig->Fill(val);
    }
    // Apply original scaling factor (from sfw2d.txt)
    h_isr3pi_orig->Scale(isr3pi_sfw);
    std::cout << "Original ISR3pi MC integral: " << h_isr3pi_orig->Integral() << std::endl;

    // ------------------------------------------------------------------
    // 4. Create signal template (correctly paired ω in MC)
    // ------------------------------------------------------------------
    int recon_indx_bdt = 0, bkg_indx = 0;
    double m3pi = 0.;
    tisr->SetBranchAddress("Br_recon_indx_bdt", &recon_indx_bdt);
    tisr->SetBranchAddress("Br_bkg_indx", &bkg_indx);
    tisr->SetBranchAddress("Br_m3pi_bdt", &m3pi);

    TH1D *h_signal_template = new TH1D("h_signal_template", "", nbins, low, high);
    h_signal_template->Sumw2();
    h_signal_template->SetDirectory(0);
    for (Long64_t i = 0; i < tisr->GetEntries(); ++i) {
        tisr->GetEntry(i);
        if (recon_indx_bdt == 2 && bkg_indx == 1)
            h_signal_template->Fill(m3pi);
    }
    // Normalise to unit area
    double sig_int = h_signal_template->Integral();
    if (sig_int > 0) h_signal_template->Scale(1.0 / sig_int);
    gSigTemplate = h_signal_template;

    // ------------------------------------------------------------------
    // 5. Define fit regions (avoid the bump 810–900 MeV)
    // ------------------------------------------------------------------
    double peak_low  = is_mev ? 740 : 0.74;
    double peak_high = is_mev ? 820 : 0.82;
    double sb_low1   = is_mev ? 680 : 0.68;
    double sb_high1  = is_mev ? 730 : 0.73;
    double sb_low2   = is_mev ? 830 : 0.83;
    double sb_high2  = is_mev ? 880 : 0.88;

    // ------------------------------------------------------------------
    // 6. Fit polynomial to data sidebands (no MC subtraction)
    // ------------------------------------------------------------------
    TH1D *h_side = (TH1D*) h_data->Clone("h_side");
    for (int bin = 1; bin <= h_side->GetNbinsX(); ++bin) {
        double x = h_side->GetBinCenter(bin);
        if (x >= peak_low && x <= peak_high) h_side->SetBinContent(bin, 0);
    }

    TF1 *bkg_poly = new TF1("bkg_poly", "pol2", sb_low1, sb_high2);
    h_side->Fit(bkg_poly, "QN", "", sb_low1, sb_high1);
    h_side->Fit(bkg_poly, "QN+", "", sb_low2, sb_high2);
    double p0 = bkg_poly->GetParameter(0);
    double p1 = bkg_poly->GetParameter(1);
    double p2 = bkg_poly->GetParameter(2);
    std::cout << "Sideband polynomial: " << p0 << " + " << p1 << "·x + " << p2 << "·x²\n";

    // ------------------------------------------------------------------
    // 7. Fit peak region: α * signal_template + polynomial (fixed)
    // ------------------------------------------------------------------
    TF1 *total_func = new TF1("total_func", template_poly, low, high, 4);
    total_func->SetParameters(1000, p0, p1, p2);
    total_func->SetParNames("alpha", "p0", "p1", "p2");
    // Fix polynomial parameters to sideband values (or allow tight limits)
    total_func->FixParameter(1, p0);
    total_func->FixParameter(2, p1);
    total_func->FixParameter(3, p2);
    // Let only alpha float
    total_func->SetParName(0, "alpha");

    h_data->Fit(total_func, "RQ", "", peak_low, peak_high);
    double alpha = total_func->GetParameter(0);
    std::cout << "Fitted signal yield α = " << alpha << std::endl;

    // ------------------------------------------------------------------
    // 8. Create signal and background histograms (for visualisation)
    // ------------------------------------------------------------------
    TH1D *h_signal = (TH1D*) h_signal_template->Clone("h_signal");
    h_signal->SetDirectory(0);
    h_signal->Scale(alpha);

    TH1D *h_background = (TH1D*) h_data->Clone("h_background");
    h_background->SetDirectory(0);
    h_background->Reset();
    for (int bin = 1; bin <= h_background->GetNbinsX(); ++bin) {
        double x = h_background->GetBinCenter(bin);
        double poly_val = p0 + p1*x + p2*x*x;
        if (poly_val < 0) poly_val = 0;
        h_background->SetBinContent(bin, poly_val);
        h_background->SetBinError(bin, TMath::Sqrt(poly_val));
    }

    // ------------------------------------------------------------------
    // 9. Compute correction weights (desired shape = h_signal)
    //    to apply to original ISR3pi MC
    // ------------------------------------------------------------------
    TH1D *h_weight = (TH1D*) h_data->Clone("h_weight");
    h_weight->SetDirectory(0);
    h_weight->Reset();
    for (int bin = 1; bin <= h_weight->GetNbinsX(); ++bin) {
        double x = h_weight->GetBinCenter(bin);
        double mc_val = h_isr3pi_orig->GetBinContent(bin);
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
    // 10. Apply correction to ISR3pi MC
    // ------------------------------------------------------------------
    TH1D *h_isr3pi_corrected = (TH1D*) h_isr3pi_orig->Clone("h_isr3pi_corrected");
    h_isr3pi_corrected->SetDirectory(0);
    for (int bin = 1; bin <= h_isr3pi_corrected->GetNbinsX(); ++bin) {
        double w = h_weight_smooth->GetBinContent(bin);
        double old = h_isr3pi_corrected->GetBinContent(bin);
        h_isr3pi_corrected->SetBinContent(bin, old * w);
        double err = h_isr3pi_corrected->GetBinError(bin);
        h_isr3pi_corrected->SetBinError(bin, err * w);
    }
    // Renormalise to keep integral in peak region unchanged
    double orig_int = h_isr3pi_orig->Integral(peak_low, peak_high);
    double new_int = h_isr3pi_corrected->Integral(peak_low, peak_high);
    if (new_int > 0 && orig_int > 0) {
        double renorm = orig_int / new_int;
        h_isr3pi_corrected->Scale(renorm);
        std::cout << "Renormalisation factor: " << renorm << std::endl;
    }

    // ------------------------------------------------------------------
    // 11. Save results (close input file first)
    // ------------------------------------------------------------------
    ftree->Close();
    delete ftree;

    TString out_path = output_path;   // from correct_omega.h
    TFile *fout = new TFile(out_path + "corrected_isr3pi_datadriven.root", "RECREATE");
    h_isr3pi_corrected->Write();
    h_signal->Write();
    h_background->Write();
    h_weight_smooth->Write();
    h_data->Write();
    fout->Close();

    // ------------------------------------------------------------------
    // 12. Visualise
    // ------------------------------------------------------------------
    TCanvas *c = new TCanvas("c", "ω peak correction (Data‑driven sideband)", 1000, 600);
    c->Divide(2,1);
    c->cd(1);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.6);
    h_data->Draw("E");
    h_isr3pi_orig->SetLineColor(kBlue);
    h_isr3pi_orig->Draw("SameHist");
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
    leg->AddEntry(h_data, "Data (BDT selected)", "lep");
    leg->AddEntry(h_isr3pi_orig, "Original ISR3pi MC", "l");
    leg->AddEntry(h_isr3pi_corrected, "Corrected ISR3pi MC", "l");
    leg->AddEntry(total_func, "Fit: α·signal_template + polynomial", "l");
    leg->AddEntry(h_signal, "Fitted signal (α·template)", "l");
    leg->AddEntry(h_background, "Fitted background (polynomial)", "l");
    leg->Draw();

    c->cd(2);
    h_weight_smooth->SetTitle("Correction weight for ISR3pi");
    h_weight_smooth->GetXaxis()->SetTitle(is_mev ? "M_{3π} [MeV]" : "M_{3π} [GeV]");
    h_weight_smooth->GetYaxis()->SetTitle("Weight");
    h_weight_smooth->Draw();
    c->SaveAs(out_path + "omega_correction_datadriven.pdf");

    // ------------------------------------------------------------------
    // 13. Clean up (fixed: no delete total_func, canvas owns it)
    // ------------------------------------------------------------------
    delete bkg_poly;
    delete h_side;
    delete c;   // canvas deletes total_func automatically

    std::cout << "\nSaved " << out_path << "corrected_isr3pi_datadriven.root and omega_correction_datadriven.pdf\n";
}
