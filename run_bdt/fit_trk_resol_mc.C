// ============================================================================
// fit_trk_resol_mc.C
//
// Fit the MC resolution vs. true dipion mass with a modified logistic function
// where the plateau is allowed to have a linear slope:
//   sigma(M) = low + (high(M) - low) / (1 + exp(-k*(M - M0)))
//   high(M) = h0 + h1 * (M - Mref)
//
// This provides a smooth, mass‑dependent resolution with a non‑constant plateau.
// Usage:
//   .x fit_logistic_mc.C()
// ============================================================================

#include <TFile.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>

// ----------------------------------------------------------------------------
// Modified logistic: plateau is linear in mass
// ----------------------------------------------------------------------------
Double_t logistic_linear(Double_t *x, Double_t *par) {
    Double_t M = x[0];
    // par[0] = low, par[1] = h0 (plateau at Mref), par[2] = h1 (slope),
    // par[3] = k, par[4] = M0, par[5] = Mref
    Double_t Mref = par[5];
    Double_t high = par[1] + par[2] * (M - Mref);
    Double_t exponent = -par[3] * (M - par[4]);
    Double_t denom = 1.0 + TMath::Exp(exponent);
    return par[0] + (high - par[0]) / denom;
}

int fit_trk_resol_mc(const TString input_file = "../trackmass_scan/pull_scan_TISR3PI_SIG_PEAK.root") {

    // ----------------------------------------------------------------------
    // Input file (output from trackmass_scan on TISR3PI_SIG_PEAK)
    // ----------------------------------------------------------------------
    //TString input_file = "../trackmass_scan/pull_scan_TISR3PI_SIG_PEAK.root";
    TFile *f = TFile::Open(input_file);
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open " << input_file << std::endl;
        return 1;
    }

    TGraphErrors *g = (TGraphErrors*)f->Get("g_resolution_vs_M");
    if (!g) {
        std::cerr << "ERROR: TGraph 'g_resolution_vs_M' not found." << std::endl;
        f->Close();
        return 1;
    }

    TGraphErrors *g_origin = (TGraphErrors*)g->Clone();
    
    // ----------------------------------------------------------------------
    // Fit with modified logistic (linear plateau)
    // ----------------------------------------------------------------------
    // Parameters: [0]=low, [1]=h0, [2]=h1, [3]=k, [4]=M0, [5]=Mref
    // We fix Mref = 500 MeV (middle of the plateau region)
    const Double_t Mref = 500.0;
    TF1 *func = new TF1("logistic_linear", logistic_linear, 200, 620, 6);
    func->SetParameters(0.4, 1.1, 0.0, 0.05, 350.0, Mref);
    func->SetParLimits(0, 0.0, 0.8);    // low
    func->SetParLimits(1, 0.8, 1.5);    // h0 (plateau at Mref)
    func->SetParLimits(2, -0.01, 0.01); // h1 (slope)
    func->SetParLimits(3, 0.01, 0.3);   // k
    //func->SetParLimits(4, 200, 420);    // M0
    // Par[5] is fixed: we keep it fixed in the fit by not setting a limit and just letting it be the given value

    // Perform fit
    g->Fit(func, "RQS");

    double low   = func->GetParameter(0);
    double low_err = func->GetParError(0);
    double h0    = func->GetParameter(1);
    double h0_err = func->GetParError(1);
    double h1    = func->GetParameter(2);
    double h1_err = func->GetParError(2);
    double k     = func->GetParameter(3);
    double k_err = func->GetParError(3);
    double M0    = func->GetParameter(4);
    double M0_err = func->GetParError(4);

    double chi2  = func->GetChisquare();
    int ndf      = func->GetNDF();
    double chi2ndf = (ndf > 0) ? chi2 / ndf : 0;

    std::cout << "\n=== Modified logistic fit results (linear plateau) ===" << std::endl;
    std::cout << "low  = " << low << " +/- " << low_err << " MeV/c^{2}" << std::endl;
    std::cout << "h0   = " << h0 << " +/- " << h0_err << " MeV/c^{2} (plateau at Mref = " << Mref << " MeV)" << std::endl;
    std::cout << "h1   = " << h1 << " +/- " << h1_err << std::endl;
    std::cout << "k    = " << k << " +/- " << k_err << " 1/MeV/c^{2}" << std::endl;
    std::cout << "M0   = " << M0 << " +/- " << M0_err << " MeV/c^{2}" << std::endl;
    std::cout << "chi2/ndf = " << chi2ndf << std::endl;

    // ----------------------------------------------------------------------
    // Also fit a standard logistic for comparison (optional)
    // ----------------------------------------------------------------------
    TF1 *logistic_std = new TF1("logistic_std", "[0]+([1]-[0])/(1+exp(-[2]*(x-[3])))", 260, 620);
    logistic_std->SetParameters(0.4, 1.1, 0.05, 350.0);
    logistic_std->SetParLimits(0, 0.0, 0.8);
    logistic_std->SetParLimits(1, 0.8, 1.5);
    logistic_std->SetParLimits(2, 0.01, 0.3);
    logistic_std->SetParLimits(3, 300, 420);
    g->Fit(logistic_std, "RQS");
    double high_std = logistic_std->GetParameter(1);
    double high_std_err = logistic_std->GetParError(1);
    double chi2std = logistic_std->GetChisquare();
    int ndfstd = logistic_std->GetNDF();
    double chi2ndf_std = (ndfstd > 0) ? chi2std / ndfstd : 0;

    std::cout << "\n=== Standard logistic fit (for comparison) ===" << std::endl;
    std::cout << "high = " << high_std << " +/- " << high_std_err << " MeV" << std::endl;
    std::cout << "chi2/ndf = " << chi2ndf_std << std::endl;

    // ----------------------------------------------------------------------
    // Save parameters to header file
    // ----------------------------------------------------------------------
    std::ofstream out("../trackmass_scan/mc_resolution_params.h");
    out << "// MC resolution parameters from modified logistic fit\n";
    out << "// sigma(M) = low + (high(M)-low)/(1+exp(-k*(M-M0)))\n";
    out << "// high(M) = h0 + h1*(M - Mref)  with Mref = " << Mref << " MeV\n";
    out << "const double SIGMA_LOW  = " << low << ";\n";
    out << "const double SIGMA_LOW_ERR  = " << low_err << ";\n";
    out << "const double H0         = " << h0 << ";\n";
    out << "const double H0_ERR     = " << h0_err << ";\n";
    out << "const double H1         = " << h1 << ";\n";
    out << "const double H1_ERR     = " << h1_err << ";\n";
    out << "const double K          = " << k << ";\n";
    out << "const double K_ERR      = " << k_err << ";\n";
    out << "const double M0         = " << M0 << ";\n";
    out << "const double M0_ERR     = " << M0_err << ";\n";
    out << "const double MREF       = " << Mref << ";\n";
    out << "const double CHI2_NDF   = " << chi2ndf << ";\n\n";
    out << "// Standard logistic plateau (constant) for comparison\n";
    out << "const double SIGMA_HIGH_CONST = " << high_std << ";\n";
    out << "const double SIGMA_HIGH_CONST_ERR = " << high_std_err << ";\n";
    out << "const double CHI2_NDF_STD = " << chi2ndf_std << ";\n\n";
    out << "// Omega width ratio (to be filled after running omega_resolution.C)\n";
    out << "// const double R_OMEGA = ...; // sigma_data_omega / sigma_mc_omega\n";
    out.close();
    std::cout << "Parameters written to ../trackmass_scan/mc_resolution_params.h" << std::endl;

    // ----------------------------------------------------------------------
    // Draw and save plot
    // ----------------------------------------------------------------------
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TCanvas *c = new TCanvas("c", "MC Resolution Fit", 900, 700);
    gPad->SetBottomMargin(0.12);
    gPad->SetLeftMargin(0.15);

    g_origin->SetMarkerStyle(20);
    g_origin->SetMarkerSize(1.2);
    g_origin->SetTitle("");
    g_origin->GetYaxis()->SetTitleSize(0.05);
    g_origin->GetYaxis()->SetTitleSize(0.05);
    g_origin->GetYaxis()->SetTitleOffset(.8);
    g_origin->GetXaxis()->SetTitle("M_{#pi#pi}^{true} [MeV/c^{2}]");
    g_origin->GetXaxis()->SetTitleSize(0.05);
    g_origin->GetXaxis()->SetTitleSize(0.05);
    g_origin->GetXaxis()->SetTitleOffset(1);
    g_origin->GetYaxis()->SetTitle("#sigma_{MC} [MeV/c^{2}]");
    g_origin->GetXaxis()->CenterTitle();
    g_origin->GetYaxis()->CenterTitle();
    g_origin->GetYaxis()->SetRangeUser(0., 2.2);
    g_origin->Draw("AP");

    func->SetLineColor(kRed);
    func->SetLineWidth(2);
    func->Draw("same");

    // Draw the standard logistic fit for comparison (dashed)
    logistic_std->SetLineColor(kGreen + 2);
    logistic_std->SetLineStyle(2);
    logistic_std->Draw("same");

    // Draw a horizontal line at the plateau value at Mref (h0)
    TF1 *plateauLine = new TF1("plateauLine", "[0]", 260, 620);
    plateauLine->SetParameter(0, h0);
    plateauLine->SetLineColor(kBlue);
    plateauLine->SetLineStyle(2);
    plateauLine->Draw("same");

    TLegend *leg = new TLegend(0.2, 0.70, 0.45, 0.90);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(g, "MC resolution", "lep");
    leg->AddEntry(func, "Logistic + linear plateau", "l");
    leg->AddEntry(logistic_std, "Standard logistic", "l");
    leg->AddEntry(plateauLine, Form("h0 (at M_ref=%.0f) = %.3f MeV/c^{2}", Mref, h0), "l");
    leg->Draw();

    TPaveText *pt = new TPaveText(0.45, 0.15, 0.7, 0.45, "NDC");
    pt->SetFillColor(0);
    pt->SetBorderSize(0);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.03);
    pt->AddText(Form("low  = %.3f #pm %.3f MeV/c^{2}", low, low_err));
    pt->AddText(Form("h0   = %.3f #pm %.3f MeV/c^{2}", h0, h0_err));
    pt->AddText(Form("h1   = %.3f #pm %.3f ", h1, h1_err));
    pt->AddText(Form("k    = %.3f #pm %.3f 1/MeV/c^{2}", k, k_err));
    pt->AddText(Form("M_{0}   = %.1f #pm %.1f MeV/c^{2}", M0, M0_err));
    pt->AddText(Form("#chi^{2}/NDF = %.2f", chi2ndf));
    pt->Draw();

    c->SaveAs("../trackmass_scan/mc_resolution_logistic_linear_fit.pdf");
    std::cout << "Plot saved to ../trackmass_scan/mc_resolution_logistic_linear_fit.pdf" << std::endl;

    f->Close();
    delete c;

    return 0;
}
