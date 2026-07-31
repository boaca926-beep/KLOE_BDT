// ============================================================================
// fit_resolution_mc.C
//
// Read the MC resolution vs. true dipion mass from the output of trackmass_scan,
// fit with sqrt((A/M)^2 + B^2), and save the parameters.
// Usage:
//   .x fit_resolution_mc.C()
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

void fit_resolution_mc() {
    // ----------------------------------------------------------------------
    // Input file from trackmass_scan run on MC (TISR3PI_SIG_PEAK)
    // ----------------------------------------------------------------------
    TString input_file = "../trackmass_scan/pull_scan_TISR3PI_SIG_PEAK.root";
    TFile *f = TFile::Open(input_file);
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open " << input_file << std::endl;
        return;
    }

    TGraphErrors *g = (TGraphErrors*)f->Get("g_resolution_vs_M");
    if (!g) {
        std::cerr << "ERROR: TGraph 'g_resolution_vs_M' not found." << std::endl;
        f->Close();
        return;
    }

    for (int i = 0; i < g->GetN(); ++i) {
      double x, y;
      g->GetPoint(i, x, y);
      std::cout << "Point " << i << ": x = " << x << ", y = " << y 
		<< ", y_err = " << g->GetErrorY(i) << std::endl;
    }
    
    // ----------------------------------------------------------------------
    // Fit function: sqrt((A/M)^2 + B^2)
    // ----------------------------------------------------------------------
    //TF1 *func = new TF1("resol_func", "sqrt(([0]/x)^2 + [1]^2)", 260, 650);
    TF1 *func = new TF1("resol_func", "sqrt(([0]/x)^2 + [1]^2)", 350, 650);
    func->SetParameters(1000.0, 0.5);   // initial guesses
    func->SetParLimits(0, 0, 5000);     // A in MeV^2
    func->SetParLimits(1, 0, 5);        // B in MeV

    // Perform fit
    g->Fit(func, "RQS");

    double A = func->GetParameter(0);
    double A_err = func->GetParError(0);
    double B = func->GetParameter(1);
    double B_err = func->GetParError(1);
    double chi2 = func->GetChisquare();
    int ndf = func->GetNDF();
    double chi2ndf = (ndf > 0) ? chi2 / ndf : 0;

    std::cout << "\n=== Fit results ===" << std::endl;
    std::cout << "A = " << A << " +/- " << A_err << " MeV^2" << std::endl;
    std::cout << "B = " << B << " +/- " << B_err << " MeV" << std::endl;
    std::cout << "chi2/ndf = " << chi2ndf << std::endl;

    // ----------------------------------------------------------------------
    // Save parameters to a header file (for inclusion in analysis)
    // ----------------------------------------------------------------------
    std::ofstream out("../trackmass_scan/mc_resolution_params.h");
    out << "// MC resolution parameters from fit to sigma vs true mass\n";
    out << "// sigma(M) = sqrt((A/M)^2 + B^2)  (M in MeV)\n";
    out << "const double A_MC = " << A << ";\n";
    out << "const double A_MC_err = " << A_err << ";\n";
    out << "const double B_MC = " << B << ";\n";
    out << "const double B_MC_err = " << B_err << ";\n";
    out << "const double chi2_ndf_MC = " << chi2ndf << ";\n";
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

    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.2);
    g->SetTitle("");
    g->GetXaxis()->SetTitle("M_{#pi#pi}^{true} (MeV/c^{2})");
    g->GetYaxis()->SetTitle("#sigma (MeV)");
    g->GetXaxis()->CenterTitle();
    g->GetYaxis()->CenterTitle();
    g->Draw("AP");

    func->SetLineColor(kRed);
    func->SetLineWidth(2);
    func->Draw("same");

    // Legend
    TLegend *leg = new TLegend(0.15, 0.70, 0.45, 0.90);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(g, "MC resolution", "lep");
    leg->AddEntry(func, Form("Fit: #sqrt{(%.0f/M)^{2}+%.3f^{2}}", A, B), "l");
    leg->Draw();

    // Text box with parameters
    TPaveText *pt = new TPaveText(0.55, 0.70, 0.90, 0.90, "NDC");
    pt->SetFillColor(0);
    pt->SetBorderSize(0);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.04);
    pt->AddText(Form("A = %.1f #pm %.1f MeV^{2}", A, A_err));
    pt->AddText(Form("B = %.3f #pm %.3f MeV", B, B_err));
    pt->AddText(Form("#chi^{2}/NDF = %.2f", chi2ndf));
    pt->Draw();

    c->SaveAs("../trackmass_scan/mc_resolution_fit.pdf");
    std::cout << "Plot saved to ../trackmass_scan/mc_resolution_fit.pdf" << std::endl;

    f->Close();
    delete c;
}
