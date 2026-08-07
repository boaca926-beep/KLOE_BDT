// ============================================================================
// fit_trk_bias_mc_linear.C
//
// Fit the MC bias (mean of residual distribution) vs. true dipion mass
// with a linear function: bias(M) = p0 + p1*M
// Excludes points with large uncertainties or failed fits.
// Usage:
//   .x fit_trk_bias_mc_linear.C()
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

int fit_trk_bias_mc(const TString input_file="../trackmass_scan/pull_scan_TISR3PI_SIG_PEAK.root") {

    // ----------------------------------------------------------------------
    // Input file (output from trackmass_scan on TISR3PI_SIG_PEAK)
    // ----------------------------------------------------------------------
    //const TString input_file = "../trackmass_scan/" + root_file + ".root";
    TFile *f = TFile::Open(input_file);
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open " << input_file << std::endl;
        return 1;
    }

    // Get the bias TGraph (g_bias_vs_M)
    TGraphErrors *g_orig = (TGraphErrors*)f->Get("g_bias_vs_M");
    if (!g_orig) {
        std::cerr << "ERROR: TGraph 'g_bias_vs_M' not found." << std::endl;
        f->Close();
        return 1;
    }

    // Also get the resolution TGraph to apply quality cuts
    TGraphErrors *g_sigma = (TGraphErrors*)f->Get("g_resolution_vs_M");
    if (!g_sigma) {
        std::cerr << "ERROR: TGraph 'g_resolution_vs_M' not found." << std::endl;
        f->Close();
        return 1;
    }

    // ----------------------------------------------------------------------
    // Apply quality cut to remove bins with unrealistic sigma or errors
    // ----------------------------------------------------------------------
    std::vector<double> x_good, y_good, ex_good, ey_good;
    for (int i = 0; i < g_orig->GetN(); ++i) {
        double x, y, ey;
        g_orig->GetPoint(i, x, y);
        ey = g_orig->GetErrorY(i);

        // Get resolution for this bin
        double sigma, sigma_err;
        g_sigma->GetPoint(i, x, sigma);
        sigma_err = g_sigma->GetErrorY(i);

        // Quality criteria:
        // - resolution error < 0.2 (fractional error < 20%)
        // - bias error < 0.1 (fractional error < 50% of bias value)
        // - sigma > 0 and sigma < 5 (physically reasonable)
        if (sigma > 0.1 && sigma < 5.0 && sigma_err > 0 && sigma_err / sigma < 0.2 &&
            ey > 0 && ey / y < 0.5) {
            x_good.push_back(x);
            y_good.push_back(y);
            ey_good.push_back(ey);
            ex_good.push_back(0);
        }
    }

    int n_good = x_good.size();
    if (n_good < 5) {
        std::cerr << "ERROR: Too few good points after quality cut." << std::endl;
        f->Close();
        return 1;
    }

    TGraphErrors *g = new TGraphErrors(n_good, x_good.data(), y_good.data(),
                                       ex_good.data(), ey_good.data());
    g->SetName("g_bias_good");

    std::cout << "Using " << n_good << " good points out of " << g_orig->GetN() << std::endl;

    // ----------------------------------------------------------------------
    // Fit with linear: bias(M) = p0 + p1*M
    // ----------------------------------------------------------------------
    const double fit_min = 300;
    const double fit_max = 600;

    cout << "Fit range = [" << fit_min << ", " << fit_max << "] MeV/c^{2}" << endl;
    
    TF1 *func = new TF1("bias_linear", "[0] + [1]*x", 300, 600);
    // Initial guesses: based on data (bias ~0 at 300, ~0.7 at 600)
    double p0 = -0.5;
    double p1 = 0.002;
    func->SetParameters(p0, p1);
    func->SetParLimits(0, -2.0, 2.0);
    func->SetParLimits(1, -0.01, 0.01);

    // Perform fit
    g->Fit(func, "RQS");

    double chi2 = func->GetChisquare();
    int ndf = func->GetNDF();
    double chi2ndf = (ndf > 0) ? chi2 / ndf : 0;

    double p0_fit = func->GetParameter(0);
    double p0_err = func->GetParError(0);
    double p1_fit = func->GetParameter(1);
    double p1_err = func->GetParError(1);

    std::cout << "\n=== Linear bias fit results ===" << std::endl;
    std::cout << "p0 = " << p0_fit << " +/- " << p0_err << " MeV" << std::endl;
    std::cout << "p1 = " << p1_fit << " +/- " << p1_err << " MeV/MeV" << std::endl;
    std::cout << "chi2/ndf = " << chi2ndf << std::endl;

    // ----------------------------------------------------------------------
    // Save parameters to header file
    // ----------------------------------------------------------------------
    std::ofstream out("../trackmass_scan/mc_bias_params.h");
    out << "// MC bias parameters from linear fit\n";
    out << "// bias(M) = p0 + p1*M   (M in MeV)\n";
    out << "const double BIAS_P0 = " << p0_fit << ";\n";
    out << "const double BIAS_P0_ERR = " << p0_err << ";\n";
    out << "const double BIAS_P1 = " << p1_fit << ";\n";
    out << "const double BIAS_P1_ERR = " << p1_err << ";\n";
    out << "const double BIAS_CHI2_NDF = " << chi2ndf << ";\n";
    out.close();
    std::cout << "Parameters written to ../trackmass_scan/mc_bias_params.h" << std::endl;

    // ----------------------------------------------------------------------
    // Draw and save plot
    // ----------------------------------------------------------------------
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    TCanvas *c = new TCanvas("c", "MC Bias Fit (Linear)", 900, 700);
    gPad->SetBottomMargin(0.12);
    gPad->SetLeftMargin(0.15);

    g_orig->SetMarkerStyle(20);
    g_orig->SetMarkerSize(1.2);
    g_orig->SetTitle("");
    g_orig->GetYaxis()->SetTitleSize(0.05);
    g_orig->GetYaxis()->SetTitleSize(0.05);
    g_orig->GetYaxis()->SetTitleOffset(.8);
    g_orig->GetYaxis()->SetTitle("Bias (Mean residual) [MeV/c^{2}]");
    g_orig->GetXaxis()->SetTitle("M_{2#pi}^{true} [MeV/c^{2}]");
    g_orig->GetXaxis()->SetTitleSize(0.05);
    g_orig->GetXaxis()->SetTitleSize(0.05);
    g_orig->GetXaxis()->SetTitleOffset(1);
    g_orig->GetXaxis()->CenterTitle();
    g_orig->GetYaxis()->CenterTitle();
    //g_orig->GetYaxis()->SetRangeUser(0.01, 2.);
    g_orig->Draw("AP");

    // Draw origal points in grey (for reference)
    g_orig->SetMarkerStyle(20);
    g_orig->SetMarkerColor(kGray);
    g_orig->Draw("AP");

    
    // Draw good points in blue
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlue);
    g->Draw("P SAME");

    // Draw fit
    func->SetLineColor(kRed);
    func->SetLineWidth(2);
    func->Draw("same");

    TLegend *leg = new TLegend(0.2, 0.70, 0.5, 0.90);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    //leg->AddEntry(g_orig, "All points", "lep");
    leg->AddEntry(g, "Good points", "lep");
    leg->AddEntry(func, Form("Linear: bias = %.3f + %.3f*M", p0_fit, p1_fit), "l");
    leg->Draw();

    TPaveText *pt = new TPaveText(0.55, 0.15, 0.8, 0.4, "NDC");
    pt->SetFillColor(0);
    pt->SetBorderSize(0);
    pt->SetTextAlign(12);
    pt->SetTextSize(0.03);
    pt->AddText(Form("a = %.4f #pm %.4f MeV/c^{2}", p0_fit, p0_err));
    pt->AddText(Form("b = %.4f #pm %.4f ", p1_fit, p1_err));
    pt->AddText(Form("#chi^{2}/NDF = %.2f", chi2ndf));
    pt->Draw();

    c->SaveAs("../trackmass_scan/mc_bias_linear_fit.pdf");
    std::cout << "Plot saved to ../trackmass_scan/mc_bias_linear_fit.pdf" << std::endl;

    f->Close();
    delete c;

    return 0;
}
