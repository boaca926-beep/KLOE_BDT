// fit_extra_smear.C
// Usage: root -l -q 'fit_extra_smear.C()'   (or with a base path argument)

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TLatex.h>
#include <iostream>
#include <vector>
#include <cmath>

// --------------------------------------------------------------
// Default base path – CHANGE THIS or pass as argument.
// --------------------------------------------------------------
const TString DEFAULT_BASE = "/home/bo/Desktop/SANITY-CHECK-PPIM-RESOLUTION/";

// --------------------------------------------------------------
void fit_extra_smear(const TString base_path = DEFAULT_BASE) {

    // List of directories to scan (add more as you run higher values)
    std::vector<TString> dir_names = {
        "bdt_tuning_post_fit_TDATA_chain_false0.000",
        "bdt_tuning_post_fit_TDATA_chain_false0.005",
        "bdt_tuning_post_fit_TDATA_chain_false0.010",
        "bdt_tuning_post_fit_TDATA_chain_false0.015",
        "bdt_tuning_post_fit_TDATA_chain_false0.020",
	"bdt_tuning_post_fit_TDATA_chain_false0.040",
	"bdt_tuning_post_fit_TDATA_chain_false0.060"
    };

    const double TARGET_RMS = 68.98;          // MeV
    const TString tree_name = "TISR3PI_SIG_PEAK";
    const TString fit_func = "pol2";         // or "pol1"
    const TString extra_branch_name = "Br_extra_smear";   // adjust to your branch name

    // --------------------------------------------------------------
    // Extend x‑axis for plotting to see the intersection.
    // Set this to a value larger than your highest scanned extra_smear.
    // --------------------------------------------------------------
    const double XMAX_PLOT = 0.10;   // now extended to 0.10

    std::vector<double> x, y, ex, ey;

    std::cout << "\n========================================" << std::endl;
    std::cout << "  Extracting RMS of Br_ppIM" << std::endl;
    std::cout << "  Base path: " << base_path << std::endl;
    std::cout << "========================================" << std::endl;

    for (size_t i = 0; i < dir_names.size(); ++i) {
        TString file_path = base_path + dir_names[i] + "/cut/tree_pre.root";
        std::cout << "Processing: " << file_path << std::endl;

        // Try to get extra_smear from the tree
        double extra_val = -1;
        TFile *f = TFile::Open(file_path);
        if (f && !f->IsZombie()) {
            TTree *tree = (TTree*)f->Get(tree_name);
            if (tree) {
                TBranch *br = tree->GetBranch(extra_branch_name);
                if (br) {
                    tree->SetBranchAddress(extra_branch_name, &extra_val);
                    tree->GetEntry(0);
                }
            }
            f->Close();
        }

        // Fallback: extract from directory name
        if (extra_val < 0) {
            std::cerr << "  WARNING: Could not read branch '" << extra_branch_name
                      << "'. Falling back to directory name parsing." << std::endl;
            TString tmp = dir_names[i];
            Int_t pos = tmp.Last('_');
            if (pos != kNPOS) {
                TString valStr = tmp(pos+1, tmp.Length() - pos - 1);
                extra_val = valStr.Atof();
            }
        }

        if (extra_val < 0) {
            std::cerr << "  ERROR: Could not determine extra_smear for " << dir_names[i] << std::endl;
            continue;
        }

        f = TFile::Open(file_path);
        if (!f || f->IsZombie()) {
            std::cerr << "  ERROR: Cannot open " << file_path << std::endl;
            continue;
        }

        TTree *tree = (TTree*)f->Get(tree_name);
        if (!tree) {
            std::cerr << "  ERROR: Cannot find tree " << tree_name << std::endl;
            f->Close();
            continue;
        }

        // Compute RMS of Br_ppIM
        TH1D *h = new TH1D("h_tmp", "", 200, 200, 650);
        tree->Draw("Br_ppIM >> h_tmp");
        double rms = h->GetRMS();
        double entries = h->GetEntries();

        if (entries < 100)
            std::cerr << "  WARNING: Low statistics (" << entries << " entries)" << std::endl;

        double rms_err = (entries > 0) ? rms / sqrt(2.0 * entries) : 0.0;

        x.push_back(extra_val);
        y.push_back(rms);
        ex.push_back(0.0);
        ey.push_back(rms_err);

        std::cout << "  extra_smear = " << extra_val
                  << "  RMS = " << rms << " +/- " << rms_err
                  << "  Entries = " << entries << std::endl;

        delete h;
        f->Close();
    }

    if (x.size() < 3) {
        std::cerr << "ERROR: Not enough valid points to fit." << std::endl;
        return;
    }

    // Build graph of data points
    TGraphErrors *gr = new TGraphErrors(x.size(), &x[0], &y[0], &ex[0], &ey[0]);
    gr->SetTitle("RMS vs extra_smear;extra_smear;RMS of Br_ppIM [MeV]");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.2);

    TGraphErrors *gr_orig = dynamic_cast<TGraphErrors*>(gr->Clone());

    double y_min = 0.9 * TMath::MinElement(y.size(), &y[0]);
    double y_max = 1.10 * TMath::Max(TMath::MaxElement(y.size(), &y[0]), TARGET_RMS);
    gr->GetYaxis()->SetRangeUser(y_min, y_max);

    // Fit over the data range
    double xmin_fit = x.front() - 0.005;
    double xmax_fit = x.back() + 0.005;
    TF1 *fit = new TF1("fit", fit_func, xmin_fit, xmax_fit);
    fit->SetLineColor(kRed);
    fit->SetLineWidth(2);

    gr->Fit(fit, "RQ");

    double chi2 = fit->GetChisquare();
    int ndf = fit->GetNDF();
    double chi2ndf = (ndf > 0) ? chi2 / ndf : 0;

    std::cout << "\n========================================" << std::endl;
    std::cout << "  Fit result" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Function: " << fit_func << std::endl;
    std::cout << "Chi2/NDF = " << chi2ndf << std::endl;
    for (int i = 0; i < fit->GetNpar(); ++i) {
        std::cout << "p" << i << " = " << fit->GetParameter(i)
                  << " +/- " << fit->GetParError(i) << std::endl;
    }

    // --------------------------------------------------------------
    // Create a TGraph with many prediction points for the full range
    // --------------------------------------------------------------
    const int Npoints = 500;
    double *x_pred = new double[Npoints];
    double *y_pred = new double[Npoints];
    double x_start = 0.0;   // start at 0
    double x_end = XMAX_PLOT;
    for (int i = 0; i < Npoints; ++i) {
        x_pred[i] = x_start + (x_end - x_start) * i / (Npoints - 1);
        y_pred[i] = fit->Eval(x_pred[i]);
    }
    TGraph *gr_pred = new TGraph(Npoints, x_pred, y_pred);
    gr_pred->SetLineColor(kRed);
    gr_pred->SetLineWidth(2);
    gr_pred->SetMarkerStyle(1);      // small circle
    gr_pred->SetMarkerSize(0.5);
    gr_pred->SetMarkerColor(kRed);

    // --------------------------------------------------------------
    // Find intersection with target over extended range
    // --------------------------------------------------------------
    const int Nscan = 1000;
    double xmin_scan = 0.0;
    double xmax_scan = XMAX_PLOT;
    double best_x = xmin_scan;
    double best_diff = 1e9;

    for (int i = 0; i < Nscan; ++i) {
        double x_test = xmin_scan + (xmax_scan - xmin_scan) * i / (Nscan - 1);
        double y_test = fit->Eval(x_test);
        double diff = std::abs(y_test - TARGET_RMS);
        if (diff < best_diff) {
            best_diff = diff;
            best_x = x_test;
        }
    }

    bool intersection_found = (best_x > xmin_scan && best_x < xmax_scan);

    std::cout << "\n----------------------------------------" << std::endl;
    std::cout << "Target RMS = " << TARGET_RMS << " MeV" << std::endl;
    if (intersection_found) {
        std::cout << "Optimal extra_smear = " << best_x << std::endl;
        std::cout << "Predicted RMS at this value = " << fit->Eval(best_x) << " MeV" << std::endl;
    } else {
        std::cout << "WARNING: Intersection not found within scanned range." << std::endl;
        std::cout << "At extra_smear = " << xmax_scan << ", RMS = " << fit->Eval(xmax_scan) << " MeV" << std::endl;
        std::cout << "Consider scanning higher extra_smear values or increasing XMAX_PLOT." << std::endl;
    }
    std::cout << "----------------------------------------" << std::endl;

    // --------------------------------------------------------------
    // Draw
    // --------------------------------------------------------------
    TCanvas *c = new TCanvas("c", "RMS vs extra_smear", 800, 600);

    // Set x-axis range
    //gr_pred->GetXaxis()->SetRangeUser(0.0, XMAX_PLOT);
    gr_pred->GetYaxis()->SetRangeUser(66.0, 72.0);

    // Draw the full prediction curve with markers
    gr_pred->Draw("AP");   // Line and Points

    // Draw data points first
    gr_orig->Draw("P same");

    
    // Horizontal target line
    TLine *line = new TLine(0.0, TARGET_RMS, XMAX_PLOT, TARGET_RMS);
    line->SetLineColor(kBlue);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");

    if (intersection_found && best_x >= 0.0 && best_x <= XMAX_PLOT) {
        // Vertical line at optimal extra_smear
        TLine *vline = new TLine(best_x, y_min, best_x, y_max);
        vline->SetLineColor(kGreen+2);
        vline->SetLineStyle(3);
        vline->SetLineWidth(2);
        vline->Draw("same");

        // Marker at intersection
        TMarker *marker = new TMarker(best_x, TARGET_RMS, 20);
        marker->SetMarkerColor(kRed);
        marker->SetMarkerSize(2);
        marker->Draw("same");

        // Label at intersection
        TLatex *label = new TLatex(best_x + 0.002, TARGET_RMS + 1.0,
                                   Form("(%.4f, %.0f)", best_x, TARGET_RMS));
        label->SetTextColor(kRed);
        label->SetTextSize(0.04);
        label->Draw("same");
    } else {
        // Message when intersection not found
        TLatex *msg = new TLatex(0.5, 0.8,
                                 Form("Intersection not reached; RMS at %.3f = %.2f MeV", XMAX_PLOT, fit->Eval(XMAX_PLOT)));
        msg->SetNDC(kTRUE);
        msg->SetTextSize(0.05);
        msg->Draw("same");
    }

    // Label for target line
    TLatex *lineLabel = new TLatex(XMAX_PLOT - 0.002, TARGET_RMS + 1.0,
                                   Form("Target RMS = %.0f MeV", TARGET_RMS));
    lineLabel->SetTextColor(kBlue);
    lineLabel->SetTextSize(0.04);
    lineLabel->Draw("same");

    // Legend
    TLegend *leg = new TLegend(0.15, 0.7, 0.45, 0.9);
    leg->AddEntry(gr, "Data points", "lep");
    leg->AddEntry(gr_pred, Form("Prediction (%s fit)", fit_func.Data()), "lp");
    leg->AddEntry(line, "Target RMS", "l");
    if (intersection_found)
        leg->AddEntry((TObject*)0, Form("Optimal = %.4f", best_x), "");
    leg->Draw();

    
    c->Update();
    c->SaveAs("extra_smear_fit.png");
    c->SaveAs("extra_smear_fit.pdf");

    std::cout << "\nPlots saved as extra_smear_fit.png/pdf" << std::endl;

    // Clean up
    delete[] x_pred;
    delete[] y_pred;

    // --------------------------------------------------------------
    // Recommendation
    // --------------------------------------------------------------
    std::cout << "\n========================================" << std::endl;
    std::cout << "  RECOMMENDATION" << std::endl;
    std::cout << "========================================" << std::endl;
    if (intersection_found) {
        std::cout << "Use extra_smear = " << best_x << std::endl;
        std::cout << "This should give RMS ≈ " << TARGET_RMS << " MeV." << std::endl;
    } else {
        std::cout << "The target RMS is not reached within the scanned range." << std::endl;
        std::cout << "You need to run additional samples with extra_smear > 0.020." << std::endl;
        std::cout << "Suggested next values: 0.025, 0.030, 0.035, 0.040, 0.045, 0.050, ..." << std::endl;
        std::cout << "Then re-run this macro with XMAX_PLOT set accordingly." << std::endl;
    }
    std::cout << "========================================\n" << std::endl;
}
