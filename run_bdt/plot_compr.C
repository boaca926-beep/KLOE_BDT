#include "../header_bdt/compr.h"
//#include "../hist.h"
#include "../header_plot/plot.h"
//#include "../fitfun.h"
#include "../header_method/method.h"

int plot_compr(){

  //gROOT->SetBatch(kTRUE);  
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  //gROOT->SetStyle("Plain");
  //gROOT->ForceStyle();
  
  TFile* intree = new TFile("../output_" + var_nm + "/hist_" + var_nm + ".root");
  
  TIter next_tree(intree -> GetListOfKeys());

  TString objnm_tree, classnm_tree;

  int i = 0;
  TKey * key;
  
  while ( (key = (TKey *) next_tree() ) ) {
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  TObjArray *Hlist = (TObjArray *)intree -> Get("Hlist");
  checkArray(Hlist);

  // Helper to retrieve histograms from the list
  auto getHist = [&](const char* name) -> TH1D* {
    TH1D* h = (TH1D*)Hlist->FindObject(name);
    if (!h) cerr << "WARNING: histogram " << name << " not found in list" << endl;
    return h;
  };

  TH1D *hist_data = getHist("hist_data");
  hist_data->Draw();
  
  /*
  TH1D *hist_data = (TH1D *) Hlist -> At(0);

  TH1D *hist_eeg_sc = (TH1D *) Hlist -> At(10);
  TH1D *hist_isr3pi_sc = (TH1D *) Hlist -> At(11);
  TH1D *hist_omegapi_sc = (TH1D *) Hlist -> At(12);
  TH1D *hist_etagam_sc = (TH1D *) Hlist -> At(13);
  TH1D *hist_ksl_sc = (TH1D *) Hlist -> At(14);
  TH1D *hist_mcrest_sc = (TH1D *) Hlist -> At(15);
  TH1D *hist_bkgsum_sc = (TH1D *) Hlist -> At(16);
  */
  
  /*
    [0]: hist_data
    [1]: hist_eeg
    [2]: hist_omegapi
    [3]: hist_ksl
    [4]: hist_kpm
    [5]: hist_rhopi
    [6]: hist_etagam
    [7]: hist_bkgrest
    [8]: hist_isr3pi
    [9]: hist_mcrest
    [10]: hist_eeg_sc
    [11]: hist_isr3pi_sc
    [12]: hist_omegapi_sc
    [13]: hist_etagam_sc
    [14]: hist_ksl_sc
    [15]: hist_mcrest_sc
    [16]: hist_bkgsum_sc
  */

  /*
  TH1D *hist_mcsum_sc = (TH1D*) hist_bkgsum_sc -> Clone();
  hist_mcsum_sc -> Add(hist_isr3pi_sc, 1.);
  hist_mcsum_sc -> SetName("hist_mcsum_sc");
  format_h(hist_mcsum_sc, 2, 2);

  // plot distri.
  const double binwidth = getbinwidth(hist_data);
  const double ymax = hist_data -> GetMaximum();
  const double residul_min = -25., residul_max = 25.;
  
  TH1D * hresidul = new TH1D("hresidul", "", binsize, var_min, var_max);
  TH1D * hresidul_distr = new TH1D("hresidul_distr", "", 200, residul_min, residul_max);

  double nb_data = 0., evnt_err = 0., Nb_data = 0., nb_mcsum = 0.;
  double residul = 0.;
  
  for (int j = 1; j <= binsize; j ++ ) {

    nb_data =  hist_data -> GetBinContent(j);
    Nb_data += nb_data;

    nb_mcsum = hist_mcsum_sc -> GetBinContent(j);  
    evnt_err = TMath::Sqrt(nb_data + nb_mcsum);
    
    // residual
    residul = (nb_data - nb_mcsum) / evnt_err;

    if (nb_data > 0. && nb_mcsum > 0.) {
      
     hresidul -> SetBinContent(j, residul);
     hresidul_distr -> Fill(residul);

    }

    //cout << j << ": nb_data = " << nb_data << " +/- " << evnt_err << ", mcsum = " << nb_mcsum << ", residul= " << residul << "\n";

    
  }

  //
  hist_data -> SetMarkerStyle(21);
  hist_data -> SetMarkerSize(0.7);

  //hist_data -> GetYaxis() -> SetNdivisions(505);
  //hist_data -> GetYaxis() -> SetLabelSize(20);
  //hist_data -> GetYaxis() -> SetTitleFont(43);
  //hist_data -> GetYaxis() -> SetLabelFont(43); // Absolute font size in pixel (precision 3)

  TPaveText *pt34 = new TPaveText(0.2, 0.75, 0.4, 0.8, "NDC");
  PteAttr(pt34); pt34 -> SetTextSize(0.1);
  pt34 -> AddText("(a)");
  //pt34 -> AddText("(b)");
  //pt34 -> AddText("(c)");
  //pt34 -> AddText("(d)");

  TCanvas *cv = new TCanvas("cv", " ", 700, 700);
  cv -> SetBottomMargin(0.15);//0.007
  cv -> SetLeftMargin(0.15);
  
  //TPad *p2 = new TPad("p2", "p2", 0., 0., 1., 0.25);
  //p2 -> Draw();
  //p2 -> SetBottomMargin(0.3);
  //p2 -> SetGrid();
  
  //TPad *p1 = new TPad("p1", "p1", 0., 0.24, 1., 1.);
  //p1 -> Draw();
  //p1 -> SetBottomMargin(0.03);//0.007
  //p1 -> cd();

  //hist_data -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  hist_data -> GetXaxis() -> SetNdivisions(505);
  hist_data -> GetXaxis() -> SetTitle(var_symb + " " + unit);
  hist_data -> GetXaxis() -> CenterTitle();
  hist_data -> GetXaxis() -> SetTitleSize(0.06);
  hist_data -> GetXaxis() -> SetTitleOffset(1.0);
  hist_data -> GetXaxis() -> SetLabelOffset(0.01);
  hist_data -> GetXaxis() -> SetLabelSize(0.05);//0.03
  hist_data -> GetXaxis() -> SetRangeUser(0., 43.); //chi2
  //hist_data -> GetXaxis() -> SetRangeUser(20., 140.); //pi0angle
  //hist_data -> GetXaxis() -> SetRangeUser(0.5, 1.);
  //hist_data -> GetXaxis() -> SetRangeUser(650., 950.); // 3pi mass omega region
  
  //hist_data -> GetYaxis() -> SetLabelFont(43); // Absolute front size in pixel (precision 3)
  //hist_data -> GetYaxis() -> SetLabelSize(0.05);
  hist_data -> GetYaxis() -> SetRangeUser(0.01, ymax * 1.5); 
  //hist_data -> GetYaxis() -> SetRangeUser(0.01, ymax * 1e2); // 3pi mass full region
  //hist_data -> GetYaxis() -> SetRangeUser(0., 5e4); // 3pi mass omega region  
  //hist_data -> GetYaxis() -> SetTitle(TString::Format("Entries/%0.2f", binwidth) + " " + unit);
  hist_data -> GetYaxis() -> SetTitle("Events");
  hist_data -> GetYaxis() -> CenterTitle();
  hist_data -> GetYaxis() -> SetTitleSize(0.06);
  hist_data -> GetYaxis() -> SetTitleOffset(1.2);
  
  //hist_data -> SetStats(0);

  format_h(hist_mcsum_sc, 1, 2);

  hist_isr3pi_sc -> SetLineStyle(2); // dashed line
  hist_etagam_sc -> SetLineStyle(3); // dotted line
  hist_ksl_sc -> SetLineStyle(4); // short-dotted-dashed line
  hist_omegapi_sc -> SetLineStyle(5); // dotted-dashed line
  hist_mcrest_sc -> SetLineStyle(6); // triple-dotted-dashed line
  hist_eeg_sc -> SetLineStyle(7); // long-dashed line
  
  hist_data -> Draw();
  hist_mcsum_sc -> Draw("SameHist");
  hist_isr3pi_sc -> Draw("SameHist");
  hist_omegapi_sc -> Draw("SameHist");
  hist_etagam_sc -> Draw("SameHist");
  hist_ksl_sc -> Draw("SameHist");
  hist_mcrest_sc -> Draw("SameHist");
  hist_eeg_sc -> Draw("SameHist");
  pt34 -> Draw("Same");
  
  //gPad->SetLogy();

  TLegend *legd_cv = new TLegend(0.65, 0.35, 0.9, 0.9);
  
  legd_cv -> SetTextFont(132);
  legd_cv -> SetFillStyle(0);
  legd_cv -> SetBorderSize(0);
  legd_cv -> SetNColumns(1);
  
  legd_cv -> AddEntry(hist_data, "Data", "lep");
  legd_cv -> AddEntry(hist_mcsum_sc, "MC sum", "l");
  legd_cv -> AddEntry(hist_isr3pi_sc, "#pi^{+}#pi^{-}#pi^{0}#gamma", "l");
  legd_cv -> AddEntry(hist_omegapi_sc, "#omega#pi^{0}", "l");
  //legd_cv -> AddEntry(hist_etagam_sc, "#varphi#rightarrow#eta#gamma#rightarrow#pi^{+}#pi^{-}#pi^{0}#gamma", "l");
  legd_cv -> AddEntry(hist_etagam_sc, "#eta#gamma", "l");
  legd_cv -> AddEntry(hist_ksl_sc, "K_{L}K_{S}", "l");
  legd_cv -> AddEntry(hist_eeg_sc, "e^{+}e^{-}#gamma", "l");
  legd_cv -> AddEntry(hist_mcrest_sc, "Others", "l");
  
  legd_cv -> Draw("Same");
  
  legtextsize(legd_cv, 0.05);
  
  // save
  
  cv -> SaveAs("./output_" + var_nm + "/cv_compr_" + var_nm + ".pdf");
  */
  
  return 0;
  
}
