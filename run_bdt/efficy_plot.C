#include "../header_plot/plot.h"


const TString TYPE_LIST[2] = {"bdt", "kloe"};
//const TString HISR_LIST[2] = {"HSIG", ""};

//const TString infile = "/home/kloe/Desktop/input_bdt_TDATA_chain/hist/hist.root";
//const TString infile = "/home/kloe/Desktop/input_kloe_TDATA_norm/hist/hist.root";

const TString sample_type = "kloe";
const int color_indx = 4;
const int mstyle_indx = 22;
const TString hefficy_type = "hefficy";
const TString cv_tit = "efficy";
const TString hgen_tit = "N^{gen}_{3#pi}";
const TString htrue_tit = "N^{true}_{3#pi}";
const TString ratio_tit = "N^{true}_{3#pi}/N^{gen}_{3#pi}";

//const TString hgen_tit = "N^{gen}_{3#pi,phok5}";
//const TString htrue_tit = "N^{true,sel}_{3#pi,phok5}";
//const TString ratio_tit = "N^{true,sel}_{3#pi,phok5}/N^{gen}_{3#pi,phok5}";
const TString y_tit = "#varepsilon_{3#pi}";
//const TString efficy_path = "efficy_output";
const double M3pi_min = 740;
const double M3pi_max = 870;

TCanvas *plot_efficy(TH1D *h1d, TH1D* h1d_1, TString cv_nm, TString cv_title) {

  
  TCanvas *cv = new TCanvas(cv_nm, cv_title, 1100, 600);

  cv -> SetBottomMargin(0.15);
  cv -> SetLeftMargin(0.15);
  cv -> SetRightMargin(0.15);

  h1d -> GetYaxis() -> SetTitle(TString::Format("Efficiency/[%0.2f MeV/c^{2}]", getbinwidth(h1d)));
  h1d -> GetYaxis() -> CenterTitle();
  //h1d -> GetYaxis() -> SetLabelSize(0.05);
  h1d -> GetYaxis() -> SetTitleSize(0.08);
  h1d -> GetYaxis() -> SetTitleFont(132);
  h1d -> GetYaxis() -> SetTitleOffset(.9);
  h1d -> GetYaxis() -> SetRangeUser(0.01, 0.05 * 1.2); 
      
  h1d -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  h1d -> GetXaxis() -> CenterTitle();
  //h1d -> GetXaxis() -> SetLabelSize(0.05);
  h1d -> GetXaxis() -> SetTitleSize(0.08); //cout << "here" << endl;
  h1d -> GetXaxis() -> SetTitleFont(132);
  h1d -> GetXaxis() -> SetRangeUser(740., 820.);
  h1d -> GetXaxis() -> SetTitleOffset(.8);

  h1d -> Draw();
  //h1d_1 -> Draw("Same");
  

  TLegend *legd_cv = new TLegend(0.4, 0.75, 0.8, 0.9);
  
  //legd_cv -> SetTextFont(42);
  legd_cv -> SetFillStyle(0);
  legd_cv -> SetBorderSize(0);
  legd_cv -> SetNColumns(2);

  //legd_cv -> AddEntry(h1d, "#varepsilon^{signal}_{3#pi}", "lep"); //ratio_tit
  //legd_cv -> AddEntry(h1d_1, "#varepsilon^{#rho#pi}_{3#pi}", "lep"); //ratio_tit
  //legd_cv -> Draw("Same");
  
  legtextsize(legd_cv, 0.1);

  return cv;
  
}

//
double binomial_err(double nb_true, double nb_gen) {
  double error = 0.;
  double ratio = 0.; 

  if (nb_gen != 0.) {
    ratio = nb_true / nb_gen;
    error = TMath::Sqrt(ratio * (1. - ratio) / nb_gen);
  }
   
  //cout << "true = " << nb_true << ", gen = " << nb_gen << ", ratio = " << ratio << ", error = " << error << endl;

  return error;
}

//
TH1D *get_efficy(TH1D *hsig_true, TH1D *hsig_gen) {



  int binsize = hsig_true -> GetNbinsX();
  
  double hmin = hsig_true -> GetXaxis() -> GetXmin();
  double hmax = hsig_true -> GetXaxis() -> GetXmax();

  //cout << "binsize = " << binsize << ", hmin = " << hmin << ", hmax = " << hmax << endl;
  
  TH1D * hefficy = new TH1D("hefficy", "", binsize, hmin, hmax);
  hefficy -> Sumw2();

  double nb_gen = 0., nb_true = 0.;
  double efficy = 0., efficy_err = 0.;
  
  for (int i = 1; i <= binsize; i ++ ) {

    nb_gen = hsig_gen -> GetBinContent(i);
    
    nb_true = hsig_true -> GetBinContent(i);
    
    efficy = nb_true / nb_gen;

    efficy_err = binomial_err(nb_true, nb_gen);
    
    if (nb_true == 0. || nb_gen == 0.) {

      efficy = 0.;
      
      efficy_err = 0.;
      
    }

    hefficy -> SetBinContent(i, efficy);
    hefficy -> SetBinError(i, efficy_err);

    //cout << "nb_true = " << nb_true << ", nb_gen = " << nb_gen << ", efficy = " << hefficy -> GetBinContent(i) << "+/-" << hefficy -> GetBinError(i) << ", efficy_err = " << efficy_err << endl;
      
    
  }

  return hefficy;
}


int efficy_plot() {

  //gROOT->SetBatch(kTRUE);
  gErrorIgnoreLevel = kError;
  TGaxis::SetMaxDigits(4);
  gStyle->SetOptStat(0);

  TString infile = "/home/kloe/Desktop/input_" +  sample_type + "_TDATA_chain/hist/hist.root";
  
  if (sample_type == TYPE_LIST[0]){
    cout << "Reading " << sample_type << " at " << infile << endl;
  }
  else if (sample_type == TYPE_LIST[1]) {
    cout << "Reading " << sample_type << " at " << infile << endl;
  }
  else {
    cout << "Wraong sample type!" << endl;
  }
  

  TKey *key;

  // get efficiency histos
  cout << "Get efficiency histos." << endl;
  //TFile *intree_histo = new TFile("./efficy_compr_output.root");

  //TH1D* hefficy_sig_apprx = (TH1D*)intree_histo -> Get("hefficy_sig_apprx");
  //format_h(hefficy_sig_apprx, 4, 2);// color 45
  //hefficy_sig_apprx -> SetMarkerStyle(21);
  //hefficy_sig_apprx -> SetMarkerSize(0.5);
  //hefficy_sig_apprx -> SetMarkerColor(4);

  // plot efficiencies
  cout << "Plot efficiencies from " << infile << endl;

  TFile *intree = new TFile(infile);
  
  TIter next_tree(intree -> GetListOfKeys());

  TString objnm_tree, classnm_tree;
  
  int i = 0;
  
  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

  TList *HSIG = (TList *) intree -> Get("HSIG");

  TH1D* htrue = nullptr;

  if (sample_type == "kloe"){
    htrue = (TH1D*)HSIG -> FindObject("hsig_true");
  }
  else if (sample_type == "bdt"){
    TList *HPeakNonReson = (TList *) intree -> Get("HPeakNonReson");
    htrue = (TH1D*)HPeakNonReson -> FindObject("h1d_IM3pi_TISR3PI_SIG_PEAK_TRUE");
  }
  format_h(htrue, color_indx, 2);
  htrue -> SetLineStyle(2);

  TH1D* hgen = (TH1D*)HSIG -> FindObject("hsig_gen");
  format_h(hgen, color_indx, 2);

  TH1D* hefficy = get_efficy(htrue, hgen);
  format_h(hefficy, color_indx, 2);

  //// Plot
  cout << "Plot efficy ... " << endl;

  //
  double binwidth = getbinwidth(htrue);
  const double xrange0 = 760., xrange1 = 820.;
  
  TCanvas *cv_ratio = new TCanvas("cv_ratio", cv_tit, 0, 0, 1000, 700);

  TPad *p2 = new TPad("p2", "p2", 0., 0., 1., 0.49);
  p2 -> Draw();
  p2 -> SetBottomMargin(0.2);
  p2 -> SetLeftMargin(0.1);
  //p2 -> SetGrid();
  
  TPad *p1 = new TPad("p1", "p1", 0., 0.47, 1., 1.);
  p1 -> Draw();
  p1 -> SetLeftMargin(0.1);
  p1 -> SetBottomMargin(0.03);//0.007
  p1 -> cd();

  //
  hgen -> SetMarkerStyle(mstyle_indx);
  hgen -> SetMarkerSize(0.7);
  hgen -> GetYaxis() -> SetTitle(TString::Format("Events/[%0.2f", binwidth) + " MeV/c^{2}]");
  hgen -> GetYaxis() -> CenterTitle();
  hgen -> GetYaxis() -> SetRangeUser(0.1, hgen -> GetMaximum() * 15);
  //hgen -> GetYaxis() -> SetRangeUser(0.1, 1200.);
  hgen -> GetYaxis() -> SetNdivisions(505);
  hgen -> GetYaxis() -> SetTitleSize(0.1);
  hgen -> GetYaxis() -> SetLabelSize(25);
  //hgen -> GetYaxis() -> SetTitleFont(43);
  hgen -> GetYaxis() -> SetTitleOffset(0.5);
  hgen -> GetYaxis() -> SetLabelFont(43); // Absolute font size in pixel (precision 3)

  //hgen -> GetXaxis() -> SetTitle();
  hgen -> GetXaxis() -> SetRangeUser(xrange0, xrange1);
  hgen -> GetXaxis() -> SetTitleOffset(1.2);
  //hgen -> GetXaxis() -> SetLabelSize(20);
  hgen -> GetXaxis() -> SetLabelOffset(15);
  hgen -> GetXaxis() -> CenterTitle();

  hgen -> SetStats(0);

  hgen -> Draw();
  htrue -> Draw("Same");
  
  gPad -> SetLogy();

  TLegend *legd_cv = new TLegend(0.15, 0.7, 0.5, 0.9);
  
  //legd_cv -> SetTextFont(132);
  legd_cv -> SetFillStyle(0);
  legd_cv -> SetBorderSize(0);
  legd_cv -> SetNColumns(2);

  legd_cv -> AddEntry(hgen, hgen_tit, "l");
  legd_cv -> AddEntry(htrue, htrue_tit, "l");  
  //legd_cv -> Draw("Same");
  
  legtextsize(legd_cv, 0.1);

  //
  hefficy -> SetMarkerStyle(mstyle_indx);
  hefficy -> SetMarkerSize(0.7);
  hefficy -> SetMarkerColor(1);
  hefficy -> SetLineColor(1);
  
  //hefficy -> GetYaxis() -> SetRangeUser(0.1, hefficy -> GetMaximum() * 1.5);
  hefficy -> GetYaxis() -> SetRangeUser(0., 0.05);
  hefficy -> GetYaxis() -> SetNdivisions(505);
  hefficy -> GetYaxis() -> SetTitleSize(0.1);
  //hefficy -> GetYaxis() -> SetTitleFont(43);
  hefficy -> GetYaxis() -> SetTitleOffset(0.5);
  //hefficy -> GetYaxis() -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hefficy -> GetYaxis() -> SetLabelSize(0.07);
  hefficy -> GetYaxis() -> CenterTitle();
  hefficy -> GetYaxis() -> SetTitle("Efficiency"); // y_tit

  //hefficy -> GetXaxis() -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hefficy -> GetXaxis() -> SetTitleOffset(.8);
  hefficy -> GetXaxis() -> SetLabelSize(0.08);
  hefficy -> GetXaxis() -> SetRangeUser(xrange0, xrange1);
  //hefficy -> GetXaxis() -> SetLabelOffset(0.03);
  hefficy -> GetXaxis() -> SetTitleSize(0.1);
  hefficy -> GetXaxis() -> CenterTitle();
  hefficy -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  
  hefficy -> SetStats(0);      // No statistics on lower plot
  
  //format_h(hefficy, color_indx, 2);

  p2 -> cd();

  hefficy -> Draw();
  //hefficy_sig_apprx -> Draw("Same");
    
  TLegend *legd_cv1 = new TLegend(0.6, 0.7, 0.8, 0.9);
  
  //legd_cv1 -> SetTextFont(42);
  legd_cv1 -> SetFillStyle(0);
  legd_cv1 -> SetBorderSize(0);
  legd_cv1 -> SetNColumns(2);

  //legd_cv1 -> AddEntry(hefficy, "#varepsilon_{3#pi}", "lep"); //ratio_tit
  //legd_cv1 -> Draw("Same");
  
  legtextsize(legd_cv1, 0.1);
  
  // save
  cv_ratio->cd();
  cv_ratio -> SaveAs("../plots_efficy/" + cv_tit + "_" + sample_type + ".pdf");

  // Clean up
  delete htrue;
  delete hgen;
  delete hefficy;
  delete cv_ratio;
  
  return 0;
  
}
