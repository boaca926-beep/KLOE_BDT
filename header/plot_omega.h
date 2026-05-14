const TString input_file = "../../input_bdt_TDATA_norm/omega_fit/omega_fit.root";
const TString file_type = "norm";

const TString outputPlot = "../plots_omega/";


TH1D *hsig; // smeared signal IM3pi MC true
TH1D *hdata;
TH1D *heeg;
TH1D *hksl;
TH1D *homegapi;
TH1D *hetagam;
TH1D *hmcrest;
TH1D *hbkgsum;
TH1D *hmcsum;

TH2D *hsmearmatr;
TH2D *hcorrmatrix;

TTree *TCRX3PI = nullptr;
TTree* TRESULT = new TTree("TRESULT", "recreate");

  
double Delta_m3pi;

TFile *f_input;

TGraphErrors *gf_n3pi_obs;
TGraphErrors *gf_n3pi_fit;
TGraphErrors *gf_n3pi_diff;

TCanvas *cv_compr;
TCanvas *cv_n3pi_fit;

const int list_size = 1000;

double M3PI_FIT[list_size], M3PI_FIT_ERR[list_size];
double N3PI_OBS_FIT[list_size], N3PI_OBS_FIT_ERR[list_size];
double N3PI_FIT[list_size], N3PI_FIT_ERR[list_size];
double N3PI_DIFF[list_size], N3PI_DIFF_ERR[list_size];
double RESIDUAL[list_size], RESIDUAL_ERR[list_size];
double ISRLUMI[list_size], CRX3PI_BW[list_size];
double EFFICY[list_size], EFFICY_ERR[list_size];
  
//TRandom *rnd=0;

double Mass_omega_fit = 0., Mass_omega_err_fit = 0.;
double Gam_omega_fit = 0., Gam_omega_err_fit = 0.;
double BB_fit = 0., BB_err_fit = 0.;
double chi2_sum = 0.;
double Lumi_int_fit = 0.;

double OMEGA_PARA[3], OMEGA_PARA_ERR[3];

int ndf = 0;

const double xrange0 = 760., xrange1 = 805.;

//
void plot_fit() {


  double residual = 0.;
  
  int bin_indx = 0;
  
  int fit_indx = 0;

  cout << TCRX3PI -> GetEntries() << endl;

  /*
  TRESULT -> Branch("Br_m3pi", &m3pi, "Br_m3pi/D");
  TRESULT -> Branch("Br_efficy", &efficy, "Br_efficy/D");
  TRESULT -> Branch("Br_efficy_err", &efficy_err, "Br_efficy_err/D");
  TRESULT -> Branch("Br_isrlumi_apprx", &isrlumi_apprx, "Br_isrlumi_apprx/D");
  TRESULT -> Branch("Br_nb_isr3pi_obs", &nb_isr3pi_obs, "Br_nb_isr3pi_obs/D");
  TRESULT -> Branch("Br_nb_isr3pi_obs_err", &nb_isr3pi_obs_err, "Br_nb_isr3pi_obs_err/D");
  */

  for (Int_t irow = 0; irow < TCRX3PI -> GetEntries(); irow++) {// start fill
    
    TCRX3PI -> GetEntry(irow); //cout << irow << endl;

    //
    bin_indx = TCRX3PI -> GetLeaf("Br_bin_indx") -> GetValue(0);
    Lumi_int_fit = TCRX3PI -> GetLeaf("Br_Lumi_int") -> GetValue(0);
    
    chi2_sum = TCRX3PI -> GetLeaf("Br_chi2_sum_crx3pi") -> GetValue(0);
    ndf = TCRX3PI -> GetLeaf("Br_ndf") -> GetValue(0);

    //
    Mass_omega_fit = TCRX3PI -> GetLeaf("Br_OMEGA_PARA") -> GetValue(0);
    Mass_omega_err_fit = TCRX3PI -> GetLeaf("Br_OMEGA_PARA_ERR") -> GetValue(0);

    Gam_omega_fit = TCRX3PI -> GetLeaf("Br_OMEGA_PARA") -> GetValue(1);
    Gam_omega_err_fit = TCRX3PI -> GetLeaf("Br_OMEGA_PARA_ERR") -> GetValue(1);
    
    BB_fit = TCRX3PI -> GetLeaf("Br_OMEGA_PARA") -> GetValue(2);
    BB_err_fit = TCRX3PI -> GetLeaf("Br_OMEGA_PARA_ERR") -> GetValue(2);

    //
    M3PI_FIT[irow] = TCRX3PI -> GetLeaf("Br_M3PI") -> GetValue(bin_indx);
    M3PI_FIT_ERR[irow] = 0.;

    EFFICY[irow] = TCRX3PI -> GetLeaf("Br_EFFICY") -> GetValue(bin_indx);
    EFFICY_ERR[irow] = TCRX3PI -> GetLeaf("Br_EFFICY_ERR") -> GetValue(bin_indx);

    N3PI_OBS_FIT[irow] = TCRX3PI -> GetLeaf("Br_N3PI_OBS") -> GetValue(bin_indx);
    N3PI_OBS_FIT_ERR[irow] = TCRX3PI -> GetLeaf("Br_N3PI_OBS_ERR") -> GetValue(bin_indx);

    N3PI_FIT[irow] = TCRX3PI -> GetLeaf("Br_N3PI_SMEAR") -> GetValue(bin_indx);
    N3PI_FIT_ERR[irow] = TCRX3PI -> GetLeaf("Br_N3PI_SMEAR_ERR") -> GetValue(bin_indx);

    N3PI_DIFF[irow] = TCRX3PI -> GetLeaf("Br_N3PI_DIFF") -> GetValue(bin_indx);
    N3PI_DIFF_ERR[irow] = TCRX3PI -> GetLeaf("Br_N3PI_DIFF_ERR") -> GetValue(bin_indx);
    //N3PI_DIFF[irow] = TCRX3PI -> GetLeaf("Br_N3PI_DIFF") -> GetValue(bin_indx) / N3PI_DIFF_ERR[irow];
    //cout << N3PI_DIFF_ERR[irow] << endl;

    // residual
    residual = N3PI_DIFF[irow] / N3PI_DIFF_ERR[irow];
    //cout << residual << endl;
    
    RESIDUAL[irow] = residual;
    RESIDUAL_ERR[irow] = 0.;

    cout << irow << ", " << bin_indx << ", "<< M3PI_FIT[irow] << ", " << N3PI_OBS_FIT[irow] << "+/-" << N3PI_OBS_FIT_ERR[irow] << ", " << N3PI_FIT[irow] << "+/-" << N3PI_FIT_ERR[irow] << ", " << RESIDUAL[irow] << endl;
    
    //
    ISRLUMI[irow] = TCRX3PI -> GetLeaf("Br_ISRLUMI") -> GetValue(bin_indx);
    CRX3PI_BW[irow] = TCRX3PI -> GetLeaf("Br_CRX3PI_BW") -> GetValue(bin_indx);
    
    fit_indx ++;

    /*
    cout << irow + 1 << ": n3pi obs. = " << N3PI_OBS_FIT[irow] << "+/-" << N3PI_OBS_FIT_ERR[irow] << ", n3pi diff. = " << N3PI_DIFF[irow] << "+/-" << N3PI_DIFF_ERR[irow] << endl;

    cout << "list indx = " << irow + 1 << "\n"
	 << "\tbin_indx = " << bin_indx << "\n"
	 << "\tm3pi = " << M3PI_FIT[irow] << " +/- " << M3PI_FIT_ERR[irow] << "\n"
	 << "\tn3pi_obs = " << N3PI_OBS_FIT[irow] << " +/- " << N3PI_OBS_FIT_ERR[irow] << "\n"
	 << "\tn3pi_fit = " << N3PI_FIT[irow] << " +/- " << N3PI_FIT_ERR[irow] << "\n"
	 << "\tn3pi_diff = " << N3PI_DIFF[irow] << " +/- " << N3PI_DIFF_ERR[irow] << "\n"
	 << "\tisr_lumi = " << ISRLUMI[irow] << "\n"
	 << "\tcrx3pi_bw = " << CRX3PI_BW[irow] << "\n"
	 << "\tefficy = " << EFFICY[irow] * 100. << " +/- " << EFFICY_ERR[irow] * 100. << " [%]\n"
	 << "\chi2_sum = " << chi2_sum << "\n\n";
    */
    
  }

  OMEGA_PARA[0] = Mass_omega_fit;
  OMEGA_PARA[1] = Gam_omega_fit;
  OMEGA_PARA[2] = BB_fit;

  OMEGA_PARA_ERR[0] = Mass_omega_err_fit;
  OMEGA_PARA_ERR[1] = Gam_omega_err_fit;
  OMEGA_PARA_ERR[2] = BB_err_fit;
  
  //TGraph *gf_n3pi_diff = new TGraphErrors(fit_indx, M3PI_FIT, N3PI_DIFF, M3PI_FIT_ERR, N3PI_DIFF_ERR);
  gf_n3pi_diff = new TGraphErrors(fit_indx, M3PI_FIT, RESIDUAL, M3PI_FIT_ERR, RESIDUAL_ERR);
  SetGFAttr(gf_n3pi_diff, "M_{3#pi} [MeV/c^{2}]", "Residual", "gf_n3pi_diff");
  
  TH1D * hdiff = new TH1D("hresidul", "", 38, 720, 820);
  int nPoints = gf_n3pi_diff -> GetN();
  cout << "nPoints = " << nPoints << endl;
  gf_n3pi_diff -> RemovePoint(0);
  
  double diff_min = 999999.;
  double diff_max = 0.;
  for (int i = 0; i < nPoints; i++) {
    
    double x, y;
    gf_n3pi_diff -> GetPoint(i, x, y);
    //cout << "i = " << i << ", x = " << x << ", y = " << y << endl;
    
    if (y > diff_max) diff_max = y;
    
    if (y < diff_min) diff_min = y;
    
    hdiff -> Fill(i+1, y);
    
  }

  /*
  cout << "diff_max = " << diff_max << "\n"
       << "diff_min = " << diff_min << endl;
  */

  //cout << N3PI_OBS_FIT_ERR[5] << endl;
  gf_n3pi_obs = new TGraphErrors(fit_indx, M3PI_FIT, N3PI_OBS_FIT, M3PI_FIT_ERR, N3PI_OBS_FIT_ERR);
  SetGFAttr(gf_n3pi_obs, "M_{3#pi} [MeV/c^{2}]", " ", "gf_n3pi_obs");
  gf_n3pi_obs -> SetMarkerColor(kBlack);
  gf_n3pi_obs -> SetMarkerSize(0.8);
  gf_n3pi_obs -> RemovePoint(0);
  
  gf_n3pi_obs -> GetYaxis() -> SetNdivisions(512);
  gf_n3pi_obs -> GetYaxis() -> SetRangeUser(hsig -> GetMinimum(), hsig -> GetMaximum() * 1.2);
  gf_n3pi_obs -> GetYaxis() -> SetTitleSize(0.07);
  gf_n3pi_obs -> GetYaxis() -> SetTitleOffset(.7);
  gf_n3pi_obs -> GetYaxis() -> SetLabelSize(0.05);
  
  gf_n3pi_obs -> GetXaxis() -> SetLabelSize(0.03);
  //gf_n3pi_obs -> SetLineColor(kBlack);
  gf_n3pi_obs -> GetYaxis() -> SetTitle(TString::Format("Events/[%0.2f MeV/c^{2}]", Delta_m3pi));
  //gf_n3pi_obs -> GetYaxis() -> SetTitle("#sigma(e^{+}e^{-}#rightarrow#pi^{+}#pi^{-}#pi^{0}) [nb]" );
  gf_n3pi_obs -> GetXaxis() -> SetRangeUser(xrange0, xrange1);
  gf_n3pi_obs -> GetXaxis() -> SetTitleOffset(1.2);
  gf_n3pi_obs -> GetXaxis() -> SetLabelSize(0);
  
  //
  gf_n3pi_fit = new TGraphErrors(fit_indx, M3PI_FIT, N3PI_FIT, M3PI_FIT_ERR, N3PI_FIT_ERR);
  SetGFAttr(gf_n3pi_fit, "M_{3#pi} [MeV/c^{2}]", " ", "gf_n3pi_fit");
  gf_n3pi_fit -> SetMarkerStyle(21);
  //gf_n3pi_fit -> SetMarkerColor(kBlue);
  gf_n3pi_fit -> SetLineColor(kBlack);
  gf_n3pi_fit -> SetLineWidth(2);
  gf_n3pi_fit -> RemovePoint(0);
  
  cv_n3pi_fit = new TCanvas("cv_n3pi_fit", "Nobs 3pi", 0, 0, 1000, 800);

  char display1[50], display2[50], display3[50], display4[50], display5[50];

  TPaveText *pt1 = new TPaveText(0.12, 0.68, 0.3, 0.70, "NDC");
  TPaveText *pt2 = new TPaveText(0.12, 0.61, 0.3, 0.63, "NDC");
  TPaveText *pt3 = new TPaveText(0.12, 0.54, 0.3, 0.57, "NDC");
  
  TPaveText *pt4 = new TPaveText(0.65, 0.72, 0.85, 0.68, "NDC");
  TPaveText *pt5 = new TPaveText(0.65, 0.60, 0.85, 0.61, "NDC");
  
  PteAttr(pt1); pt1 -> SetTextSize(0.05); 
  PteAttr(pt2); pt2 -> SetTextSize(0.05);
  PteAttr(pt3); pt3 -> SetTextSize(0.05);
  PteAttr(pt4); pt4 -> SetTextSize(0.05);
  PteAttr(pt5); pt5 -> SetTextSize(0.05);

  TLatex t(0.1, 0.2,"M_{#omega}");
  
  //sprintf(display1,"#chi^{2}/ndf=%0.2f", chi2_sum / ndf);
  pt1 -> AddText(Form("#Gamma_{#omega}=%0.2f#pm%0.2f MeV", Gam_omega_fit, Gam_omega_err_fit)); 
  pt2 -> AddText(Form("B_{ee}B_{3#pi}=(%0.2f#pm%0.2f)#times10^{-5}", BB_fit * 1e5, BB_err_fit * 1e5));
  pt3 -> AddText(Form("M_{#omega}=%0.2f#pm%0.2f MeV/c^{2}", Mass_omega_fit, Mass_omega_err_fit)); 
  pt4 -> AddText(Form("#chi^{2}/ndf=%0.1f/%0.2d", chi2_sum, ndf)); 
  //pt4 -> AddText(Form("#chi^{2}/ndf=%0.2f", chi2_sum / ndf)); 
  pt5 -> AddText(Form("L_{int}=%0.2f fb^{-1}", Lumi_int_fit * 1e-6)); 
  
  TPad *p2 = new TPad("p2", "p2", 0., 0., 1., 0.25);
  p2 -> Draw();
  p2 -> SetBottomMargin(0.4);
  p2 -> SetLeftMargin(0.1);
  p2 -> SetGrid();
  
  TPad *p1 = new TPad("p1", "p1", 0., 0.25, 1., 1.);
  p1 -> Draw();
  p1 -> SetBottomMargin(0.02);//0.007
  p1 -> SetLeftMargin(0.1);
  p1 -> cd();
  
  gf_n3pi_obs -> Draw("PA");
  gf_n3pi_fit -> Draw("cSame"); //P
  //hsig -> Draw("HistSame");
  pt1 -> Draw("Same");
  pt2 -> Draw("Same");
  pt3 -> Draw("Same");
  pt4 -> Draw("Same");
  pt5 -> Draw("Same");
  
  //
  TLegend * legd_cv = new TLegend(0.12, 0.75, 0.5, 0.9);
  
  SetLegend(legd_cv);
  legd_cv -> SetNColumns(1);
  
  //
  //legd_cv -> AddEntry(gf_n3pi_obs, "N^{obs}_{3#pi}", "lep");
  legd_cv -> AddEntry(gf_n3pi_obs, "Data", "lep");
  legd_cv -> AddEntry(gf_n3pi_fit, "Breit-Wigner", "l");
  //legd_cv -> AddEntry(hsig,  "#pi^{+}#pi^{-}#pi^{0}", "lep");
  
  legd_cv -> Draw("Same");
  
  legtextsize(legd_cv, 0.06);
  
  p2 -> cd();
  
  gf_n3pi_diff -> SetLineColor(0);
  gf_n3pi_diff -> GetYaxis() -> SetNdivisions(505);
  gf_n3pi_diff -> GetYaxis() -> SetTitleSize(33);
  gf_n3pi_diff -> GetYaxis() -> SetTitleFont(43); //43, 132
  gf_n3pi_diff -> GetYaxis() -> SetTitleOffset(1.6);
  gf_n3pi_diff -> GetYaxis() -> SetLabelSize(25);
  gf_n3pi_diff -> GetYaxis() -> SetLabelFont(43); // Absolute front size in pixel (precision 3)
  gf_n3pi_diff -> GetYaxis() -> CenterTitle();
  gf_n3pi_diff -> GetYaxis() -> SetTitle("Residual");
  //gf_n3pi_diff -> GetYaxis() -> SetRangeUser(diff_min * 3., diff_max * 3.);
  gf_n3pi_diff -> GetYaxis() -> SetRangeUser(-5., 5.);
  //gf_n3pi_diff -> GetYaxis() -> SetTitle("(N_{d} - N_{mc})/#sqrt{N_{d}}");
  //gf_n3pi_diff -> GetXaxis() -> SetLabelFont(43); // Absolute font size in pixel (precision 3)

  gf_n3pi_diff -> GetXaxis() -> SetRangeUser(xrange0, xrange1);
  gf_n3pi_diff -> GetXaxis() -> SetTitleOffset(0.85);
  //gf_n3pi_diff -> GetXaxis() -> SetTitleSize(45);
  gf_n3pi_diff -> GetXaxis() -> SetLabelSize(0.15);
  //gf_n3pi_diff -> GetXaxis() -> SetLabelOffset(0.2);
  gf_n3pi_diff -> GetXaxis() -> SetTitleSize(0.2);
  gf_n3pi_diff -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  gf_n3pi_diff -> GetXaxis() -> CenterTitle();
  gf_n3pi_diff -> SetMarkerStyle(21);
  gf_n3pi_diff -> SetMarkerSize(0.7);
  
  //TMultiGraph *mg = new TMultiGraph();
  //mg -> SetTitle("Exclusion graphs");
  
  //mg -> Add(gf_n3pi_diff);
  gf_n3pi_diff -> SetFillColor(1);
  gf_n3pi_diff -> SetLineStyle(9);
  gf_n3pi_diff -> SetFillStyle(3001);
  
  //gf_n3pi_diff -> Draw("a3");
  gf_n3pi_diff -> Draw();

  //mg -> Draw("PC3");

  TLine *line = new TLine(740., 0., 845., 0.);
  line -> SetLineColor(kRed);
  line -> SetLineWidth(2);

  //line -> Draw("Same");
  
}

//
void plot_IM3pi() {

  // attributes
  format_h(hdata, 1, 2);
  format_h(heeg, 6, 2);
  format_h(homegapi, 7, 2);
  format_h(hetagam, 3, 2);
  format_h(hksl, 28, 2);
  format_h(hsig, 4, 2);
  format_h(hmcrest, 37, 2);
  format_h(hmcsum, 2, 2);
  format_h(hbkgsum, 37, 2);

  cv_compr = new TCanvas("cv_compr", "cv_compr", 1000, 800);
  cv_compr -> SetLeftMargin(0.1);

  hdata -> GetYaxis() -> SetNdivisions(512);
  hdata -> GetYaxis() -> SetRangeUser(0., hsig -> GetMaximum() * 2);
  hdata -> GetYaxis() -> SetTitle(TString::Format("Events/[%0.2f MeV/c^{2}]", Delta_m3pi));
  hdata -> GetYaxis() -> SetTitleSize(0.05);
  hdata -> GetYaxis() -> SetTitleOffset(1.);
  hdata -> GetYaxis() -> SetLabelSize(0.035);
  hdata -> GetYaxis() -> CenterTitle();
  
  hdata -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
  hdata -> GetXaxis() -> SetTitleOffset(.8);
  hdata -> GetXaxis() -> SetTitleSize(0.05);
  hdata -> GetXaxis() -> SetLabelSize(0.04);
  hdata -> GetXaxis() -> SetRangeUser(650., 950.);
  //hdata -> GetXaxis() -> SetRangeUser(760., 800.);
  hdata -> GetXaxis() -> CenterTitle();
  
  hdata -> SetLineColor(1);
  hsig -> SetLineWidth(2);

  hdata -> SetName("hdata");
  hmcsum -> SetName("hmcsum");
  hbkgsum -> SetName("hbkgsum");
  hsig -> SetName("hsig");
  
  hdata -> Draw();
  hmcsum -> Draw("samehist");
  hbkgsum -> Draw("samehist");
  hsig -> Draw("samehist");

  //
  TLegend * legd_cv_compr = new TLegend(0.55, 0.6, 0.9, 0.85);
  
  SetLegend(legd_cv_compr);
  legd_cv_compr -> SetTextSize(0.04);
  legd_cv_compr -> SetNColumns(1);
  
  legd_cv_compr -> AddEntry(hdata, "Data", "lep");
  legd_cv_compr -> AddEntry(hsig,  "Signal", "l");
  legd_cv_compr -> AddEntry(hmcsum, "MC sum", "l");
  legd_cv_compr -> AddEntry(hbkgsum, "Backgound sum", "l");
  
  legd_cv_compr -> Draw("Same");
  
}

//
void plot_matr() {

  cout << "Smearing matrix ..." << endl;

  // check smearing matrix reference
  TCanvas * cv = new TCanvas("cv", "Corrlation and smearing matrices", 0, 0, 1500, 800);

  cv -> Divide(2,1);

  cv -> cd(1);
  hsmearmatr -> Draw("COLZ");
  gPad -> SetLogz();
  gPad -> SetLeftMargin(0.12);
  gPad -> SetRightMargin(0.12);
  gPad -> SetBottomMargin(0.12);
  
  
  cv -> cd(2);
  hcorrmatrix -> Draw("COLZ");
  gPad -> SetLogz();
  gPad -> SetLeftMargin(0.12);
  gPad -> SetRightMargin(0.12);
  gPad -> SetBottomMargin(0.12);
  
  
}

