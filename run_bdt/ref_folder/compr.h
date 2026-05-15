//const TString infile_nm = "/home/bo/Desktop/analysis/crx3pi/output_chains_norm/tree_cut0.root";
//const TString infile_nm = "/home/bo/Desktop/analysis/crx3pi/output_chi2/tree_cut0.root";
//const TString infile1_nm = "/home/bo/Desktop/analysis/sfw1d/sfw1d.root";

//const TString var_nm = "lagvalue_min_7C";
//const TString unit = "";
//const TString var_symb = "#chi^{2}_{7C}";

const TString var_nm = "lagvalue_min_7C";
const TString unit = "";
const TString var_symb = "#chi^{2}_{7C}";

const int binsize = 200;
const double var_min = 0;
const double var_max = 100;

const double IM3pi_min = 720; //760 720
const double IM3pi_max = 820; //800 620

//bkg scaling factors
const double sfw2d_eeg =1.66807;// 1.69177;//1.18501;
const double sfw2d_isr3pi = 0.049098; //0.0491606;//0.0521485;
const double sfw2d_omegapi = 1.24699; //1.26532;//1.25099;
const double sfw2d_etagam = 1.01702; //1.01622;//1.01886;
const double sfw2d_ksl = 1.03706; //1.01702; //1.04951;//1.11415;
const double sfw2d_mcrest = 4.81468; //4.94752;//3.3577;

//signal scaling factor, /home/bo/Desktop/analysis/sf_2g/plot_scan.C

const double sfw1d_isr3pi = 4.60022e-02;
