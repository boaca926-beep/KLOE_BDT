//void Proceed(const TString input_str = "Analys_class(\"../../path_small/ksl_path1\",\"../../ksl_pre\")"){
//void Proceed(const TString input_str = "Analys_class(\"../path_chain/SIG_sum\",\"../../sig_sum\")"){
//void Proceed(const TString input_str = "Analys_class(\"../../path_norm/sig_path_sum\",\"../../sig_sum\")"){
//void Proceed(const TString input_str = "Analys_class(\"../../path_small/sig_path\",\"../../sig\")"){
//void Proceed(const TString input_str = "Analys_class(\"../../path_small/sig_tmp\",\"../../sig_tmp1\")"){
//void Proceed(const TString input_str = "Analys_class(\"../path_chain/exp_path\",\"../../exp\")"){
//void Proceed(const TString input_str = "Analys_class(\"../path_chain/sig_path_tmp\",\"../../sig\")"){
//void Proceed(const TString input_str = "Analys_class(\"../path_chain/sig_path_chain\",\"../../sig\")"){
void Proceed(const TString input_str = "Analys_class(\"../path_chain/sig_path_small\",\"../../sig\")"){
//void Proceed(const TString input_str = "Analys_class(\"../path_chain/eeg_path_temp\",\"../../sig\")"){

  cout << "READING INPUT FILES !" << endl;
  //gROOT -> ProcessLine( "gErrorIgnoreLevel = 1001;")
  //gROOT -> ProcessLine(".L ../run/MyClass.C");
  //gROOT -> ProcessLine(".L ../run/Analys_class.C");

  // BDT
  gROOT -> ProcessLine(".L ../run_bdt/MyClass.C");
  gROOT -> ProcessLine(".L ../run_bdt/Analys_class.C");

  // further bkg rejection, seems crashed
  //gROOT -> ProcessLine(".L ../run_vertex_bkgrej/MyClass.C");
  //gROOT -> ProcessLine(".L ../run/Analys_class.C");

  // ML training, pi0 mass, smearing pi0 decay photons
  //gROOT -> ProcessLine(".L ../run_vertex/MyClass.C");
  //gROOT -> ProcessLine(".L ../run/Analys_class.C");


  // ufo sample
  //gROOT -> ProcessLine(".L ../run_ufo/MyClass.C");
  //gROOT -> ProcessLine(".L ../run/Analys_class.C");
  
  gROOT -> ProcessLine(input_str);
}
