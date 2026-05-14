void run_script() {
  gROOT->ProcessLine(".L ../run/MyClass.C");
  gROOT->ProcessLine(".L ../run_bdt/Analys_class.C");
  gROOT->ProcessLine("Analys_class(rootFile, sampleFile)");
}
