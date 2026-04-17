#include <iostream>
void plot_script() {
  gROOT->ProcessLine(".L ../run_bdt/plot_omega.C");
  gROOT->ProcessLine("plot_omega()");
}
