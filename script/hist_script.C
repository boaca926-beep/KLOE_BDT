#include <iostream>
void hist_script() {
gROOT->ProcessLine(".L ../run_bdt/gethist.C");
gROOT->ProcessLine("gethist()");
}
