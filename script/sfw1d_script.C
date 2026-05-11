#include <iostream>
void sfw1d_script() {
gROOT->ProcessLine(".L ../run_bdt/sfw1d.C");
gROOT->ProcessLine("sfw1d()");
}
