#include <iostream>
void sfw2d_script() {
gROOT->ProcessLine(".L ../run_bdt/sfw2d.C");
gROOT->ProcessLine("sfw2d()");
}
