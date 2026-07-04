#include <iostream>
void hist_script() {
gROOT->ProcessLine(".L ../run/gethist.C");
gROOT->ProcessLine("gethist()");
}
