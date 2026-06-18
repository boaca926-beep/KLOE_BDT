#include <iostream>
void sfw1d_script() {
gROOT->ProcessLine(".L ../run/sfw1d.C");
gROOT->ProcessLine("sfw1d()");
}
