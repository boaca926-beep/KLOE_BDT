#include <iostream>
void omega_fit_script() {
gROOT->ProcessLine(".L ../run/omega_fit.C");
gROOT->ProcessLine("omega_fit()");
}
