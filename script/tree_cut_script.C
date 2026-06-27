#include <iostream>
void tree_cut_script() {
gROOT->ProcessLine(".L ../run_bdt/tree_cut_bdt_scaled.C");
gROOT->ProcessLine("tree_cut_bdt_scaled()");
}
