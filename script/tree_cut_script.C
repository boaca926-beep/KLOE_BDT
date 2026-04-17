#include <iostream>
void tree_cut_script() {
gROOT->ProcessLine(".L ../run_bdt/tree_cut.C");
gROOT->ProcessLine("tree_cut()");
}
