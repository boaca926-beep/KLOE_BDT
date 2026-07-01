#include <iostream>
void tree_cut_script() {
gROOT->ProcessLine(".L ../run/tree_cut_raw.C");
gROOT->ProcessLine("tree_cut_raw()");
}
