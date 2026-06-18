#include <iostream>
void tree_cut_script() {
gROOT->ProcessLine(".L ../run/tree_cut.C");
gROOT->ProcessLine("tree_cut()");
}
