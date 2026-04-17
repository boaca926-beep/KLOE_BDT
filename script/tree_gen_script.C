#include <iostream>
void tree_gen_script() {
gROOT->ProcessLine(".L ../run_bdt/tree_gen.C");
gROOT->ProcessLine("tree_gen()");
}
