#include <iostream>
void tree_gen_script() {
gROOT->ProcessLine(".L ../run/tree_gen.C");
gROOT->ProcessLine("tree_gen()");
}
