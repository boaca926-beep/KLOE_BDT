#include <iostream>
void tree_cut_script() {
gROOT->ProcessLine(".L ../run_bdt/tree_cut_bdt_tuning_post_fit.C");
gROOT->ProcessLine("tree_cut_bdt_tuning_post_fit()");
}
