TFile *f = TFile::Open("tree_pre.root");
TParameter<Long64_t> *param = (TParameter<Long64_t>*)f->Get("nb_pre_TDATA");
if (param) {
    Long64_t nb = param->GetVal();
    cout << "TDATA pre-selected = " << nb << endl;
}
// Similarly for other channels
