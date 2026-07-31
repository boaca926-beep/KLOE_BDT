void cut_flow() {
  // Declare filename first!
  const char *filename = "/home/bo/Desktop/bdt_raw_TDATA_chain_false/cut/tree_pre.root";
  
  TFile *f = TFile::Open(filename);
  
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: Cannot open " << filename << std::endl;
    return;
  }
  
  TList *keys = f->GetListOfKeys();
  TIter next(keys);
  TObject *obj;
  std::cout << "\n=== Cut Flow (entries per tree) ===\n";
  
  while ((obj = next())) {
    TString name = obj->GetName();
    TTree *tree = dynamic_cast<TTree*>(f->Get(name));
    if (tree) {
      TNamed *nb = (TNamed*)tree->GetUserInfo()->FindObject("nb_pre");
      Long64_t nb_pre = 0;
      if (nb) nb_pre = TString(nb->GetTitle()).Atoll();
      Long64_t n_final = tree->GetEntries();
      std::cout << Form("%-15s : pre‑cut = %lld, final = %lld", name.Data(), nb_pre, n_final) << std::endl;
    }
  }
  
  f->Close();
}
