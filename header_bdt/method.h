TRandom *rnd=0;

void getObj(TFile * f){

  TIter next_tree(f -> GetListOfKeys());

  TString objnm_tree, classnm_tree;
  
  int i = 0;
  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop

    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();

    cout << "classnm = " << classnm_tree << ", objnm = " << objnm_tree << ", " << key -> GetSeekKey() << endl;

  }

  //f -> Close();
  
}

//
void checkFile(TFile *f_input){

  TIter next_tree(f_input -> GetListOfKeys());

  TString objnm_tree, classnm_tree;

  int i = 0;
  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {
    
    i ++;
    
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();
    
    cout << "tree" << i << ": classnm = " << classnm_tree << ", objnm = " << objnm_tree << endl;
    
  }

}

//
double GetFBeta(double a1_temp, double b1_temp, double c1_temp, double m2pi_temp) {
  m2pi_temp = m2pi_temp / 1000.;
  double fbeta = a1_temp + 1. / (exp((m2pi_temp - c1_temp) / b1_temp) - 1.);
  /*cout << "a1 = " << a1 << ", a2 = " << a2 << "\n"
    << "b1 = " << b1 << ", b2 = " << b2 << "\n"
    << "c1 = " << c1 << ", c2 = " << c2 << "\n\n";*/
  //cout << "fbeta = " << fbeta << endl;
  return fbeta;
}

//
void checkList(TList *list_tmp){

  TIter next(list_tmp);
  TObject* object = 0;
  int obj_indx = 0;
  while ((object = next()))
    {

      cout << "[" << obj_indx << "]: " << object -> GetName() << endl;
      obj_indx ++;
      
    }
  
  

}

//
void gethist(TFile *f, TList *list_tmp, int list_key){

  TIter next_tree(f -> GetListOfKeys());

  TString objnm_tree, classnm_tree;
  
  int i = 0;
  TKey *key;
  
  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop
    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();

    cout << "classnm = " << classnm_tree << ", objnm = " << objnm_tree << ", " << key -> GetSeekKey() << endl;

    
    TList *list_loop = (TList*)f -> Get(objnm_tree);

    // loop over list objects
    TIter next(list_loop);
    TObject* object = 0;
    int obj_indx = 0;
    while ((object = next())){
      if (list_key == key -> GetSeekKey()) {
	list_tmp -> Add(object);
	cout << "[" << obj_indx << "]: " << object -> GetName() << endl;
      }
      obj_indx ++;
    }
    //break;
  }
  //cout << "!!" << endl;
  //f -> Close();
  
}

//
double DetectorEvent_fcn(double m, double *para_tmp) {
  // m[0]: mTrue
  // para[0] = frac
  // para[1] = smallBias
  // para[2] = smallSigma
  // para[3] = wideBias
  // para[4] = wideSigma

  // smear by double-gaussian
  if(rnd->Rndm()>para_tmp[0]) {
    return rnd->Gaus(m+para_tmp[1],para_tmp[2]);
  } else {
    return rnd->Gaus(m+para_tmp[3],para_tmp[4]);
  }
}

//
double DetectorEvent(double mTrue) {
  // smear by double-gaussian
  if(rnd->Rndm()>frac) {
    return rnd->Gaus(mTrue+smallBias,smallSigma);
  } else {
    return rnd->Gaus(mTrue+wideBias,wideSigma);
  }
}


