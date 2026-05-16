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
void checkArray(TObjArray *array){

  // Create a TIter object for the TObjArray
  TIter next(array);

  TObject* object = 0;
  int obj_indx = 0;
  while ((object = next()))
    {

      cout << "[" << obj_indx << "]: " << object -> GetName() << endl;
      obj_indx ++;
      
    }
  
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
