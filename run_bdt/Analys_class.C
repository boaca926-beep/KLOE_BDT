#include <TString.h>
#include <TChain.h>
#include <iostream>
#include <fstream>
#include <string>

#include "../header_bdt/path.h"

using namespace std;

int Analys_class(TString rootFile, TString sampleFile) {


  // char *rootFile = new char[strlen(argv[2]) + 1];
  // strcpy(rootFile, argv[2]);
  Int_t RUNNR=0;
  TChain *tree = new TChain("ETAPPG/h1");
  // reading list of ROOT files from the file list
  // files can be commented with '!'
  string line;
  ifstream filelist(rootFile); // = delete
  //  ifstream filelist("filelist.txt");
  if (filelist.is_open()) {
    while (!filelist.eof()) {
      if (getline(filelist, line, '\n'))
        if (line[0] != '!') {
          tree->Add(line.data());
          cout << "Adding file: " << line << " to the chain of files" << endl;
          tree->SetBranchAddress("run_nr",&RUNNR); 
          /*for (Int_t irow=0;irow<tree->GetEntries();irow++) {
      	 	tree->GetEntry(irow);      			
      	 }*/
      	 //cout<<RUNNR<<endl;
        }
    }
    filelist.close();
  } else {
    cout << "Unable to open filelist" << endl;
    return 0;
  }
  
  

 //  //MC del, kommentera bort för att köra data
 //  //open histrogram file here

  TFile *myfile;
  myfile = new TFile(sampleFile + ".root","RECREATE");
  myfile->cd();
  
  MyClass *analysis = new MyClass(tree);
  analysis->Main();

  //saving and closing histogram file
  myfile->Write();
  myfile->Close();

  delete analysis;

  return 0;
}
