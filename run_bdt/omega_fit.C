#include "../header_bdt/sm_para.h"
#include "../header_bdt/path.h"
#include "../header_bdt/method.h"
#include "../header_bdt/omega_fit.h"
#include "../header_bdt/cut_para.h"

int omega_fit(){

  gErrorIgnoreLevel = kError;
  
  cout << "Extract omega parameters ..." << endl;

  // ------------------------------------------------------------------
  // 1. MC-Data tuning in omega region, extract combinatorical background
  // ------------------------------------------------------------------
  // using gethist.C output
  
  // efficiency
  get_efficy();
  //hsig_true->Draw();
  hefficy->Draw();
  
  return 0;
  
}
