// input root files
const TString inpath_data = "/home/bo/Desktop/analysis/chains_102_norm"; // Data with nominal conditions
const TString inpath_ufo_norm = "/home/bo/Desktop/analysis/chains_ufo_norm"; // UFO with nominal conditions
double sqrtS, Lumi_int;
double frac, smallBias, smallSigma;
double wideBias, wideSigma;
double pi;
double crx_scale;

TRandom *rnd;

TKey * key;

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
void checkList(TObjArray *list_tmp){

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
void getObj(TFile *f){

  TIter next_tree(f -> GetListOfKeys());

  TString objnm_tree, classnm_tree;
  
  int i = 0;

  while ( (key = (TKey *) next_tree() ) ) {// start tree while lop

    objnm_tree   =  key -> GetName();
    classnm_tree = key -> GetClassName();
    key -> GetSeekKey();

    cout << "classnm = " << classnm_tree << ", objnm = " << objnm_tree << ", " << key -> GetSeekKey() << endl;

  }

  //f -> Close();
  
}

//
void gethist(TFile *f, TList *list_tmp, int list_key){

  TIter next_tree(f -> GetListOfKeys());

  TString objnm_tree, classnm_tree;
  
  int i = 0;

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
/*
double binomial_err(double nb_true, double nb_gen) {
  double error = 0.;
  double ratio = 0.; 

  if (nb_gen != 0.) {
    ratio = nb_true / nb_gen;
    error = TMath::Sqrt(ratio * (1. - ratio) / nb_gen);
  }
   
  //cout << "true = " << nb_true << ", gen = " << nb_gen << ", ratio = " << ratio << ", error = " << error << endl;

  return error;
}
*/

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
double GetISRLumi_apprx(double m3pi, double m3pi_lower, double m3pi_upper, double W0_full) {

  double isrlumi = 0.;
  double Delta_m3pi = (m3pi_upper - m3pi_lower) * 1e-3;
  
  m3pi = m3pi * 1e-3;

  		 
  isrlumi = W0_full * (2. * m3pi / TMath::Power(sqrtS, 2.)) * Lumi_int * Delta_m3pi;
  //isrlumi = W0_full * (2. * m3pi / TMath::Power(sqrtS, 2.)) * 38794.1 * Delta_m3pi;
    
  //cout << "m3pi = " << m3pi << " [GeV/c^2], width = " << Delta_m3pi << " [GeV/c^2], W0_full = " << W0_full << ", isrlumi = " << isrlumi << ", Lumi_int = " << Lumi_int << "\n";

  return isrlumi;
  
}

//
double GetBkgErr(double NData, double nMC, double NMC, double f, double f_err) {

  double r = nMC / NMC;
  double r_err = 0.; 
  double r_ratio = 0.;
  
  if (nMC != 0.) {
    r = nMC / NMC;
    r_err = binomial_err(nMC, NMC);
    r_ratio = r_err / r;
  }
  
  double f_ratio = f_err / f;
  
  double n_scaled = NData * r * f;

  //cout << "n_scaled = " << n_scaled << endl;
  
  double n_scaled_err = n_scaled * TMath::Sqrt(r_ratio * r_ratio + f_ratio * f_ratio);
  
  return n_scaled_err;
  
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




//
double getloglh(double n_d, double mu) {
  double value = 0.;

  value = n_d * TMath::Log(mu) - mu + n_d - n_d * TMath::Log(n_d);
  //value = n_d * TMath::Log(mu) - mu;
  
  return value;
}

double GetBB_new(double sigma_max, double Mass_V) {

  double BB = sigma_max * 1e-6 / crx_scale * Mass_V * Mass_V / 12. / pi;

  return BB;

}

double crystalball(double x, double alpha, double n, double sigma, double mean) {
  //cout << "sigma = " << sigma << endl;
  
  if (sigma < 0.) return 0.;
  double z = (x - mean)/sigma;
  if (alpha < 0) z = -z;
  double abs_alpha = TMath::Abs(alpha);
  if (z  > - abs_alpha)
    return TMath::Exp(- 0.5 * z * z);
  else {
    double nDivAlpha = n/abs_alpha;
    double AA = TMath::Exp(-0.5*abs_alpha*abs_alpha);
    double B = nDivAlpha -abs_alpha;
    double arg = nDivAlpha/(B-z);
    return AA * TMath::Power(arg,n);
  }

}
