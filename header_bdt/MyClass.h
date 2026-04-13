//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb 17 11:06:32 2020 by ROOT version 5.34/36
// from TTree h1/etappg
// found on file: mcksl30300.root
//////////////////////////////////////////////////////////

#ifndef MyClass_h
#define MyClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TEnv.h>
#include <TCutG.h>
#include <TMatrixD.h>
#include <TRandom.h>
#include <TVector2.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include "TVectorD.h"
#include "TFitter.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class MyClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run_nr;
   Int_t           ev_nr;
   Int_t           mcflag;
   Int_t           phid;
   Int_t           atyp[3];
   Int_t           btyp[3];
   Int_t           t3down;
   Int_t           t3flag;
   Int_t           ecltrgw;
   Int_t           eclfilfo;
   Int_t           necls;
   Int_t           eclword[8];   //[necls]
   Int_t           eclstream[8];   //[necls]
   Int_t           ecltagnum[8];   //[necls]
   Int_t           kineid;
   Int_t           bit_select;
   Int_t           trgtype;
   Int_t           bphi;
   Int_t           ephi;
   Int_t           wphi;
   Int_t           t1d;
   Int_t           t2d;
   Float_t         bx;
   Float_t         by;
   Float_t         bz;
   Float_t         bpx;
   Float_t         bpy;
   Float_t         bpz;
   Float_t         bene;
   Int_t           nv;
   Int_t           vtx[20];   //[nv]
   Float_t         xv[20];   //[nv]
   Float_t         yv[20];   //[nv]
   Float_t         zv[20];   //[nv]
   Int_t           ntv;
   Int_t           iv[30];   //[ntv]
   Int_t           trknumv[30];   //[ntv]
   Float_t         curv[30];   //[ntv]
   Float_t         phiv[30];   //[ntv]
   Float_t         cotv[30];   //[ntv]
   Float_t         pxtv[30];   //[ntv]
   Float_t         pytv[30];   //[ntv]
   Float_t         pztv[30];   //[ntv]
   Float_t         pmodv[30];   //[ntv]
   Float_t         lenvtx[30];   //[ntv]
   Float_t         vtxchi2[30];   //[ntv]
   Int_t           pidtv[30];   //[ntv]
   Int_t           nt;
   Int_t           trkind[100];   //[nt]
   Int_t           trkver[100];   //[nt]
   Float_t         cur[100];   //[nt]
   Float_t         phi[100];   //[nt]
   Float_t         cot[100];   //[nt]
   Float_t         pxt[100];   //[nt]
   Float_t         pyt[100];   //[nt]
   Float_t         pzt[100];   //[nt]
   Float_t         pmod[100];   //[nt]
   Float_t         xfirst[100];   //[nt]
   Float_t         yfirst[100];   //[nt]
   Float_t         zfirst[100];   //[nt]
   Float_t         length[100];   //[nt]
   Float_t         curla[100];   //[nt]
   Float_t         phila[100];   //[nt]
   Float_t         cotla[100];   //[nt]
   Float_t         pxtla[100];   //[nt]
   Float_t         pytla[100];   //[nt]
   Float_t         pztla[100];   //[nt]
   Float_t         pmodla[100];   //[nt]
   Float_t         xlast[100];   //[nt]
   Float_t         ylast[100];   //[nt]
   Float_t         zlast[100];   //[nt]
   Float_t         xpca[100];   //[nt]
   Float_t         ypca[100];   //[nt]
   Float_t         zpca[100];   //[nt]
   Float_t         curpca[100];   //[nt]
   Float_t         cotpca[100];   //[nt]
   Float_t         phipca[100];   //[nt]
   Int_t           nhit[100];   //[nt]
   Int_t           nfithit[100];   //[nt]
   Float_t         sigcurv[100];   //[nt]
   Float_t         sigcot[100];   //[nt]
   Float_t         sigphi[100];   //[nt]
   Float_t         xp_dcw[100][7][2];   //[nt]
   Float_t         xp_bp[100][7][2];   //[nt]
   Float_t         xp_ip[100][7][2];   //[nt]
   Int_t           n_seg[100][2];   //[nt]
   Float_t         len_seg[100][3][2];   //[nt]
   Float_t         dist_ip[100][2];   //[nt]
   Int_t           good_extrap[100][2];   //[nt]
   Int_t           mc_corr[100];   //[nt]
   Int_t           ntfmc;
   Int_t           ncontr[100];   //[ntfmc]
   Int_t           trkine1[100];   //[ntfmc]
   Int_t           trtype1[100];   //[ntfmc]
   Int_t           trhits1[100];   //[ntfmc]
   Int_t           trkine2[100];   //[ntfmc]
   Int_t           trtype2[100];   //[ntfmc]
   Int_t           trhits2[100];   //[ntfmc]
   Int_t           trkine3[100];   //[ntfmc]
   Int_t           trtype3[100];   //[ntfmc]
   Int_t           trhits3[100];   //[ntfmc]
   Int_t           nclu;
   Float_t         enecl[100];   //[nclu]
   Float_t         tcl[100];   //[nclu]
   Float_t         xcl[100];   //[nclu]
   Float_t         ycl[100];   //[nclu]
   Float_t         zcl[100];   //[nclu]
   Float_t         mmclu[100];   //[nclu]
   Int_t           charged_old[100];   //[nclu]
   Int_t           charged_new[100];   //[nclu]
   Float_t         ene_plane[100][5];   //[nclu]
   Int_t           intime[100];   //[nclu]
   Int_t           nclumc;
   Int_t           Npar[100];   //[nclumc]
   Int_t           Pnum1[100];   //[nclumc]
   Int_t           Pid1[100];   //[nclumc]
   Int_t           Pnum2[100];   //[nclumc]
   Int_t           Pid2[100];   //[nclumc]
   Int_t           Pnum3[100];   //[nclumc]
   Int_t           Pid3[100];   //[nclumc]
   Int_t           ntmc;
   Int_t           kine[50];   //[ntmc]
   Int_t           pidmc[50];   //[ntmc]
   Int_t           virmom[50];   //[ntmc]
   Int_t           indv[50];   //[ntmc]
   Float_t         pxmc[50];   //[ntmc]
   Float_t         pymc[50];   //[ntmc]
   Float_t         pzmc[50];   //[ntmc]
   Int_t           nvtxmc;
   Int_t           kinmom[30];   //[nvtxmc]
   Int_t           mother[30];   //[nvtxmc]
   Float_t         tofvmc[30];   //[nvtxmc]
   Float_t         xvmc[30];   //[nvtxmc]
   Float_t         yvmc[30];   //[nvtxmc]
   Float_t         zvmc[30];   //[nvtxmc]
   Float_t         trkvtxmc[30];   //[nvtxmc]
   Int_t           ntclo;
   Float_t         charge[100];   //[ntclo]
   Int_t           trk_quality[100];   //[ntclo]
   Int_t           extr_zone[100];   //[ntclo]
   Float_t         extrapolation[100][9];   //[ntclo]
   Int_t           best_clu[100];   //[ntclo]
   Int_t           associated[100];   //[ntclo]
   Int_t           nassclu[100];   //[ntclo]
   Float_t         trk_distalong[100];   //[ntclo]
   Float_t         trk_dista_tra[100];   //[ntclo]
   Int_t           n_parts[100];   //[ntclo]
   Float_t         leng_parts[100][3];   //[ntclo]
   Float_t         mome_parts[100][3];   //[ntclo]
   Float_t         mome_partsend[100][3];   //[ntclo]
   Float_t         beta_parts[100][3];   //[ntclo]

   // List of branches
   TBranch        *b_run_nr;   //!
   TBranch        *b_ev_nr;   //!
   TBranch        *b_mcflag;   //!
   TBranch        *b_phid;   //!
   TBranch        *b_atyp;   //!
   TBranch        *b_btyp;   //!
   TBranch        *b_t3down;   //!
   TBranch        *b_t3flag;   //!
   TBranch        *b_ecltrgw;   //!
   TBranch        *b_eclfilfo;   //!
   TBranch        *b_necls;   //!
   TBranch        *b_eclword;   //!
   TBranch        *b_eclstream;   //!
   TBranch        *b_ecltagnum;   //!
   TBranch        *b_kineid;   //!
   TBranch        *b_bit_select;   //!
   TBranch        *b_trgtype;   //!
   TBranch        *b_bphi;   //!
   TBranch        *b_ephi;   //!
   TBranch        *b_wphi;   //!
   TBranch        *b_t1d;   //!
   TBranch        *b_t2d;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_by;   //!
   TBranch        *b_bz;   //!
   TBranch        *b_bpx;   //!
   TBranch        *b_bpy;   //!
   TBranch        *b_bpz;   //!
   TBranch        *b_bene;   //!
   TBranch        *b_nv;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_xv;   //!
   TBranch        *b_yv;   //!
   TBranch        *b_zv;   //!
   TBranch        *b_ntv;   //!
   TBranch        *b_iv;   //!
   TBranch        *b_trknumv;   //!
   TBranch        *b_curv;   //!
   TBranch        *b_phiv;   //!
   TBranch        *b_cotv;   //!
   TBranch        *b_pxtv;   //!
   TBranch        *b_pytv;   //!
   TBranch        *b_pztv;   //!
   TBranch        *b_pmodv;   //!
   TBranch        *b_lenvtx;   //!
   TBranch        *b_vtxchi2;   //!
   TBranch        *b_pidtv;   //!
   TBranch        *b_nt;   //!
   TBranch        *b_trkind;   //!
   TBranch        *b_trkver;   //!
   TBranch        *b_cur;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_cot;   //!
   TBranch        *b_pxt;   //!
   TBranch        *b_pyt;   //!
   TBranch        *b_pzt;   //!
   TBranch        *b_pmod;   //!
   TBranch        *b_xfirst;   //!
   TBranch        *b_yfirst;   //!
   TBranch        *b_zfirst;   //!
   TBranch        *b_length;   //!
   TBranch        *b_curla;   //!
   TBranch        *b_phila;   //!
   TBranch        *b_cotla;   //!
   TBranch        *b_pxtla;   //!
   TBranch        *b_pytla;   //!
   TBranch        *b_pztla;   //!
   TBranch        *b_pmodla;   //!
   TBranch        *b_xlast;   //!
   TBranch        *b_ylast;   //!
   TBranch        *b_zlast;   //!
   TBranch        *b_xpca;   //!
   TBranch        *b_ypca;   //!
   TBranch        *b_zpca;   //!
   TBranch        *b_curpca;   //!
   TBranch        *b_cotpca;   //!
   TBranch        *b_phipca;   //!
   TBranch        *b_nhit;   //!
   TBranch        *b_nfithit;   //!
   TBranch        *b_sigcurv;   //!
   TBranch        *b_sigcot;   //!
   TBranch        *b_sigphi;   //!
   TBranch        *b_xp_dcw;   //!
   TBranch        *b_xp_bp;   //!
   TBranch        *b_xp_ip;   //!
   TBranch        *b_n_seg;   //!
   TBranch        *b_len_seg;   //!
   TBranch        *b_dist_ip;   //!
   TBranch        *b_good_extrap;   //!
   TBranch        *b_mc_corr;   //!
   TBranch        *b_ntfmc;   //!
   TBranch        *b_ncontr;   //!
   TBranch        *b_trkine1;   //!
   TBranch        *b_trtype1;   //!
   TBranch        *b_trhits1;   //!
   TBranch        *b_trkine2;   //!
   TBranch        *b_trtype2;   //!
   TBranch        *b_trhits2;   //!
   TBranch        *b_trkine3;   //!
   TBranch        *b_trtype3;   //!
   TBranch        *b_trhits3;   //!
   TBranch        *b_nclu;   //!
   TBranch        *b_enecl;   //!
   TBranch        *b_tcl;   //!
   TBranch        *b_xcl;   //!
   TBranch        *b_ycl;   //!
   TBranch        *b_zcl;   //!
   TBranch        *b_mmclu;   //!
   TBranch        *b_charged_old;   //!
   TBranch        *b_charged_new;   //!
   TBranch        *b_ene_plane;   //!
   TBranch        *b_intime;   //!
   TBranch        *b_nclumc;   //!
   TBranch        *b_Npar;   //!
   TBranch        *b_Pnum1;   //!
   TBranch        *b_Pid1;   //!
   TBranch        *b_Pnum2;   //!
   TBranch        *b_Pid2;   //!
   TBranch        *b_Pnum3;   //!
   TBranch        *b_Pid3;   //!
   TBranch        *b_ntmc;   //!
   TBranch        *b_kine;   //!
   TBranch        *b_pidmc;   //!
   TBranch        *b_virmom;   //!
   TBranch        *b_indv;   //!
   TBranch        *b_pxmc;   //!
   TBranch        *b_pymc;   //!
   TBranch        *b_pzmc;   //!
   TBranch        *b_nvtxmc;   //!
   TBranch        *b_kinmom;   //!
   TBranch        *b_mother;   //!
   TBranch        *b_tofvmc;   //!
   TBranch        *b_xvmc;   //!
   TBranch        *b_yvmc;   //!
   TBranch        *b_zvmc;   //!
   TBranch        *b_trkvtxmc;   //!
   TBranch        *b_ntclo;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_trk_quality;   //!
   TBranch        *b_extr_zone;   //!
   TBranch        *b_extrapolation;   //!
   TBranch        *b_best_clu;   //!
   TBranch        *b_associated;   //!
   TBranch        *b_nassclu;   //!
   TBranch        *b_trk_distalong;   //!
   TBranch        *b_trk_dista_tra;   //!
   TBranch        *b_n_parts;   //!
   TBranch        *b_leng_parts;   //!
   TBranch        *b_mome_parts;   //!
   TBranch        *b_mome_partsend;   //!
   TBranch        *b_beta_parts;   //!

   int pstrnb;

   double speedc;
   double massomega;
   double masschpion;
   double massneupion;

   double egammamin;
   double Zvmax;
   double Rhovmax;   
   double nb_sigma_T_clust;

   //for smearing the MC tracks to be more like data
   TRandom *generator;
   
   MyClass(TTree *tree=0);
   virtual ~MyClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Main();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   bool IfStreamed(int pstrnb);
   bool IfTriggered();
   bool IfTriggered_EMC();
   bool IfTriggered_DC();
   bool IfTriggered_EMCandDC();
   bool IfFilfoed();
   bool IfFilfoed_28();
   bool IfBroken(int idx1, int idx2);
   
   int pho_indx[100], clu_indx[100], pid_indx[100];
   
   double E_list[100], E_true_list[100], X_list[100], Y_list[100], Z_list[100], T_list[100];
   double Sigma2_E_list[100], Sigma2_X_list[100], Sigma2_Y_list[100], Sigma2_Z_list[100], Sigma2_T_list[100];
   //double etakinfit_min_7C[15];   
   double Trkmass(TLorentzVector bestppl, TLorentzVector bestpmi);
   double DeltaE(TLorentzVector bestppl, TLorentzVector bestpmi);

   
   TVectorD Getpiontrnb();
   TVectorD FillInputvector_7C(int size, int index1, int index2, int index3, int indtr1, int indtr2);
   TVectorD FillSigma2vector_7C(int size, int index1, int index2, int index3, int indtr1, int indtr2);
   TVectorD Gfunc_7C(TVectorD etavector, int rownb, TLorentzVector beam, int indtr1, int indtr2);
   TVectorD Lambdavector(int size, TMatrixD inveSmatrix, TMatrixD dgmatrix, TVectorD constr, TVectorD diffv);
   TVectorD GetTVectorD(int size);
   TVectorD etakinfitfunc(TMatrixD Vmatrix, TMatrixD dgmatrixTrans, TVectorD lambdavector, TVectorD etatilde, int rownb);
   TVectorD Fillsigma2vectorkinfit(int size, TMatrixD matrix, int matrixcol, int matrixrow);
   TVectorD Fillpullsvector(int size, TVectorD sigma2vector_old, TVectorD sigma2vector_new, TVectorD inputvector_old, TVectorD inputvector_new);
   TVectorD Fillpermutvector(int size, TVectorD input, int index1, int index2, int index3);
   
   TLorentzVector Beam;
   TLorentzVector GetLorentzVector(TVector3 vector, double mass);
   //TLorentzVector Gettrack4vectorkinfit(double curv, double cotan, double phi);
   TLorentzVector Gettrack4vectorkinfit(int index);
   TLorentzVector Getphoton4vector(double E, double x, double y, double z);
   
   TMatrixD Getafunc_7C(TVectorD etavector, int rownb, int colnb);
   TMatrixD CovMatrix(TVectorD vector, int row, int col);
   TMatrixD Mtrans(TMatrixD ma);
   TMatrixD MInvert(TMatrixD ma);

   TH1D *hmpi0_iden, *hmpi0_sel, *hrecon_typeI, *hIM3pi, *hIM3pi_bkg_indx_rest, *hIM3pi_bkg_indx0, *hprompt_distr, *hstr_distr, *hTof_clust, *hRhov, *hZv;
   TH1D *h_iv_ip, *h_iv_indx;
   TH1D *h_iv_ip_1, *h_iv_indx_1;
  
  TH2D *hsmearmatr_clust, *hsmearmatr_trk, *hsmearmatr_typeI, *hsmearmatr_typeII, *herror_type, *h_nv_ip, *h_nv_ip_1, *h_xpca, *h_ypca, *h_zpca;
};

#endif

#ifdef MyClass_cxx
MyClass::MyClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  /*if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("mcksl30300.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("mcksl30300.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("mcksl30300.root:/ETAPPG");
      dir->GetObject("h1",tree);

      }*/
   Init(tree);

   // initializing vaiables
   pstrnb = 2;
   generator = new TRandom();
   speedc = 29.9792458;
   massomega = 782.65;
   masschpion = 139.57; // MeV
   massneupion = 134.977; // MeV

   egammamin = 15;
   Rhovmax = 4;
   Zvmax = 10;
   nb_sigma_T_clust = 3;

   
   for (int i = 0; i < 100; i ++) {
     E_true_list[i] = 0.;
     E_list[i] = 0.;
     X_list[i] = 0.; 
     Y_list[i] = 0.;
     Z_list[i] = 0.; 
     T_list[i] = 0.;

     Sigma2_E_list[i] = 0.;
     Sigma2_X_list[i] = 0.;
     Sigma2_Y_list[i] = 0.;
     Sigma2_Z_list[i] = 0.;
     Sigma2_T_list[i] = 0.;

     pho_indx[i] = 0;
     clu_indx[i] = 0;
     pid_indx[i] = 0;
     
   }

   //for (int i = 0; i < 15; i ++) {
   //etakinfit_min_7C[i] = 0.;
   //}

   hprompt_distr = new TH1D("prompt_distr", "number of prompt photon distribution", 10, 0, 10);
   hstr_distr = new TH1D("str_distr", "stream distribution", 10, 0, 10);
   hmpi0_sel = new TH1D("hmpi0_sel", "hmpi0_sel", 200, 0., 1000.);
   hmpi0_iden = new TH1D("hmpi0_iden", "hmpi0_iden", 200, 0., 1000.);
   
   hTof_clust = new TH1D("hTof_clust", "hTof_clust", 100, -1., 1.);
   hRhov = new TH1D("hRhov", "hRhov", 200, 0., 250.);
   hZv = new TH1D("hZv", "hZv", 200, 0., 250.);
   h_iv_ip = new TH1D("h_iv_ip", "h_iv_ip", 10, 0, 10);
   h_iv_indx = new TH1D("h_iv_indx", "h_iv_indx", 10, 0, 10);
   h_iv_ip_1 = new TH1D("h_iv_ip_1", "h_iv_ip_1", 10, 0, 10);
   h_iv_indx_1 = new TH1D("h_iv_indx_1", "h_iv_indx_1", 10, 0, 10);
   hIM3pi = new TH1D("hIM3pi", "hIM3pi", 100, 400., 1000.);
   hIM3pi_bkg_indx_rest = new TH1D("hIM3pi_bkg_indx_rest", "hIM3pi_bkg_indx_rest", 100, 400., 1000.);
   hIM3pi_bkg_indx0 = new TH1D("hIM3pi_bkg_indx0", "hIM3pi_bkg_indx0", 100, 0., 1000.);
   hrecon_typeI = new TH1D("hrecon_typeI", "hrecon_typeI", 100, 400., 1000.);
   
   hsmearmatr_clust = new TH2D("hsmearmatr_clust", "hsmearmatr_clust", 100, 400., 1000., 100, 400., 1000.);
   hsmearmatr_typeI = new TH2D("hsmearmatr_typeI", "IM3pi true v.s. recon box", 100, 400., 1000., 100, 400., 1000.);
   hsmearmatr_typeII = new TH2D("hsmearmatr_typeII", "IM3pi true v.s. recon typeII", 100, 400., 1000., 100, 400., 1000.);
   hsmearmatr_trk = new TH2D("hsmearmatr_trk", "hsmearmatr_trk", 200, 0., 1000., 200, 0., 1000.);
   herror_type = new TH2D("herror_type", "herror_type", 3, 0, 3, 3, 0, 3);
   h_nv_ip = new TH2D("h_nv_ip", "h_nv_ip", 10, 0, 10, 10, 0, 10);
   h_nv_ip_1 = new TH2D("h_nv_ip_1", "h_nv_ip_1", 10, 0, 10, 10, 0, 10);
   h_xpca = new TH2D("h_xpca", "h_xpca", 10, 0, 10, 10, 0, 10);
   h_ypca = new TH2D("h_ypca", "h_ypca", 10, 0, 10, 10, 0, 10);
   h_zpca = new TH2D("h_zpca", "h_zpca", 10, 0, 10, 10, 0, 10);
   
}
   

MyClass::~MyClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run_nr", &run_nr, &b_run_nr);
   fChain->SetBranchAddress("ev_nr", &ev_nr, &b_ev_nr);
   fChain->SetBranchAddress("mcflag", &mcflag, &b_mcflag);
   fChain->SetBranchAddress("phid", &phid, &b_phid);
   fChain->SetBranchAddress("atyp", atyp, &b_atyp);
   fChain->SetBranchAddress("btyp", btyp, &b_btyp);
   fChain->SetBranchAddress("t3down", &t3down, &b_t3down);
   fChain->SetBranchAddress("t3flag", &t3flag, &b_t3flag);
   fChain->SetBranchAddress("ecltrgw", &ecltrgw, &b_ecltrgw);
   fChain->SetBranchAddress("eclfilfo", &eclfilfo, &b_eclfilfo);
   fChain->SetBranchAddress("necls", &necls, &b_necls);
   fChain->SetBranchAddress("eclword", eclword, &b_eclword);
   fChain->SetBranchAddress("eclstream", eclstream, &b_eclstream);
   fChain->SetBranchAddress("ecltagnum", ecltagnum, &b_ecltagnum);
   fChain->SetBranchAddress("kineid", &kineid, &b_kineid);
   fChain->SetBranchAddress("bit_select", &bit_select, &b_bit_select);
   fChain->SetBranchAddress("trgtype", &trgtype, &b_trgtype);
   fChain->SetBranchAddress("bphi", &bphi, &b_bphi);
   fChain->SetBranchAddress("ephi", &ephi, &b_ephi);
   fChain->SetBranchAddress("wphi", &wphi, &b_wphi);
   fChain->SetBranchAddress("t1d", &t1d, &b_t1d);
   fChain->SetBranchAddress("t2d", &t2d, &b_t2d);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("by", &by, &b_by);
   fChain->SetBranchAddress("bz", &bz, &b_bz);
   fChain->SetBranchAddress("bpx", &bpx, &b_bpx);
   fChain->SetBranchAddress("bpy", &bpy, &b_bpy);
   fChain->SetBranchAddress("bpz", &bpz, &b_bpz);
   fChain->SetBranchAddress("bene", &bene, &b_bene);
   fChain->SetBranchAddress("nv", &nv, &b_nv);
   fChain->SetBranchAddress("vtx", vtx, &b_vtx);
   fChain->SetBranchAddress("xv", xv, &b_xv);
   fChain->SetBranchAddress("yv", yv, &b_yv);
   fChain->SetBranchAddress("zv", zv, &b_zv);
   fChain->SetBranchAddress("ntv", &ntv, &b_ntv);
   fChain->SetBranchAddress("iv", iv, &b_iv);
   fChain->SetBranchAddress("trknumv", trknumv, &b_trknumv);
   fChain->SetBranchAddress("curv", curv, &b_curv);
   fChain->SetBranchAddress("phiv", phiv, &b_phiv);
   fChain->SetBranchAddress("cotv", cotv, &b_cotv);
   fChain->SetBranchAddress("pxtv", pxtv, &b_pxtv);
   fChain->SetBranchAddress("pytv", pytv, &b_pytv);
   fChain->SetBranchAddress("pztv", pztv, &b_pztv);
   fChain->SetBranchAddress("pmodv", pmodv, &b_pmodv);
   fChain->SetBranchAddress("lenvtx", lenvtx, &b_lenvtx);
   fChain->SetBranchAddress("vtxchi2", vtxchi2, &b_vtxchi2);
   fChain->SetBranchAddress("pidtv", pidtv, &b_pidtv);
   fChain->SetBranchAddress("nt", &nt, &b_nt);
   fChain->SetBranchAddress("trkind", trkind, &b_trkind);
   fChain->SetBranchAddress("trkver", trkver, &b_trkver);
   fChain->SetBranchAddress("cur", cur, &b_cur);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("cot", cot, &b_cot);
   fChain->SetBranchAddress("pxt", pxt, &b_pxt);
   fChain->SetBranchAddress("pyt", pyt, &b_pyt);
   fChain->SetBranchAddress("pzt", pzt, &b_pzt);
   fChain->SetBranchAddress("pmod", pmod, &b_pmod);
   fChain->SetBranchAddress("xfirst", xfirst, &b_xfirst);
   fChain->SetBranchAddress("yfirst", yfirst, &b_yfirst);
   fChain->SetBranchAddress("zfirst", zfirst, &b_zfirst);
   fChain->SetBranchAddress("length", length, &b_length);
   fChain->SetBranchAddress("curla", curla, &b_curla);
   fChain->SetBranchAddress("phila", phila, &b_phila);
   fChain->SetBranchAddress("cotla", cotla, &b_cotla);
   fChain->SetBranchAddress("pxtla", pxtla, &b_pxtla);
   fChain->SetBranchAddress("pytla", pytla, &b_pytla);
   fChain->SetBranchAddress("pztla", pztla, &b_pztla);
   fChain->SetBranchAddress("pmodla", pmodla, &b_pmodla);
   fChain->SetBranchAddress("xlast", xlast, &b_xlast);
   fChain->SetBranchAddress("ylast", ylast, &b_ylast);
   fChain->SetBranchAddress("zlast", zlast, &b_zlast);
   fChain->SetBranchAddress("xpca", xpca, &b_xpca);
   fChain->SetBranchAddress("ypca", ypca, &b_ypca);
   fChain->SetBranchAddress("zpca", zpca, &b_zpca);
   fChain->SetBranchAddress("curpca", curpca, &b_curpca);
   fChain->SetBranchAddress("cotpca", cotpca, &b_cotpca);
   fChain->SetBranchAddress("phipca", phipca, &b_phipca);
   fChain->SetBranchAddress("nhit", nhit, &b_nhit);
   fChain->SetBranchAddress("nfithit", nfithit, &b_nfithit);
   fChain->SetBranchAddress("sigcurv", sigcurv, &b_sigcurv);
   fChain->SetBranchAddress("sigcot", sigcot, &b_sigcot);
   fChain->SetBranchAddress("sigphi", sigphi, &b_sigphi);
   fChain->SetBranchAddress("xp_dcw", xp_dcw, &b_xp_dcw);
   fChain->SetBranchAddress("xp_bp", xp_bp, &b_xp_bp);
   fChain->SetBranchAddress("xp_ip", xp_ip, &b_xp_ip);
   fChain->SetBranchAddress("n_seg", n_seg, &b_n_seg);
   fChain->SetBranchAddress("len_seg", len_seg, &b_len_seg);
   fChain->SetBranchAddress("dist_ip", dist_ip, &b_dist_ip);
   fChain->SetBranchAddress("good_extrap", good_extrap, &b_good_extrap);
   fChain->SetBranchAddress("mc_corr", mc_corr, &b_mc_corr);
   fChain->SetBranchAddress("ntfmc", &ntfmc, &b_ntfmc);
   fChain->SetBranchAddress("ncontr", ncontr, &b_ncontr);
   fChain->SetBranchAddress("trkine1", trkine1, &b_trkine1);
   fChain->SetBranchAddress("trtype1", trtype1, &b_trtype1);
   fChain->SetBranchAddress("trhits1", trhits1, &b_trhits1);
   fChain->SetBranchAddress("trkine2", trkine2, &b_trkine2);
   fChain->SetBranchAddress("trtype2", trtype2, &b_trtype2);
   fChain->SetBranchAddress("trhits2", trhits2, &b_trhits2);
   fChain->SetBranchAddress("trkine3", trkine3, &b_trkine3);
   fChain->SetBranchAddress("trtype3", trtype3, &b_trtype3);
   fChain->SetBranchAddress("trhits3", trhits3, &b_trhits3);
   fChain->SetBranchAddress("nclu", &nclu, &b_nclu);
   fChain->SetBranchAddress("enecl", enecl, &b_enecl);
   fChain->SetBranchAddress("tcl", tcl, &b_tcl);
   fChain->SetBranchAddress("xcl", xcl, &b_xcl);
   fChain->SetBranchAddress("ycl", ycl, &b_ycl);
   fChain->SetBranchAddress("zcl", zcl, &b_zcl);
   fChain->SetBranchAddress("mmclu", mmclu, &b_mmclu);
   fChain->SetBranchAddress("charged_old", charged_old, &b_charged_old);
   fChain->SetBranchAddress("charged_new", charged_new, &b_charged_new);
   fChain->SetBranchAddress("ene_plane", ene_plane, &b_ene_plane);
   fChain->SetBranchAddress("intime", intime, &b_intime);
   fChain->SetBranchAddress("nclumc", &nclumc, &b_nclumc);
   fChain->SetBranchAddress("Npar", Npar, &b_Npar);
   fChain->SetBranchAddress("Pnum1", Pnum1, &b_Pnum1);
   fChain->SetBranchAddress("Pid1", Pid1, &b_Pid1);
   fChain->SetBranchAddress("Pnum2", Pnum2, &b_Pnum2);
   fChain->SetBranchAddress("Pid2", Pid2, &b_Pid2);
   fChain->SetBranchAddress("Pnum3", Pnum3, &b_Pnum3);
   fChain->SetBranchAddress("Pid3", Pid3, &b_Pid3);
   fChain->SetBranchAddress("ntmc", &ntmc, &b_ntmc);
   fChain->SetBranchAddress("kine", kine, &b_kine);
   fChain->SetBranchAddress("pidmc", pidmc, &b_pidmc);
   fChain->SetBranchAddress("virmom", virmom, &b_virmom);
   fChain->SetBranchAddress("indv", indv, &b_indv);
   fChain->SetBranchAddress("pxmc", pxmc, &b_pxmc);
   fChain->SetBranchAddress("pymc", pymc, &b_pymc);
   fChain->SetBranchAddress("pzmc", pzmc, &b_pzmc);
   fChain->SetBranchAddress("nvtxmc", &nvtxmc, &b_nvtxmc);
   fChain->SetBranchAddress("kinmom", kinmom, &b_kinmom);
   fChain->SetBranchAddress("mother", mother, &b_mother);
   fChain->SetBranchAddress("tofvmc", tofvmc, &b_tofvmc);
   fChain->SetBranchAddress("xvmc", xvmc, &b_xvmc);
   fChain->SetBranchAddress("yvmc", yvmc, &b_yvmc);
   fChain->SetBranchAddress("zvmc", zvmc, &b_zvmc);
   fChain->SetBranchAddress("trkvtxmc", trkvtxmc, &b_trkvtxmc);
   fChain->SetBranchAddress("ntclo", &ntclo, &b_ntclo);
   fChain->SetBranchAddress("charge", charge, &b_charge);
   fChain->SetBranchAddress("trk_quality", trk_quality, &b_trk_quality);
   fChain->SetBranchAddress("extr_zone", extr_zone, &b_extr_zone);
   fChain->SetBranchAddress("extrapolation", extrapolation, &b_extrapolation);
   fChain->SetBranchAddress("best_clu", best_clu, &b_best_clu);
   fChain->SetBranchAddress("associated", associated, &b_associated);
   fChain->SetBranchAddress("nassclu", nassclu, &b_nassclu);
   fChain->SetBranchAddress("trk_distalong", trk_distalong, &b_trk_distalong);
   fChain->SetBranchAddress("trk_dista_tra", trk_dista_tra, &b_trk_dista_tra);
   fChain->SetBranchAddress("n_parts", n_parts, &b_n_parts);
   fChain->SetBranchAddress("leng_parts", leng_parts, &b_leng_parts);
   fChain->SetBranchAddress("mome_parts", mome_parts, &b_mome_parts);
   fChain->SetBranchAddress("mome_partsend", mome_partsend, &b_mome_partsend);
   fChain->SetBranchAddress("beta_parts", beta_parts, &b_beta_parts);
   Notify();
}

Bool_t MyClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyClass_cxx
