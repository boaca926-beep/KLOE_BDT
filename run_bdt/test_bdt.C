#include <TMVA/RBDT.hxx>
#include <TMVA/RTensor.hxx>
#include <ROOT/RDataFrame.hxx>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <iostream>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TMath.h>

#include "../header_bdt/helper.h"
#include "../header_bdt/cut_para.h"

using namespace TMVA::Experimental;

/*
const double chi2_cut = 43;
const double angle_cut = 138;
const double deltaE_cut = -150;
const double beta_cut = 1.98;
const double c0 = 0.11;
const double c1 = 0.8;
*/

//const TString phys_ch[2] = {"TETAGAM", "#eta#gamma"};
//const TString phys_ch[2] = {"TISR3PI_SIG", "3#pi#gamma"};
const TString phys_ch[2] = {"TDATA", "Data"};

const TString ch_nm = phys_ch[0];
const TString ch_type = phys_ch[1];

	
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

void test_bdt(const char* model_filename = "/home/kloe/Desktop/KLOE_BDT/models/bdt_pi0_TCOMB.root",
                const char* data_filename = "/home/kloe/Desktop/KLOE_BDT/dataset/kloe_bdt.root"){

    gErrorIgnoreLevel = kError;
    TGaxis::SetMaxDigits(4);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetFitFormat("6.4g");
  
    
    // Manually load libraries
    gSystem->Load("libTMVA");
    gSystem->Load("libTMVAUtils");
    
    std::cout << "Testing model_KLOE ..." << std::endl;  // Added std::

    // 1. Check if model file exists
    if (gSystem->AccessPathName(model_filename)) {
        std::cout << "" << model_filename << " does not exists!" << std::endl;  // Added std::
        return;
    }

    // 2. Load the BDT model
    std::cout << "Loading model from " << model_filename << std::endl;  // Added std::
    TMVA::Experimental::RBDT bdt("BDT_pi0", model_filename);  // ← Declaration is HERE

    std::cout << "✓ Model loaded successfully!" << std::endl;  // Added this line

    //std::cout << "Number of input features: " << bdt.GetNInputDim() << std::endl;

    // Define histos
    TH1D* he1 = new TH1D("he1", "", 200, 0, 500);
    TH1D* he2 = new TH1D("he2", "", 200, 0, 500);
    TH1D* he3 = new TH1D("he3", "", 200, 0, 500);

    TH1D* hm_gg = new TH1D("hm_gg", "", 200, 0, 1000);
    TH1D* hcos_theta = new TH1D("hcos_theta", "", 200, -1, 1);
    TH1D* hopen_angle = new TH1D("hopen_angle", "", 200, 0, pi);
    TH1D* hE_asym = new TH1D("hE_asym", "", 200, 0, 1);
    TH1D* he_min_x_angle = new TH1D("he_min_x_angle", "", 200, 0, 1000);
    TH1D* hE_diff = new TH1D("hE_diff", "", 200, 0, 500);
    TH1D* hasym_x_angle = new TH1D("hasym_x_angle", "", 200, 0, pi);

    TH1D* hM3pi_BDT = new TH1D("hM3pi_BDT", "", 200, 400, 1000); // BDT selection
    TH1D* hM3pi_BDT_good = new TH1D("hM3pi_BDT_good", "", 200, 400, 1000);
    TH1D* hM3pi_BDT_bad = new TH1D("hM3pi_BDT_bad", "", 200, 400, 1000);
    TH1D* hM3pi_BDT_best = new TH1D("hM3pi_BDT_best", "", 200, 400, 1000);
    
    TH1D* hM_gg_BDT = new TH1D("hM_gg_BDT", "", 200, 50, 200); 
    TH1D* hM_gg_BDT_good = new TH1D("hM_gg_BDT_good", "", 200, 50, 200);
    TH1D* hM_gg_BDT_bad = new TH1D("hM_gg_BDT_bad", "", 200, 50, 200);
    
    TH1D* hE1_BDT_good = new TH1D("hE1_BDT_good", "", 200, 0, 500); 
    TH1D* hE1_BDT_bad = new TH1D("hE1_BDT_bad", "", 200, 0, 500); 
    TH1D* hE1_BDT_best = new TH1D("hE1_BDT_best", "", 200, 0, 500); 
    
    TH1D* hE2_BDT_good = new TH1D("hE2_BDT_good", "", 200, 0, 500); 
    TH1D* hE2_BDT_bad = new TH1D("hE2_BDT_bad", "", 200, 0, 500); 

    TH1D* hE1 = new TH1D("hE1", "", 200, 0, 500); // KLOE selection
    TH1D* hE1_good = new TH1D("hE1_good", "", 200, 0, 500); 
    TH1D* hE1_bad = new TH1D("hE1_bad", "", 200, 0, 500); 

    TH1D* hE2 = new TH1D("hE2", "", 200, 0, 500); 
    TH1D* hE2_good = new TH1D("hE2_good", "", 200, 0, 500); 
    TH1D* hE2_bad = new TH1D("hE2_bad", "", 200, 0, 500); 

    TH1D* hM_gg = new TH1D("hM_gg", "", 200, 50, 200); 
    TH1D* hM_gg_good = new TH1D("hM_gg_good", "", 200, 50, 200);
    TH1D* hM_gg_bad = new TH1D("hM_gg_bad", "", 200, 50, 200);

    TH1D* hM3pi = new TH1D("hM3pi", "", 200, 400, 1000);
    TH1D* hM3pi_good = new TH1D("hM3pi_good", "", 200, 400, 1000);
    TH1D* hM3pi_bad = new TH1D("hM3pi_bad", "", 200, 400, 1000);
    
    int evnt_KLOE = 0;
    int evnt_good = 0;
    int evnt_bad = 0;

    int bdt_indx = -999;
    int kloe_indx = -999;

    // If data file exists, process it with RDataFame
    if(!gSystem -> AccessPathName(data_filename)){
    
        cout << "\nProcessing data file: " << data_filename << endl;

        // Open the root file
        TFile* file = TFile::Open(data_filename);
        if (!file || file -> IsZombie())
        {
            cout << "Error: Cannot open file " << data_filename << endl;
            return;
        }

        // Get the tree
	cout << ch_nm << endl;
	
        //TTree* tree = (TTree*)file -> Get("TISR3PI_SIG");
	TTree* tree = (TTree*)file -> Get(ch_nm);
	//TTree* tree = (TTree*)file -> Get("TDATA");
        
        if (!tree) {
            cout << "Error: Cannot find 'tree' in file" << endl;
            file -> Close();
            return;
        }

        int nentries = tree -> GetEntries();
        cout << "Tree has " << nentries << " entries" << endl;

        // Set branch addres for input features
        double lagvalue_min_7C = 0., deltaE = 0., betapi0 = 0.,  angle_pi0gam12 = 0.;
	double ppIM = 0.;
        double E1 = 0., px1 = 0., py1 = 0., pz1 = 0.;
        double E2 = 0., px2 = 0., py2 = 0., pz2 = 0.;
        double E3 = 0., px3 = 0., py3 = 0., pz3 = 0.;
	double ppl_E = 0., ppl_px = 0., ppl_py = 0., ppl_pz = 0.;
        double pmi_E = 0., pmi_px = 0., pmi_py = 0., pmi_pz = 0.;
        
        int bkg_indx, recon_indx;

	tree -> SetBranchAddress("Br_deltaE", &deltaE);
	tree -> SetBranchAddress("Br_angle_pi0gam12", &angle_pi0gam12);
	tree -> SetBranchAddress("Br_betapi0", &betapi0);
    
        tree -> SetBranchAddress("Br_lagvalue_min_7C", &lagvalue_min_7C);
	tree -> SetBranchAddress("Br_ppIM", &ppIM);
	
        tree -> SetBranchAddress("Br_bkg_indx", &bkg_indx);
        tree -> SetBranchAddress("Br_recon_indx", &recon_indx);

	tree -> SetBranchAddress("Br_ppl_E", &ppl_E);
        tree -> SetBranchAddress("Br_ppl_px", &ppl_px);
        tree -> SetBranchAddress("Br_ppl_py", &ppl_py);
        tree -> SetBranchAddress("Br_ppl_pz", &ppl_pz);

	tree -> SetBranchAddress("Br_pmi_E", &pmi_E);
        tree -> SetBranchAddress("Br_pmi_px", &pmi_px);
        tree -> SetBranchAddress("Br_pmi_py", &pmi_py);
        tree -> SetBranchAddress("Br_pmi_pz", &pmi_pz);

        tree -> SetBranchAddress("Br_E1", &E1);
        tree -> SetBranchAddress("Br_px1", &px1);
        tree -> SetBranchAddress("Br_py1", &py1);
        tree -> SetBranchAddress("Br_pz1", &pz1);

        tree -> SetBranchAddress("Br_E2", &E2);
        tree -> SetBranchAddress("Br_px2", &px2);
        tree -> SetBranchAddress("Br_py2", &py2);
        tree -> SetBranchAddress("Br_pz2", &pz2);
       
        tree -> SetBranchAddress("Br_E3", &E3);
        tree -> SetBranchAddress("Br_px3", &px3);
        tree -> SetBranchAddress("Br_py3", &py3);
        tree -> SetBranchAddress("Br_pz3", &pz3);

	//tree -> SetBranchAddress("Br_event_id", &event_id);

	
        // Create output file and tree
        TFile* outfile = TFile::Open("output_with_bdt.root", "RECREATE");
        TTree* outtree = new TTree("new_tree", "Tree with BDT response");

        // Output variables
        int out_event;
        int pi0_pho1_idx, pi0_pho2_idx; // indices of bdt pi0 photons
        int prompt_pho_idx; // prompt photon index
	double bdt_score;
        
        outtree -> Branch("event", &out_event);

        // Copy original branches if needed
        //outtree -> Branch("E1", &E1);

        const double energy_threshold = 5.0;

	int n_found = 0;
	
        double m_gg_bdt = 0, m_gg = 0;
	double m3pi = 0, m3pi_bdt = 0;

	double e1_bdt = 0, e2_bdt = 0, e3_bdt = 0;
	double opening_angle = 0, cos_theta = 0;
        double E_asym = 0, e_min_x_angle = 0;
        double asym_x_angle = 0, E_diff = 0;
		
	
	// Loop over entries
        for (int i = 0; i < nentries; i++) {
            tree -> GetEntry(i);

            // Cuts
            //cout << lagvalue_min_7C << endl;
            if (lagvalue_min_7C > chi2_cut) continue;
	    else if (deltaE > deltaE_cut) continue;
	    else if (angle_pi0gam12 > angle_cut) continue;
	    else if (betapi0 > GetFBeta(beta_cut, c0, c1, ppIM)) continue;

	    //cout << "angle_pi0gam12 = " << angle_pi0gam12 << endl;
	    /*
	    cout << "angle_cut = " << angle_cut << "\n"
		 << "lagvalue_min_7C = " << chi2_cut << "\n"
		 << "delta_cut = " << deltaE_cut << "\n";
	    */
	    
	    // Clean data
            if (TMath::IsNaN(E1) || TMath::IsNaN(E2) || TMath::IsNaN(E3)) continue;
            //if (!TMath::IsNaN(px1)) continue;

	    //if (TMath::IsNaN(angle_pi0gam12)) {
	    //  cout << angle_pi0gam12 << endl;
	    //}

	    //cout << betapi0 << endl;
	    //cout << deltaE << endl;
	    //cout << ppIM << endl;
	    
            
	    // Store tracks
	    double trk[2][4] = {
	      {ppl_E, ppl_px, ppl_py, ppl_pz},
	      {pmi_E, pmi_px, pmi_py, pmi_pz}
	    };

	    //cout << "(ppl_E, ppl_px, ppl_py, ppl_z) = (" << ppl_E << ", " << ppl_px << ", " << ppl_py << ", " << ppl_pz <<  ")" << endl;
              
	  
        
            // Store photons
            double photons[3][4] = {
                {E1, px1, py1, pz1},
                {E2, px2, py2, pz2},
                {E3, px3, py3, pz3}
            };

            //cout << "(E1, px1, py1, pz1) = (" << E1 << ", " << px1 << ", " << py1 << ", " << pz1 <<  ")" << endl;
            //cout << "(E2, px2, py2, pz2) = (" << E2 << ", " << px2 << ", " << py2 << ", " << pz2 <<  ")" << endl;
            //cout << "(E3, px3, py3, pz3) = (" << E3 << ", " << px3 << ", " << py3 << ", " << pz3 <<  ")" << endl;
            
            // All 3 possible pairs
            int pair_indicies[3][2] = {{0, 1}, {2, 0}, {1, 2}};

            /*
            cout << pair_indicies[0][0] << ", " << pair_indicies[0][1] << "\n"
                << pair_indicies[1][0] << ", " << pair_indicies[1][1] << "\n"
                << pair_indicies[2][0] << ", " << pair_indicies[2][1] << "\n\n";  
            */

            // Stor scores and pi0 masses for each pair
            double scores[3] = {0., 0., 0.};
            double masses[3] = {0., 0., 0.};

            // Calculate BDT score for each pair and store them
            for (int p = 0; p < 3; p++){

	      // Local variables for THIS pair only

	      int i_idx = pair_indicies[p][0];
	      int j_idx = pair_indicies[p][1];
	      //cout << "(i, j) = (" << i_idx << ", " << j_idx << ")\n"; 
	      
	      // Paired photon energies
	      double e1 = photons[i_idx][0];
	      double e2 = photons[j_idx][0];
	      //cout << "e1 = " << e1 << ", e2 = " << e2 << endl;
	      
	      // Found unpaired photon energy
	      int unpaired_idx = -1;
	      for (int k = 0; k < 3; k++) {
		if (k != i_idx && k != j_idx) {
		  unpaired_idx = k;
		  break;
		}
	      }
	      double e3 = photons[unpaired_idx][0];

	      // Fill histos
	      he1 -> Fill(e1);
	      he2 -> Fill(e2);
	      he3 -> Fill(e3);
	      //cout << "i_idx = " << i_idx << ", j_idx = " << j_idx << ", unpaired_idx = " << unpaired_idx << endl;
	      
	      // Calculate features
	      double m_gg_local = 0.0;
	      double opening_angle_local = 0.0;
	      double cos_theta_local = 0.0;
	      double E_asym_local = 0.0;
	      double e_min_x_angle_local = 0.0;
	      double E_diff_local = 0.0;
	      double asym_x_angle_local = 0.0;
	      
	      if (e1 >= energy_threshold && e2 >= energy_threshold) {
		// Invariant mass
		m_gg_local = inv_mass_4vector(i_idx, j_idx, photons);
		masses[p] = m_gg_local;
		
		//cout << m3pi << endl;
		
		//hm_gg -> Fill(m_gg_local);
		
		// cos_theta
		cos_theta_local = get_cos_theta(i_idx, j_idx, photons);
		hcos_theta -> Fill(cos_theta_local);
		
		// opening angle
		opening_angle_local = TMath::ACos(cos_theta_local);
		hopen_angle -> Fill(opening_angle_local);
		
		//cout << opening_angle_local << endl;
		
		// energy features
		E_asym_local = TMath::Abs(e1 - e2) / (e1 + e2 + 1e-10);
		E_asym_local = TMath::Max(0.0, TMath::Min(1.0, E_asym_local));
		E_diff_local = TMath::Abs(e1 - e2);
		e_min_x_angle_local = TMath::Min(e1, e2) * opening_angle_local;
		asym_x_angle_local = E_asym_local * opening_angle_local;
		
		//cout << E_asym_local << endl;
		//cout << e_min_x_angle_local << endl;
		hE_asym -> Fill(E_asym_local);
		he_min_x_angle -> Fill(e_min_x_angle_local);
		hE_diff -> Fill(E_diff_local);
		hasym_x_angle -> Fill(asym_x_angle_local);
		
	      }
	      
	      // Prepare features for BDT
	      std::vector<float> features = {
		(float)m_gg_local, (float)opening_angle_local, (float)cos_theta_local, 
		(float)E_asym_local, (float)e_min_x_angle_local, 
		(float)asym_x_angle_local, (float)E_diff_local,
		(float)e1, (float)e2, (float)e3	
	      };
	      
	      // Get BDT score
	      // Convert vector to tensor (1 event, n_features)
	      //TMVA::Experimental::RTensor<float> input(dummy.data(), {1, (size_t)n_features});
	      
	      TMVA::Experimental::RTensor<float> input_tensor(features.data(), {1, features.size()});
	      auto result = bdt.Compute(input_tensor);
	      scores[p] = result(0, 0);
	      //cout << "p: " << p << ", score: " << scores[p] << endl;
	      
            }
	    
            // Find the best pair (highest BDT score)
            int best_pair = 0;
            if(scores[1] > scores[best_pair]) best_pair = 1;
            if(scores[2] > scores[best_pair]) best_pair = 2;

            // Get the indices for the best pair
            int best_i = pair_indicies[best_pair][0];
            int best_j = pair_indicies[best_pair][1];

            // Find prompt photon (the one not in the best pair)
            int prompt_idx = -1;
            for (size_t k = 0; k < 3; k++)
            {
                if (k != best_i && k != best_j)
                {
                    prompt_idx = k;
                    break;
                }
                
            }
            
            //cout << "best pair indices: (" << best_i << ", " << best_j << "), prompt index: " << prompt_idx << endl;

            // Calculate pi0 4-vector
            pi0_pho1_idx = best_i;
            pi0_pho2_idx = best_j;
            prompt_pho_idx = prompt_idx;

            e1_bdt = photons[pi0_pho1_idx][0];
	    e2_bdt = photons[pi0_pho2_idx][0];
	    
            m_gg_bdt = inv_mass_4vector(pi0_pho1_idx, pi0_pho2_idx, photons);
	    //m_gg = inv_mass_4vector(0, 1, photons);
	    m_gg = masses[0];
	    //cout << m_gg << ", " << masses[0] << endl;
	    m3pi = inv_3pimass_4vector(0, 1, photons, trk);
            m3pi_bdt = inv_3pimass_4vector(pi0_pho1_idx, pi0_pho2_idx, photons, trk);
            	 
            //m_gg_bdt = inv_mass_4vector(0, 1, photons);
            hM_gg -> Fill(m_gg);
	    hM3pi -> Fill(m3pi);
            hE1 -> Fill(photons[0][0]);
            hE2 -> Fill(photons[1][0]);
            
	    bdt_score = scores[best_pair];
            
            //cout << "m_gg_bdt = " << m_gg_bdt << endl;

	    //out_event = event_id;
	    //cout << out_event << endl;
            outtree -> Fill();
            
	    // BDT selection
	    if (scores[best_pair] > 0.5) {
	      n_found ++;
	      hE1_BDT_good -> Fill(e1_bdt);
	      hE1_BDT_best -> Fill(e1_bdt);
	      
	      hE2_BDT_good -> Fill(e2_bdt);
	      
	      hM_gg_BDT_good -> Fill(m_gg_bdt);
	      hM3pi_BDT_good -> Fill(m3pi_bdt);
	      hM3pi_BDT_best -> Fill(m3pi_bdt);

	      bdt_indx = 1;
            }
	    else {
	      hE1_BDT_bad -> Fill(e1_bdt);
	      hE2_BDT_bad -> Fill(e2_bdt);
	      hM_gg_BDT_bad -> Fill(m_gg_bdt);
	      hM3pi_BDT_bad -> Fill(m3pi_bdt);

	      bdt_indx = 0;
	    }
            
            // KLOE selection
            if (recon_indx == 2 && bkg_indx == 1){//  true pi0 gg
	        hE1_good -> Fill(photons[0][0]);
		hE2_good -> Fill(photons[1][0]);
		hM_gg_good -> Fill(m_gg);
		hM3pi_good -> Fill(m3pi);
            
		evnt_good += 1;
		kloe_indx = 1;

		//cout << recon_indx << endl;
            }
            else{// false pi0 gg
                hE1_bad -> Fill(photons[0][0]);
		hE2_bad -> Fill(photons[1][0]);
		hM_gg_bad -> Fill(m_gg);
		hM3pi_bad -> Fill(m3pi);
            
		evnt_bad += 1;
		kloe_indx = 0;

		//cout << recon_indx << endl;
            } 

	    evnt_KLOE += 1;
    
            //if (i > 10) break; // Fisrt 10 events



        }


	hM3pi -> Write();
	hM3pi_good -> Write();
	hM3pi_bad -> Write();
	hM3pi_BDT_good -> Write();
	hM3pi_BDT_bad -> Write();
	hM3pi_BDT_best -> Write();
	
	outfile -> Write();
	outfile -> Close();
	file -> Close();

	cout << "\n pi0 finding complete!" << endl;
    }    

    if (evnt_KLOE > 0) {
      TCanvas *cv0 = new TCanvas("c1", "BDT Selection (" + ch_nm + ")", 1200, 600);
      cv0 -> SetLeftMargin(0.1);
      cv0 -> SetBottomMargin(0.1);//0.007

      cv0 -> Divide(2, 1);  // [0] columns, [1] rows
      cv0 -> cd(1);

      Double_t ymax_e1 = hE1 -> GetBinContent(hE1 -> GetMaximumBin());
      //cout << ymax_e1 << endl;
    
      TPaveText *pt1 = new TPaveText(0.11, 0.87, 0.80, 0.89, "NDC");
    
      PteAttr(pt1); pt1 -> SetTextSize(0.03); pt1 -> SetTextColor(kBlack);
    
      pt1 -> AddText(Form("Events=%d, BDT Selected=%d, Discarded=%d", evnt_KLOE, evnt_good, evnt_bad));
    
    
      format_h(hE1, 1, 2);
      format_h(hE1_good, 4, 1);
      format_h(hE1_bad, 2, 1);

      formatfill_h(hE1_BDT_good, 3, 3001);
      formatfill_h(hE1_BDT_bad, 2, 3001);
      format_h(hE1_BDT_best, 3, 2);
 
      hE1 -> GetYaxis() -> SetNdivisions(505);
      hE1 -> GetYaxis() -> SetRangeUser(0.1, ymax_e1 * 1.2); 
      hE1 -> GetXaxis() -> SetTitle("E_{1} [MeV]");
      hE1 -> GetXaxis() -> CenterTitle();
      hE1 -> GetXaxis() -> SetTitleSize(0.04);
      //hE1 -> GetXaxis() -> SetTitleOffset(1.0);
      //hE1 -> GetXaxis() -> SetLabelOffset(0.01);
      //hE1 -> GetXaxis() -> SetLabelSize(0.05);//0.03
      
      
      hE1 -> Draw();
      hE1_good -> Draw("Same");
      hE1_bad -> Draw("Same");
      hE1_BDT_good -> Draw("Same");
      hE1_BDT_bad -> Draw("Same");
      //hE1_BDT_best -> Draw("Same");
      
      pt1 -> Draw("Same");
      
      TLegend *legd_cv = new TLegend(0.5, 0.5, 0.9, 0.85);
      
      legd_cv -> SetTextFont(132);
      legd_cv -> SetFillStyle(0);
      legd_cv -> SetBorderSize(0);
      legd_cv -> SetNColumns(1);
      
      legd_cv -> AddEntry(hE1, "#chi^{2}_{m_{#gamma#gamma}} Sum", "l");
      legd_cv -> AddEntry(hE1_good, "#chi^{2}_{m_{#gamma#gamma}} Selected", "l");
      legd_cv -> AddEntry(hE1_bad, "#chi^{2}_{m_{#gamma#gamma}} Discarded", "l");
      legd_cv -> AddEntry(hE1_BDT_good, "BDT Selected", "f");
      legd_cv -> AddEntry(hE1_BDT_bad, "BDT Discarded", "f");
      
      legd_cv -> Draw("Same");
      
      legtextsize(legd_cv, 0.04);
      
      //
      cv0 -> cd(2);
      
      format_h(hE2, 1, 2);
      format_h(hE2_good, 4, 1);
      format_h(hE2_bad, 2, 1);
      
      formatfill_h(hE2_BDT_good, 3, 3001);
      formatfill_h(hE2_BDT_bad, 2, 3001);
      //format_h(hE2_BDT_best, 3, 2);
      
      hE2 -> GetYaxis() -> SetNdivisions(505);
      hE2 -> GetYaxis() -> SetRangeUser(0.1, ymax_e1 * 1.2); 
      hE2 -> GetXaxis() -> SetTitle("E_{2} [MeV]");
      hE2 -> GetXaxis() -> CenterTitle();
      hE2 -> GetXaxis() -> SetTitleSize(0.04);
      
      
      hE2 -> Draw();
      hE2_good -> Draw("Same");
      hE2_bad -> Draw("Same");
      hE2_BDT_good -> Draw("Same");
      hE2_BDT_bad -> Draw("Same");
      
      TCanvas *cv01 = new TCanvas("cv01", "BDT Selection (" + ch_nm + ")", 1200, 600);
      cv01 -> SetLeftMargin(0.1);
      cv01 -> SetBottomMargin(0.1);//0.007
      
      cv01 -> Divide(2, 1);  // [0] columns, [1] rows
      cv01 -> cd(1);
      
      double ymax_m_gg = hM_gg -> GetBinContent(hM_gg -> GetMaximumBin());
      
      format_h(hM_gg, 1, 2);
      format_h(hM_gg_good, 4, 1);
      format_h(hM_gg_bad, 2, 1);
      
      format_h(hM_gg_BDT_good, 3, 2);
      formatfill_h(hM_gg_BDT_bad, 2, 3001);
      
      hM_gg -> GetYaxis() -> SetNdivisions(505);
      hM_gg -> GetYaxis() -> SetRangeUser(0.1, ymax_m_gg * 1.2); 
      hM_gg -> GetXaxis() -> CenterTitle();
      hM_gg -> GetXaxis() -> SetTitleSize(0.04);
      hM_gg -> GetXaxis() -> SetTitle("M(#gamma_{1}#gamma_{2}) [MeV/c^{2}]");
      hM_gg -> GetXaxis() -> CenterTitle();
      hM_gg -> GetXaxis() -> SetTitleSize(0.04);
      
      hM_gg -> Draw();
      hM_gg_good -> Draw("Same");
      hM_gg_bad -> Draw("Same");
      hM_gg_BDT_good -> Draw("Same");
      hM_gg_BDT_bad -> Draw("Same");
      gPad->SetLogy(1); 

      //
      cv01 -> cd(2);
      
      double ymax_m3pi = hM3pi -> GetBinContent(hM3pi -> GetMaximumBin());
    
      format_h(hM3pi, 1, 2);
      format_h(hM3pi_good, 4, 1);
      format_h(hM3pi_bad, 2, 1);
      
      formatfill_h(hM3pi_BDT_good, 3, 3001);
      formatfill_h(hM3pi_BDT_bad, 2, 3001);
      //format_h(hM3pi_BDT_best, 3, 2);
      
      hM3pi -> GetYaxis() -> SetNdivisions(505);
      hM3pi -> GetYaxis() -> SetRangeUser(0.1, ymax_m3pi * 1.5); 
      hM3pi -> GetXaxis() -> SetTitle("M_{3#pi} [MeV/c^{2}]");
      hM3pi -> GetXaxis() -> CenterTitle();
      hM3pi -> GetXaxis() -> SetTitleSize(0.04);
      
      hM3pi -> Draw();
      hM3pi_good -> Draw("Same");
      hM3pi_bad -> Draw("Same");
      hM3pi_BDT_good -> Draw("Same");
      hM3pi_BDT_bad -> Draw("Same");
      gPad->SetLogy(1); 

      TLegend *legd1 = new TLegend(0.5, 0.6, 0.9, 0.9);
      
      legd1 -> SetTextFont(132);
      legd1 -> SetFillStyle(0);
      legd1 -> SetBorderSize(0);
      legd1 -> SetNColumns(1);
      
      legd1 -> AddEntry(hM3pi, "#chi^{2}_{m_{#gamma#gamma}} Sum", "l");
      legd1 -> AddEntry(hM3pi_good, "#chi^{2}_{m_{#gamma#gamma}} Selected", "l");
      legd1 -> AddEntry(hM3pi_bad, "#chi^{2}_{m_{#gamma#gamma}} Discarded", "l");
      legd1 -> AddEntry(hM3pi_BDT_good, "BDT Selected", "f");
      legd1 -> AddEntry(hM3pi_BDT_bad, "BDT Discarded", "f");
      
      legd1 -> Draw("Same");

      legtextsize(legd1, 0.03);

      cv0 -> SaveAs("../plots_bdt/energy_spectra_comparison.pdf");
      cv01 -> SaveAs("../plots_bdt/invariant_mass_comparison.pdf");
      
      delete he1; delete he2; delete he3;
      delete hm_gg; delete hcos_theta; delete hopen_angle;
      delete hE_asym; delete he_min_x_angle; delete hE_diff; delete hasym_x_angle;
      delete hM3pi_BDT; delete hM3pi_BDT_good; delete hM3pi_BDT_bad; delete hM3pi_BDT_best;
      delete hM_gg_BDT; delete hM_gg_BDT_good; delete hM_gg_BDT_bad;
      delete hE1_BDT_good; delete hE1_BDT_bad; delete hE1_BDT_best;
      delete hE2_BDT_good; delete hE2_BDT_bad;
      delete hE1; delete hE1_good; delete hE1_bad;
      delete hE2; delete hE2_good; delete hE2_bad;
      delete hM_gg; delete hM_gg_good; delete hM_gg_bad;
      delete hM3pi; delete hM3pi_good; delete hM3pi_bad;

      delete cv0;
      delete cv01;
      
    }
}
