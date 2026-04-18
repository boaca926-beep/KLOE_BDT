double pi = TMath::Pi();

void legtextsize(TLegend* l, Double_t size) {
  for(int i = 0 ; i < l -> GetListOfPrimitives() -> GetSize() ; i++) {
    TLegendEntry *header = (TLegendEntry*)l->GetListOfPrimitives()->At(i);
    header->SetTextSize(size);
  }
}

//
void PteAttr(TPaveText *pt) {

  pt -> SetTextSize(0.04);
  pt -> SetFillColor(0);
  pt -> SetTextAlign(12);
  pt -> SetBorderSize(0);
}

void format_h(TH1D* h, Int_t linecolor, Int_t width) {
  h->SetLineColor(linecolor);
  h->SetMarkerColor(linecolor);
  //cout << "histo format" << endl;
  h->SetLineWidth(width);
}

void formatfill_h(TH1D* h, Int_t fillcolor, Int_t fillstyle) {
  h -> SetFillStyle(fillstyle);
  h -> SetFillColor(fillcolor);
  h -> SetLineColor(0);
}

double get_cos_theta(int i_idx, int j_idx, double photons[3][4]){

    double p1_mag = TMath::Sqrt(TMath::Max(0.0, 
        photons[i_idx][1] * photons[i_idx][1] + 
        photons[i_idx][2] * photons[i_idx][2] +
        photons[i_idx][3] * photons[i_idx][3]) 
    );

    double p2_mag = TMath::Sqrt(TMath::Max(0.0,
        photons[j_idx][1] * photons[j_idx][1] + 
        photons[j_idx][2] * photons[j_idx][2] +
        photons[j_idx][3] * photons[j_idx][3])
    );

    double dot = photons[i_idx][1] * photons[j_idx][1] + 
    photons[i_idx][2] * photons[j_idx][2] +
    photons[i_idx][3] * photons[j_idx][3];

    double cos_theta = dot / (p1_mag * p2_mag + 1e-10);
    cos_theta = TMath::Max(-1.0, TMath::Min(1.0, cos_theta));

    return cos_theta;
}

//
double inv_mass_4vector(int i_idx, int j_idx, double photons[3][4]){
    double inv_mass = 0;

    // Calculate total energy and momentum
    double e = photons[i_idx][0] + photons[j_idx][0];
    
    double px = photons[i_idx][1] + photons[j_idx][1];
    double py = photons[i_idx][2] + photons[j_idx][2];
    double pz = photons[i_idx][3] + photons[j_idx][3];

    double mass_square = e*e - (px*px + py*py + pz*pz);
    
    if (mass_square < 0 && mass_square > -1e-10) {
        return 0;
    }
    
    return (mass_square > 0) ? TMath::Sqrt(mass_square) : 0;

}

//
double inv_3pimass_4vector(int i_idx, int j_idx, double photons[3][4], double trk[2][4]){
    double inv_mass = 0;

    // Calculate total energy and momentum
    double e = photons[i_idx][0] + photons[j_idx][0] + trk[0][0] + trk[1][0];
    
    double px = photons[i_idx][1] + photons[j_idx][1] + trk[0][1] + trk[1][1];
    double py = photons[i_idx][2] + photons[j_idx][2] + trk[0][2] + trk[1][2];
    double pz = photons[i_idx][3] + photons[j_idx][3] + trk[0][3] + trk[1][3];

    double mass_square = e*e - (px*px + py*py + pz*pz);
    
    if (mass_square < 0 && mass_square > -1e-8) {
        return 0;
    }
    
    return (mass_square > 0) ? TMath::Sqrt(mass_square) : 0;

}


