

const double rho_c = 550; 
const double rho_i = 917; 

double double_exp(double *x, double * par) 
{

  double rho_s = par[0]; 
  double L1 = par[1]; 
  double L2 = par[2]; 

  double z_c = L1 * log( (rho_i - rho_s) / (rho_i - rho_c)); 

  double z = *x; 
  if (z <=z_c) 
    return rho_i - (rho_i - rho_s) *exp(-z/L1); 
  else 
    return rho_i - (rho_i - rho_c) * exp(-(z-z_c)/L2); 
}


void setErrors(TGraph * g, double ex, double ey) 
{
  for (int i= 0; i < g->GetN(); i++)
  {
    g->GetEY()[i] = ey; 
    g->GetEX()[i] = ex; 
  }

}

void firns(bool use_hawley_06 = false) 
{

  gStyle->SetOptFit(1); 
  TGraph * hawley_06 = new TGraphErrors("data/hawley06.txt"); 
  hawley_06->SetTitle("Hawley'06 Neutron Data ( ??? )"); 
  setErrors(hawley_06, 0.2, 20);

  TGraph * hawley_neutron = new TGraphErrors("data/hawley08_neutron.txt"); 
  setErrors(hawley_neutron, 0.2, 20);
  TGraph * hawley_icecore = new TGraphErrors("data/hawley08_icecore.txt"); 
  setErrors(hawley_icecore, 0.1, 10);
  TGraph * hawley_snowpit = new TGraphErrors("data/hawley08_snowpit.txt"); 
  setErrors(hawley_snowpit, 0.1, 5);

  hawley_neutron->SetTitle("Hawley'08 Neutron Data"); 
  hawley_icecore->SetTitle("Hawley'08 IceCore Data"); 
  hawley_snowpit->SetTitle("Hawley'08 Snowpit Data"); 

  TGraph * alley = new TGraphErrors("data/alley_koci.txt"); 
  setErrors(alley, 0.5, 10); 
  alley->SetTitle("Alley & Koci '88, IceCore Data"); 

  TGraph * gisp2 = new TGraphErrors("data/gisp2.txt"); 
  setErrors(gisp2,0.1,10); 
  gisp2->SetTitle("GISP2"); 

  TMultiGraph *  mg = new TMultiGraph; 

  if (use_hawley_06) 
    mg->Add(hawley_06); 

  mg->Add(hawley_neutron); 
  mg->Add(hawley_icecore); 
  mg->Add(hawley_snowpit); 
  mg->Add(alley); 
  mg->Add(gisp2); 


  mg->Draw("a pmc plc"); 
  TF1 * arthern = new TF1("arthern", double_exp,0,500,3); 
  arthern->SetTitle("Arthern Nominal");
  arthern->SetParameters(280,27,42); 
  arthern->SetLineColor(3); 

  arthern->Draw("lsame"); 

  TGraph * critical = new TGraph(2); 
  critical->SetPoint(0, 0, 550); 
  critical->SetPoint(1, alley->GetX()[alley->GetN()-1], 550); 
  critical->SetTitle("#rho_{c}"); 
  critical->SetLineStyle(2); 
  critical->Draw("lsame"); 


  TF1 * fit = new TF1("fit", double_exp,0,500,3); 
  fit->SetTitle("Best Fit"); 
  fit->SetParameters(280,27,42); 
  fit->SetLineColor(2); 
  fit->SetParName(0,"Best Fit #rho_{s}"); 
  fit->SetParName(1,"Best Fit L_{1}"); 
  fit->SetParName(2,"Best Fit L_{2}"); 

  TF1 * neutron_only = new TF1("neutron_only", double_exp,0,500,3); 
  neutron_only->SetTitle(use_hawley_06 ? "Fit to hawley'06" : "Fit to hawley'08 "); 
  neutron_only->SetLineColor(11); 
  neutron_only->SetParameters(280,27,42); 

  TF1 * icecore_only = new TF1("icecore_only", double_exp,0,500,3); 
  icecore_only->SetTitle("Fit to alley'88 "); 
  icecore_only->SetLineColor(12); 
  icecore_only->SetParameters(280,27,42); 


  mg->GetXaxis()->SetTitle("z (m)"); 
  mg->GetYaxis()->SetTitle("#rho (kg/m^{3})"); 
  mg->GetXaxis()->SetRangeUser(0,120); 
  mg->Fit(fit,""); 
  fit->Draw("lsame"); 

  (use_hawley_06 ? hawley_06 : hawley_neutron)->Fit(neutron_only,"N"); 
  neutron_only->Draw("lsame"); 
  alley->Fit(icecore_only,"N"); 
  icecore_only->Draw("lsame"); 
 
  fit->Print(); 
  gPad->BuildLegend(0.5,0.15,0.8,0.45,"","lp"); 

}
