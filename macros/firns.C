

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
  TCanvas * c = new TCanvas("firns","Firns", 800,1200); 
  c->Divide(1,3); 

  gStyle->SetOptFit(11); 
  gStyle->SetStatX(0.49); 
  gStyle->SetStatY(0.88); 
  gStyle->SetTitleH(0.1); 
  gStyle->SetPalette(kRainBow); 
  gStyle->SetStatStyle(0); 
  gStyle->SetLegendBorderSize(0); 
  gStyle->SetStatBorderSize(0); 
  gStyle->SetLineScalePS(1); 
  gStyle->SetFitFormat(".2f"); 
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

  hawley_neutron->SetMarkerStyle(2); 
  hawley_icecore->SetMarkerStyle(3); 
  hawley_snowpit->SetMarkerStyle(5); 
  alley->SetMarkerStyle(29); 
  gisp2->SetMarkerStyle(34); 
  mg->Add(hawley_neutron,"PLX"); 
  mg->Add(hawley_icecore,"PLX"); 
  mg->Add(hawley_snowpit,"PLX"); 
  mg->Add(alley,"PLX"); 
  mg->Add(gisp2,"PLX"); 
  mg->SetTitle("Data and Fits"); 
  gStyle->SetMarkerSize(2); 

  c->cd(1); 
  mg->Draw("a pmc plc"); 
  mg->GetXaxis()->SetLabelSize(0.045); 
  mg->GetYaxis()->SetLabelSize(0.045); 
  mg->GetXaxis()->SetTitleSize(0.06); 
  mg->GetYaxis()->SetTitleSize(0.06); 
  mg->GetXaxis()->SetTitleOffset(0.7); 
  mg->GetYaxis()->SetTitleOffset(0.6); 


  TF1 * arthern = new TF1("arthern", double_exp,0,500,3); 
  arthern->SetTitle("Arthern'13 Model");
  arthern->SetParameters(280,27,42); 
  arthern->SetLineColor(3); 
  arthern->SetMarkerColor(3); 

  arthern->Draw("lsame"); 

  TGraph * critical = new TGraph(2); 
  critical->SetPoint(0, 0, 550); 
  critical->SetPoint(1, 120, 550); 
  critical->SetTitle("#rho_{c}"); 
  critical->SetLineStyle(2); 
  critical->Draw("lsame"); 


  TF1 * fit = new TF1("fit", double_exp,0,500,3); 
  fit->SetTitle("Best Fit"); 
  fit->SetParameters(280,27,42); 
  fit->SetLineColor(2); 
  fit->SetMarkerColor(2); 
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
//  neutron_only->Draw("lsame"); 
  alley->Fit(icecore_only,"N"); 
//  icecore_only->Draw("lsame"); 
 
  gPad->BuildLegend(0.6,0.15,0.85,0.45,"","lp"); 
  gPad->SetRightMargin(0.02); 

  gPad->SetGridx(); 
  gPad->SetGridy(); 
  gPad->SetTickx(); 
  gPad->SetTicky(); 

  c->cd(2); 

  iceprop::MultiDatasetFit multi_fit; 
  iceprop::Firn * Hawley =  new iceprop::DensityTableFirn ("data/hawley08_neutron.txt", &multi_fit); 
  iceprop::Firn * Alley =  new iceprop::DensityTableFirn ("data/alley_koci.txt", &multi_fit); 
  TGraph *Gh = Hawley->makeGraph(1000); 
  Gh->SetLineColor(31); 
  Gh->SetTitle("Hawley-derived Model"); 
  Gh->Draw("alp"); 
  Gh->GetXaxis()->SetRangeUser(0,120); 
  Gh->GetXaxis()->SetLabelSize(0.049); 
  Gh->GetYaxis()->SetLabelSize(0.049); 
  Gh->GetXaxis()->SetTitleSize(0.06); 
  Gh->GetYaxis()->SetTitleSize(0.06); 
  Gh->GetXaxis()->SetTitleOffset(0.7); 
  Gh->GetYaxis()->SetTitleOffset(0.6); 

  gPad->SetRightMargin(0.02); 


  gPad->SetGridx(); 
  gPad->SetGridy(); 
  gPad->SetTickx(); 
  gPad->SetTicky(); 
  critical->Draw("lsame"); 
  c->cd(3); 
  TGraph *Ga = Alley->makeGraph(1000); 
  Ga->SetTitle("Alley-derived Model"); 
  Ga->SetLineColor(32); 
  Ga->Draw("alp"); 
  Ga->GetXaxis()->SetLabelSize(0.045); 
  Ga->GetYaxis()->SetLabelSize(0.045); 
  Ga->GetXaxis()->SetTitleSize(0.06); 
  Ga->GetYaxis()->SetTitleSize(0.06); 
  Ga->GetXaxis()->SetTitleOffset(0.7); 
  Ga->GetYaxis()->SetTitleOffset(0.6); 


  Ga->GetXaxis()->SetRangeUser(0,120); 

  gPad->SetGridx(); 
  gPad->SetGridy(); 
  gPad->SetRightMargin(0.02); 
  gPad->SetTickx(); 
  gPad->SetTicky(); 
  critical->Draw("lsame"); 

  fit->Print(); 
  c->SaveAs("firn_fits.pdf"); 

}
