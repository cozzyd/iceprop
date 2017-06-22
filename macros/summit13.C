/* simulation of summit13 code */ 


const double ft_to_mtr = 0.3048; 
const double receiver_depth = -ft_to_mtr; 



void summit13(int depth_in_ft = 100, double f = 0.25, double w = 0.2) 
{

  /* just make the output dir */ 
  TString dir; 
  dir.Form("summit13_%d", depth_in_ft); 

  TString cmd; 
  cmd.Form("mkdir -p %s", dir.Data()); 
  system(cmd.Data()); 


  iceprop::ArthernFirn frn ;

  iceprop::SimGeometry g; 
  g.max_depth =depth_in_ft * ft_to_mtr * 1.5; 
  g.sky_height=20; 
  g.max_r = 400; 

  iceprop::GaussianPulseSource s(0, -depth_in_ft*ft_to_mtr, f,w); 

  iceprop::Sim sim(&frn, &g, &s); 

//  sim.makeHist(meep::Dielectric)->Draw("colz"); 


  /* for timing */ 
  sim.addTimeDomainMeasurement( meep::Ez, 0, -depth_in_ft*ft_to_mtr); 
  sim.addTimeDomainMeasurement( meep::Ez, 200 * ft_to_mtr, receiver_depth); 
  sim.addTimeDomainMeasurement( meep::Ez, 480 * ft_to_mtr,receiver_depth ); 
  sim.addTimeDomainMeasurement( meep::Ez, 1000*ft_to_mtr, receiver_depth); 
  sim.addTimeDomainMeasurement( meep::Ez, 1050*ft_to_mtr, receiver_depth); 

  iceprop::StepOutput o; 
  o.out_dir = dir.Data(); 
  o.format = iceprop::O_PNG; 
  sim.addStepOutput(o); 
  sim.run(2000); 

  TString fname;
  fname.Form("%s/time_graphs.root", dir.Data());
  TFile fout(fname.Data(),"RECREATE"); 

  int n = (int)  sim.getMeasurements().size(); 
  TCanvas * ms = new TCanvas("measurements","Measurements",1000,800); 
  ms->Divide(3,2); 
  for (int i = 0; i< n; i++)
  {
    ms->cd(i+1); 
    TGraph * greal = sim.getMeasurements()[i].makeGraph(iceprop::Real); 
    TGraph * gimag = sim.getMeasurements()[i].makeGraph(iceprop::Imag); 
    TGraph * gmag = sim.getMeasurements()[i].makeGraph(iceprop::Mag); 
    greal->Draw("al"); 
    gimag->SetLineColor(2); 
    gimag->Draw("lsame");
    gmag->SetLineColor(3); 
    gmag->Draw("lsame");

    greal->Write(); 
    gimag->Write(); 
    gmag->Write(); 
  }

  fname.Form("%s/time_graphs.pdf", dir.Data()); 
  ms->SaveAs(fname.Data()); 


}
