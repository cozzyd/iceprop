/* simulation of summit'13 measurements
 *
 *
 * */ 


const double ft_to_mtr = 0.3048; 
const double receiver_depth = -ft_to_mtr; 



void summit13(int depth_in_ft = 200, double f = 0.35, double w = 0.2, bool movie = true) 
{

  /* define and make the output dir */ 
  TString dir; dir.Form("out/summit13_%d", depth_in_ft); 
  TString cmd; cmd.Form("mkdir -p %s", dir.Data()); system(cmd.Data()); 

  /* use "arthern model" */
  iceprop::ArthernFirn frn ;

  /* Define simulation geometry */ 
  iceprop::SimGeometry g; 
  g.max_depth =(depth_in_ft+50) * ft_to_mtr; 
  g.sky_height=20; 
  g.max_r = 400; 


  /* Define the source */ 
  iceprop::GaussianPulseSource s(0, -depth_in_ft*ft_to_mtr, f,w); 

  /* Initialize the simulation */ 
  iceprop::Sim sim(&frn, &g, &s); 

  /* add time domain measurements */ 
  sim.addTimeDomainMeasurement( meep::Ez, 1*ft_to_mtr, -depth_in_ft*ft_to_mtr); 
  sim.addTimeDomainMeasurement( meep::Ez, 200 * ft_to_mtr, receiver_depth); 
  sim.addTimeDomainMeasurement( meep::Ez, 480 * ft_to_mtr,receiver_depth ); 
  sim.addTimeDomainMeasurement( meep::Ez, 1000*ft_to_mtr, receiver_depth); 
  sim.addTimeDomainMeasurement( meep::Ez, 1050*ft_to_mtr, receiver_depth); 


  /* If movie, add step output*/ 
  if (movie) 
  {
    iceprop::StepOutput o; 
    o.out_dir = dir.Data(); 
    o.format = iceprop::O_ROOT | iceprop::O_PNG; 
    sim.addStepOutput(o); 
  }

  /* uncomment to track global maximum of Ez (very slow!) */ 
//   sim.trackGlobalMaximum(meep::Ez); 

  /* run the simulation for 2000 ns*/ 
  sim.run(2000); 

  /* that's it for running the simulation, the rest plots and saves some of the stuff we made*/ 
  
  TString fname; fname.Form("%s/time_graphs.root", dir.Data());
  TFile fout(fname.Data(),"RECREATE"); 

  /* how many measurements did I give? */ 
  int n = (int)  sim.getMeasurements().size(); 

  TCanvas * ms = new TCanvas("measurements","Measurements",1000,800); 
   ms->SetLogy(); 
  for (int i = 0; i< n; i++)
  {
    TGraph * greal = sim.getMeasurements()[i].makeGraph(iceprop::Real); 
    TGraph * gimag = sim.getMeasurements()[i].makeGraph(iceprop::Imag); 
    TGraph * gmag = sim.getMeasurements()[i].makeGraph(iceprop::Mag); 
    greal->Write(); 
    gimag->Write(); 
    gmag->Write(); 

    gmag->SetLineColor(i+1); 
    gmag->Draw(i == 0 ? "al" : "lsame");
  }

  sim.getMaximums()[0].max->Write(); 
  sim.getMaximums()[0].tmax->Write(); 
  fname.Form("%s/time_graphs.pdf", dir.Data()); 
  ms->SaveAs(fname.Data()); 


}
