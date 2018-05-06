/* simulation of summit'13 measurements
 *
 *
 * */ 


const double ft_to_mtr = 0.3048; 
const double receiver_depth = -ft_to_mtr; 



void layers(int depth_in_ft = 20, int firn = 0, double f = 0.35, double w = 0.2, bool movie = true) 
{

  /* define and make the output dir */ 
  TString dir; dir.Form("layers/depth_%d_f_%g_w_%g_firn_%d", depth_in_ft,f,w, firn); 
  TString cmd; cmd.Form("mkdir -p %s", dir.Data()); system(cmd.Data()); 

  /* use "arthern model" */
  iceprop::ArthernFirn arthern ;
  iceprop::Firn * frn = 0;

  if (firn == 0) frn = &arthern; 
  if (firn == 1 || firn == 2) 
  {
    double z[] = {-3,-6,-9}; 
    double A1[] = {150,100,50}; 
    double A2[] = {-150,-100,50} ; 
    double sigma[] = {0.3,0.2,0.1}; 
    frn = new iceprop::PerturbedFirn(arthern, 3,z,firn == 1 ?A1 : A2,sigma); 
  }

  /* Define simulation geometry */ 
  iceprop::SimGeometry g; 
  g.max_depth =TMath::Max(3*depth_in_ft,30); 
  g.sky_height=20; 
  g.output_skip_factor=10; 
  g.resolution=20; 
  g.max_r = 100; 


  /* Define the source */ 
  iceprop::GaussianPulseSource s(0, -depth_in_ft*ft_to_mtr, f,w); 

  /* Initialize the simulation */ 
  iceprop::Sim sim(frn, &g, &s); 


  /* If movie, add step output*/ 
  if (movie) 
  {
    iceprop::StepOutput o; 
    o.out_dir = dir.Data(); 
    o.format = iceprop::O_ROOT | iceprop::O_PNG; 
    sim.addStepOutput(o); 
  }

  /* uncomment to track global maximum of Ez (very slow!) */ 
  sim.trackGlobalMaximum(meep::Ez); 
  sim.trackGlobalIntegral(meep::Ez); 

  /* run the simulation for 2000 ns*/ 
  sim.run(300); 

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

  sim.getMaximums()[0].max->SetDirectory(&fout);
  sim.getMaximums()[0].max->Write(); 
  sim.getMaximums()[0].tmax->SetDirectory(&fout); 
  sim.getMaximums()[0].tmax->Write(); 
  sim.getIntegrals()[0].integ->SetDirectory(&fout);
  sim.getIntegrals()[0].integ->Write(); 
  sim.getIntegrals()[0].tfirst->SetDirectory(&fout); 
  sim.getIntegrals()[0].tfirst->Write(); 
  fname.Form("%s/time_graphs.pdf", dir.Data()); 
  ms->SaveAs(fname.Data()); 


}
