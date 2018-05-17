/* simulation of summit'13 measurements
 *
 *
 * */ 


const double ft_to_mtr = 0.3048; 
const double receiver_depth = -ft_to_mtr; 



void summit13(int depth_in_ft = 1, int firn = 0, double f = 0.35, double w = 0.15, bool movie = true) 
{

  /* define and make the output dir */ 
  TString dir; dir.Form("summit13_out/depth_%d_f_%g_w_%g_firn_%d", depth_in_ft,f,w, firn); 
  TString cmd; cmd.Form("mkdir -p %s", dir.Data()); system(cmd.Data()); 

  /* use "arthern model" */
  iceprop::ArthernFirn arthern ;
  iceprop::Firn * frn = 0;

  if (firn == 0) frn = &arthern; 
  if (abs(firn) == 1) 
  {

    double z = -3; 
    double A = firn >0 ? 100 : -100; 
    double sigma = 0.2; 
    frn = new iceprop::PerturbedFirn(arthern, 1,&z,&A,&sigma); 
  }
  if (abs(firn) == 2) 
  {
    double zs[] = {-2,-4}; 
    double As[] = {50,50}; 

    if (firn < 0)
    {
      for (int i = 0; i < sizeof(As)/sizeof(*As); i++) As[i] = -As[i]; 
    }

    double sigmas[] = {0.2,0.2}; 
    frn = new iceprop::PerturbedFirn(arthern, 2,zs,As,sigmas); 
  }

  if (abs(firn) == 3) 
  {
    double zs[] = {-2,-4,-10}; 
    double As[] = {50,50,25}; 
    if (firn < 0)
    {
      for (int i = 0; i < sizeof(As)/sizeof(*As); i++) As[i] = -As[i]; 
    }
    double sigmas[] = {0.2,0.2,0.2}; 
    frn = new iceprop::PerturbedFirn(arthern, 3,zs,As,sigmas); 
  }



  /* Define simulation geometry */ 
  iceprop::SimGeometry g; 
  g.max_depth =(TMath::Max(depth_in_ft+100,400) * ft_to_mtr); 
  g.sky_height=40; 
  g.resolution=20; // 5cm, enough for 500 MHz 
  g.courant_factor = 0.5; 
  g.pml_size = 10; 
  g.output_skip_factor=10; 
  g.max_r = 400; 

  /* Define the source */ 
  iceprop::Source * s = 0 ; 
  if ( f> 0) 
  {
    s = new iceprop::GaussianPulseSource(0, -depth_in_ft*ft_to_mtr, f,2*w); 
  }
  else
  {
    s = new iceprop::ButterworthSource(0, -depth_in_ft * ft_to_mtr,  -f,w); 
  }

  /* Initialize the simulation */ 
  iceprop::Sim sim(frn, &g, s); 

  /* add time domain measurements */ 
  sim.addTimeDomainMeasurement( meep::Ez, 0, -depth_in_ft*ft_to_mtr); 
  sim.addTimeDomainMeasurement( meep::Ez, 1, -depth_in_ft*ft_to_mtr); 
  sim.addTimeDomainMeasurement( meep::Ez, 100 * ft_to_mtr, receiver_depth); 
  sim.addTimeDomainMeasurement( meep::Ez, 200 * ft_to_mtr, receiver_depth); 
  sim.addTimeDomainMeasurement( meep::Ez, 400 * ft_to_mtr, receiver_depth); 
  sim.addTimeDomainMeasurement( meep::Ez, 480 * ft_to_mtr,receiver_depth ); 
  sim.addTimeDomainMeasurement( meep::Ez, 1000*ft_to_mtr, receiver_depth); 
  sim.addTimeDomainMeasurement( meep::Ez, 1050*ft_to_mtr, receiver_depth); 
  sim.addTimeDomainMeasurement( meep::Ez, 50*ft_to_mtr, -depth_in_ft); 
  sim.addTimeDomainMeasurement( meep::Ez, 500*ft_to_mtr, -depth_in_ft); 
  sim.addTimeDomainMeasurement( meep::Ez, 1000*ft_to_mtr, -depth_in_ft); 
  sim.addTimeDomainMeasurement( meep::Ez, 1200*ft_to_mtr, -depth_in_ft); 


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

  /* run the simulation for 3.5 us*/ 
  sim.run(3500); 

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
