/* simulation of summit'13 measurements
 *
 *
 * */ 

const double ft_to_mtr = 0.3048; 
const double transmitter_depth = -3*ft_to_mtr; 

const bool use_butterworth = true; 
const double f= 0.25; 
const double w = 0.05; 

const bool movie = true; 


iceprop::ArthernFirn arthern ;
iceprop::MultiDatasetFit multi_fit; 


iceprop::Source * s = new iceprop::ButterworthSource(0, transmitter_depth,  f,w); 

const char * firn_descs[] = 
{
  "Arthern", 
  "Multi Fit", 
  "Hawley08 + Multi Fit", 
  "Alley88", 
  "" 
}; 


iceprop::Firn * getFirn(int firn) 
{

  if (firn == 0) 
    return &arthern; 

  if (firn == 1) 
    return &multi_fit; 

  if (firn == 2) 
    return new iceprop::DensityTableFirn ("data/hawley08_neutron.txt", &multi_fit); 

  if (firn == 3) 
  {
    return new iceprop::DensityTableFirn("data/alley_koci.txt", &multi_fit); 
  }

  return 0; 
}


void summit13_single_transmit(int firn = 0, bool vpol = true) 
{

  /* define and make the output dir */ 
  TString dir; dir.Form("summit13_single_low_f/firn_%d_%s/", firn, vpol ? "vpol" : "hpol"); 
  TString cmd; cmd.Form("mkdir -p %s", dir.Data()); system(cmd.Data()); 
  cmd.Form("echo %s  > %s/firn.txt", firn_descs[firn], dir.Data()); system(cmd.Data()); 

  /* use "arthern model" */
  iceprop::Firn * frn = getFirn(firn); 

  /* Define simulation geometry */ 
  iceprop::SimGeometry g; 
  g.resolution=20; // 5cm, enough for 300 MHz 
  g.courant_factor = 0.5; 
  g.output_skip_factor=10; 
  g.pml_size = 30; 
  g.sky_height=g.pml_size+10; 
  g.max_depth =(620 * ft_to_mtr) + g.pml_size; 
  g.max_r = 360; 

  /* Initialize the simulation */ 
  iceprop::Sim sim(frn, &g, s); 


  meep::component c = vpol ? meep::Ez: meep::Er;

  /* add time domain measurements */ 
  sim.addTimeDomainMeasurement( c, 0, transmitter_depth); 
  sim.addTimeDomainMeasurement( c, 1, transmitter_depth); 
  sim.addTimeDomainMeasurement( c, 2, transmitter_depth); 
  sim.addTimeDomainMeasurement( c, 4, transmitter_depth); 
  sim.addTimeDomainMeasurement( c, 8, transmitter_depth); 
  sim.addTimeDomainMeasurement( c, 25, transmitter_depth); 
  sim.addTimeDomainMeasurement( c, 16, transmitter_depth); 
  sim.addTimeDomainMeasurement( c, 32, transmitter_depth); 
  sim.addTimeDomainMeasurement( c, 50, transmitter_depth); 
  sim.addTimeDomainMeasurement( c, 64, transmitter_depth); 
  sim.addTimeDomainMeasurement( c, 100, transmitter_depth); 
  sim.addTimeDomainMeasurement( c, 128, transmitter_depth); 
  sim.addTimeDomainMeasurement( c, 200, transmitter_depth); 
  sim.addTimeDomainMeasurement( c, 256, transmitter_depth); 


  sim.addTimeDomainMeasurement( c, 200 * ft_to_mtr,-150*ft_to_mtr ); 
  sim.addTimeDomainMeasurement( c, 200 * ft_to_mtr,-200*ft_to_mtr ); 
  sim.addTimeDomainMeasurement( c, 200 * ft_to_mtr,-400*ft_to_mtr ); 
  sim.addTimeDomainMeasurement( c, 200 * ft_to_mtr,-580*ft_to_mtr ); 
  sim.addTimeDomainMeasurement( c, 480 * ft_to_mtr,-200*ft_to_mtr ); 
  sim.addTimeDomainMeasurement( c, 480 * ft_to_mtr,-400*ft_to_mtr ); 
  sim.addTimeDomainMeasurement( c, 480 * ft_to_mtr,-580*ft_to_mtr ); 
  sim.addTimeDomainMeasurement( c, 1000 * ft_to_mtr,-200*ft_to_mtr ); 
  sim.addTimeDomainMeasurement( c, 1000 * ft_to_mtr,-400*ft_to_mtr ); 
  sim.addTimeDomainMeasurement( c, 1050 * ft_to_mtr,-600*ft_to_mtr ); 


  /* If movie, add step output*/ 
  if (movie) 
  {
    iceprop::StepOutput o; 
    o.out_dir = dir.Data(); 
    o.what = c; 
    o.type = iceprop::Mag; 
    o.format = iceprop::O_ROOT; 
    sim.addStepOutput(o); 
  }
  TString fname; 
  fname.Form("%s/intermediate.root", dir.Data()); 
  sim.saveIntermediateGlobals(fname.Data(),500); 
  /* uncomment to track global maximum of Ez (very slow!) */ 
  sim.trackGlobalMaximum(c); 
  sim.trackGlobalIntegral(c); 

  /* run the simulation for 2.5 us*/ 
  sim.run(500); 

  /* that's it for running the simulation, the rest plots and saves some of the stuff we made*/ 
  fname.Form("%s/measurements.root", dir.Data());
  TFile fout(fname.Data(),"RECREATE"); 

  /* how many measurements did I give? */ 
  int n = (int)  sim.getMeasurements().size(); 

  for (int i = 0; i< n; i++)
  {
    TGraph * greal = sim.getMeasurements()[i].makeGraph(iceprop::Real); 
    TGraph * gimag = sim.getMeasurements()[i].makeGraph(iceprop::Imag); 
    TGraph * gmag = sim.getMeasurements()[i].makeGraph(iceprop::Mag); 
    greal->Write(); 
    gimag->Write(); 
    gmag->Write(); 
  }

  sim.getMaximums()[0].max->SetDirectory(&fout);
  sim.getMaximums()[0].max->Write(); 
  sim.getMaximums()[0].tmax->SetDirectory(&fout); 
  sim.getMaximums()[0].tmax->Write(); 
  sim.getIntegrals()[0].integ->SetDirectory(&fout);
  sim.getIntegrals()[0].integ->Write(); 
  sim.getIntegrals()[0].tfirst->SetDirectory(&fout); 
  sim.getIntegrals()[0].tfirst->Write(); 


}
