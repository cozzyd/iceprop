
/** Test macro to check air propagation works as expected */ 


void air2air() 
{
  iceprop::mpi::init init(0,0); 
  system("mkdir -p air2air"); 

  iceprop::ArthernFirn frn ;

  iceprop::SimGeometry g; 
  g.max_depth =0; 
  g.sky_height=100; 
  g.max_r = 200; 

  iceprop::GaussianPulseSource s(0, 50, 0.2, 0.2); 

  iceprop::Sim sim(&frn, &g, &s); 

  sim.addTimeDomainMeasurement( meep::Ez, 0, 50); 
  sim.addTimeDomainMeasurement( meep::Ez, 10, 50); 
  sim.addTimeDomainMeasurement( meep::Ez, 20, 50); 
  sim.addTimeDomainMeasurement( meep::Ez, 40, 50); 

  iceprop::StepOutput o; 
  o.out_dir = "air2air"; 
  o.format = iceprop::O_PNG; 
  sim.addStepOutput(o); 
  sim.run(500); 

  if (iceprop::mpi::am_master()) 
  {

    TCanvas * ms = new TCanvas("measurements","Measurements",1000,800); 

    ms->Divide(2,2); 
    ms->cd(1); 
    sim.getMeasurements()[0].makeGraph()->Draw("al"); 
    ms->cd(2); 
    sim.getMeasurements()[1].makeGraph()->Draw("al"); 
    ms->cd(3); 
    sim.getMeasurements()[2].makeGraph()->Draw("al"); 
    ms->cd(4); 
    sim.getMeasurements()[3].makeGraph()->Draw("al"); 

    ms->SaveAs("air2air/air2air.pdf"); 
  }

}
