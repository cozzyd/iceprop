

void snapshot() 
{
  iceprop::mpi::init init(0,0); 
  system("mkdir -p snapshot"); 

  iceprop::MultiDatasetFit fit ;
  iceprop::DensityTableFirn frn("data/hawley08_neutron.txt",&fit);

  iceprop::SimGeometry g; 
  g.max_depth =200; 
  g.resolution=12; 
  g.output_skip_factor=20; 
  g.sky_height=50; 
  g.max_r = 2000; 

  iceprop::ButterworthSource s(0, -100, 0.12, 0.04); 

  iceprop::Sim sim(&frn, &g, &s); 
  //do this first, otherwise the outputs can get screwed up if you add them before loading the snapshot (because trees will be recreated instead of updated and images will get the wrong index). 
  sim.loadSnapshot("snapshot/snaps",-1); 

  sim.addTimeDomainMeasurement( meep::Ez, 0, 50); 
  sim.addTimeDomainMeasurement( meep::Ez, 10, 50); 
  sim.addTimeDomainMeasurement( meep::Ez, 20, 50); 
  sim.addTimeDomainMeasurement( meep::Ez, 40, 50); 

  iceprop::StepOutput o; 
  o.out_dir = "snapshot"; 
  o.format = iceprop::O_ROOT; 
  o.type = iceprop::Mag; 
  sim.addStepOutput(o); 


  //try to load from a snapshot (but nothign will happen if it's not available). 
  sim.enableSnapshots("snapshot/snaps",1000); 
  sim.run(8000); 

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

  }

}
