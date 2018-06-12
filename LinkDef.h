#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link off all namespaces;

//Firn.h
#pragma link C++ namespace iceprop; 
#pragma link C++ class iceprop::Firn; 
#pragma link C++ class iceprop::PerturbedFirn; 
#pragma link C++ class iceprop::DoubleExponentialDensityFirn; 
#pragma link C++ class iceprop::ArthernFirn; 
#pragma link C++ class iceprop::MultiDatasetFit; 
#pragma link C++ class iceprop::DensityTableFirn; 
#pragma link C++ class iceprop::ConstantFirn; 

//Source.h 
#pragma link C++ class iceprop::Source; 
#pragma link C++ class iceprop::GaussianPulseSource; 
#pragma link C++ class iceprop::GraphSource; 
#pragma link C++ class iceprop::ButterworthSource; 

//Sim.h 

#pragma link C++ class iceprop::SimGeometry; 
#pragma link C++ class iceprop::Sim; 
#pragma link C++ class iceprop::StepOutput; 
#pragma link C++ class iceprop::GlobalMaximum; 
#pragma link C++ class iceprop::TimeDomainMeasurement; 


#endif 
