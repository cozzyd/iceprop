{
  gSystem->Load("lib/libiceprop.so"); 
  gSystem->Load("libRootFftwWrapper.so"); 

  //who knows what this will do in ROOT 5
  gStyle->SetPalette(kLightTemperature); 
  gStyle->SetNumberContours(255); 
}
