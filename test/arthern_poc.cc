/** Quick and extremely dirty proof of concept test based on one of the examples. 
 *  
 *
 *  Cosmin Deaconu
 *  <cozzyd@kicp.uchicago.edu>
 *
 */ 

#include <meep.hpp> 
#include "TGraph.h" 
#include "TStyle.h" 
#include "TFile.h" 
#include "TCanvas.h" 
#include "TH2.h" 
using namespace meep; 


const int height =300; //ft ?
const int width = 1200; //ft? 

const int surface=250;// ft; 
const int antenna_depth = 200; //ft; 

const double graph_zs[] = {200,200,200, 100,100,100,1,1,1}; 
const double graph_xs[] = {200,600,1000, 200,600,1000,200,600,1000}; 
const int ngraphs = sizeof(graph_zs)/sizeof(*graph_zs); 

double density(double z) 
{


  double rho_surf = 280; 
  double rho_c = 550; 
  double rho_deep = 917; 

  double scale_1 = 27; 
  double scale_2 = 42; 

  double z_crit = -1* scale_1 * log((rho_deep - rho_surf) / (rho_deep - rho_c)); 

  if (z < z_crit) return (rho_deep - (rho_deep - rho_c) * exp((z-z_crit)/scale_2)); 
  else return (rho_deep - (rho_deep - rho_surf)* exp(z/scale_1)); 

}

double eps(const vec&p)
{
  double k = 0.845 * 1e-3; 
  if (p.z() > surface) return 1; 
  return pow(1.+ k * density(0.3048*(p.z()-surface)),2); 

}


void fillHist( TH2 * h, const fields * f) 
{
  for (int i = 1; i <=  h->GetNbinsX(); i++)
  {
    for (int j = 1; j <= h->GetNbinsY(); j++) 
    {
      double x= h->GetXaxis()->GetBinCenter(i);  
      double y= h->GetYaxis()->GetBinCenter(i)+surface;  

      h->SetBinContent(i,j, f->get_field(Ey,veccyl(x,y)).real()); 

    }
  }

}

int main(int argc, char **argv) {

  initialize mpi(argc, argv); // do this even for non-MPI Meep, apparently



  double resolution = 3; // pixels per distance
  grid_volume v = volcyl(width,height, resolution); 
  structure s(v, eps, pml(5.0));// pml is a ficitions perfectly matched layer, I'm setting it to a size of 5 ft to try to dampen reflections from the edge. 
  fields f(&s);
  
  f.output_hdf5(Dielectric, v.surroundings());
  
  double freq = 0.3, fwidth = 0.1;

  gaussian_src_time src(freq, fwidth);

  f.add_point_source(Ez, src, veccyl(0, surface-antenna_depth));

  double last= -100; 

  TFile fout("arthern.root","RECREATE"); 


  std::vector<TGraph*> graphs; 
  for (int i =0; i < ngraphs; i++)
  {
    TGraph * g = new TGraph; 
    g->SetTitle(TString::Format("x=%g ft, depth = %g ft", graph_xs[i], graph_zs[i])); 
    graphs.push_back(g); 
  }


  TGraph gsurface(2);
  gsurface.SetPoint(0,0,0); 
  gsurface.SetPoint(1,width,0); 

  while (f.time() < 2000) {
    f.step();
    if (f.time() - last >= 5 && f.time() == int(f.time())) 
    {
//      TCanvas c("cey","Ey",width*resolution,height*resolution); 
//      TH2F  ey(TString::Format("ey_%g",f.time()),TString::Format("Ey t =%g",f.time()),width*resolution, 0,width,height*resolution, -surface, -surface+height); 
//      ey.SetStats(0); 
//      fillHist(&ey,&f); 
//      ey.DrawCopy("colz"); 
//      gsurface.Draw("lsame"); 
//      c.SaveAs(TString::Format("rootout/ey_%g.png", f.time())); 
      f.output_hdf5(Ez, v.surroundings(),0,false,true);
      last = f.time(); 
//      ey.Write(); 
    }
    for (int ig = 0; ig < ngraphs; ig++)
    {
      graphs[ig]->SetPoint(graphs[ig]->GetN(),f.time(), f.get_field(Ez, veccyl(graph_xs[ig], -graph_zs[ig]+surface)).real()); 
    }
  }

  for (int i = 0; i < ngraphs; i++) graphs[i]->Write(TString::Format("g_%g_%g", graph_xs[i], graph_zs[i])); 
  
  return 0;
} 

