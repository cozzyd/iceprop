 /************************************************************************************
 *
 * Source.cc 
 * Implements field sources
 * Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 * 
 *
 * This file is part of iceprop. 
 *
 * iceprop is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * iceprop is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with iceprop.  If not, see <http://www.gnu.org/licenses/>.
 *
 *************************************************************************************/ 



#include "iceprop/Source.h" 
#include "iceprop/Units.h" 
#include "TGraph.h" 
#include "DigitalFilter.h" 
#include "FFTtools.h" 


iceprop::GaussianPulseSource::GaussianPulseSource(double r_, double z_, double f, double w, meep::component c) 
  :  src( f * GHz_to_meep, w * GHz_to_meep)
{
  r = r_; 
  z = z_; 
  component = c; 
  //that's it for now
} 



static std::complex<double> graph_eval_meep(double tmeep, void * data)
{
  double t = tmeep * iceprop::meep_to_ns; 
  TGraph ** g = (TGraph**) data; 
  if ( t > g[0]->GetX()[g[0]->GetN()-1] ) return std::complex<double>(0,0); 
  return std::complex<double>(g[0]->Eval(t), g[1]->Eval(t)); 

}

iceprop::GraphSource::~GraphSource()
{
  delete g[0]; 
  delete g[1]; 
}

iceprop::GraphSource::GraphSource(double r_, double z_, const TGraph*  gr, TGraph * gi, meep::component c)
: src(graph_eval_meep, (void*) &(g[0]) )
{
  src.is_integrated= true; 
  r = r_; 
  z = z_; 
  component = c; 
  g[0] = new TGraph(gr->GetN(), gr->GetX(), gr->GetY()); 
  g[0]->SetBit(TGraph::kIsSortedX); 
  g[1] = gi ? new TGraph(gi->GetN(), gi->GetX(), gi->GetY()) :  FFTtools::getHilbertTransform(g[0]); 
  g[0]->SetBit(TGraph::kIsSortedX); 
}


void iceprop::GraphSource::offset(double t) 
{
  for (int i = 0; i < g[0]->GetN(); i++) 
  {
    g[0]->GetX()[i]+=t; 
    g[1]->GetX()[i]+=t; 
  }
}


iceprop::ButterworthSource::~ButterworthSource() 
{
  delete graph_source; 
}

iceprop::ButterworthSource::ButterworthSource(double r_, double z_, double fc, double w, double fnyq, meep::component c, int order, int os, int max_length) 
{
  r = r_; 
  z = z_; 
  component = c; 
  FFTtools::ButterworthFilter filt(FFTtools::BANDPASS, order, fc/fnyq, w/fnyq); 
  TGraph * g = filt.impulseGraph(max_length, 0.5/fnyq,fnyq/(fc-w)); 
  TGraph * gos = FFTtools::getInterpolatedGraphFreqDom(g,0.5/(fnyq*os)); 

  // shift to zero 
  double shift = gos->GetX()[0]; 
  if (shift) 
  {
    for (int i = 0; i < gos->GetN(); i++) gos->GetX()[i] -=shift; 
  }
  
  TGraph * gos_i =FFTtools::getHilbertTransform(gos); 
  long max= TMath::LocMax(gos->GetN(), gos->GetY()); 
  for (int i = 0; i < gos->GetN(); i++) 
  {
    gos->GetY()[i] *= TMath::Gaus(gos->GetX()[i],gos->GetX()[max], fnyq/(fc-w)); 
    gos_i->GetY()[i] *= TMath::Gaus(gos->GetX()[i],gos->GetX()[max], fnyq/(fc-w)); 
  }
  graph_source = new GraphSource(r,z,gos,gos_i,component) ; 
  delete g; 
  delete gos; 
}


