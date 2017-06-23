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
  return std::complex<double>(g[0]->Eval(t), g[1]->Eval(t)); 

}

iceprop::GraphSource::~GraphSource()
{
  delete g[0]; 
  delete g[1]; 
}

iceprop::GraphSource::GraphSource(double r_, double z_, const TGraph*  g_ , meep::component c)
: src(graph_eval_meep, (void*) &(g[0]) )
{
  src.is_integrated= false; 
  r = r_; 
  z = z_; 
  component = c; 
  g[0] = new TGraph(g_->GetN(), g_->GetX(), g_->GetY()); 
  g[0]->SetBit(TGraph::kIsSortedX); 
  g[1] = FFTtools::getHilbertTransform(g[0]); 
  g[0]->SetBit(TGraph::kIsSortedX); 
}

iceprop::ButterworthSource::~ButterworthSource() 
{
  delete graph_source; 
}

iceprop::ButterworthSource::ButterworthSource(double r_, double z_, double fmin, double fmax, meep::component c, int order, int os, int max_length) 
{
  r = r_; 
  z = z_; 
  component = c; 
  FFTtools::ButterworthFilter filt(FFTtools::BANDPASS, order, (fmin/fmax + 1) / (2.*os), (1-fmin/fmax)/ (2.*os)); 
  TGraph * g = filt.impulseGraph(max_length, 1./(2*fmax*os), os/fmax); 
  graph_source = new GraphSource(r,z,g,component); 
  delete g; 
}


