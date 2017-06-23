/**************************************************************************************
 * Source.h 
 * 
 * Describes sources 
 *
 * Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
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
#ifndef ICEPROP_SOURCE_H
#define ICEPROP_SOURCE_H

#include "meep.hpp" 
class TGraph; 
namespace FFTtools
{

}

namespace iceprop 
{

  /** The iceprop source is a wrapper around the meep source. 
   * Right now that seems very silly... but maybe in the future it'll be interesting. 
   *
   * TODO: 
   * Right now it only goes in one component, but in general we should be able to have a weighted
   * component list for like a neutrino or something. 
   *
   */ 

  class Source
  {

    public: 
      virtual const meep::src_time & getSource() const = 0; 
      double getR() const { return r; } 
      double getZ() const { return z; } 
      meep::component getComponent() const { return component; } 
      virtual ~Source() {; } 

    protected: 
      double r, z; 
      meep::component component;
  }; 


 /* dumb wrapper around meep gaussian pulse... does take care of the frequency conversion from our units though */ 
 class GaussianPulseSource  : public Source
 {
   public: 
     GaussianPulseSource(double r, double z, double f=1, double w=0.5, meep::component c = meep::Ez); 
     virtual const meep::src_time & getSource() const { return src; } 
   private: 
     meep::gaussian_src_time src; 
 };


 /* This uses a graph as the model. It takes the hilbert transform for the imaginary part. I'm having trouble getting it to do something
  * reasonale, I don't think it likes discontinuities :( */ 
 class GraphSource  : public Source
 {
   public: 
     GraphSource(double r, double z, const TGraph * g, meep::component c = meep::Ez); 
     virtual const meep::src_time & getSource() const { return src; } 
     virtual ~GraphSource(); 
     TGraph * getReal() { return g[0]; } 
     TGraph * getImag() { return g[1]; } 

   private:
     meep::custom_src_time src; 
     TGraph *g[2]; 
 }; 

 /* Attempt to use the impulse response of a digital butterworth filter a source. Not working out so great for me ... */ 
 class ButterworthSource : public Source
 {
   public: 
     ButterworthSource(double r, double z, double fmin = 0.2, double fmax = 0.7,
                        meep::component c = meep::Ez, int order = 2, int oversample = 16, int max_length = 1000);
     virtual const meep::src_time & getSource() const { return graph_source->getSource(); } 
     GraphSource * getGraphSource() { return graph_source; } 
     virtual ~ButterworthSource(); 
   private:
     GraphSource * graph_source; 
 }; 

} 


#endif




