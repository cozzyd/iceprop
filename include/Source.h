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

namespace iceprop 
{

  /** The iceprop source is a wrapper around the meep source. 
   * Right now that seems very silly... but maybe in the future it'll be interesting. 
   */ 
  class Source
  {

    public: 
      virtual const meep::source * source() const = 0; 

      double getR() const { return r; } 
      double getZ() const { return z; } 
    protected: 
      double r, z; 
  }; 


 /* dumb wrapper around meep gaussian pulse... does take care of the frequency conversion from our units though */ 
 class GaussianPulseSource 
 {
   public: 
     GaussianPulseSource(double r, double z, double f=0.5, double w=0.2); 
     virtual const meep::source & source() const { return src }; 
   private: 
     meep::gaussian_src_time src; 
 };

} 


#endif




