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



#include "Source.h" 
#include "Units.h" 


iceprop::GaussianPulseSource::GaussianPulseSource(double r_, double z_, double f, double w, meep::component c) 
  :  src( f * GHz_to_meep, w * GHz_to_meep)
{
  r = r_; 
  z = z_; 
  component = c; 
  //that's it for now
} 


