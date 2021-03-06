/**************************************************************************************
 * Units.h 
 * 
 * Constants to convert between our unit system and Meep's internal unit system.  
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


#ifndef ICEPROP_UNITS_H
#define ICEPROP_UNITS_H


namespace iceprop 
{
  //TODO: I'm not so good at algebra so we should fix this 
  static const double meep_to_ns = 3.3354095; 
  static const double ns_to_meep = 1./meep_to_ns; 
  static const double GHz_to_meep = meep_to_ns; 
}

#endif
