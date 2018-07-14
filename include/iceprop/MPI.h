#ifndef _iceprop_mpi_h
#define _iceprop_mpi_h

/***************************************************************************
 * MPI.h 
 *
 * Just a few helper things 
 * 
 *
 * Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 *
 * This file is part of iceprop. 
 *
 * iceprop is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * iceprop is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * iceprop.  If not, see <http://www.gnu.org/licenses/>.
 *
 *************************************************************************/ 

#include <meep.hpp>
#include <stdlib.h>

#ifdef ENABLE_MPI
#include <mpi.h> 
#endif

namespace iceprop 
{
  namespace mpi
  {
#ifdef ENABLE_MPI
    inline bool enabled() { return true; } 

    inline bool am_master() 
    { 
      int rank; 
      MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
      return rank ==0 ; 
    }

    inline int barrier() { return MPI_Barrier(MPI_COMM_WORLD); } 
    inline void abort(int code) { MPI_Abort(MPI_COMM_WORLD, code); } 

#else
    inline bool enabled() { return false; }
    inline bool am_master() { return true; }
    inline int barrier() { return 0; }
    inline void abort(int code) { exit(code); }
#endif

    //this is just a wrapper around the meep class
    class init
    {
      public:
        init (int nargs, char ** args); 
      private: 
        meep::initialize minit; 
    };
  }
}

#endif
