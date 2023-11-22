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

    inline int rank() 
    { 
      int rank; 
      MPI_Comm_rank(MPI_COMM_WORLD,&rank); 

      return rank ; 
    }
    inline int max_rank() 
    {
      int size; 
      MPI_Comm_size(MPI_COMM_WORLD,&size); 
      return size; 
    }

    inline int am_master() { return rank() == 0; } 

    inline int barrier() { return MPI_Barrier(MPI_COMM_WORLD); } 
    inline void abort(int code) { MPI_Abort(MPI_COMM_WORLD, code); } 

    inline int broadcast(double * data, int N) { return MPI_Bcast(data, N, MPI_DOUBLE, 0, MPI_COMM_WORLD); }

#else
    inline bool enabled() { return false; }
    inline bool am_master() { return true; }
    int broadcast(double * data, int N) { return 0; }  
    inline int barrier() { return 0; }
    inline int rank() { return 0; } 
    inline int max_rank() { return 1; } 
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
