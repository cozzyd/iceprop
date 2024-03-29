/***************************************************************************
 * Sim.h 
 *
 *
 * Primary simulation driver for iceprop. Dragons be here. 
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
 *************************************************************************/ 

#ifndef ICEPROP_SIM_H
#define ICEPROP_SIM_H

#include <vector> 
#include "TGraph.h" 
#include <meep.hpp> 

class TFile; 
class TTree; 
class TObject; 
class TH2; 
class TPad; 

namespace iceprop 
{
  /** output formats, should probably be a enum but I'm too lazy to make the necessary overloads */
   const int O_PNG =1; 
   const int O_PDF =2;
   const int O_ROOT=4; 
   const int O_HDF5=8; 
   typedef int OutputFormat_t; 


  /** This defines the geometry of the simulation */ 
  struct SimGeometry
  {
    int max_depth = 200; /// The maximum  depth (in m) of the simulation (for now, integer)
    int max_r = 500; /// The maximum radius (or extent) (in m) of the simulation (for now, integer)
    int sky_height = 10; /// The amount of sky to simulate (in m)  (for now, integer)
    int  resolution = 10; ///The ``resolution" (number of elements per m) 
    int output_skip_factor = 1; ///record every nth bin in the output (to avoid having to calculate the field everywhere) 

    double courant_factor = 0.5; /// The courant factor controls the scaling between delta_t and delta_x. 
    double pml_size = 5; ///The size of the "perfectly matched layer" "free space" boundary condition. See Meep documentation. 
    bool is_cylindrical = true; 

    SimGeometry() 
    {
    }
  }; 



  /* Since many things are complex quantities, we need to decide what we want to output */ 
  enum ScalarType
  {
    Real, 
    Imag, 
    Mag, 
    Mag2, 
    Phase
  }; 
  

  /* This is the output of a time domain measurement added with addTimeDomainMeasurement. Not intended to be filled in by user. */ 
  struct TimeDomainMeasurement
  {
    double r, z; 
    meep::component what; 
    std::vector<double> t; 
    std::vector<std::complex<double> >  val; 
    int nskip; 
    const char * name; 

    /** make a graph of this measurement. Allocates new memory*/ 
    TGraph * makeGraph(ScalarType type = Real) const; 
    void makeName(TString & str, ScalarType type = Real) const; 
  }; 

  /* This is the output of the result of trackMaximum. Not intended to be filled in by the user. */ 
  struct GlobalMaximum
  {
    meep::component what; 
    ScalarType type; 
    TH2 * max; 
    TH2 * tmax; 
  }; 

  /* This is the output of the result of trackMaximum. Not intended to be filled in by the user. */ 
  struct GlobalIntegral
  {
    meep::component what; 
    ScalarType type; 
    TH2 * integ; 
    TH2 * tfirst; 
  }; 


  /* This is used to define outputs in the middle of the simulation (typically not at every step unless you like waiting forever :) ) */ 
  struct StepOutput
  {
    meep::component what; //default is Ez
    ScalarType type;  //default is Real
    int skip_factor;  //default is 50 
    int skip_offset; //first step to write out. default is 0. 
    const char * out_dir;  // default is current_dir 

    //For image output, the maximum scale. If zero (default), auto-scaled. For phase, if non-zero, 0-2 pi is scaled to these units. 
    double output_scale; 




    OutputFormat_t format;// output format, default is PNG | ROOT

   /* output prefix. Default of 0 makes it  autogenerated from what and type. 
    * the outputs will look like 
    *   out_dir/output_prefix.root (this will be a tree of TH2's, where the tree name will be also be the output_prefix. The branch name will be "hist." 
    *   TODO: I thought about putting everything in the same tree, but that doesn't work well with different output cadences. I will revisit this, there is probably a smarter way to organize the output (for example, allowing a stepOutput to have multiple things output). 
    *   or 
    *   out_dir/output_prefix.%05d.[png|pdf|hf5] where %05d is the index (i.e. the ith output file, NOT the time or step number). This format makes it easy to use with ffmpeg. 
    *   Note that having multiple outputs with the same output prefix (auto generated or not) may confuse each other. 
    *
    */ 
    const char * output_prefix; 

    // dimensions of canvas (important for image output) . Default of 0 corresponds to  110% of histogram size. 
    int cw,ch; 
    bool double_precision; // true for double precision, otherwise single is used. Default is false. 


    TFile * outf; //default is 0. this will be used for the ROOT output file. You can overwrite it I suppose (e.g. to put multiple trees in the same file) but you should know what you're doing. 
    int index; //This is used for output naming for images. It is zero by default, but I guess you could set it to something else if you want to start at a different number. 
    TPad * c; // you're welcome to pass your own pad... otherwise a canvas of size cw and ch is created.

    TH2* h; //no touchy... will be overwritten
    TTree* t;// no touchy... will be overwritten 


    StepOutput() 
    {
      what = meep::Ez; 
      type = Real; 
      skip_factor=50; 
      skip_offset = 0; 
      out_dir = "." ; 
      output_scale = 0; 
      format = O_PNG | O_ROOT; 
      output_prefix = 0;
      cw = 0; 
      ch =0; 
      double_precision = false; 
      outf = 0; 
      index = 0; 
      h = 0; 
      t = 0; 
      c= 0; 
    }

    
  };

  class Firn; 
  class Source; 
  class Sim
  {

    public: 
      Sim(const Firn * firn, const SimGeometry * geom, const Source * source = 0); 

      void addSource(const Source * source); 


      /** To work around job time limits on clusters,
       *  snapshots allow taking a snapshot of a simulation 
       *  Note that the snapshot only saves the fields. It does NOT save the
       *  dielectrics or geometry (since that is handled by firn/geom). 
       *
       *  If you use a different geometry, it will probably crash. If you use a different firn... funny things might happen. 
       *
       *
       *  A snapshot directory can hold any number of snapshots, indexed by i. Choosing i=-1 will load the latest snapshot. 
       */ 


      //Load a snaphsot from the path. Returns 0 if it worked. 
      //you should do this BEFORE adding outputs / globals otherwise any trees might get overwritten. 
      //and indices for image /hdf5 outputs won't be correct. 
      int loadSnapshot(const char * path, int i = -1); 
      int snapshot(const char *path, int i = -1); 
      void enableSnapshots(const char * path, int nsteps = 5000 , int start_index = -1); 

      /* Adds a time domain measuremnt at the given position 
       *
       * skipfactor: If you don't want to record every time step, set skipfactor to something bigger than 1 
       *
       **/ 
      void addTimeDomainMeasurement(meep::component what, double r, double z, const char * name = 0,  int skipfactor = 1);  

      /* Enables tracking the maximum of this quantity and the time at which the maximum occurs. */ 
      void trackGlobalMaximum(meep::component what, ScalarType type = Mag); 
      void trackGlobalIntegral(meep::component what, ScalarType type = Mag); 

      int  saveIntermediateGlobals(const char * file, int skip_factor = 100) ; 
      
      void addStepOutput(StepOutput output);  
      
      /** Runs the simulation for some time steps. There is no way to go back right now. Sorry.*/ 
      void run(double time, bool relative = false); 

      /* Returns the current simulation time */ 
      double getCurrentTime() const; 

      /** Make a histogram based on the current state. fillme will be used if not null, otherwise a new one will be allocated.
       * For now, it must be the same size as the simulation or it will crash. This will likely be addressed later. 
       **/ 
      TH2 * makeHist(meep::component what = meep::Ez, ScalarType type = Real,  TH2 * fillme = 0, bool double_precision = true) const; 

      /* This makes a histogram, DrawCopy's it, and deletes it */ 
      void draw(meep::component what=meep::Ez, ScalarType typ = Real) const; 

      const std::vector<TimeDomainMeasurement> & getMeasurements() const  { return measurements; } 

      const std::vector<GlobalMaximum> & getMaximums() const  { return maxima; } 
      const std::vector<GlobalIntegral> & getIntegrals() const  { return integrals; } 

      /** Low-level meep accessors. May need to hide from ROOT 5. */ 
      meep::fields * getFields() { return f; } 
      meep::grid_volume & getVolume() { return gv; } 
      meep::structure * getStructure() { return s; } 
      virtual ~Sim(); 

    private: 

      /* our stuff */ 
      const Firn * firn; 
      const SimGeometry * geom; 

      std::vector<TimeDomainMeasurement> measurements; 
      std::vector<GlobalMaximum> maxima; 
      std::vector<GlobalIntegral> integrals; 
      std::vector<StepOutput*> outputs; 
      std::vector<TObject *> delete_list; 
      TGraph gsurf; 

      /** meeps stuff. May need to hide from ROOT 5 at some point... */ 
      meep::grid_volume gv; 
      meep::fields *f; 
      meep::structure *s; 

      TFile * intermediate_globals_file; 
      TTree * intermediate_globals_tree; 
      int intermediate_globals_interval; 

      int snapshot_nsteps; 
      const char * snapshot_path; 
      int snapshot_i; 
      bool snapshot_loaded; 
      double time_offset; 
      int step; 
      TString loaded_measurements_path; 
      TFile * loaded_measurements; 
  }; 
}


#endif 
