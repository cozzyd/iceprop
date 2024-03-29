/**************************************************************************************
 * Firn.h 
 * 
 * Describes firn properties. 
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
#ifndef ICEPROP_FIRN_H
#define ICEPROP_FIRN_H

#include "TGraph.h" 
#include <vector> 

namespace ROOT
{
  namespace Math
  {
    class Interpolator; 
  }
}

namespace iceprop
{

  const double k = 0.845*1e-3; 

  /**
   *  The Firn class encompasses the ice properties. 
   *  All that matters for us right now is the index of refraction, but usually that is derived from the density. 
   *
   *
   */
  class Firn
  {
    public: 
      /** Returns n, the index of refraction. By default, this is based on density, but can be overriden. 
       *  z should be negative (relative to surface) if you want ice and in m. 
       *
       **/
      virtual double getIndexOfRefraction(double z) const;   

      /** Returns the density in kg/m^3. This is just used to get the index of refraction, so if you override that you don' tneed to override this. 
       *  z should be negative (relative to surface)  if you want ice and in m. 
       *
       **/
      virtual double getDensity(double z) const { return 0; } // kg/m^3. Will be used for refraction unless overriden./  

      TGraph * makeGraph(bool density = true, int npoints = 200, double min = 0, double max = 120) const; 

      /** Returns the dielectric constant, epsilon. By default this is derived from the index of refraction (assuming no relative permeability).
       *  z should be negative (relative to surface) if you want ice and in m
       **/ 
      virtual double eps(double z) const; 


      /** Key-based firn retrieval, for use with config file 
       *
       * Possible values right now: 
       *
       *   DoubleExponentialDensityFirn shallow_scale deep_scale surface_dens crit_dens deep_dens; 
       *   ArthernFirn
       *   DensityTableFirn file
       **/ 
      static Firn * getFirn(const char * key); 


      static double rho2n(double rho) 
      {
        return 1 + k * rho; 

      }

      static double n2rho(double n) 
      {
        return (n-1)/k; 
      }

      virtual ~Firn() { ; } 
  }; 

  /** Firn with a ``double exponential model'' */
  class DoubleExponentialDensityFirn : public Firn 
  {
    public: 
      DoubleExponentialDensityFirn (double shallow_length_scale, double deep_length_scale, double surface_density, double critical_density = 550, double deep_density = 917); 
      virtual double getDensity(double z) const; 

    private: 
      double scale_shallow; 
      double scale_deep; 
      double rho_surf; 
      double rho_c; 
      double rho_deep; 
      double z_c; 
  };  

  class ExponentialFirn : public Firn 
  {
    public:
      ExponentialFirn(double scale_m, double surface_density, double deep_density = 917, double max_z=2800, double n_under =2) 
       : rho_surf(surface_density), scale(scale_m), rho_deep(deep_density) , max_depth(max_z), rho_under(n2rho(n_under))
      {}
    private: 
      virtual double getDensity(double z) const; 
      double rho_surf; 
      double scale; 
      double rho_deep; 
      double max_depth; 
      double rho_under; 
  }; 

  /** I'm sure Peter loves that this is named after him in the "literature" (https://icecube.wisc.edu/~mnewcomb/radio/#icedensity) */ 
  class GorhamFirn 
    : public ExponentialFirn
  {
    public: 
      GorhamFirn() 
        : ExponentialFirn (-0.01535, 374.15) 
      {
      }
  };

  /** This uses the attenuation depth model used in Arthern et al */ 
  class ArthernFirn  : public DoubleExponentialDensityFirn 
  {
    public:
      ArthernFirn(); 
  }; 

  /** Fit to a bunch of density data */
  class MultiDatasetFit : public DoubleExponentialDensityFirn
  {
    public: 
      MultiDatasetFit(); 
  }; 

  /** Density data from Hawley and Morris 06 , with fit outside of bounds */ 


  /** This takes a firn model and adds gaussian to it */ 
  class PerturbedFirn : public Firn
  {
    public: 
      /** Perturbed firn with fixed density perturbations */ 
      PerturbedFirn(const Firn & base, int nperturbs, const double * zs, const double * A, const double * sigmas); 
      virtual double getDensity(double z) const; 

    private: 

      const Firn& base; 
      std::vector<double> zs; 
      std::vector<double> As; 
      std::vector<double> sigmas; 


  }; 

  /** Firn based on interpolating a density table */ 
  class DensityTableFirn : public Firn 
  {

    public: 
      /** Pass the points in a density table for the firn.  Note that
       * the table should have the opposite sign for z (i.e. in units of depth)
       * to match most references */ 
      DensityTableFirn(int Npoints, const double * z, const double * rho, const Firn * outside = 0); 
      /** Pass a file. It should be a whitespaec separated depth and density */ 
      DensityTableFirn(const char * file, const Firn * outside = 0); 

      enum InterpolationType
      {
        LINEAR, 
        CSPLINE,
        AKIMA 
      }; 

      void setInterpolationType(InterpolationType typ); 

      virtual double getDensity(double z) const; 

      /* Retrieve the graph. While non-const to allow plotting,
       * be careful what you do with it. The graph is assumed to be sorted
       * in z, so if mess that up, it will not evaluate properly. Also, if
       * you modify after calling getDensity already and you are using spline interpolation,
       * the spline will not be updated (unless you toggle interpolation types...)!*/ 
      TGraph * getGraph() { return &g; } 

      virtual ~DensityTableFirn(); 

    private: 
      TGraph g;
      ROOT::Math::Interpolator * ipl; 
      const Firn * outside; 
      InterpolationType ipl_type;  
  }; 

  class ConstantFirn : public Firn
  {

    public:
      ConstantFirn(double dens) : rho(dens) { ; } 

      virtual double getDensity(double z) const { return rho; } 
    
    private: 
      double rho; 

  };
}

#endif
