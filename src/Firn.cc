/************************************************************************************
 *
 * Firn.cc 
 * Implements firn models. 
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


#include "Firn.h" 
#include <cmath>
#include <stdio.h> 
#include <string.h> 
#include "TAxis.h" 
#include "TSpline.h" 


static iceprop::ArthernFirn arthern;

const double k = 0.845e-3; 

double iceprop::Firn::getIndexOfRefraction(double z) const 
{
  return 1 + k * getDensity(z); 
}

double iceprop::Firn::eps(double z) const 
{
  double n = getIndexOfRefraction(z); 
  return n*n; 
}


iceprop::Firn * iceprop::Firn::getFirn(const char * key) 
{

  if (strstr(key,"DoubleExponentialDensityFirn")) 
  {
    double shallow,deep,surface_rho,crit_rho,deep_rho =917; 
    sscanf(key,"DoubleExponentialDensityFirn %lf %lf %lf %lf %lf" , &shallow, &deep, &surface_rho,&crit_rho,&deep_rho); 
    return new DoubleExponentialDensityFirn(shallow,deep,surface_rho,crit_rho,deep_rho); 
  }

  else if (strstr(key,"ArthernFirn"))
  {
    return &arthern; 
  }

  else if (strstr(key,"DensityTableFirn"))
  {
    char file[513]; 
    sscanf(key,"DensityTableFirn %512s", file); 
    return new DensityTableFirn(file); 
  }
  else
  {
    fprintf(stderr,"Sorry, cannot understand %s\n", key); 
  }
  return 0; 
} 

iceprop::DoubleExponentialDensityFirn::DoubleExponentialDensityFirn(double shallow_length_scale,
                                                                    double deep_length_scale,
                                                                    double surface_density, 
                                                                    double critical_density,
                                                                    double deep_density)
{
  scale_shallow = shallow_length_scale; 
  scale_deep = deep_length_scale; 
  rho_surf = surface_density; 
  rho_c = critical_density; 
  rho_deep = deep_density; 
  z_c= -1* scale_shallow * log((rho_deep - rho_surf) / (rho_deep - rho_c)); 
}

double iceprop::DoubleExponentialDensityFirn::getDensity(double z) const 
{
  if ( z > 0) return 0; 

  if (z < z_c) return (rho_deep - (rho_deep - rho_c) * exp((z-z_c)/scale_deep)); 
  else return (rho_deep - (rho_deep - rho_surf)* exp(z/scale_shallow)); 
}


iceprop::ArthernFirn::ArthernFirn() 
  : DoubleExponentialDensityFirn(27,42,280,550,917) 
{
}




iceprop::DensityTableFirn::DensityTableFirn(int Npts, const double * z, const double * rho) 
  : g(Npts,z,rho) 
{
  g.SetBit(TGraph::kIsSortedX); 
  g.GetXaxis()->SetTitle("z (m)") ; 
  g.GetYaxis()->SetTitle("#rho (kg/m^3)") ; 
  ipl = LINEAR; 
  spl = 0;
}

iceprop::DensityTableFirn::DensityTableFirn(const char * f) 
  : g(f) 
{
  g.SetBit(TGraph::kIsSortedX); 
  g.GetXaxis()->SetTitle("z (m)") ; 
  g.GetYaxis()->SetTitle("#rho (kg/m^3)") ; 
  ipl = LINEAR; 
  spl = 0;
}

static int spline_counter = 0; 

double iceprop::DensityTableFirn::getDensity(double z) const
{

  if (ipl == SPLINE3) 
  {
    if (!spl) 
    {
      TString str; str.Form("firnspline%d",spline_counter++); 
      spl = new TSpline3(str.Data(), g.GetX(), g.GetY(), g.GetN()); 
    }

    return spl->Eval(z); 
  }

  else
  {
    return g.Eval(z); 
  }

}

iceprop::DensityTableFirn::~DensityTableFirn() 
{
  if (spl) delete spl; 
}





