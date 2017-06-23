/************************************************************************************
 *
 * Sim.cc 
 * Implements bulk of simulation
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

#include "iceprop.h"
#include <map>

#include "TCanvas.h" 
#include "TFile.h" 
#include "TTree.h" 
#include "TH2.h" 
#include "TGraph.h" 
#include "TROOT.h" 



/** necessary for epsilon with no void */
static iceprop::Firn * the_firn; 
static iceprop::SimGeometry * the_geom;

static double eps(const meep::vec & v) 
{

  double val = 0;
  // check if in sky 
  if (v.z() > the_geom->max_depth) val =1.; 
  else val = the_firn->eps(v.z()-the_geom->max_depth); 
//  printf("%g %g %g\n", v.z(), v.r(), val); 
  return val; 

}

__attribute__((pure)) 
static double get_scalar_val(iceprop::ScalarType type, std::complex<double> c) 
{
  return 
    type == iceprop::Real  ? c.real()     : 
    type == iceprop::Imag  ? c.imag()     : 
    type == iceprop::Mag   ? std::abs(c)  : 
    type == iceprop::Mag2  ? std::norm(c) : 
    type == iceprop::Phase ? std::arg(c)  : 
    0; 
}

__attribute__((pure)) 
static const char * get_scalar_name(iceprop::ScalarType type) 
{ 
  return 
    type == iceprop::Real     ?  "real"   : 
    type == iceprop::Imag     ?  "imag"   : 
    type == iceprop::Mag      ?  "mag"    : 
    type == iceprop::Mag2     ?  "mag2"   : 
    type == iceprop::Phase    ?  "phase"  : 
    0; 
}

__attribute__((pure)) 
static const char * get_meep_component_name(meep::component type) 
{ 
  return 
    type == meep::Er          ?   "Er"    : 
    type == meep::Ez          ?   "Ez"    : 
    type == meep::Ep          ?   "Ephi"  : 
    type == meep::Ex          ?   "Ex"    : 
    type == meep::Ey          ?   "Ey"    : 
    type == meep::Hr          ?   "Hr"    : 
    type == meep::Hz          ?   "Hz"    : 
    type == meep::Hp          ?   "Hphi"  : 
    type == meep::Hx          ?   "Hx"    : 
    type == meep::Hy          ?   "Hy"    : 
    type == meep::Dielectric  ?   "eps"   : 
    0; 
}



static TH2 * make_hist(const char * name, const char * title, const iceprop::SimGeometry * geom, bool double_precision = true) 
{
  TH2 * h = 
  double_precision ? 
      (TH2*) new TH2D(name, title, geom->max_r*geom->resolution, 0, geom->max_r, 
                    (geom->max_depth+geom->sky_height)*geom->resolution, -geom->max_depth, geom->sky_height) 
     :
      (TH2*) new TH2F(name, title, geom->max_r*geom->resolution, 0, geom->max_r, 
                   (geom->max_depth+geom->sky_height)*geom->resolution, -geom->max_depth, geom->sky_height); 
 
  h->GetXaxis()->SetTitle("r (m)"); 
  h->GetYaxis()->SetTitle("z (m)"); 
  h->SetStats(0); 

  return h; 

}


iceprop::Sim::~Sim()
{
  delete f; 
  delete s; 

  /* make sure we write the trees */ 

  for (size_t o = 0; o < outputs.size(); o++)
  {
    if (outputs[o].t) 
    {
      outputs[o].outf->cd(); 
      outputs[o].t->Write(); 
    }
  }
 


  for (size_t i =0; i < delete_list.size(); i++) 
  {
    if (delete_list[i]) delete delete_list[i]; 
  }
}

iceprop::Sim::Sim(const Firn * firn, const SimGeometry * geom, const Source * source) 
:  firn(firn), geom(geom), gv(meep::volcyl(geom->max_r, geom->max_depth + geom->sky_height, geom->resolution))
{

  //evil thread-unsafe code here , will probably need to fx
  the_firn = (Firn*)firn; 
  the_geom = (SimGeometry*)geom; 
  s = new meep::structure(gv, eps, meep::pml(geom->pml_size)); 
  f = new meep::fields(s); 
  if (source)  addSource(source); 
}


void iceprop::Sim::addSource(const Source * s) 
{
  f->add_point_source(s->getComponent(), s->getSource(), meep::veccyl(s->getR(), geom->max_depth + s->getZ())); 
}


void iceprop::Sim::addTimeDomainMeasurement(meep::component what, double r, double z, int skip_factor) 
{
  TimeDomainMeasurement m;
  m.r = r; 
  m.z = z; 
  m.what = what; 
  m.nskip = skip_factor; 
  if (m.nskip < 1) m.nskip = 1; 
  measurements.push_back(m); 
}

void iceprop::Sim::trackGlobalMaximum(meep::component what, ScalarType type )
{

  GlobalMaximum gm; 
  gm.what = what; 
  gm.type = type; 

  TString title; 
  title.Form("Maximum %s (%s)", 
              get_meep_component_name(what), get_scalar_name(type)); 
  TString name; 
  name.Form("max_%s_%s)", 
              get_meep_component_name(what), get_scalar_name(type)); 

  gm.max = make_hist(name.Data(), title.Data(), geom); 

  title.Form("Maximum Time %s (%s)", get_meep_component_name(what), get_scalar_name(type)); 
  name.Form("max_t_%s_%s)", get_meep_component_name(what), get_scalar_name(type)); 

  gm.tmax = make_hist(name.Data(), title.Data(), geom); 


  maxima.push_back(gm); 
}


/* dumb bookkeeping of output prefixes to avoid the never-will happen scenario
* of memory pressure from storing too many identical ones..*/ 

static const char * get_output_prefix(meep::component what, iceprop::ScalarType type)
{
    static std::map<std::pair<meep::component, iceprop::ScalarType>, char * > default_output_prefixes; 
    std::pair<meep::component,iceprop::ScalarType> pair(what,type);

    if (default_output_prefixes.count(pair)) 
      return default_output_prefixes[pair];

    char * str; 
    asprintf(&str,"%s_%s", get_meep_component_name(what), get_scalar_name(type)); 

    default_output_prefixes[pair] = str; 

    return str; 
}


static int count_the_canvases = 0; 

void iceprop::Sim::addStepOutput(StepOutput out)
{


  /* initialize stuff */ 
  if (!out.cw) out.cw = geom->max_r*geom->resolution * 1.1; 
  if (!out.ch) out.ch =(geom->max_depth + geom->sky_height)*geom->resolution * 1.1; 


  if (!out.output_prefix)
  {
    out.output_prefix = get_output_prefix(out.what, out.type); 
  }


  if (!out.c && ( out.format & (O_PNG | O_PDF)))
  {
    TString cname;
    cname.Form("c%d_%s", count_the_canvases++, out.output_prefix); 
    out.c = new TCanvas(cname.Data(), cname.Data(), out.cw, out.ch); 
    delete_list.push_back(out.c); 
  }

  if (!out.outf && (out.format & O_ROOT)) 
  {
    TString path;
    path.Form("%s/%s.root", out.out_dir, out.output_prefix); 
    out.outf = new TFile(path.Data(),"RECREATE"); //TODO: think if this what we want... 
    delete_list.push_back(out.outf); 
    out.outf->cd(); 
  }

  /* store a histogram to avoid constant reallocation */ 
  if (out.format & (O_ROOT | O_PNG | O_PDF)) 
  {
    out.h = make_hist(out.output_prefix, out.output_prefix, geom, out.double_precision); 
    if (!out.format & (O_ROOT))
    {
      delete_list.push_back(out.h); 
    }
  }

  if (out.format & O_ROOT) 
  {
    out.t = new TTree(out.output_prefix,out.output_prefix);  
    out.t->Branch("hist",out.h);   
    gROOT->cd(); 
  }


  outputs.push_back(out); 

}

double iceprop::Sim::getCurrentTime() const 
{
  return f->time() * meep_to_ns; 
}


struct store_calculation
{

  std::vector<std::complex<double> > v; 

  void fill( meep::component what, const meep::fields  * f) 
  {
    for (int j = 0; j < nbinsy; j++)
    {
      for (int i = 0; i < nbinsx; i++)
      {
        v[j*nbinsx+i] = f->get_field(what, meep::veccyl(i*dx,j*dx)); 
      }
    }
  }


  void updateHist(TH2 * h,  iceprop::ScalarType type)
  {
    double max = 0; 
    double min = 0; 
    for (int j = 0; j < nbinsy; j++)
    {
      for (int i = 0; i < nbinsx; i++ )
      {
        double v = get_scalar_val(type, get(i,j)); 
        if (v > max) max = v; 
        if (v < min) min = v; 
//        printf("%d %d %g\n",i,j,v);
        h->SetBinContent(i+1,j+1,v); 
      }
    }
    double scale = TMath::Max(fabs(max),fabs(min));
    h->SetMaximum(scale); 
    h->SetMinimum(-scale); 
  }

  void updateMax( TH2 * m, TH2 * t, iceprop::ScalarType type, double now) 
  {
    for (int j = 0; j < nbinsy; j++)
    {
      for (int i = 0; i < nbinsx; i++ )
      {
        double v = get_scalar_val(type, get(i,j)); 
        if (v > m->GetBinContent(i+1,j+1))
        {
          m->SetBinContent(i+1,j+1,v); 
          t->SetBinContent(i+1,j+1,now); 
        }
      }
    }
  }


  int nbinsx; 
  int nbinsy; 
  double dx; 

  std::complex<double> get(int i, int j) 
  {
    return v[ j * nbinsx + i]; 
  }

  store_calculation(const iceprop::SimGeometry * g) 
  {
    nbinsx = g->max_r * g->resolution; 
    nbinsy = (g->max_depth+g->sky_height) * g->resolution;  
    dx = 1./g->resolution;
    v.resize(nbinsx*nbinsy); 
  }


}; 



static double meep_scalar_function( const std::complex<double> * vals, const meep::vec & loc, void * type) 
{
  (void) loc; 
  return get_scalar_val( *((iceprop::ScalarType*)type), *vals); 
}


void iceprop::Sim::run(double time) 
{
  double meep_time = time * ns_to_meep;  
  double t0 = f->time(); 

  /* this is just the line for the surface in output plots */ 
  TGraph gsurf(2); 
  gsurf.SetPoint(0, 0,0); 
  gsurf.SetPoint(1, geom->max_r,0); 
 
  int i = 0; 

  /* bookkeeping has a lot of double letters, which is good, because we have a lot of doubles! */ 
  std::vector<store_calculation * > store(meep::NUM_FIELD_COMPONENTS,0);
  std::vector<meep::component> always_need_to_calculate;
  std::vector<bool> am_i_always_calculated(meep::NUM_FIELD_COMPONENTS,false); 
  always_need_to_calculate.reserve(maxima.size()); 


  /* If we look for maxima, we always have to calculate the entire thing everywhere. It blows */ 
  for (size_t m = 0; m < maxima.size(); m++)
  {
    if (!am_i_always_calculated[maxima[m].what])
    {
      always_need_to_calculate.push_back(maxima[m].what); 
      am_i_always_calculated[maxima[m].what] = true; 
    }
  }

  while (f->time() < t0 + meep_time) 
  {
    double now = getCurrentTime(); 
    std::vector<meep::component> need_to_calculate = always_need_to_calculate; 
    std::vector<bool> am_i_calculated = am_i_always_calculated; 

    /** figure out what additional things we need to fully calculate at this step */ 
    for (int o = 0; o < outputs.size(); o++)
    {
      if ( (i - outputs[o].skip_offset > 0 )  && ( (i - outputs[o].skip_offset) % outputs[o].skip_factor) == 0)
      {
         outputs[o].index++; //we will have to subtract 1. Otherwise we have to loop over again and increment later. 

        /* Use meep's internals for hdf5 since I think it's smarter than we are and I don't want to
         * figure this stuff out. */
        if (outputs[o].format & O_HDF5 ) 
        {
          TString path; 
          path.Form("%s/%s-%d.h5", outputs[o].out_dir, outputs[o].output_prefix, outputs[o].index-1); 
          meep::h5file h5(path.Data(),meep::h5file::WRITE); 

          f->output_hdf5( get_output_prefix(outputs[o].what, outputs[o].type),
                         1, &(outputs[o].what), &meep_scalar_function, (void*) &(outputs[o].type),
                         gv.surroundings(), 
                         &h5, false, !outputs[o].double_precision);
        }


        //we only need to calculate it if we need it for something other than O_HDF5 

        if ( (outputs[o].format & ( O_PNG | O_PDF | O_ROOT )) &&  ! am_i_calculated[outputs[o].what])
        {
          need_to_calculate.push_back(outputs[o].what); 
          am_i_calculated[outputs[o].what] = true; 
        }
      }
    }

    /* Compute the massive things we need to calculate
     *
     *  This currently doesn't do anything smart in the MPI case 
     *  (and honestly, all the ROOT stuff probably doesn't work in the MPI case anyway)
     *  If this becomes a problem I'll need to figure out how to use all that loop over chunks business
     *  or use an intermediate O_HDF5 file even if not outputting to one. 
     *
     **/ 
    for (size_t calc = 0;  calc < need_to_calculate.size(); calc++) 
    {
      meep::component what = need_to_calculate[calc]; 

      printf("Need to calculate %s at time %g\n", get_meep_component_name(what), getCurrentTime()); 

      /* make sure we have a store for the solution */ 
      if (! store[what])
      {
         printf("Making store for %s\n",get_meep_component_name(what)); 
         store[what] = new store_calculation(geom); 
      }

      store[what]->fill(what,f) ; 
    }

    /* Fill in maxima. This probably could be done while filling in the store for a bit of extra efficiency.
     * Will investigate if necessary */ 
         
    for (size_t m = 0; m < maxima.size(); m++) 
    {
      store[maxima[m].what]->updateMax(maxima[m].max, maxima[m].tmax, maxima[m].type,now); 
    }

    /* Output the non-O_HDF5 stuff. Remember that index should be one minus what we have */ 
    for (size_t o = 0; o < outputs.size(); o++) 
    {
      //not this step 
      if ( ! ((i - outputs[o].skip_offset > 0 )  && ( (i - outputs[o].skip_offset) % outputs[o].skip_factor) == 0)) 
        continue; 

      // only hdf5
      if ( ! (outputs[o].format & (O_PDF | O_PNG | O_ROOT)))
        continue; 

      //fill the histogram 
      store[outputs[o].what]->updateHist(outputs[o].h, outputs[o].type); 

      if (outputs[o].output_scale)
      {
        outputs[o].h->SetMaximum(outputs[o].output_scale); 
        outputs[o].h->SetMinimum(-outputs[o].output_scale); 
      }

      //change the title to what it is and when it is 

      if (outputs[o].format & (O_PNG | O_PDF) ) 
      {

        TString titl; 
        titl.Form("%s (t=%g, step=%d, index=%d)", get_output_prefix(outputs[o].what,outputs[o].type), 
                                                  getCurrentTime(), i, outputs[o].index-1); 
        outputs[o].h->SetTitle(titl.Data()); 
        outputs[o].c->cd(); 



        outputs[o].h->Draw("colz"); 

        /** add a line for the surface. We can do something smarter later probably */ 
        gsurf.Draw("lsame"); 


        TString path;  

        if (outputs[o].format & O_PNG) 
        {
          path.Form("%s/%s.%d.png", outputs[o].out_dir, outputs[o].output_prefix, outputs[o].index-1); 
          outputs[o].c->SaveAs(path); 
        }

        if (outputs[o].format & O_PDF) 
        {
          path.Form("%s/%s.%d.pdf", outputs[o].out_dir, outputs[o].output_prefix, outputs[o].index-1); 
          outputs[o].c->SaveAs(path); 
        }
      }

      if (outputs[o].format & O_ROOT)  
      {
        // this is cargo culting at this point for me, so I don't know if this is necessary or not 
        outputs[o].outf->cd(); 
        outputs[o].t->Fill(); 
        gROOT->cd(); 
      }
    }
          
    /** Fill in the time domain measurements */ 


    for (size_t m = 0; m < measurements.size(); m++) 
    {
       if (i % measurements[m].nskip == 0)
       {
         measurements[m].t.push_back(now); 
         measurements[m].val.push_back(f->get_field( measurements[m].what, meep::veccyl( measurements[m].r, measurements[m].z + geom->max_depth))); 
       }
    }

    f->step(); 
    i++; 
  }
} 


TGraph * iceprop::TimeDomainMeasurement::makeGraph(ScalarType type) const
{
  TGraph * g = new TGraph(t.size()); 

  for (int i =0; i < g->GetN(); i++) 
  {
    g->GetX()[i] = t[i]; 
    g->GetY()[i] = get_scalar_val(type,val[i]);
  }

  TString title; 
  title.Form("%s (%s) at r=%g m z=%g m", get_meep_component_name(what), get_scalar_name(type), r,z);  
  g->SetTitle(title.Data()); 
  title.Form("g_%s_%g_%g", get_output_prefix(what,type),r,z);
  g->SetName(title.Data()); 
  title.Form("%s (%s)", get_meep_component_name(what), get_scalar_name(type)); 
  g->GetYaxis()->SetTitle(title.Data()); 
  g->GetXaxis()->SetTitle("t (ns)"); 
  g->SetBit(TGraph::kIsSortedX) ; 
  g->SetBit(TGraph::kNotEditable); 

  return g; 
}



TH2 * iceprop::Sim::makeHist(meep::component what, ScalarType type, TH2 * h, bool double_precision) const 
{
  if (!h)
  {
    TString titl;
    titl.Form("%s t=%g", get_output_prefix(what,type), getCurrentTime()) ; 
    h = make_hist(get_output_prefix(what,type), titl.Data(), geom, double_precision); 
  }

  store_calculation s(geom); 
  s.fill(what,f); 
  s.updateHist(h,type); 
  return h; 
}

void  iceprop::Sim::draw(meep::component what, ScalarType typ ) const
{
  TH2 * h = makeHist(what,typ); 
  h->DrawCopy("colz"); 
  delete h;
}






