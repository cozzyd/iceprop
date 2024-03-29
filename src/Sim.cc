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
#include <boost/container/flat_map.hpp>

#include "TCanvas.h" 
#include "TFile.h" 
#include "TTree.h" 
#include "TH2.h" 
#include "TMutex.h" 
#include "TROOT.h" 
#include <dirent.h> 
#include "TSystem.h" 
#include "TGraph2D.h" 
#include "TChain.h" 
#include <sys/types.h> 
#include <unistd.h> 





/** necessary for epsilon with no void */
static iceprop::Firn * the_firn; 
static iceprop::SimGeometry * the_geom;

static int pid = 0; 

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
    type == meep::Dr          ?   "Dr"    : 
    type == meep::Dz          ?   "Dz"    : 
    type == meep::Dp          ?   "Dphi"  : 
    type == meep::Dx          ?   "Dx"    : 
    type == meep::Dy          ?   "Dy"    : 
    type == meep::Br          ?   "Br"    : 
    type == meep::Bz          ?   "Bz"    : 
    type == meep::Bp          ?   "Bphi"  : 
    type == meep::Bx          ?   "Bx"    : 
    type == meep::By          ?   "By"    : 
 
    type == meep::Dielectric  ?   "eps"   : 
    0; 
}



static TH2 * make_hist(const char * name, const char * title, const iceprop::SimGeometry * geom, bool double_precision = true) 
{
  TH2 * h = 
  double_precision ? 
      (TH2*) new TH2D(name, title, geom->max_r*geom->resolution/geom->output_skip_factor, 0, geom->max_r, 
                    (geom->max_depth+geom->sky_height)*geom->resolution/geom->output_skip_factor, -geom->max_depth, geom->sky_height) 
     :
      (TH2*) new TH2F(name, title, geom->max_r*geom->resolution/geom->output_skip_factor, 0, geom->max_r, 
                   (geom->max_depth+geom->sky_height)*geom->resolution/geom->output_skip_factor, -geom->max_depth, geom->sky_height); 
 
  h->GetXaxis()->SetTitle(geom->is_cylindrical ? "r (m)" : "x (m)" ); 
  h->GetYaxis()->SetTitle("z (m)"); 
  h->SetStats(0); 

  return h; 

}


iceprop::Sim::~Sim()
{
  delete f; 
  delete s; 

  if (loaded_measurements) delete loaded_measurements; 

  /* make sure we write the trees */ 

  for (size_t o = 0; o < outputs.size(); o++)
  {
    if (outputs[o]->t) 
    {
      outputs[o]->outf->cd(); 
      outputs[o]->t->Write(); 
    }
    delete outputs[o]; 
  }
 

  if (intermediate_globals_file) 
  {
    intermediate_globals_file->cd(); 
    intermediate_globals_tree->Write(); 
    delete intermediate_globals_file; 
  }


  for (size_t i =0; i < delete_list.size(); i++) 
  {
    if (delete_list[i]) delete delete_list[i]; 
  }
}


iceprop::Sim::Sim(const Firn * firn, const SimGeometry * geom, const Source * source) 
:  firn(firn), geom(geom), gsurf(2), 
   gv(geom->is_cylindrical ? meep::volcyl(geom->max_r, geom->max_depth + geom->sky_height, geom->resolution) 
                           : meep::vol2d(geom->max_r, geom->max_depth + geom->sky_height, geom->resolution))
{

  static TMutex construction_lock; 
  TLockGuard l (&construction_lock); 
  //evil thread-unsafe code here , will probably need to fx
  the_firn = (Firn*)firn; 
  the_geom = (SimGeometry*)geom; 
  s = new meep::structure(gv, eps, meep::pml(geom->pml_size), 
                          meep::identity(), 0,geom->courant_factor); 
  loaded_measurements = 0; 
  f = new meep::fields(s); 
  if (source)  addSource(source); 
  gsurf.SetPoint(0, 0,0); 
  gsurf.SetLineWidth(4); 
  gsurf.SetPoint(1, geom->max_r,0); 
  intermediate_globals_file = 0; 
  intermediate_globals_tree = 0; 
  intermediate_globals_interval = 100; 

  if (!pid) pid = getpid(); 
  snapshot_nsteps = 0; 
  snapshot_path = 0; 
  snapshot_i = 0; 
  step =  0; 
  time_offset = 0; 
  snapshot_loaded = false; 
}


int iceprop::Sim::saveIntermediateGlobals(const char * file, int skip_factor)
{
  if (intermediate_globals_file) 
  {
    fprintf(stderr,"Already saving intermediate globals to %s\n", intermediate_globals_file->GetName()); 
    return 1; 
  }

  if (!mpi::am_master()) return 0; 


  intermediate_globals_file = new TFile(file,snapshot_loaded ? "UPDATE" : "RECREATE"); 
  intermediate_globals_tree = new TTree("globals","globals"); 
  intermediate_globals_tree->SetBranchAddress("pid",&pid); 
  intermediate_globals_tree->SetAutoSave(10); 
  gROOT->cd(); 

  for (size_t i = 0; i < maxima.size(); i++)
  {

    intermediate_globals_tree->Branch(maxima[i].max->GetName(), &maxima[i].max); 
    intermediate_globals_tree->Branch(maxima[i].tmax->GetName(), &maxima[i].tmax); 
  }

  for (size_t i = 0; i < integrals.size(); i++)
  {

    intermediate_globals_tree->Branch(integrals[i].integ->GetName(), &integrals[i].integ); 
    intermediate_globals_tree->Branch(integrals[i].tfirst->GetName(), &integrals[i].tfirst); 
  }

  return 0; 
}


void iceprop::Sim::addSource(const Source * s) 
{
  f->add_point_source(s->getComponent(), s->getSource(), 
      geom->is_cylindrical ? meep::veccyl(s->getR(), geom->max_depth + s->getZ())
      : meep::vec(s->getR(), geom->max_depth + s->getZ())); 
}


void iceprop::Sim::addTimeDomainMeasurement(meep::component what, double r, double z, const char * name, int skip_factor) 
{
  TimeDomainMeasurement mm;
  measurements.emplace_back(mm); 
  TimeDomainMeasurement & m = measurements.back(); 
  m.r = r; 
  m.z = z; 
  m.what = what; 
  m.nskip = skip_factor; 
  m.name = name; 
  if (m.nskip < 1) m.nskip = 1; 


  if (snapshot_loaded) 
  {
    if (!loaded_measurements) loaded_measurements = new TFile(loaded_measurements_path.Data()); 
    gROOT->cd(); 

    TString key; 
    m.makeName(key,Real); 
    TGraph * greal = (TGraph*) loaded_measurements->Get(key.Data()); 
    m.makeName(key,Imag); 
    TGraph * gimag = (TGraph*) loaded_measurements->Get(key.Data()); 

    m.t.reserve(greal->GetN()); 
    m.val.reserve(greal->GetN()); 
    for (int i = 0; i < greal->GetN(); i++) 
    {
      m.t.push_back(greal->GetX()[i]); 
      m.val.push_back(std::complex<double>(greal->GetY()[i], gimag->GetY()[i])); 
    }

    delete greal;
    delete gimag; 
  }

}

void iceprop::Sim::trackGlobalIntegral(meep::component what, ScalarType type )
{

  GlobalIntegral gm; 
  gm.what = what; 
  gm.type = type; 

  TString title; 
  title.Form("Integral %s (%s)", 
              get_meep_component_name(what), get_scalar_name(type)); 
  TString name; 
  name.Form("integ_%s_%s", 
              get_meep_component_name(what), get_scalar_name(type)); 

  gm.integ = make_hist(name.Data(), title.Data(), geom); 

  title.Form("First Time %s (%s)", get_meep_component_name(what), get_scalar_name(type)); 
  name.Form("first_%s_%s", get_meep_component_name(what), get_scalar_name(type)); 

  gm.tfirst = make_hist(name.Data(), title.Data(), geom); 
  integrals.push_back(gm); 

  if (intermediate_globals_tree && mpi::am_master()) 
  {
    intermediate_globals_tree->Branch(gm.integ->GetName(), &integrals[integrals.size()-1].integ); 
    intermediate_globals_tree->Branch(gm.tfirst->GetName(), &integrals[integrals.size()-1].tfirst); 
  }

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
  name.Form("max_%s_%s", 
              get_meep_component_name(what), get_scalar_name(type)); 

  gm.max = make_hist(name.Data(), title.Data(), geom); 

  title.Form("Maximum Time %s (%s)", get_meep_component_name(what), get_scalar_name(type)); 
  name.Form("max_t_%s_%s", get_meep_component_name(what), get_scalar_name(type)); 

  gm.tmax = make_hist(name.Data(), title.Data(), geom); 

  maxima.push_back(gm); 
  if (intermediate_globals_tree && mpi::am_master()) 
  {
    intermediate_globals_tree->Branch(gm.max->GetName(), &maxima[maxima.size()-1].max); 
    intermediate_globals_tree->Branch(gm.tmax->GetName(), &maxima[maxima.size()-1].tmax); 
  }



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



void iceprop::Sim::addStepOutput(StepOutput o)
{

  StepOutput * out = new StepOutput(o); 
  outputs.push_back(out); 

  out->index = 0; 
  /* initialize stuff */ 
  if (!out->cw) out->cw = geom->max_r*geom->resolution * 1.1; 
  if (!out->ch) out->ch =(geom->max_depth + geom->sky_height)*geom->resolution * 1.1; 


  if (!out->output_prefix && mpi::am_master())
  {
    out->output_prefix = get_output_prefix(out->what, out->type); 
  }

  if (!out->c && ( out->format & (O_PNG | O_PDF))  && mpi::am_master())
  {
    TString cname;
    cname.Form("c%d_%s", count_the_canvases++, out->output_prefix); 
    out->c = new TCanvas(cname.Data(), cname.Data(), out->cw, out->ch); 
    delete_list.push_back(out->c); 
  }

  TString newpath; 
  bool found_old = false; 
  if (!out->outf && (out->format & O_ROOT) && mpi::am_master()) 
  {
    TString path;
    path.Form("%s/%s.root", out->out_dir, out->output_prefix); 

    // Since I can't get the TTree to update properly, I'll cheat and move the old file out of the way. 

    if (!access(path.Data(),R_OK)) found_old = true; 
    if (snapshot_loaded && found_old) 
    {
      newpath.Form("%s/%s.%d.root", out->out_dir, out->output_prefix, pid); 
      gSystem->Rename(path.Data(), newpath.Data()); 
    }
    out->outf = new TFile(path.Data(),"RECREATE"); 
    delete_list.push_back(out->outf); 
  }

  /* store a histogram to avoid constant reallocation */ 
  if ((out->format & (O_ROOT | O_PNG | O_PDF)) && mpi::am_master()) 
  {
    if (out->outf) out->outf->cd(); 

    out->h =  make_hist(out->output_prefix, out->output_prefix, geom, out->double_precision); 

    if (!out->outf) delete_list.push_back(out->h); 
  }

  if ((out->format & O_ROOT) && mpi::am_master()) 
  {
    out->outf->cd(); 
    out->t =  new TTree(out->output_prefix,out->output_prefix);  

    out->t->Branch("hist",out->h);   
    out->t->Branch("pid",&pid);   
    out->t->Branch("index",&out->index); 
    out->t->Branch("step",&step);   
 
    if (snapshot_loaded && found_old) //copy over the old stuff
    {
      printf("Copying tree!\n"); 
      int old_step = step; 
      int old_pid = pid; 
      int old_index = out->index; 

      TFile old(newpath); 
      TTree * ot = (TTree*) old.Get(out->output_prefix); 
      ot->SetBranchAddress("hist",&out->h);   
      ot->SetBranchAddress("pid",&pid);   
      ot->SetBranchAddress("index",&out->index); 
      ot->SetBranchAddress("step",&step);   

      for (int i = 0; i < ot->GetEntries(); i++) 
      {
          ot->GetEntry(i); 
          out->outf->cd(); 
          out->t->Fill(); 
      }

      pid = old_pid; 
      out->index = old_index; 
      step = old_step; 
    }
  }

  gROOT->cd(); 

  if (snapshot_loaded) 
  {
    out->index = step < out->skip_offset ? 0 : (step-out->skip_offset) /out->skip_factor; 
  } 
  else out->index = 0; 



}

double iceprop::Sim::getCurrentTime() const 
{
  return f->time() * meep_to_ns + time_offset; 
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
        v[j*nbinsx+i] = f->get_field(what, cylindrical ? meep::veccyl(i*dx,j*dx)
                                                       : meep::vec(i*dx,j*dx)); 
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
    h->SetMaximum(scale*1.1); 
    h->SetMinimum(-scale*1.1); 
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
  void updateInteg( TH2 * m, TH2 * t, iceprop::ScalarType type, double now) 
  {
    for (int j = 0; j < nbinsy; j++)
    {
      for (int i = 0; i < nbinsx; i++ )
      {
        double v = get_scalar_val(type, get(i,j)); 

        double content = m->GetBinContent(i+1,j+1); 
        if (v && !content)
        {
          t->SetBinContent(i+1,j+1,now); 
        }

        m->SetBinContent(i+1,j+1,content+v); 

      }
    }
  }



  int nbinsx; 
  int nbinsy; 
  double dx; 
  bool cylindrical; 

  std::complex<double> get(int i, int j) 
  {
    return v[ j * nbinsx + i]; 
  }

  store_calculation(const iceprop::SimGeometry * g) 
  {
    int nskip = g->output_skip_factor; 
    nbinsx = (g->max_r * g->resolution) / nskip; 
    nbinsy = ((g->max_depth+g->sky_height) * g->resolution) / nskip;  
    dx = 1./g->resolution * nskip;
    cylindrical = g->is_cylindrical; 
    v.resize(nbinsx*nbinsy); 
  }

}; 



static double meep_scalar_function( const std::complex<double> * vals, const meep::vec & loc, void * type) 
{
  (void) loc; 
  return get_scalar_val( *((iceprop::ScalarType*)type), *vals); 
}


static std::complex<double> meep_component_function ( const std::complex<double> * vals, const meep::vec & loc, void * type) 
{
  (void) loc; 
  (void) type; 
  return *vals; 
}


void iceprop::Sim::run(double time, bool relative) 
{

  bool first = true; 
  double t0 = relative ? getCurrentTime() : 0; 

  /* this is just the line for the surface in output plots */ 


  /* bookkeeping has a lot of double letters, which is good, because we have a lot of doubles! */ 
  std::vector<store_calculation * > store(meep::NUM_FIELD_COMPONENTS,0);
  std::vector<meep::component> always_need_to_calculate;
  std::vector<bool> am_i_always_calculated(meep::NUM_FIELD_COMPONENTS,false); 
  always_need_to_calculate.reserve(maxima.size()); 



  /* if we look for maxima, we always have to calculate the entire thing everywhere. it blows */ 
  for (size_t m = 0; m < maxima.size(); m++)
  {
    if (!am_i_always_calculated[maxima[m].what])
    {
      always_need_to_calculate.push_back(maxima[m].what); 
      am_i_always_calculated[maxima[m].what] = true; 
    }
  }

  /* ditto for integral */ 
  for (size_t m = 0; m < integrals.size(); m++)
  {
    if (!am_i_always_calculated[integrals[m].what])
    {
      always_need_to_calculate.push_back(integrals[m].what); 
      am_i_always_calculated[integrals[m].what] = true; 
    }
  }


  while (getCurrentTime() < t0 + time) 
  {
    double now = getCurrentTime(); 
    std::vector<meep::component> need_to_calculate = always_need_to_calculate; 
    std::vector<bool> am_i_calculated = am_i_always_calculated; 

    /** figure out what additional things we need to fully calculate at this step */ 
    for (int o = 0; o < outputs.size(); o++)
    {
      if ( (step - outputs[o]->skip_offset >= 0 )  && ( (step - outputs[o]->skip_offset) % outputs[o]->skip_factor) == 0)
      {
         outputs[o]->index++; //we will have to subtract 1. Otherwise we have to loop over again and increment later. 

        /* Use meep's internals for hdf5 since I think it's smarter than we are and I don't want to
         * figure this stuff out. */
        if (outputs[o]->format & O_HDF5 ) 
        {
          TString path; 
          path.Form("%s/%s-%d.h5", outputs[o]->out_dir, outputs[o]->output_prefix, outputs[o]->index-1); 
          meep::h5file h5(path.Data(),meep::h5file::WRITE); 

          f->output_hdf5( get_output_prefix(outputs[o]->what, outputs[o]->type),
                         1, &(outputs[o]->what), &meep_scalar_function, (void*) &(outputs[o]->type),
                         gv.surroundings(), 
                         &h5, false, !outputs[o]->double_precision);
        }


        //we only need to calculate it if we need it for something other than O_HDF5 

        if ( (outputs[o]->format & ( O_PNG | O_PDF | O_ROOT )) &&  ! am_i_calculated[outputs[o]->what])
        {
          need_to_calculate.push_back(outputs[o]->what); 
          am_i_calculated[outputs[o]->what] = true; 
        }
      }
    }

    for (size_t calc = 0;  calc < need_to_calculate.size(); calc++) 
    {
      meep::component what = need_to_calculate[calc]; 

      if (mpi::am_master())
      {
        printf("Need to calculate %s at time %g\n", get_meep_component_name(what), getCurrentTime()); 
      }

      /* make sure we have a store for the solution */ 
      if (! store[what])
      {
        if (mpi::am_master())
        {
           printf("Making store for %s\n",get_meep_component_name(what)); 
        }
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

    for (size_t m = 0; m < integrals.size(); m++) 
    {
      store[integrals[m].what]->updateInteg(integrals[m].integ, integrals[m].tfirst, maxima[m].type,now); 
    }


    if (intermediate_globals_tree && (step % intermediate_globals_interval == 0) && mpi::am_master())
    {
      intermediate_globals_file->cd(); 
      intermediate_globals_tree->Fill(); 
      gROOT->cd(); 
    }

    /* Output the non-O_HDF5 stuff. Remember that index should be one minus what we have */ 
    for (size_t o = 0; o < outputs.size(); o++) 
    {
      if (!mpi::am_master()) break; 

      //not this step 
      if ( ! ((step - outputs[o]->skip_offset >= 0 )  && ( (step - outputs[o]->skip_offset) % outputs[o]->skip_factor) == 0)) 
        continue; 

      // only hdf5
      if ( ! (outputs[o]->format & (O_PDF | O_PNG | O_ROOT)))
        continue; 

      //fill the histogram 
      store[outputs[o]->what]->updateHist(outputs[o]->h, outputs[o]->type); 

      if (outputs[o]->output_scale)
      {
        outputs[o]->h->SetMaximum(outputs[o]->output_scale); 
        outputs[o]->h->SetMinimum(-outputs[o]->output_scale); 
      }

      //change the title to what it is and when it is 
      if (outputs[o]->format & (O_PNG | O_PDF | O_ROOT) ) 
      {
        TString titl; 
        titl.Form("%s (t=%g, step=%d, index=%d)", get_output_prefix(outputs[o]->what,outputs[o]->type), 
                                                  getCurrentTime(), step, outputs[o]->index-1); 
      
        outputs[o]->h->SetTitle(titl.Data()); 

      }
 
      if (outputs[o]->format & (O_PNG | O_PDF ) ) 
      {

        outputs[o]->c->cd(); 



        outputs[o]->h->Draw("colz"); 

        /** add a line for the surface. We can do something smarter later probably */ 
        gsurf.Draw("lsame"); 


        TString path;  

        if (outputs[o]->format & O_PNG) 
        {
          path.Form("%s/%s.%d.png", outputs[o]->out_dir, outputs[o]->output_prefix, outputs[o]->index-1); 
          outputs[o]->c->SaveAs(path); 
        }

        if (outputs[o]->format & O_PDF) 
        {
          path.Form("%s/%s.%d.pdf", outputs[o]->out_dir, outputs[o]->output_prefix, outputs[o]->index-1); 
          outputs[o]->c->SaveAs(path); 
        }
      }

      if (outputs[o]->format & O_ROOT)  
      {
       
        // this is cargo culting at this point for me, so I don't know if this is necessary or not 
       
//        bool hist_status = TH1::AddDirectoryStatus(); 
//        TH1::AddDirectory(false); 
        
        outputs[o]->outf->cd(); 
        outputs[o]->index--; 
        outputs[o]->t->Fill(); 
        outputs[o]->index++; 
        gROOT->cd(); 
 //       TH1::AddDirectory(hist_status); 
      }
    }
          
    /** Fill in the time domain measurements */ 


    for (size_t m = 0; m < measurements.size(); m++) 
    {
       if (step % measurements[m].nskip == 0)
       {
         measurements[m].t.push_back(now); 
         measurements[m].val.push_back(f->get_field( measurements[m].what,
                                         geom->is_cylindrical ? meep::veccyl( measurements[m].r, measurements[m].z + geom->max_depth)
                                                              : meep::vec( measurements[m].r, measurements[m].z + geom->max_depth)
                                         )); 
       }
    }


    if (snapshot_nsteps && (step % snapshot_nsteps == 0) && !first) snapshot(snapshot_path, snapshot_i++); 

    f->step(); 
    step++;
    first = false; 
  }
} 



void iceprop::TimeDomainMeasurement::makeName(TString & str, ScalarType type) const
{
  if (name) 
  {
    str.Form("%s_%s",name,get_output_prefix(what,type)); 
  }
  else
  {
    str.Form("g_%s_%g_%g", get_output_prefix(what,type),r,z);
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
  title.Form("%s %s (%s) at r=%g m z=%g m", name? name:"", get_meep_component_name(what), get_scalar_name(type), r,z);  
  g->SetTitle(title.Data()); 
  makeName(title,type); 
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


#define snapshot_index_format "%06d" 
static int getNewestSnapshotIndex(const char * path) 
{
  struct dirent * de; 
  DIR * dr = opendir(path); 

  if (dr==NULL) 
  {
    fprintf(stderr,"Could not open %s\n", path); 
    return -1; 
  }

  int max_i = -1; 
  while ((de = readdir(dr)) != NULL) 
  {
    int i; 
    if ( 1 == sscanf(de->d_name, snapshot_index_format,&i)) 
    {
      char testbuf[128]; 
      snprintf(testbuf,128,snapshot_index_format, i); 
      if (!strcmp(testbuf,de->d_name))
      {
        if (i > max_i) 
        {
          max_i = i; 
        }
      }
    }
  }

  return max_i; 
}



static meep::component cyl_components[] = {meep::Dr, meep::Dz, meep::Dp,meep::Br, meep::Bz, meep::Bp}; 
static meep::component twod_components[] = {meep::Ex, meep::Ey, meep::Ez, meep::Hz, meep::Hy, meep::Hz}; 

struct chunk_data
{
  bool is_cylindrical; 
}; 

static void chunk_output_tree(meep::fields* f, chunk_data *data) 
{
  chunk_data * cd = (chunk_data*) data; 

  int nfields = cd->is_cylindrical ? sizeof(cyl_components)/ sizeof(*cyl_components)
                                    : sizeof(twod_components)/ sizeof(*twod_components); 
  const meep::component * components = cd->is_cylindrical ? cyl_components : twod_components; 
           
  for (int ifield = 0; ifield < nfields; ifield++) 
  {
    meep::component c = components[ifield]; 

    TTree * t = new TTree( get_meep_component_name(c), get_meep_component_name(c)); 
    t->SetAutoSave(0); 
    t->SetAutoFlush(0); 
    double val[2]; 
    double x[2]; 

    t->Branch("val_real",val); 
    t->Branch("val_imag",val+1); 
    t->Branch(cd->is_cylindrical ? "r" : "x" ,x); 
    t->Branch(cd->is_cylindrical ? "z" : "y" ,x+1); 

    using namespace meep;
    for (int ichunk=0; ichunk < f->num_chunks; ichunk++)
    {
      meep::fields_chunk * chunk = f->chunks[ichunk]; 
      if (!chunk->f[c][0])
      {
          printf("No %s in chunk %d\n", get_meep_component_name(c), ichunk); 
          break; 
      }
 
      LOOP_OVER_VOL(chunk->gv, c, i)
      {
        IVEC_LOOP_LOC(chunk->gv,here); 
        x[0] = cd->is_cylindrical ? here.r() : here.x(); 
        x[1] = cd->is_cylindrical ? here.z() : here.y(); 
        val[0] = chunk->f[c][0][i]; 
        val[1] = chunk->f[c][1][i]; 
        t->Fill(); 
      }
      t->Write(0, TObject::kOverwrite); 
    }
  }
}

int iceprop::Sim::snapshot(const char * path, int i) 
{
  //make the path if it doesn't exist 
  gSystem->mkdir(path,true); 
  
  //If i < 0, see if there is any older snapshot 
  if (i < 0) 
  {
    i = getNewestSnapshotIndex(path); 
    if (i < 0)
    {
      i = 0; 
    }
    else i++; 
  }

  TString str; 
  str.Form("%s/" snapshot_index_format , path, i); 

  if (mpi::am_master()) 
  {
    printf("Making snapshot at %s\n", str.Data()); 
   gSystem->mkdir(str.Data(),true); 

    //save the current simulation time
    //this could be saved more smartlyfile, but it's nice to have it more easily human readable
    str.Form("%s/" snapshot_index_format"/time",path,i); 
    FILE * fi = fopen(str.Data(),"w"); 
    fprintf(fi,"step: %d\n", step); 
    fprintf(fi,"t_ns: %g\n", getCurrentTime()); 
    fprintf(fi,"do not edit this file if you don't know what you're doing... the parser is stupid\n"); 
    fclose(fi); 

  }

  //now we save the fields 
  
  int nfields = geom->is_cylindrical ? sizeof(cyl_components)/ sizeof(*cyl_components)
                                    : sizeof(twod_components)/ sizeof(*twod_components); 
  const meep::component * components = geom->is_cylindrical ? cyl_components : twod_components; 
 
  mpi::barrier(); 

  str.Form("%s/" snapshot_index_format "/fields%d.root",path,i, mpi::rank()); 
  //we will write this in the most inefficient form possible :) 
  {
    TFile fi(str.Data(),"CREATE"); 

    for (int ifield = 0; ifield < nfields; ifield++) 
    {
      meep::component c = components[ifield]; 

      TTree * t = new TTree( get_meep_component_name(c), get_meep_component_name(c)); 
      t->SetAutoSave(0); 
      t->SetAutoFlush(0); 
      double val[2]; 
      double x[2]; 

      t->Branch("val_real",val); 
      t->Branch("val_imag",val+1); 
      t->Branch(geom->is_cylindrical ? "r" : "x" ,x); 
      t->Branch(geom->is_cylindrical ? "z" : "y" ,x+1); 

      using namespace meep;
      for (int ichunk=0; ichunk < f->num_chunks; ichunk++)
      {
        meep::fields_chunk * chunk = f->chunks[ichunk]; 
        if (!chunk->f[c][0])
        {
//            printf("No %s in chunk %d\n", get_meep_component_name(c), ichunk); 
            continue; 
        }
   
        LOOP_OVER_VOL(chunk->gv, c, i)
        {
          IVEC_LOOP_LOC(chunk->gv,here); 
          x[0] = geom->is_cylindrical ? here.r() : here.x(); 
          x[1] = geom->is_cylindrical ? here.z() : here.y(); 
          val[0] = chunk->f[c][0][i]; 
          val[1] = chunk->f[c][1][i]; 
          t->Fill(); 
        }

//        t->BuildIndex( geom->is_cylindrical ? " int(r * 1000)" : "int(x * 1000)" , 
//                       geom->is_cylindrical ? " int(z * 1000)" : "int(y * 1000)"  ); 
        t->Write(0, TObject::kOverwrite); 
      }
    }
  }

  TString cmd; 
  cmd.Form("gzip %s", str.Data()); 
  system(cmd.Data()); 


  if (mpi::am_master()) 
  {
    //flush the outputs tree 
    for (size_t o =0; o< outputs.size(); o++) 
    {
      if (outputs[o]->format & O_ROOT)
      {
        printf("Flushing %s\n", outputs[o]->outf->GetName()); 
        outputs[o]->outf->cd(); 
        outputs[o]->t->FlushBaskets(); 
        outputs[o]->t->AutoSave(); 
        outputs[o]->outf->SaveSelf(); 
        gROOT->cd(); 
      }
    }

    //and save the time domain measurements! 
    str.Form("%s/" snapshot_index_format "/measurements.root",path,i); 
    TFile of(str.Data(),"RECREATE"); 

    for (unsigned m = 0; m < getMeasurements().size(); m++) 
    {
      TGraph * gr = getMeasurements()[m].makeGraph(Real); 
      TGraph * gi = getMeasurements()[m].makeGraph(Imag); 
      gr->Write(); 
      gi->Write(); 
      delete gr; 
      delete gi; 
    }
  }
  return 0; 
}

static bool comp( const std::pair<std::pair<double,double>, std::complex<double> > & a, 
                    const std::pair<std::pair<double,double>, std::complex<double> > & b) 
                {
                  return a.first < b.first; 
                }; 


int iceprop::Sim::loadSnapshot(const char * path, int i) 
{
  if (i < 0) 
  {
    i = getNewestSnapshotIndex(path); 
    if (i < 0)
    {
      fprintf(stderr,"No snapshots in %s!\n", path); 
      return 1; 
    }
  }

  //check if we can access the directory 
  TString str; 
  str.Form("%s/" snapshot_index_format , path, i); 
  if (access(str.Data(), R_OK | X_OK))
  {
    fprintf(stderr,"Could not access %s\n", str.Data()); 
    return 1; 
  }

  // grab the time 
  str.Form("%s/" snapshot_index_format "/time",path,i); 
  FILE * fi = fopen(str.Data(),"r"); 
  if (!fi) 
  {
    fprintf(stderr,"Could not open %s\n", str.Data()); 
    return 1; 

  }

  fscanf(fi,"step: %d\n", &step); 
  f->t = step; 
//  fscanf(fi,"t_ns: %lg", &time_offset); 
  fclose(fi); 

  TString cmd; 
  if (mpi::am_master())
  {
    printf("LOADING SNAPSHOT FROM %s\n", path); 
    cmd.Form("cd  %s/" snapshot_index_format "; for i in *.root.gz; do echo ${i} ->${i%%.gz};  gunzip -c ${i} > ${i%%.gz}; done", path,i); 
    printf("%s\n", cmd.Data()); 
    system(cmd.Data()); 
  }


  mpi::barrier(); 

  //now we load the fields 
  {

    int nfields = geom->is_cylindrical ? sizeof(cyl_components)/ sizeof(*cyl_components)
                                       : sizeof(twod_components)/ sizeof(*twod_components); 
    
    const meep::component * components = geom->is_cylindrical ? cyl_components : twod_components; 

    //this is probably super super super slow. Oh well. 
    for (int ifield = 0; ifield < nfields; ifield++) 
    {

      // this will have min,max x0, min,max x1 to a file name
      static std::vector<std::pair<std::array<double,4>, TString> > _tree_map; 

      _tree_map.clear(); 

      //loop over all the trees 
      int ich = 0; 
      while(true) 
      {
        str.Form("%s/" snapshot_index_format "/fields%d.root",path,i,ich); 

        if (access(str.Data(),W_OK)) break; // we are out of files! 

        std::array<double,4>  X;


        if (mpi::am_master()) 
        {

          //TODO: this would go substantially faster if each process took care of a tree, but that's a bit more complicated 
          TFile ft(str,mpi::am_master ? "UPDATE" : "READ"); 
          TTree * tr = (TTree*) ft.Get(get_meep_component_name(components[ifield])); 
          X[0] = tr->GetMinimum(geom->is_cylindrical ? "r" : "x"); 
          X[1] = tr->GetMaximum(geom->is_cylindrical ? "r" : "x"); 
          X[2] = tr->GetMinimum(geom->is_cylindrical ? "z" : "y"); 
          X[3] = tr->GetMaximum(geom->is_cylindrical ? "z" : "y"); 

          printf("%s %s %d  %d %g %g %g %g\n", str.Data(),get_meep_component_name(components[ifield]),ich, tr->GetEntries(),  X[0], X[1], X[2], X[3] ); 
          if (!tr->GetTreeIndex()) 
          {
            printf("BUILDING INDEX FOR %s in  %s\n", get_meep_component_name(components[ifield]),str.Data()); 
            tr->BuildIndex( geom->is_cylindrical ? " int(r * 1000)" : "int(x * 1000)" , 
                            geom->is_cylindrical ? " int(z * 1000)" : "int(y * 1000)"  ); 

            tr->Write(); 

          }
        }

        mpi::broadcast(&X[0], 4); 
        _tree_map.push_back(std::pair<std::array<double,4>, TString> (X, str)); 
        ich++; 
      }



      double val[2]; 
      double x[2]; 

      static bool _is_cyl; 
      static meep::component _c; 
      _c = components[ifield]; 

      _is_cyl = geom->is_cylindrical; 

      static TFile * _cur_file = 0; 
      static TTree * _cur_tree = 0; 
      static double  _vreal = 0; 
      static double  _vimag = 0; 
      static double  _x0 = 0; 
      static double  _x1 = 0; 
      static double _cur_min0 =0; 
      static double _cur_min1 =0; 
      static double _cur_max0 =0; 
      static double _cur_max1 =0; 



      f->initialize_field(components[ifield], 
           [] (const meep::vec & v) -> std::complex<double> 
           {

              
              double x0 = _is_cyl ? v.r() : v.x(); 
              double x1 = _is_cyl ? v.z() : v.y(); 

              if (!_cur_file || x0 < _cur_min0 || x0 > _cur_max0 || x1 < _cur_min1 || x1 > _cur_max1)
              {
                for (auto it : _tree_map) 
                {
                  std::array<double,4> X = it.first; 

                  if ( x0 >= X[0] && x0 <= X[1] && x1 >= X[2] && x1<=X[3])
                  {
                    _cur_min0 = X[0]; 
                    _cur_min1 = X[2]; 
                    _cur_max0 = X[1]; 
                    _cur_max1 = X[3]; 
                    if (_cur_file) delete _cur_file; 

                    // do this one process at a time otherwise it's too slow
                    for (int rank = 0; rank < mpi::max_rank(); rank++) 
                    {
                      if (mpi::rank() == rank) 
                      {
                        printf("Rank %d switching to file %s\n", rank, it.second.Data()); 
                        _cur_file = new TFile(it.second.Data()); 
                        _cur_tree = (TTree*) _cur_file->Get(get_meep_component_name(_c)); 
                      }
                    }

                    _cur_tree->SetBranchAddress(_is_cyl ? "r" : "x", &_x0); 
                    _cur_tree->SetBranchAddress(_is_cyl ? "z" : "y", &_x1); 
                    _cur_tree->SetBranchAddress("val_real",&_vreal); 
                    _cur_tree->SetBranchAddress("val_imag",&_vimag); 
                  }
                }
              }
      
              _cur_tree->GetEntryWithIndex((int) (x0*1000), (int) (x1*1000)); 

              if (x0!= _x0 || x1!=_x1) printf("uh oh!!: %g %g %g %g\n", x0,_x0,x1,_x1);  

              return std::complex<double> (_vreal,_vimag); 

           }); 


      delete _cur_file; 
      _cur_file = 0; 

    }

  }

  loaded_measurements_path.Form("%s/" snapshot_index_format "/measurements.root",path,i); 

  if (mpi::am_master())
  {
   cmd.Form("rm %s/" snapshot_index_format "/fields*.root", path,i); 
   system(cmd.Data()); 
  }


  snapshot_loaded = true; 
//  snapshot(path); //to compare
  return 0; 
}


void iceprop::Sim::enableSnapshots(const char * path, int nsteps, int i) 
{
  gSystem->mkdir(path,true); 
  if (i < 0) 
  {
    i = getNewestSnapshotIndex(path); 
    if (i < 0) i = 0; 
    else i++; 
  }
  snapshot_i = i; 
  snapshot_path = path; 
  snapshot_nsteps = nsteps; 
}




