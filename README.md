#iceprop  
---

Cosmin Deaconu <cozzyd@kicp.uchicago.edu>

My attempt at using meep (http://ab-initio.mit.edu/wiki/index.php/Meep) to
simulate radio propagation in firn This package has a quite narrow purpose. If
you want to do something way more complicated, it's probably best just to use
meep directly. The meep writers are much better coders than I am. 


---
### Requirements (for now): 
  - meep (which requires guile, hdf5 and libctl.) 
  - ROOT  (tested with ROOT 6, but ROOT 5 might work) 

---
### Usage:

 Currently You can either use it a library (in e.g. a ROOT macro or your own
 compiled code). There should be examples under macros (or, there will be soon,
 I wrote this README before writing the macros...). Eventually, I'd like to
 have a input file  that would be read by an iceprop binary so that simulations
 could be run by people without C++ knowledge. Some day...  

---
### Units: 
  - Meep is dimensionless, so we are free to use our own units. 
  - Distance units are in meters. 
  - Cylindrical coordinates are used, where the transmitter is always at r=0. Negative z is below the ice surface, positive z is above. 
  - Time units are in ns, therefore frequency is in GHz. 
  - Electric field units are in Volts / m  (that implies current is in units of 1 / 377 A, I think). 
---
### Features: 
  - Cylindrically symmetric simulation of a transmitter antenna in the firn
  - Ability to track electric field vs. time at any point
  - Ability to track maximum and maximum time of electric field within the whole volume (but slow and not currently parallelized).
  - Output to ROOT tree of histograms (undoubtedly slower than creading HDF5 files and converting to ROOT later) 
  - Output to png / pdf through ROOT (undoubtedly slower than creating an HDF5 file and then converting to png/pdf later). 

---

### Notes/Caveats: 

  - meep supports MPI, including parallel output to hdf5 files. iceprop does
    nothing special for MPI, but I suppose if you compile meep with MPI and
    only output hdf5 files (and don't do maximum tracking), then it's plausible
    it'll mostly work in parallel.  Or it might crash. I don't know yet. If you
    use maximum tracking and none-HDF5 output (and it manages not to crash), I
    think you'll still get parallelized timesteps, but filling in outputs /
    doing maximum tracking might become a serious bottleneck. This can probably
    be addressed, but I'm not sure yet if meep exposes enough of its API to do
    the parallel whole-field evaluation without an intermediate HDF5 file. I
    guess I could probably trick it with a virtual HDF5 file and then use that
    for output in a separate process. Honestly, that sounds like a much better
    design and I'll probably switch to that at some point if this package
    survives more than a few weeks. 
  - Normal caveats about the perfectly matched boundary layers apply (they are
    slightly reflective due to numerical issues). This can be addressed by
    making it bigger.
  - Because it's cylindrically symmetric, if you put 

---

### ToDos / Planned features
  - Useful documentation 
  - Flux measurements
  - More firn models , including hoar frost layers, real density measurements, whatever. 
  - More source models. Would be great to incorporate a realistic neutrino model! 
  - Driver program with config file 
  - Python bindings? (maybe already works with pyROOT?) 
  - Parallelize better (see Notes/Caveats), perhaps via an intermediate hdf5 file. 
  - Noise simulation 
  - Geometries other than cylindrically symmetric 
  - Make output prettier.

---

### License: 
  GPLv3 (see COPYING).


  


 

