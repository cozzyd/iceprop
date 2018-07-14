### iceprop Makefile ### 


### Configurable options: 
## turn on to enable MPI support (But will be disabled anyway if meep wasn't built with MPI support) 
USE_MPI=yes



### End configurable options 

CXXFLAGS=`pkg-config --cflags meep` `root-config --cflags` -g -Os  -fPIC -march=native
LDFLAGS=`pkg-config --libs meep` `root-config --libs` -lRootFftwWrapper

MEEP_MPI=`meep --version | grep "Using MPI"` 

ifeq ($(MEEP_MPI),) 
	USE_MPI=no
endif 

MPI_MESSAGE="MPI is disabled" 
ifeq ($(USE_MPI),yes)
	CXXFLAGS+= -DENABLE_MPI  
	CXX=mpiCC  
	MPI_MESSAGE="MPI is enabled" 
endif

.PHONY: clean all 

LIBDIR=lib
BUILDDIR=build
INCLUDEDIR=include
BINDIR=bin

OBJS := $(addprefix $(BUILDDIR)/, Firn.o Sim.o Source.o MPI.o icepropDict.o )

BINARIES = $(addprefix $(BINDIR)/, iceprop )
INCLUDES = $(addprefix $(INCLUDEDIR)/, iceprop.h iceprop/Sim.h iceprop/Firn.h iceprop/Source.h iceprop/Units.h iceprop/MPI.h)

LINKLIBNAME=iceprop

LIBNAME = $(LIBDIR)/lib$(LINKLIBNAME).so 

all: $(LIBNAME) $(BINARIES) 
	@echo $(MPI_MESSAGE)

$(LIBNAME): $(OBJS) | $(LIBDIR)
	@echo Building shared library $@
	@$(CXX) $(LDFLAGS) $(OBJS) -shared -o $@ $(LIBS)
	cp $(BUILDDIR)/*.pcm $(LIBDIR) 

$(BUILDDIR): 
	mkdir -p $(BUILDDIR)

$(BINDIR): 
	mkdir -p $(BINDIR)

$(LIBDIR): 
	mkdir -p $(LIBDIR)


$(BUILDDIR)/%.o: src/%.cc $(INCLUDES) Makefile | $(BUILDDIR)
	@echo Compiling  $< 
	@$(CXX)  -I./include $(CXXFLAGS) -o $@ -c $< 

$(BUILDDIR)/%.o: build/%.cc $(INCLUDES) Makefile | $(BUILDDIR) 
	@echo Compiling  $< 
	@$(CXX)  -I../include -I./ $(CXXFLAGS) -o $@ -c $< 

$(BINDIR)/%: src/%.cc $(INCLUDES) Makefile $(LIBNAME) | $(BINDIR)
	@echo Compiling $<
	@$(CXX)  -I./include -I./ $(CXXFLAGS) -o $@ -L./$(LIBDIR) $(LDFLAGS) -l$(LINKLIBNAME)  $<  $(LIBS) 

$(BUILDDIR)/icepropDict.cc: $(INCLUDES) LinkDef.h | $(BUILDDIR)
	@echo Running rootcint
	@rootcint  -f $@ -c -p -I$(ANITA_UTIL_INSTALL_DIR)/include $(INCLUDES) LinkDef.h

clean: 
	rm -rf build
	rm -rf bin
	rm -rf lib
	









