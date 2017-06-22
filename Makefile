#someday I'll probably convert this to cmake... but not today 

CXXFLAGS=`pkg-config --cflags meep` `root-config --cflags` -g -O2 
LDFLAGS=`pkg-config --libs meep` `root-config --libs` 

.PHONY: clean all 

LIBDIR=lib
BUILDDIR=build
INCLUDEDIR=include
BINDIR=bin

OBJS := $(addprefix $(BUILDDIR)/, Firn.o Sim.o Source.o icepropDict.o )

BINARIES := $(addprefix $(BINDIR)/, iceprop )
INCLUDES := $(addprefix $(INCLUDEDIR)/, $(shell ls $(INCLUDEDIR)))
LINKLIBNAME=firnprop 
LIBNAME = $(LIBDIR)/lib$(LINKLIBNAME).so 

all: $(LIBNAME) $(BINARIES) 

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
	$(CXX)  -I../include -I./ $(CXXFLAGS) -o $@ -c $< 

$(BINDIR)/%: src/%.cc $(INCLUDES) Makefile $(LIBNAME) | $(BINDIR)
	@echo Compiling $<
	@$(CXX)  -I./include -I./ $(CXXFLAGS) -o $@ -L./$(LIBDIR) $(LDFLAGS) -l$(LINKLIBNAME)  $<  $(LIBS) 

$(BUILDDIR)/icepropDict.cc: $(INCLUDES) LinkDef.h | $(BUILDDIR)
	@echo Running rootcint
	rootcint  -f $@ -c -p -I$(ANITA_UTIL_INSTALL_DIR)/include $(INCLUDES) LinkDef.h

clean: 
	rm -rf build
	rm -rf bin
	rm -rf lib
	










