
CPPFLAGS += -lstdc++ -Df2cFortran
###
MYTEMP =  -L/opt/gcc-5.3.0/lib/gcc/x86_64-unknown-linux-gnu/5.3.0 -L/opt/gcc-5.3.0/lib64 -lgfortran
EXTRALIBS += $(MYTEMP)

CERN_ROOT = /home/t2k/tatsuya1/cern/cernlib/2006
CERNINC = $(CERN_ROOT)/include
CERNLIB = -L$(CERN_ROOT)/lib 
CERNLIBS = -lmathlib -lpacklib -lgraflib -lgrafX11 -lpacklib -lkernlib -lmathlib
CPPFLAGS += -I$(CERNINC)
EXTRALIBS += $(CERNLIBS)

DATALIBDIR = /gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/lib
CPPFLAGS += -I$(DATALIBDIR)

MYLIB = -lm $(DATALIBDIR)/EVENTSUMMARY.so $(DATALIBDIR)/HitSummary.so $(DATALIBDIR)/SimHitSummary.so $(DATALIBDIR)/SimVertexSummary.so $(DATALIBDIR)/SimParticleSummary.so $(DATALIBDIR)/BeamInfoSummary.so $(DATALIBDIR)/FirstReducSummary.so $(DATALIBDIR)/BasicReconSummary.so $(DATALIBDIR)/TrackSummary.so $(DATALIBDIR)/NeutInfoSummary.so $(DATALIBDIR)/TwoDimReconSummary.so $(DATALIBDIR)/ThreeDimReconSummary.so
EXTRALIBS += -L$(MYLIBDIR) $(MYLIB)

###
ROOTCINT=YES
ROOTLIBS =  $(shell $(ROOTSYS)/bin/root-config --glibs)
CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
EXTRALIBS += $(ROOTLIBS)

CPPFLAGS += -I$(G4INCLUDE)

###
name := b2mc
G4TARGET := $(name)
G4EXLIB := true

.PHONY: all
all: bin

include /gpfs/fs03/t2k/beam/work/rtam/software/GEANT/4.10.03.01/geant4.10.03.p01-install/share/Geant4-10.3.1/geant4make/config/binmake.gmk

LDFLAGS += $(CERNLIB) 

CXXFLAGS_WITHOUT_O := $(filter-out -O% , $(CXXFLAGS)) 
CXXFLAGS_WITHOUT_O := $(filter-out +O% , $(CXXFLAGS_WITHOUT_O))
