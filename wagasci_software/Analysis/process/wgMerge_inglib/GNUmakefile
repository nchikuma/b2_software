TARGET = wgMerge_inglib

INCLUDE = -I ../../include
OBJDIR = ../../obj
SRCDIR = ../../src
BINDIR = ../../bin
LIBDIR = /gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/lib
INCLUDE += -I$(LIBDIR)

SRCS = $(shell ls $(SRCDIR)/*.cc)
OBJS = $(subst $(SRCDIR),$(OBJDIR),$(SRCS:.cc=.o))	

SRCS += $(TARGET).cc
OBJS += $(TARGET).o

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)


INGLIBS = -lm $(LIBDIR)/EVENTSUMMARY.so $(LIBDIR)/HitSummary.so \
          $(LIBDIR)/SimHitSummary.so $(LIBDIR)/SimVertexSummary.so \
          $(LIBDIR)/SimParticleSummary.so $(LIBDIR)/BeamInfoSummary.so \
          $(LIBDIR)/BasicReconSummary.so $(LIBDIR)/FirstReducSummary.so \
          $(LIBDIR)/TrackSummary.so $(LIBDIR)/NeutInfoSummary.so \
          $(LIBDIR)/TwoDimReconSummary.so $(LIBDIR)/ThreeDimReconSummary.so 

CXXFLAGS   = -g -MMD -MP -Wall -fPIC 
CXXFLAGS   += $(INCLUDE) $(ROOTCFLAGS)
CXXLIBS    = -lSpectrum
CXXLIBS    += $(INCLUDE) $(ROOTLIBS)
CXXLIBS    += $(INGLIBS)


DEPENDS = $(OBJS:.o=.d)
CC = g++

all: $(BINDIR)/$(TARGET)

$(BINDIR)/$(TARGET): $(OBJS)
	$(CC) $(CXXLIBS) -o $@ $^

$(OBJDIR)/%.o:$(SRCDIR)/%.cc
	@if [ ! -d $(OBJDIR) ]; \
		then echo "mkdir -p $(OBJDIR)"; mkdir -p $(OBJDIR); \
		fi
	$(CC) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

clean:
	$(RM) -f $(BINDIR)/$(TARGET) $(OBJS) $(DEPENDS)

-include $(DEPENDS)

.PHONY: all clean

