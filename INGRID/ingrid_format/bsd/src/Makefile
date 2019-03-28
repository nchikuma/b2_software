CXX=g++
AR=ar
LD_SHARED=g++

TOP_DIR=$(shell cd ..; pwd)

CXXFLAGS= -Wall -g -fPIC
SOFLAGS= -shared -Wl,-export-dynamic
LDFLAGS=
ARFLAGS= rcsv


ROOTLIBS   = -Wl,-rpath,$(shell root-config --libdir) $(shell root-config --libs) -lThread 
ROOTGLIBS  = -Wl,-rpath,$(shell root-config --libdir) $(shell root-config --glibs)
ROOTCFLAGS = $(shell root-config --cflags)
CXXFLAGS += $(ROOTCFLAGS) -I$(TOP_DIR)

LIBS = $(TOP_DIR)/lib/libBeamData.a
LIBS_SO = $(TOP_DIR)/lib/so/libBeamData.so


LIB_OBJECTS=BeamData.o

all: mklibdir $(LIBS) $(LIBS_SO)


%.o:%.cc
	$(CXX) $(CXXFLAGS) -c $<


mklibdir:
	mkdir -p $(TOP_DIR)/lib/so


$(LIBS): $(LIB_OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

$(LIBS_SO): $(LIBS) $(LIB_OBJECTS)
	$(LD_SHARED) $(SOFLAGS) -o $@ $(LIB_OBJECTS)


clean:
	rm -vf *.o  $(LIBS) $(LIBS_SO)
