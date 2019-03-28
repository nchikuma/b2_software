CPP 		= g++
CXXFLAGS	= -g -O3 -Wall -fPIC -D_REENTRANT -Wno-deprecated

ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTLIBS     	:= $(shell root-config --libs)
ROOTGLIBS    	:= $(shell root-config --glibs)
CXXFLAGS	+= $(ROOTCFLAGS)

WORKDIR=/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software
LIB  = $(WORKDIR)/lib
BSD  = $(WORKDIR)/INGRID/ingrid_format/bsd/lib/so
LIBS = $(ROOTGLIBS) \
          -lm $(LIB)/EVENTSUMMARY.so $(LIB)/HitSummary.so \
          $(LIB)/SimHitSummary.so $(LIB)/SimVertexSummary.so \
          $(LIB)/SimParticleSummary.so $(LIB)/BeamInfoSummary.so \
          $(LIB)/BasicReconSummary.so $(LIB)/FirstReducSummary.so \
          $(LIB)/TrackSummary.so $(LIB)/NeutInfoSummary.so \
          $(LIB)/TwoDimReconSummary.so $(LIB)/ThreeDimReconSummary.so \
          $(BSD)/libBeamData.so


SRC = ../src
INC = ../inc

CXXFLAGS += -I$(LIB) -I$(SRC) -I$(INC)

TARGET = IngCalib_new add_spilldb IngAddBSD

all: $(TARGET)

IngCalib_new: IngCalib_new.o
add_spilldb: add_spilldb.o
IngAddBSD: IngAddBSD.o



%: %.o
	@echo "Now make $@"
	@$(CPP) -o $@ $< $(CXXFLAGS) $(LIBS)
	@echo "..Compile done! "

%.o: %.cc
	@echo "Start Compiling $<"
	@$(CPP) $(CXXFLAGS) -c $<
	@echo ".. Compiling Object Files $<   --> done"
	@echo ""

clean: 
	@echo "Now Clean Up"
	rm -f $(TARGET) *~ *.o *.o~ core
