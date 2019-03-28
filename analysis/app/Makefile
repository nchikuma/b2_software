CPP 		= g++
CXXFLAGS	= -g -O3 -Wall -fPIC -D_REENTRANT -Wno-deprecated

ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTLIBS     	:= $(shell root-config --libs)
ROOTGLIBS    	:= $(shell root-config --glibs)
CXXFLAGS	+= $(ROOTCFLAGS)

LIB = /gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/lib

LIBS 	= $(ROOTGLIBS) \
          -lm $(LIB)/EVENTSUMMARY.so $(LIB)/HitSummary.so \
          $(LIB)/SimHitSummary.so $(LIB)/SimVertexSummary.so \
          $(LIB)/SimParticleSummary.so $(LIB)/BeamInfoSummary.so \
          $(LIB)/BasicReconSummary.so $(LIB)/FirstReducSummary.so \
          $(LIB)/TrackSummary.so $(LIB)/NeutInfoSummary.so \
          $(LIB)/TwoDimReconSummary.so $(LIB)/ThreeDimReconSummary.so \
	  $(LIB)/FluxTuning.so

SRC = ../src
INC = ../inc
OBJ = ../obj

CXXFLAGS += -I$(LIB) -I$(SRC) -I$(INC) -I$(LIB)

TARGET = B2Disp_cyc B2Plot B2Plot_BG TwoDimRecon ThreeDimRecon B2gain HitEfficiency HitEfficiency_2 TrackEfficiency TrackEfficiency2 TrackEfficiency_sim B2pe B2pe_Scinti DetEff HitEff2 Selection FluxWeight

all: $(TARGET)

HitEfficiency: HitEfficiency.o DetectorDimension.o INGRID_BadCh_mapping.o

HitEfficiency_2: HitEfficiency_2.o DetectorDimension.o INGRID_BadCh_mapping.o

TrackEfficiency: TrackEfficiency.o DetectorDimension.o INGRID_BadCh_mapping.o
TrackEfficiency2: TrackEfficiency2.o DetectorDimension.o INGRID_BadCh_mapping.o
HitEff2: HitEff2.o DetectorDimension.o INGRID_BadCh_mapping.o

TrackEfficiency_sim: TrackEfficiency_sim.o DetectorDimension.o INGRID_BadCh_mapping.o

B2Disp_cyc: B2Disp_cyc.o DetectorDimension.o INGRID_BadCh_mapping.o

B2Plot: B2Plot.o DetectorDimension.o INGRID_BadCh_mapping.o
Selection: Selection.o DetectorDimension.o INGRID_BadCh_mapping.o

DetEff: DetEff.o DetectorDimension.o INGRID_BadCh_mapping.o


B2Plot_BG: B2Plot_BG.o

TwoDimRecon: TwoDimRecon.o DetectorDimension.o INGRID_BadCh_mapping.o

ThreeDimRecon: ThreeDimRecon.o DetectorDimension.o INGRID_BadCh_mapping.o

B2gain: B2gain.o INGRID_BadCh_mapping.o

FluxWeight: FluxWeight.o

B2pe: B2pe.o DetectorDimension.o INGRID_BadCh_mapping.o

B2pe_Scinti: B2pe_Scinti.o DetectorDimension.o INGRID_BadCh_mapping.o

% : %.o
	@echo " ===================================================== "
	@echo "    Now make $@"
	@echo " ===================================================== "
	@echo " $(CPP) -o $@ $^ $(CXXFLAGS) $(LIBS)"
	@$(CPP) -o $@ $^ $(CXXFLAGS) $(LIBS)
	@echo " ===================================================== "
	@echo "    Compile ...done! $@"
	@echo " ===================================================== "
	@echo ""

%.o: %.cc
	@echo " ----------------------------------------------------- "
	@echo "    Start Compiling Object File: $<"
	@echo " ----------------------------------------------------- "
	@$(CPP) $(CXXFLAGS) -o $@ -c $<
	@echo " ----------------------------------------------------- "
	@echo ".. Compiling Object Files $<   --> done"
	@echo " ----------------------------------------------------- "
	@echo ""

%.o: $(SRC)/%.cc
	@echo " ----------------------------------------------------- "
	@echo "    Start Compiling Object File: $<"
	@echo " ----------------------------------------------------- "
	@$(CPP) $(CXXFLAGS) -o $@ -c $<
	@echo " ----------------------------------------------------- "
	@echo ".. Compiling Object Files $<   --> done"
	@echo " ----------------------------------------------------- "
	@echo ""

clean: 
	@echo "Now Clean Up"
	rm -f $(TARGET) *~ *.o *.o~ core
