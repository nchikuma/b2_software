ROOTINC = -I`root-config --incdir`
CLASSDIR = .
CLASSINC= -I${CLASSDIR}
INCDIRS = ${ROOTINC}${CLASSINC}

CXXFLAGS  = `root-config --cflags`
COPTFLAGS = -g

SOBJS = NeutInfoSummary.so EVENTSUMMARY.so HitSummary.so SimHitSummary.so SimVertexSummary.so SimParticleSummary.so BeamInfoSummary.so BasicReconSummary.so TwoDimReconSummary.so ThreeDimReconSummary.so FirstReducSummary.so TrackSummary.so FluxTuning.so

OBJS  =  ${SOBJS} 

.SUFFIXES:      .so

GENROOTSO = env COPTFLAGS="${COPTFLAGS}" INCDIRS="${INCDIRS}" \
                        ./gen_root_so.sh


.cc.o:
	$(CXX) -c $(COPTFLAGS) $(CXXFLAGS) $(INCDIRS) -o $@ $<

.cc.so:
	$(GENROOTSO) $(basename $<)

.F.o:
	$(FC) -c $(FOPTFLAGS) $(FINCDIRS) -o $@ $<

all: $(OBJS)

clean:
	$(RM) -f *.o *~ *LinkDef.h *Dict.* ${OBJS} chkfort
