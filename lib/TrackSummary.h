#ifndef TRACKSUMMARY_H
#define TRACKSUMMARY_H
#include <iostream>
#include <vector>
#include "TObject.h"

#include "TRefArray.h"
#include "TRef.h"

#include "HitSummary.h"
#include "SimParticleSummary.h"

#define TRACK_MAXHITS 1000
#define TRACK_MAXSIMPARTICLES 10


//......................................................................

class TrackSummary : public TObject
{
public:
    TrackSummary();
    TrackSummary(const TrackSummary& trk);
    virtual ~TrackSummary();
    
    void Clear   (Option_t* option="");
    void Print();

public:
  
    float vtxi[3];
    float vtxf[3];
    float length;
    float ekin;
    float tx;
    float ty;
    float etx;
    float ety;
    float ex0;
    float ey0;
    float covx;
    float covy;
    float chi2x;
    float chi2y;
    float btheta;
    float bphi;
    int mrdhitid[2];
    float mucl;
    float vpe;
    int view;
    
    int NHits() const {return nhits;}
    int NSimParticles() const {return nsimparticles;}
    int NSimEMShowers() const {return nsimemshowers;}

    void AddHit(HitSummary* hit);
    HitSummary* GetHit(int i) const;

    void AddSimParticle(SimParticleSummary* part);
    SimParticleSummary* GetSimParticle(int i) const;



private:

    int nhits;
    int nsimparticles;
    int nsimemshowers;

    TRef fHit[TRACK_MAXHITS];
    TRef fSimParticle[TRACK_MAXSIMPARTICLES];



    ClassDef(TrackSummary, 4) // SciBar Track Summary
        };

#endif // SBSCIBARTRACKSUMMARY_H
////////////////////////////////////////////////////////////////////////

