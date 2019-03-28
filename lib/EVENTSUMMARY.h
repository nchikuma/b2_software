#ifndef EVENTSUMMARY_H
#define EVENTSUMMARY_H
#include <iostream>
#include "TObject.h"
#include "TClonesArray.h"
#include "TRef.h"
#include "TTimeStamp.h"

#include <stdlib.h>
#include <stdio.h>
#include <cstdio>
#include <vector>
#include <algorithm>
using namespace std;
#include "SimHitSummary.h"
#include "HitSummary.h"
#include "SimVertexSummary.h"
#include "SimParticleSummary.h"
#include "BeamInfoSummary.h"
#include "BasicReconSummary.h"
#include "TrackSummary.h"
#include "FirstReducSummary.h"
#include "NeutInfoSummary.h"
#include "TwoDimReconSummary.h"
#include "ThreeDimReconSummary.h"

//......................................................................

class EventSummary : public TObject
{
  public:
    EventSummary();
    EventSummary(const EventSummary& evt);
    virtual ~EventSummary();

    void NewTClonesArrays();
    void Clear     (Option_t* option="");

    void Print();

    SimHitSummary*        AddSimHit         (SimHitSummary* hit);
    SimHitSummary*        GetSimHit         (int i) const;
    HitSummary*           AddHit            (HitSummary* hit);
    HitSummary*           GetHit            (int i) const;
    HitSummary*           AddModHit         (HitSummary* hit, int nmod, int ncyc);
    HitSummary*           GetModHit         (int i, int nmod, int ncyc) const;

    SimVertexSummary*     AddSimVertex      (SimVertexSummary* hit);
    SimVertexSummary*     GetSimVertex      (int i) const;
    SimParticleSummary*   AddSimParticle    (SimParticleSummary* hit);
    SimParticleSummary*   GetSimParticle    (int i) const;

    BeamInfoSummary*      AddBeamSummary    (BeamInfoSummary*  beamsummary);
    BeamInfoSummary*      GetBeamSummary    (int i) const;

    BasicReconSummary*    AddBasicRecon     (BasicReconSummary* basicrecon);
    BasicReconSummary*    GetBasicRecon     (int i) const;
    TrackSummary*         AddTrack          (TrackSummary* trk);
    TrackSummary*         GetTrack          (int i) const;


    FirstReducSummary*    AddFirstReduc     (FirstReducSummary* reduc); 
    FirstReducSummary*    GetFirstReduc     (int i) const;

    NeutInfoSummary*      AddNeut           (NeutInfoSummary* neut); 
    NeutInfoSummary*      GetNeut           (int i) const;


    TwoDimReconSummary*   AddTwoDimRecon    (TwoDimReconSummary* twodrecon);
    TwoDimReconSummary*   GetTwoDimRecon    (int i) const;
    TwoDimReconSummary*   AddModTwoDimRecon (TwoDimReconSummary* hit, int nmod, int ncyc, int view);
    TwoDimReconSummary*   GetModTwoDimRecon (int i, int nmod, int ncyc, int view) const;
    ThreeDimReconSummary* AddThreeDimRecon  (ThreeDimReconSummary* threedrecon);
    ThreeDimReconSummary* GetThreeDimRecon  (int i) const;


  public:

    unsigned int      run;        // run number
    unsigned int    event;        // event number
    int           runmode;        // run mode
    int             trgid;        // trgger id(1:Beam, 2:Calib. , 128:Cosmic)
    unsigned char version;        // data structure version
    int              date;
    int              time;
    int           trgtime;
    int       nd280nspill;        // spill # at nd280 
                                  // = ( spill # at beam line ) 0xffff  + 1
    bool     bunch_flag[23];

    int NSimHits         ()    const {return nsimhits;}
    int NHits            ()    const {return nhits;}  
    int NSimVertexes     ()    const {return nsimvertexes;}
    int NSimParticles    ()    const {return nsimparticles;}
    int NNeutInfo        ()    const {return nneutinfos;}    
    int NBeamSummarys    ()    const {return nbeamsummarys;}
    int NBasicRecons     ()    const {return nbasicrecons;}
    int NTwoDimRecons        ()    const {return ntwodrecons;}  
    int NThreeDimRecons        ()    const {return nthreedrecons;}  
    int NTracks          ()    const {return ntracks;}  
    int NFirstReduc        ()    const {return nFirstreducs;}
    int NModHits (int nmod, int ncyc ) const {return nmodhits[nmod][ncyc];}  
    int NModTwoDimRecons   (int nmod, int ncyc , int view) const {return nmodtwodrecons[nmod][ncyc][view];}  
    //int NTracks  () const {return ntracks;}

  private:
    int nsimhits;         // number of SimHitSummaries  in this event
    int nhits;            // number of HitSummaries     in this event
    int nsimparticles;    // number of SimParticles     in this event
    int nsimvertexes;     // number of SimVertiexes     in this event
    int nbeamsummarys;    // number of BeamInfoSummarys       in this event  
    int nbasicrecons;     // number of BasicReconSummarys  in this event    
    int ntwodrecons;      // number of TwoDimReconSummarys  in this event    
    int nthreedrecons;    // number of ThreeDimReconSummarys  in this event    
    int ntracks;          // number of TrackSummarys  in this event    
    int nFirstreducs;     // number of FirstReducs  in this event    
    int nneutinfos;       // number of FirstReducs  in this event    

    int nmodhits[30][23];
    vector<int> nidmodhits[30*23];

    int nmodtwodrecons[30][23][2];
    vector<int> nidmodtwodrecons[30*23*2];


    TClonesArray* fSimHit; 
    TClonesArray* fHit;
    TClonesArray* fSimParticle;
    TClonesArray* fSimVertex;
    TClonesArray* fBeamSummary;
    TClonesArray* fBasicRecon;
    TClonesArray* fTwoDimRecon;
    TClonesArray* fThreeDimRecon;
    TClonesArray* fFirstReduc;
    TClonesArray* fTrack;
    TClonesArray* fNeutInfo;

    ClassDef(EventSummary, 5) // DST Summary of event

};

#endif // EVENTSUMMARY_H
////////////////////////////////////////////////////////////////////////

