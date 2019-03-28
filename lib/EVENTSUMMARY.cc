#include "./EVENTSUMMARY.h"

#include <vector>
#include <algorithm>
using namespace std;

//......................................................................

EventSummary::EventSummary():
  fSimHit(0),fHit(0), 
  fSimParticle(0), fSimVertex(0),
  fBeamSummary(0), fBasicRecon(0), fTwoDimRecon(0), fThreeDimRecon(0),
  fTrack(0),fFirstReduc(0),
  fNeutInfo(0)
{ 
  this->NewTClonesArrays();
  this->Clear("C");
}

//......................................................................

EventSummary::EventSummary(const EventSummary& evt) 
{ 

  this->NewTClonesArrays();
  this->Clear("C");

  run      = evt.run;
  event    = evt.event;
  runmode  = evt.runmode;
  trgid    = evt.trgid;
  version  = evt.version;
  date     = evt.date;
  time     = evt.time;
  trgtime  = evt.trgtime;

  for(int i=0; i<23;i++)bunch_flag[i] = evt.bunch_flag[i];

  nsimhits = 0;
  nhits    = 0;
  nbasicrecons   = 0;
  ntwodrecons      = 0;
  nthreedrecons      = 0;
  nbeamsummarys  = 0;
  nsimvertexes   = 0;
  nsimparticles  = 0;
  nFirstreducs     = 0;
  nneutinfos     = 0;

  for(int mod=0; mod<30; mod++){
    for(int cyc=0; cyc<23; cyc++){
      nidmodhits    [mod * 23 + cyc].clear();
      nmodhits[mod][cyc] = 0;
    }
  }

  for(int mod=0; mod<30; mod++){
    for(int cyc=0; cyc<23; cyc++){
      for(int view=0; view<2; view++){
        nidmodtwodrecons    [(mod * 23 + cyc )*2 + view].clear();
        nmodtwodrecons[mod][cyc][view] = 0;
      }
    }
  }

  //ntracks  = 0;

  // Copy all arrays of detector summaries
  int i;
  for (i=0; i < evt.nsimhits;    ++i) {
    AddSimHit(evt.GetSimHit(i));
  }
  for (i=0; i < evt.nhits;       ++i) {
    HitSummary*  hitsum = evt.GetHit(i);
    int mod = hitsum -> mod;
    int cyc = hitsum -> cyc;
    AddModHit   ( evt.GetHit(i), mod, cyc);
  }
  for (i=0; i < evt.nbeamsummarys;    ++i) {
    AddBeamSummary(evt.GetBeamSummary(i));
  }
  for (i=0; i < evt.nbasicrecons;    ++i) {
    AddBasicRecon(evt.GetBasicRecon(i));
  }
  for (i=0; i < evt.ntwodrecons;    ++i) {
    TwoDimReconSummary*  twodreconsum = evt.GetTwoDimRecon(i);
    int  mod  = twodreconsum -> hitmod;
    int  cyc  = twodreconsum -> hitcyc;
    int  view = twodreconsum -> view;
    AddModTwoDimRecon   ( evt.GetTwoDimRecon(i), mod, cyc, view);
  }
  for (i=0; i < evt.nthreedrecons;    ++i) {
    AddThreeDimRecon(evt.GetThreeDimRecon(i));
  }
  for (i=0; i < evt.ntracks;    ++i) {
    AddTrack(evt.GetTrack(i));
  }
  for (i=0; i < evt.nFirstreducs;    ++i) {
    AddFirstReduc(evt.GetFirstReduc(i));
  }
  for (i=0; i < evt.nneutinfos; ++i) {
    AddNeut(evt.GetNeut(i));
  }

}

//......................................................................

EventSummary::~EventSummary() 
{ 
  this->Clear("C");

  if (fSimHit) {
    delete fSimHit;
    fSimHit = 0;
  }

  if (fHit) {
    delete fHit;
    fHit = 0;
  }

  if (fSimParticle) {
    delete fSimParticle;
    fSimParticle = 0;
  }

  if (fSimVertex) {
    delete fSimVertex;
    fSimVertex = 0;
  }

  if (fBeamSummary) {
    delete fBeamSummary;
    fBeamSummary = 0;
  }
  if (fBasicRecon) {
    delete fBasicRecon;
    fBasicRecon = 0;
  }
  if (fTwoDimRecon) {
    delete fTwoDimRecon;
    fTwoDimRecon = 0;
  }
  if (fThreeDimRecon) {
    delete fThreeDimRecon;
    fThreeDimRecon = 0;
  }
  if (fTrack) {
    delete fTrack;
    fTrack = 0;
  }
  if (fFirstReduc) {
    delete fFirstReduc;
    fFirstReduc = 0;
  }
  if (fNeutInfo) {
    delete fNeutInfo;
    nneutinfos = 0;
  }

}

//......................................................................

void EventSummary::NewTClonesArrays()  
{ 

  fSimHit        = new TClonesArray("SimHitSummary"       , 1000);
  fHit           = new TClonesArray("HitSummary"          , 1000);
  fSimVertex     = new TClonesArray("SimVertexSummary"    , 1000);
  fSimParticle   = new TClonesArray("SimParticleSummary"  , 1000);
  fBeamSummary   = new TClonesArray("BeamInfoSummary"     , 1000);
  fBasicRecon    = new TClonesArray("BasicReconSummary"   , 1000);
  fTwoDimRecon   = new TClonesArray("TwoDimReconSummary"  , 1000);
  fThreeDimRecon = new TClonesArray("ThreeDimReconSummary", 1000);
  fTrack         = new TClonesArray("TrackSummary"        , 1000);
  fFirstReduc    = new TClonesArray("FirstReducSummary"   , 1000);
  fNeutInfo      = new TClonesArray("NeutInfoSummary"     , 1000);

}

//......................................................................

void EventSummary::Clear(Option_t* option) 
{ 

  if (fSimHit)       fSimHit       -> Clear(option);
  if (fHit)          fHit          -> Clear(option);
  if (fSimParticle)  fSimParticle  -> Clear(option);
  if (fSimVertex)    fSimVertex    -> Clear(option);
  if (fBeamSummary)  fBeamSummary  -> Clear(option);
  if (fBasicRecon)   fBasicRecon   -> Clear(option);
  if (fTwoDimRecon)      fTwoDimRecon      -> Clear(option);
  if (fThreeDimRecon)      fThreeDimRecon      -> Clear(option);    
  if (fTrack)        fTrack        -> Clear(option);
  if (fFirstReduc)     fFirstReduc     -> Clear(option);
  if (fNeutInfo)     fNeutInfo     -> Clear(option);

  nsimhits        = 0;
  nhits           = 0;
  nsimparticles   = 0;
  nsimvertexes    = 0;
  nbeamsummarys   = 0;
  nbasicrecons    = 0;
  ntwodrecons       = 0;
  nthreedrecons       = 0;
  ntracks         = 0;
  nFirstreducs      = 0;
  nneutinfos      = 0;
  for(int mod=0; mod<30; mod++){
    for(int cyc=0; cyc<23; cyc++){
      nidmodhits[mod*23+cyc].clear();
      nmodhits[mod][cyc] = 0;
    }
  }
  for(int mod=0; mod<30; mod++){
    for(int cyc=0; cyc<23; cyc++){
      for(int view=0; view<2; view++){
        nidmodtwodrecons    [(mod * 23 + cyc )*2 + view].clear();
        nmodtwodrecons[mod][cyc][view] = 0;
      }
    }
  }


  // initialize rest of data members
  run     = 0;
  event   = 0;
  runmode = -10;
  trgid   = -10;
  version = ' ';
  date = -1;
  time = -1;
  trgtime = -1;
  for(int i=0;i<23;i++)bunch_flag[i]=false;
}

//......................................................................

void EventSummary::Print()
{
  std::cout << "Event summary: " << std::endl;

}

//......................................................................


SimHitSummary* EventSummary::AddSimHit(SimHitSummary* hitsum) 
{
  TClonesArray &simhit_s = *fSimHit;
  new(simhit_s[nsimhits++]) SimHitSummary(*hitsum);
  return (SimHitSummary*)(fSimHit->At(nsimhits-1));
}


//......................................................................

SimHitSummary* EventSummary::GetSimHit(int i) const
{ 
  if (i < nsimhits && i>=0 ) return (SimHitSummary*)(fSimHit->At(i));
  return 0;
}


//......................................................................

HitSummary* EventSummary::AddHit(HitSummary* hitsum) 
{
  TClonesArray &hit_s = *fHit;
  new(hit_s[nhits++]) HitSummary(*hitsum);
  return (HitSummary*)(fHit->At(nhits-1));
}



//......................................................................

HitSummary* EventSummary::GetHit(int i) const
{ 
  if (i < nhits && i>=0 ) return (HitSummary*)(fHit->At(i));
  return 0;
}



HitSummary* EventSummary::AddModHit(HitSummary* hitsum, int nmod, int ncyc) 
{
  TClonesArray &hit_s = *fHit;
  nidmodhits    [ nmod * 23 + ncyc ].push_back( nhits );
  nmodhits[ nmod ][ ncyc ]++;
  new(hit_s[nhits++]) HitSummary(*hitsum);
  return (HitSummary*)(fHit->At(nhits-1));
}


HitSummary* EventSummary::GetModHit(int i, int nmod, int ncyc ) const
{ 
  if (i < nmodhits[nmod][ncyc] && i>=0 ) return (HitSummary*)(fHit->At( nidmodhits[nmod * 23 + ncyc][i] ));
  return 0;
}

SimVertexSummary* EventSummary::AddSimVertex(SimVertexSummary* simvertex) 
{
  TClonesArray &simvertex_s = *fSimVertex;
  new(simvertex_s[nsimvertexes++]) SimVertexSummary(*simvertex);
  return (SimVertexSummary*)(fSimVertex->At(nsimvertexes-1));
}


//......................................................................

SimVertexSummary* EventSummary::GetSimVertex(int i) const
{ 
  if (i < nsimvertexes && i>=0 ) return (SimVertexSummary*)(fSimVertex->At(i));
  return 0;
}


SimParticleSummary* EventSummary::AddSimParticle(SimParticleSummary* simparticle) 
{
  TClonesArray &simparticle_s = *fSimParticle;
  new(simparticle_s[nsimparticles++]) SimParticleSummary(*simparticle);
  return (SimParticleSummary*)(fSimParticle->At(nsimparticles-1));
}


//......................................................................

SimParticleSummary* EventSummary::GetSimParticle(int i) const
{ 
  if (i < nsimparticles && i>=0 ) return (SimParticleSummary*)(fSimParticle->At(i));
  return 0;
}
//......................................................................
BeamInfoSummary* EventSummary::AddBeamSummary(BeamInfoSummary* beamsummary) 
{
  TClonesArray &beamsummary_s = *fBeamSummary;
  new(beamsummary_s[nbeamsummarys++]) BeamInfoSummary(*beamsummary);
  return (BeamInfoSummary*)(fBeamSummary->At(nbeamsummarys-1));
}

//......................................................................

BeamInfoSummary* EventSummary::GetBeamSummary(int i) const
{ 
  if (i < nbeamsummarys && i>=0 ) return (BeamInfoSummary*)(fBeamSummary->At(i));
  return 0;
}

//......................................................................
BasicReconSummary* EventSummary::AddBasicRecon(BasicReconSummary* basicrecon) 
{
  TClonesArray &basicrecon_s = *fBasicRecon;
  new(basicrecon_s[nbasicrecons++]) BasicReconSummary(*basicrecon);
  return (BasicReconSummary*)(fBasicRecon->At(nbasicrecons-1));
}


//......................................................................

BasicReconSummary* EventSummary::GetBasicRecon(int i) const
{ 
  if (i < nbasicrecons && i>=0 ) return (BasicReconSummary*)(fBasicRecon->At(i));
  return 0;
}




//......................................................................
TwoDimReconSummary* EventSummary::AddTwoDimRecon(TwoDimReconSummary* twodrecon) 
{

  TClonesArray &twodrecon_s = *fTwoDimRecon;
  new(twodrecon_s[ntwodrecons++]) TwoDimReconSummary(*twodrecon);
  return (TwoDimReconSummary*)(fTwoDimRecon->At(ntwodrecons-1));
}


//......................................................................

TwoDimReconSummary* EventSummary::GetTwoDimRecon(int i) const
{ 
  if (i < ntwodrecons && i>=0 ) return (TwoDimReconSummary*)(fTwoDimRecon->At(i));
  return 0;
}




TwoDimReconSummary* EventSummary::AddModTwoDimRecon(TwoDimReconSummary* twodreconsum, int nmod, int ncyc, int nview) 
{
  TClonesArray &twodrecon_s = *fTwoDimRecon;
  //##### push back nhits to hit # ID vector ######
  nidmodtwodrecons    [ (nmod * 23 + ncyc )*2 + nview ].push_back( ntwodrecons );
  nmodtwodrecons[ nmod ][ ncyc ][ nview ]++;
  new(twodrecon_s[ntwodrecons++]) TwoDimReconSummary(*twodreconsum);
  return (TwoDimReconSummary*)(fTwoDimRecon->At(ntwodrecons-1));
}


TwoDimReconSummary* EventSummary::GetModTwoDimRecon(int i, int nmod, int ncyc, int nview ) const
{ 
  if (i < nmodtwodrecons[nmod][ncyc][nview] && i>=0 ) return (TwoDimReconSummary*)(fTwoDimRecon->At( nidmodtwodrecons[(nmod * 23 + ncyc)*2 + nview][i] ));
  return 0;
}




ThreeDimReconSummary* EventSummary::AddThreeDimRecon(ThreeDimReconSummary* threedrecon) 
{
  TClonesArray &threedrecon_s = *fThreeDimRecon;
  new(threedrecon_s[nthreedrecons++]) ThreeDimReconSummary(*threedrecon);
  return (ThreeDimReconSummary*)(fThreeDimRecon->At(nthreedrecons-1));
}

//......................................................................

ThreeDimReconSummary* EventSummary::GetThreeDimRecon(int i) const
{ 
  if (i < nthreedrecons && i>=0 ) return (ThreeDimReconSummary*)(fThreeDimRecon->At(i));
  return 0;
}



//......................................................................
FirstReducSummary* EventSummary::AddFirstReduc(FirstReducSummary* reduc) 
{
  TClonesArray &reduc_s = *fFirstReduc;
  new(reduc_s[nFirstreducs++]) FirstReducSummary(*reduc);
  return (FirstReducSummary*)(fFirstReduc->At(nFirstreducs-1));
}


//......................................................................

FirstReducSummary* EventSummary::GetFirstReduc(int i) const
{ 
  if (i < nFirstreducs && i>=0 ) return (FirstReducSummary*)(fFirstReduc->At(i));
  return 0;
}


//......................................................................
NeutInfoSummary* EventSummary::AddNeut(NeutInfoSummary* neut) 
{
  TClonesArray &neutinfo_s = *fNeutInfo;
  new(neutinfo_s[nneutinfos++]) NeutInfoSummary(*neut);
  return (NeutInfoSummary*)(fNeutInfo->At(nneutinfos-1));
}


//......................................................................

NeutInfoSummary* EventSummary::GetNeut(int i) const
{ 
  if (i < nneutinfos && i>=0 ) return (NeutInfoSummary*)(fNeutInfo->At(i));
  return 0;
}



TrackSummary* EventSummary::AddTrack(TrackSummary* trk)
{
  TClonesArray &track_s = *fTrack;
  new(track_s[ntracks++]) TrackSummary(*trk);
  return (TrackSummary*)(fTrack->At(ntracks-1));
}


//......................................................................

TrackSummary* EventSummary::GetTrack(int i) const
{ 
  if (i < ntracks && i>=0 ) return (TrackSummary*)(fTrack->At(i));
  return 0;
}



ClassImp(EventSummary)


  ////////////////////////////////////////////////////////////////////////
