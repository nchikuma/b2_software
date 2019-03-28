#include "ThreeDimReconSummary.h"

//......................................................................

ThreeDimReconSummary::ThreeDimReconSummary():
  nhits(0), ntracks(0)
{ 
  Ntrack      = -1; 
  Ningtrack   = -1; 
  clstime     = -1.e5;  
  clstimecorr = -1.e5;
  exptime     = -1.e5;
  hitcyc      = -1;
  ontime      =  false;
  vetowtracking=  false;
  edgewtracking=  false;
  for(int i=0;i<100;i++) vact[i]=-1;

  startmod.clear();
  stopmodx.clear();
  stopmody.clear();

  x.clear();
  y.clear();
  z.clear();
  zx.clear();
  zy.clear();
  startxpln.clear();
  startypln.clear();
  startxch.clear();
  startych.clear();
  endxpln.clear();
  endypln.clear();
  endxch.clear();
  endych.clear();
  thetax.clear();
  thetay.clear();
  angle.clear();
  ing_startmod.clear();
  ing_endmod.clear();
  ing_startpln.clear();
  ing_endpln.clear();
  ing_trk.clear();
  pm_stop.clear();
  ing_stop.clear();
  sci_range.clear();
  iron_range.clear();
  iron_pene.clear();
  veto.clear();
  edge.clear();
  pdg.clear();
  mucl.clear();
  trkpe.clear();
  oneview.clear();
  diff_posx.clear();
  diff_posy.clear();
  diff_angx.clear();
  diff_angy.clear();
  diff_timex.clear();
  diff_timey.clear();
  diff_posx2.clear();
  diff_posy2.clear();
  diff_angx2.clear();
  diff_angy2.clear();
  diff_timex2.clear();
  diff_timey2.clear();
  diff_posx3.clear();
  diff_posy3.clear();
  diff_angx3.clear();
  diff_angy3.clear();
  diff_timex3.clear();
  diff_timey3.clear();

  for(int i=0;i<RECON_MAXTRACKS;i++) nhitTs[i]=0;

}

//......................................................................
ThreeDimReconSummary::ThreeDimReconSummary(const ThreeDimReconSummary& basicsum) :
  nhits(0), ntracks(0)
{
  Ntrack      = basicsum. Ntrack ;
  Ningtrack   = basicsum. Ningtrack ;
  clstime     = basicsum. clstime ;
  clstimecorr = basicsum. clstimecorr;
  exptime     = basicsum. exptime;
  hitcyc      = basicsum. hitcyc;
  ontime      = basicsum. ontime;
  vetowtracking  = basicsum. vetowtracking ;
  edgewtracking  = basicsum. edgewtracking ;
  for(int i=0;i<100;i++) vact[i]=basicsum.vact[i];

  startmod     = basicsum. startmod;
  stopmodx     = basicsum. stopmodx;
  stopmody     = basicsum. stopmody;

  x            = basicsum. x;
  y            = basicsum. y;
  z            = basicsum. z;
  zx           = basicsum. zx;
  zy           = basicsum. zy;
  startxpln    = basicsum. startxpln;
  startypln    = basicsum. startypln;
  startxch     = basicsum. startxch;
  startych     = basicsum. startych;
  endxpln      = basicsum. endxpln;
  endypln      = basicsum. endypln;
  endxch       = basicsum. endxch;
  endych       = basicsum. endych;
  thetax       = basicsum. thetax;
  thetay       = basicsum. thetay;
  angle        = basicsum. angle;
  ing_startmod = basicsum. ing_startmod;
  ing_endmod   = basicsum. ing_endmod;
  ing_startpln = basicsum. ing_startpln;
  ing_endpln   = basicsum. ing_endpln;
  ing_trk      = basicsum. ing_trk;
  pm_stop      = basicsum. pm_stop;
  ing_stop     = basicsum. ing_stop;
  sci_range    = basicsum. sci_range;
  iron_range   = basicsum. iron_range;
  iron_pene    = basicsum. iron_pene;
  veto         = basicsum. veto;
  edge         = basicsum. edge;
  pdg          = basicsum. pdg;
  mucl         = basicsum. mucl;
  trkpe        = basicsum. trkpe;
  oneview      = basicsum. oneview;
  diff_posx    = basicsum. diff_posx;
  diff_posy    = basicsum. diff_posy;
  diff_angx    = basicsum. diff_angx;
  diff_angy    = basicsum. diff_angy;
  diff_timex   = basicsum. diff_timex;
  diff_timey   = basicsum. diff_timey;
  diff_posx2   = basicsum. diff_posx2;
  diff_posy2   = basicsum. diff_posy2;
  diff_angx2   = basicsum. diff_angx2;
  diff_angy2   = basicsum. diff_angy2;
  diff_timex2  = basicsum. diff_timex2;
  diff_timey2  = basicsum. diff_timey2;
  diff_posx3   = basicsum. diff_posx3;
  diff_posy3   = basicsum. diff_posy3;
  diff_angx3   = basicsum. diff_angx3;
  diff_angy3   = basicsum. diff_angy3;
  diff_timex3  = basicsum. diff_timex3;
  diff_timey3  = basicsum. diff_timey3;

  ntracks = 0;
  nhits = 0; 
  for(int i=0;i<RECON_MAXTRACKS;i++) nhitTs[i] = 0; 

  for (int i=0; i<HIT_MAXHITS; ++i) {
    fHit[i] = TRef(NULL);
  }
  for (int i=0; i < basicsum.nhits; ++i) 
    AddHit(basicsum.GetHit(i));

  for (int i=0; i<RECON_MAXTRACKS; ++i) {
    fTrack[i] = TRef(NULL);
  }
  for (int i=0; i < basicsum.ntracks; ++i){
    AddTrack(basicsum.GetTrack(i)); 
  }
  for (int j=0; j< HIT_MAXHITS; ++j) {
    for (int i=0; i < RECON_MAXTRACKS; ++i) {
      fHitTrk[i][j] = TRef(NULL);
    }
  }
  for (int j=0; j< basicsum.Ntrack; ++j) {
    for (int i=0; i < basicsum.nhitTs[j]; ++i) {
      AddHitTrk(basicsum.GetHitTrk(i,j),j);
    }
  }

}

//......................................................................



ThreeDimReconSummary::~ThreeDimReconSummary() 
{

}

//......................................................................


void ThreeDimReconSummary::Clear(Option_t* option)
{
  for (int i=0; i<HIT_MAXHITS; ++i)
    fHit[i] = TRef(NULL);
  nhits = 0;
  
  for (int i=0; i<RECON_MAXTRACKS; ++i)
    fTrack[i] = TRef(NULL);
  ntracks = 0;

  for(int i=0; i<RECON_MAXTRACKS;i++){
     for(int j=0;j<HIT_MAXHITS;j++){
	     fHitTrk[j][i] = TRef(NULL);
     }
     nhitTs[i]=0;
  }
}

//......................................................................

void ThreeDimReconSummary::Print()
{

}
//......................................................................
void ThreeDimReconSummary::AddHit(HitSummary* hit) 
{
  if (nhits < HIT_MAXHITS) {
    fHit[nhits] = TRef((HitSummary*) hit);
    ++nhits;
  }
}
//......................................................................


HitSummary* ThreeDimReconSummary::GetHit(int i) const
{ 
  return (HitSummary*)fHit[i].GetObject();
}

//......................................................................
void ThreeDimReconSummary::AddHitTrk(HitSummary* hit, int trk)
{

  if (trk < RECON_MAXTRACKS) {
    if (nhitTs[trk] < HIT_MAXHITS) {
      fHitTrk[nhitTs[trk]][trk] = TRef((HitSummary*) hit);
      fHit[nhits] = TRef((HitSummary*) hit);
      ++nhitTs[trk];
      ++nhits;
    }
  }  

}

//......................................................................

HitSummary* ThreeDimReconSummary::GetHitTrk(int i, int j) const
{
  return (HitSummary*)fHitTrk[i][j].GetObject();
}

//......................................................................

void ThreeDimReconSummary::AddTrack(TrackSummary* trk) 
{
  if (ntracks < RECON_MAXTRACKS) {
    fTrack[ntracks] = TRef((TrackSummary*) trk);
    ++ntracks;
  }
}
//......................................................................


TrackSummary* ThreeDimReconSummary::GetTrack(int i) const
{ 

  return (TrackSummary*)fTrack[i].GetObject();
}

//......................................................................


ClassImp(ThreeDimReconSummary)

////////////////////////////////////////////////////////////////////////
