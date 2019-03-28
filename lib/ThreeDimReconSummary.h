#ifndef ThreeDimRECONSUMMARY_H
#define ThreeDimRECONSUMMARY_H
#include <iostream>
#include "TObject.h"
#include "TClonesArray.h"

#include "TRefArray.h"
#include "TRef.h"
#include "vector"
using namespace std;

#include "HitSummary.h"
#include "TrackSummary.h"


#define HIT_MAXHITS 1000  //## Temporary
#define RECON_MAXTRACKS 10 //## Temporary

//......................................................................
class ThreeDimReconSummary : public TObject{
 public:

  int       Ntrack;
  int       Ningtrack;
  float     clstime;      // time of cluster defined by most high p.e.
  float     clstimecorr;  // time after correction by measurement time of CT
  float     exptime;      // diff. from expected time
  int       hitcyc;       //
  bool      ontime;       //
  bool      vetowtracking; // Upstream VETO
  bool      edgewtracking; // Fiducial CUT
  float     vact[100];

  vector<Float_t>  x;
  vector<Float_t>  y;
  vector<Float_t>  z;
  vector<Float_t>  zx;
  vector<Float_t>  zy;

  vector<Int_t>    startmod;
  vector<Int_t>    stopmodx;
  vector<Int_t>    stopmody;

  vector<Int_t>    startxpln;
  vector<Int_t>    startypln;
  vector<Float_t>  startxch;
  vector<Float_t>  startych;
  vector<Int_t>    endxpln;
  vector<Int_t>    endypln;
  vector<Float_t>  endxch;
  vector<Float_t>  endych;
  vector<Float_t>  thetax;
  vector<Float_t>  thetay;
  vector<Float_t>  angle;
  vector<Int_t>    ing_startmod;
  vector<Int_t>    ing_endmod;
  vector<Int_t>    ing_startpln;
  vector<Int_t>    ing_endpln;
  vector<Bool_t>   ing_trk;
  vector<Bool_t>   pm_stop;
  vector<Bool_t>   ing_stop;
  vector<Float_t>  sci_range;
  vector<Float_t>  iron_range;
  vector<Int_t>    iron_pene;
  vector<Bool_t>   veto; // Upstream VETO
  vector<Bool_t>   edge; // Fiducial CUT
  vector<Int_t>    pdg;
  vector<Float_t>  mucl;
  vector<Float_t>  trkpe;
  vector<Int_t>    oneview;
  vector<Float_t>  diff_posx; //WM-ING
  vector<Float_t>  diff_posy;
  vector<Float_t>  diff_angx;
  vector<Float_t>  diff_angy;
  vector<Float_t>  diff_timex;
  vector<Float_t>  diff_timey;
  vector<Float_t>  diff_posx2; //PM-WM
  vector<Float_t>  diff_posy2;
  vector<Float_t>  diff_angx2;
  vector<Float_t>  diff_angy2;
  vector<Float_t>  diff_timex2;
  vector<Float_t>  diff_timey2;
  vector<Float_t>  diff_posx3; //PM-ING
  vector<Float_t>  diff_posy3;
  vector<Float_t>  diff_angx3;
  vector<Float_t>  diff_angy3;
  vector<Float_t>  diff_timex3;
  vector<Float_t>  diff_timey3;

 
  
  //###########################################
  //###########################################
  ThreeDimReconSummary();
  ThreeDimReconSummary(const ThreeDimReconSummary& basicsum);
  virtual ~ThreeDimReconSummary();
  void Clear   (Option_t* option="");
  void Print();
  void AddHit(HitSummary* hit);
  HitSummary* GetHit(int i) const;
  void AddHitTrk(HitSummary* hit, int track);
  HitSummary* GetHitTrk(int i, int track) const;
  void AddTrack(TrackSummary* trk);
  TrackSummary* GetTrack(int i) const;
  int nhits;
  int nhitTs[RECON_MAXTRACKS];
  
 private:
  
  TRef fHit[HIT_MAXHITS];
  TRef fTrack[RECON_MAXTRACKS];
  TRef fHitTrk[HIT_MAXHITS][RECON_MAXTRACKS];
 public:
  int ntracks;
  int Nhits(){return nhits;}
  int NhitTs(int trk) {return nhitTs[trk];}
  int Ntracks(){return ntracks;}
  ClassDef(ThreeDimReconSummary, 4)
    };

#endif
////////////////////////////////////////////////////////////////////////
