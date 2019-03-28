#ifndef __THREEDIMCLASS_HXX__
#define __THREEDIMCLASS_HXX__

#include "TwoDimRecon.hxx"

class Hits;
class TrackIng;
class TrackWM;
class TrackPM;
class Trk;
class AnaTrack;

vector<AnaTrack> analyzed_trk;
vector<TrackWM>  hwmtrack;
vector<TrackWM>  vwmtrack;
vector<TrackPM>  hpmtrack;
vector<TrackPM>  vpmtrack;
vector<TrackIng> hingtrack;
vector<TrackIng> vingtrack;

float nonrechits     [Cmod][Cview][Cpln][Cch];
float nonrechits_lope[Cmod][Cview][Cpln][Cch];
float nonrechits_pdg [Cmod][Cview][Cpln][Cch];
int   nonrechits_id  [Cmod][Cview][Cpln][Cch];

Int_t INGMODNUM_start = -1;
Int_t INGMODNUM_end   = -1;
Int_t INGMODNUM_mid   = -1;
Int_t WMMODNUM        = -1;
Int_t PMMODNUM        = -1;

Int_t MODE_DET = -1;

class Hits{
public:

  Int_t   mod;
  Int_t   view;
  Int_t   pln;
  Int_t   ch;
  Int_t   pdg;
  Int_t   cyc;
  Float_t pe;
  Float_t lope;
  Bool_t  isohit;
  Int_t   recon_id;
  Int_t   hit_id;

  void clear(){

    mod      = -1;
    view     = -1;
    pln      = -1;
    ch       = -1;
    pdg      = -1;
    cyc      = -1;
    pe       = -1.e-5;
    lope     = -1.e-5;
    isohit   = false;
    recon_id = -1;
    hit_id   = -1;
  }
};


class TrackIng{
public:
  Int_t   mod;
  Int_t   view;
  Int_t   ipln;
  Int_t   fpln;
  Float_t ixy;
  Float_t fxy;
  Float_t iz;
  Float_t fz;
  Float_t slope;
  Float_t intcpt;
  Float_t ang;
  Float_t clstime;
  Bool_t  veto;
  Bool_t  edge;
  Bool_t  stop;
  Float_t vetodist;

  vector<Hits> hit;  

  void clear(){
    mod      = -1;
    view     = -1;
    ipln     = -1;
    fpln     = -1;
    ixy      = -1.e-5;
    fxy      = -1.e-5;
    iz       = -1.e-5;
    fz       = -1.e-5;
    slope    = -1.e-5;
    intcpt   = -1.e-5;
    ang      = -1.e-5;
    clstime  = -1.e-5;
    veto     = false;
    edge     = false;
    stop     = false;
    vetodist = -1.e-5;
    hit.clear();
  }
};

class TrackWM{
public:

  Int_t   view;
  Int_t   ipln;
  Int_t   fpln;
  Float_t ixy;
  Float_t fxy;
  Float_t iz;
  Float_t fz;
  Float_t slope;
  Float_t intcpt;
  Float_t ang;
  Float_t clstime;
  Bool_t  veto;
  Bool_t  edge;
  Bool_t  stop;  
  Int_t   ing_imod;
  Int_t   ing_fmod;
  Int_t   ing_ipln;
  Int_t   ing_fpln;
  Bool_t  ing_trk;
  Int_t   ing_num;
  Bool_t  pm_stop;
  Bool_t  ing_stop;
  Bool_t  pm_match;
  Int_t   iron_pene;
  vector<Hits> hit;  

  //For joint matching study 
  Float_t diff_pos;
  Float_t diff_time;
  Float_t diff_ang;

  void clear(){
    view      = -1;
    ipln      = -1;
    fpln      = -1;
    ixy       = -1.e-5;
    fxy       = -1.e-5;
    iz        = -1.e-5;
    fz        = -1.e-5;
    slope     = -1.e-5;
    intcpt    = -1.e-5;
    ang       = -1.e-5;
    clstime   = -1.e-5;
    veto      = false;
    edge      = false;
    stop      = false;
    ing_imod  = -1;
    ing_fmod  = -1;
    ing_ipln  = -1;
    ing_fpln  = -1;
    ing_trk   = false;
    ing_num   = -1;
    pm_stop   = false;
    ing_stop  = false;
    pm_match  = false;
    iron_pene = -1;

    //For joint matching study 
    diff_pos  = -1000;
    diff_time = -1000;
    diff_ang  = -1000;

    hit.clear();
  }
};


class TrackPM{
public:

  Int_t   view;
  Int_t   ipln;
  Int_t   fpln;
  Float_t ixy;
  Float_t fxy;
  Float_t iz;
  Float_t fz;
  Float_t slope;
  Float_t intcpt;
  Float_t ang;
  Float_t clstime;
  Bool_t  veto;
  Bool_t  edge;
  Bool_t  stop;  
  Int_t   ing_imod;
  Int_t   ing_fmod;
  Int_t   ing_ipln;
  Int_t   ing_fpln;
  Int_t   wm_ipln;
  Int_t   wm_fpln;
  Bool_t  ing_trk;
  Int_t   ing_num;
  Bool_t  wm_trk;
  Int_t   wm_num;
  Bool_t  pm_stop;
  Bool_t  wm_stop;
  Bool_t  ing_stop;
  Int_t   iron_pene;
  vector<Hits> hit;  

  //For joint matching study 
  Float_t diff_pos;
  Float_t diff_time;
  Float_t diff_ang;


  void clear(){
    view      = -1;
    ipln      = -1;
    fpln      = -1;
    ixy       = -1.e-5;
    fxy       = -1.e-5;
    iz        = -1.e-5;
    fz        = -1.e-5;
    slope     = -1.e-5;
    intcpt    = -1.e-5;
    ang       = -1.e-5;
    clstime   = -1.e-5;
    veto      = false;
    edge      = false;
    stop      = false;
    ing_imod  = -1;
    ing_fmod  = -1;
    ing_ipln  = -1;
    ing_fpln  = -1;
    wm_ipln   = -1;
    wm_fpln   = -1;
    ing_trk   = false;
    ing_num   = -1;
    wm_trk    = false;
    wm_num    = -1;
    pm_stop   = false;
    wm_stop   = false;
    ing_stop  = false;
    iron_pene = -1;

    //For joint matching study 
    diff_pos  = -1000;
    diff_time = -1000;
    diff_ang  = -1000;

    hit.clear();
  }
};

class Trk{
public:

  Int_t   startmod;
  Int_t   stopmodx;
  Int_t   stopmody;

  Float_t x;
  Float_t y;
  Float_t z;
  Float_t zx;
  Float_t zy;

  Int_t   startxpln;
  Int_t   startypln;
  Float_t startxch;
  Float_t startych;
  Int_t   endxpln;
  Int_t   endypln;
  Float_t endxch;
  Float_t endych;
  Float_t thetax;
  Float_t thetay;
  Float_t angle;
  Int_t   ing_startmod;
  Int_t   ing_endmod;
  Int_t   ing_startpln;
  Int_t   ing_endpln;
  Bool_t  ing_trk;
  Bool_t  pm_stop;
  Bool_t  ing_stop;
  Float_t sci_range;
  Float_t iron_range;
  Int_t   iron_pene;
  Bool_t  vetowtracking; //Upstream VETO
  Bool_t  edgewtracking; //Fiducial CUT
  Int_t   pdg;
  Float_t trkpe;
  Float_t slopex;
  Float_t slopey;
  Float_t intcptx;
  Float_t intcpty;
  Float_t mucl;
  Int_t   oneview;

  Int_t hnum;
  Int_t vnum;


  //For joint matching study 
  //WM-ING
  Float_t diff_posx;
  Float_t diff_timex;
  Float_t diff_angx;
  Float_t diff_posy;
  Float_t diff_timey;
  Float_t diff_angy;

  //PM-WM
  Float_t diff_posx2;
  Float_t diff_timex2;
  Float_t diff_angx2;
  Float_t diff_posy2;
  Float_t diff_timey2;
  Float_t diff_angy2;

  //PM-ING
  Float_t diff_posx3;
  Float_t diff_timex3;
  Float_t diff_angx3;
  Float_t diff_posy3;
  Float_t diff_timey3;
  Float_t diff_angy3;

  void clear(){

    startmod      = -1;
    stopmodx      = -1;
    stopmody      = -1;

    x             = -1.e-5;
    y             = -1.e-5;
    z             = -1.e-5;
    zx            = -1.e-5;
    zy            = -1.e-5;
    startxpln     = -1;
    startypln     = -1;
    startxch      = -1.e-5;
    startych      = -1.e-5;
    endxpln       = -1;
    endypln       = -1;
    endxch        = -1.e-5;
    endych        = -1.e-5;
    thetax        = -1.e-5;
    thetay        = -1.e-5;
    angle         = -1.e-5;
    ing_startmod  = -1;
    ing_endmod    = -1;
    ing_startpln  = -1;
    ing_endpln    = -1;
    ing_trk       = false;
    pm_stop       = false;
    ing_stop      = false;
    sci_range     = -1.e-5;
    iron_range    = -1.e-5;
    iron_pene     = -1;
    vetowtracking = false; //Upstream VETO
    edgewtracking = false; //Fiducial CUT
    pdg           = -1;
    trkpe         = -1.e-5;
    slopex        = -1.e-5;
    slopey        = -1.e-5;
    intcptx       = -1.e-5;
    intcpty       = -1.e-5;
    mucl          = -1.e-5;
    oneview       = -1;

    hnum =  -1;
    vnum =  -1;

    //For joint matching study 
    diff_posx  = -1000;
    diff_timex = -1000;
    diff_angx  = -1000;
    diff_posy  = -1000;
    diff_timey = -1000;
    diff_angy  = -1000;

    diff_posx2  = -1000;
    diff_timex2 = -1000;
    diff_angx2  = -1000;
    diff_posy2  = -1000;
    diff_timey2 = -1000;
    diff_angy2  = -1000;

    diff_posx3  = -1000;
    diff_timex3 = -1000;
    diff_angx3  = -1000;
    diff_posy3  = -1000;
    diff_timey3 = -1000;
    diff_angy3  = -1000;

  }
};


class AnaTrack{
public:

  Int_t   Ntrack;
  Int_t   Ningtrack;
  Float_t clstime;
  Bool_t  vetowtracking; //Upstream VETO
  Bool_t  edgewtracking; //Fiducial CUT

  vector<Trk> trk;

  void clear(){
    Ntrack        = -1;
    Ningtrack     = -1;
    clstime       = -1.e-5;
    vetowtracking = false; //Upstream VETO
    edgewtracking = false; //Fiducial CUT
    trk.clear();
  }
};
#endif
