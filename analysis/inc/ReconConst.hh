#ifndef __RECON_CONST__
#define __RECON_CONST__
#include<iostream>
#include<sstream>
using namespace std;

//====================================
//========== Data directory ==========
//====================================
#define dst_ingrid     "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_dst"
#define dst_ingrid2    "/home/t2k/nchikuma/b2_data2/data_dst"
#define dst_wagasci    "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_dst_wagasci"
#define cosmic_data    "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_cosmic"
#define cosmic_wagasci "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_cosmic_wagasci"
#define mc_sandmu      "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/mc_sandmu"
#define mc_cosmic      "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/mc_cosmic"
#define mc_neut        "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/mc_neut"
#define mc_muon        "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/mc_muon_test"


const static int nMod    = 17;
const static int nModIng = 16;
const static int nView   =  2;
const static int nTPL    = 11;
const static int nVETO   =  4;
const static int nPln    = nTPL+nVETO;
const static int nCh     = 24;
const static int nCyc    = 23;
const static int RVETO   = 11;
const static int LVETO   = 12;
const static int BVETO   = 13;
const static int UVETO   = 14;
const static int FromX   =  0;
const static int FromY   =  1;
const static int X       =  0;
const static int Y       =  1;
const static int Z       =  2;

const static float CTtimeBase   = 1.08e-6; //Base time of CT [sec]
const static int   GapbwBunch   =     581; //[nsec]
const static int   TDCOffset    =     300; //[nsec]
const static int   bunch1st_cyc =       4;
const static int   fIntTime     =     480; //[nsec]
const static int   fRstTime     =     100; //[nsec]
const static int   fExpBeamTime =     200; //[nsec]


//========================================================
//========== Threshold for Basic Reconstruction ==========
//========================================================
const static int   cut_nactpln =   1;
const static float cut_layerpe = 6.5;
const static int   beamontime  = 100; //Ontime


//================================
//========== IngEvtDisp ==========
//================================
const static Int_t  NIron       = 10; //Actual #of iron: 9
const static Int_t  NumPln      = 11;
const static Int_t  Nscinti     = 24;
const static Int_t  NscintiVETO = 22;
const static Int_t  doffset     =  4; //*ScintiWidth

const static float  IronDensity      =         7.87; //g/cm^3
const static double IronThick        =          6.5; //cm
const static double ScintiWidth      =          5.0; //cm
const static double ScintiThick      =          1.0; //cm
const static double PlnThick         =          4.2; //cm
const static double VetoOffsetZX     =         -2.0; //cm
const static double VetoStartZ       =         +1.0; //cm
const static double VetoOffsetZY     =         -1.0; //cm
const static double VetoOffsetRight  = -10.575-2.5 ; //cm
const static double VetoOffsetLeft   = 130.9  -2.5 ; //cm
const static double VetoOffsetBottom =         -5.9; //cm
const static double VetoOffsetUp     =        128.4; //cm
const static double PlnIronThick     = PlnThick + IronThick;
const static double BWModule         =         30.0;

const static float DISTANCE_MODULE = 150; //cm
const static float StartModuleZ    =   0;
const static float EndModuleZ      = IronThick * NIron + PlnThick * NumPln;
const static float StartVetoZ      = 1.5;
const static float EndVetoZ        = 1.5 + ScintiWidth * (NscintiVETO-1);
const static float StartTPLXY      =   0;
const static float EndTPLXY        = ScintiWidth * (Nscinti-1);

const static float  MIPdepEIron = 1.45; //MeV/cm
const static float  MIPdepE     = 1.99; //MeV/cm
const static float  Mp          =  938; //MeV
const static float  Mn          =  940; //MeV
const static float  Mmu         =  106; //MeV
const static float  NuclearPotE =    0; //MeV | 0 is temporary value


const static double VetoZY[15] =
  {
    0, //0
    0, //1
    0, //2
    0, //3
    0, //4
    0, //5
    0, //6
    0, //7
    0, //8
    0, //9
    0, //10
    -10.575-2.5, //Right 
    130.9  -2.5, //Left
    -5.9,        //Bottom
    128.4        //Up
  };

const static double PlnZ[11] =
  {
    0  * ( PlnIronThick ),
    1  * ( PlnIronThick ),
    2  * ( PlnIronThick ),
    3  * ( PlnIronThick ),
    4  * ( PlnIronThick ),
    5  * ( PlnIronThick ),
    6  * ( PlnIronThick ),
    7  * ( PlnIronThick ),
    8  * ( PlnIronThick ),
    9  * ( PlnIronThick ),
    10 * ( PlnIronThick )
  };

float fVetoXY(int mod, int pln){
  float list[7] = 
    {
      -13.075, //Right
        128.4, //Left
         -8.4, //Bottom
        125.9, //Up
        -21.6, //Double used Right
        -24.1, //Double used Bottom
        -24.1  //Bottom of mod#7

      //### 2010/5/19???
      //### IT SHOULD BE UPDATED
      //-10.575, //Right
      //  130.9, //Left
      //   -5.9, //Bottom
      //  128.4, //Up
      //  -19.1, //Double used Right
      //  -21.6, //Double used Bottom
      //  -21.6  //Bottom of mod#7
    };

  int idcase=-1;
 
  if(mod==0){
    if(pln == RVETO) idcase = 0;
    if(pln == LVETO) idcase = 1;
    if(pln == BVETO) idcase = 2;
    if(pln == UVETO) idcase = 3;
  }
  else if(1<=mod && mod<=6){
    if(pln == RVETO) idcase = 4;
    if(pln == LVETO) idcase = 1;
    if(pln == BVETO) idcase = 2;
    if(pln == UVETO) idcase = 3;
  }
  else if(mod==7){
    if(pln == RVETO) idcase = 0;
    if(pln == LVETO) idcase = 1;
    if(pln == BVETO) idcase = 6;
    if(pln == UVETO) idcase = 3;
  }
  //else if(8<=mod && mod<=13){
  else if(8<=mod && mod<=14){
    if(pln == RVETO) idcase = 0;
    if(pln == LVETO) idcase = 1;
    if(pln == BVETO) idcase = 5;
    if(pln == UVETO) idcase = 3;
  }

  if(idcase==-1) return -1.e-5;
  else           return list[idcase];
}

bool isVeto(int pln){
  if(0<=pln && pln<=10 )      return false;
  else if(pln>=11 && pln<=14) return true;

  return false;
}

bool isPM(int rmm, int tfb){
  if(rmm==4 && tfb>=26) return true;
  else                  return false;
}

bool isPM(int mod){
  if(mod==16) return true;
  else        return false;
}
  
#endif
