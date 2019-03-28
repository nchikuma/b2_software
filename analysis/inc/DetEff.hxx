#ifndef ____B2_Plot___
#define ____B2_Plot___

#include "FluxTuning.h"
#include "TwoDimRecon.hxx"
#include "ThreeDimRecon.hxx"
#include <TArrayF.h>

const int num_reweight = 147;

const float hittime_off[3] = {490.,470.,480.}; //ING, PM, WM

const int   IronPeneCut0     = 2; //pln ING
const int   IronPeneCut1     = 1; //pln PM
const int   IronPeneCut2     = 2; //pln WM
const float BeamTimingCut    = 100.; //ns
const float BeamTimingCut_MC = 300.; //ns

const int   VetoCutPln0 = 1; //ING
const int   VetoCutPln1 = 1; //PM
const int   VetoCutPln2 = 4; //WM  //TopView (view:1,pln:1)

const float FVcutX0     = 350.; //mm half  //ING
const float FVcutY0     = 350.; //mm half
const int   FVcutPlnZ0  = 8   ; //SideView
const float FVcutX1     = 350.; //mm half  //PM
const float FVcutY1     = 350.; //mm half
const int   FVcutPlnZ1  = 15  ; //SideView
const float FVcutX2     = 350.; //mm half  //WM
const float FVcutY2     = 350.; //mm half
const int   FVcutPlnZ2  = 15  ; // SideView (view:0,pln:5)

const int   INGStopPlnZ =    9;  //pln
const float INGStopXY   = 500.; //mm

const int   CoT[3][3] = {
  {C_B2MotherPosX+C_B2INGPosX,C_B2MotherPosY+C_B2INGPosY,C_B2MotherPosZ+C_B2INGPosZ},
  {C_B2MotherPosX+C_B2CHPosX ,C_B2MotherPosY+C_B2CHPosY ,C_B2MotherPosZ+C_B2CHPosZ },
  {C_B2MotherPosX+C_B2WMPosX ,C_B2MotherPosY+C_B2WMPosY ,C_B2MotherPosZ+C_B2WMPosZ }
};

const float density[3] = {
  //6.087,  //7.860,  //iron         //g/cm3
  7.8494,  //7.860,  //iron         //g/cm3
  1.01366,//1.032,  //scintillator //g/cm3
  1.00075 //1.000   //water        //g/cm3
};
const float density_wallbg = 2.2; //g/cm3

//INGRID
//const float prob_fe_ch = 1./(1.+C_INGTotMassIron/C_INGTotMassScinti);
//const float thick_ch   = C_INGScintiThick*2*11/10.;
//const float thick_fe   = C_INGIronThick    *9 /10.;

const float Avogadro = 6.02214e+23;
const float thickness[3] = {
  //C_INGScintiThick*2*11/10.+C_INGIronThick*9/10.,   //cm ING, thickness of iron
  C_INGIronThick*9/10.,                             //cm ING, thickness of iron
  C_INGScintiThick*2/10.+C_PMScintiThick*2*17/10.,  //cm PM, thickness of matter
  C_WMWaterTargetSizeZ/10.                          //cm WM, thickness of water
};
const float thickness_wallbg   = 600.; //cm
const float thickness_ceiling  = 904.3; //cm
const float thickness_floor    = 904.3; //cm
const float thickness_pillar_r = 400.; //cm
const float thickness_pillar_l = 713.3; //cm


const float area0[3] = {
  C_INGScintiLength   /10. * C_INGScintiLength   /10.,  //cm2 ING, MC flux area
  C_PMScintiLength    /10. * C_PMScintiLength    /10.,  //cm2 PM
  C_WMWaterTargetSizeX/10. * C_WMWaterTargetSizeY/10.   //cm2 WM
};
const float area[3] = {
  FVcutX0/10.*2. * FVcutY0/10.*2.,  //cm2 ING, MC flux area
  FVcutX1/10.*2. * FVcutY1/10.*2.,  //cm2 PM
  FVcutX2/10.*2. * FVcutY2/10.*2.   //cm2 WM
};


const int NumHist  = 18;
const int NumHist2 =  9;

void SetHist(int i,
    int &nbin, double &min, double &max,
    string &title, string &xtitle, string &ytitle,
    string &name
    )
{
  if(i<NumHist){
    name = Form("h1_%02d",i);
    if(i<=2){
      nbin   =  50 ; 
      min    =   0.;  //GeV/c
      max    =   5.;  //GeV/c
      if(i==0) title = "Initial Momentum : ING ";
      if(i==1) title = "Initial Momentum : PM ";
      if(i==2) title = "Initial Momentum : WM ";
      xtitle = "Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=5){
      nbin   =    180; 
      min    =     0.;  //deg
      max    =   180.;  //deg
      if(i==3) title = "Initial Angle: ING ";
      if(i==4) title = "Initial Angle: PM ";
      if(i==5) title = "Initial Angle: WM ";
      xtitle = "Angle [deg]";
      ytitle = "";
    }
    else if(i<=8){
      nbin   =  50 ; 
      min    =   0.;  //GeV/c
      max    =   5.;  //GeV/c
      if(i==6) title = "Initial Momentum (Detected Event): ING ";
      if(i==7) title = "Initial Momentum (Detected Event): PM ";
      if(i==8) title = "Initial Momentum (Detected Event): WM ";
      xtitle = "Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=11){
      nbin   =    180; 
      min    =     0.;  //deg
      max    =   180.;  //deg
      if(i== 9) title = "Initial Angle (Detected Event): ING ";
      if(i==10) title = "Initial Angle (Detected Event): PM ";
      if(i==11) title = "Initial Angle (Detected Event): WM ";
      xtitle = "Angle [deg]";
      ytitle = "";
    }
    else if(i<=14){
      nbin   =  50 ; 
      min    =   0.;  //GeV/c
      max    =   5.;  //GeV/c
      if(i==12) title = "Initial Momentum (Initial Particle Detected): ING ";
      if(i==13) title = "Initial Momentum (Initial Particle Detected): PM ";
      if(i==14) title = "Initial Momentum (Initial Particle Detected): WM ";
      xtitle = "Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=17){
      nbin   =    180; 
      min    =     0.;  //deg
      max    =   180.;  //deg
      if(i==15) title = "Initial Angle (Initial Particle Detected): ING ";
      if(i==16) title = "Initial Angle (Initial Particle Detected): PM ";
      if(i==17) title = "Initial Angle (Initial Particle Detected): WM ";
      xtitle = "Angle [deg]";
      ytitle = "";
    }
    else{
      cout << "Unexpected ID for 1D histogram." << endl;
      cout << " i=" << i << endl;
      exit(1);
    }
  }
  else{
    cout << "Unexpected ID for 1D histogram." << endl;
    cout << " i=" << i << endl;
    exit(1);
  }
}

void SetHist2(int i,
    int& nbin1, double& min1, double& max1,
    int& nbin2, double& min2, double& max2,
    string &title, string &xtitle, string &ytitle,
    string &name
    )
{
  if(i<NumHist2){
    name = Form("h2_%02d",i);
    if(i<=8){
      nbin1   =  50; 
      min1    =   0.; // GeV/c
      max1    =   5.; // GeV/c
      nbin2   = 180; 
      min2    =   0.; // deg
      max2    = 180.; // deg 
      title   = "";
      if(i== 0) title += "Initial Kinematics : ING";
      if(i== 1) title += "Initial Kinematics : PM ";
      if(i== 2) title += "Initial Kinematics : WM ";
      if(i== 3) title += "Detected Events : ING";
      if(i== 4) title += "Detected Events : PM ";
      if(i== 5) title += "Detected Events : WM ";
      if(i== 6) title += "Initial Particle Detected : ING";
      if(i== 7) title += "Initial Particle Detected : PM ";
      if(i== 8) title += "Initial Particle Detected : WM ";
      xtitle = "Momentum [GeV/c]";
      ytitle = "Angle [deg]";
    }
    else{
      cout << "Unexpected ID for 2D histogram." << endl;
      cout << " i=" << i << endl;
      exit(1);
    }
  }
  else{
    cout << "Unexpected ID for 2D histogram." << endl;
    cout << " i=" << i << endl;
    exit(1);
  }
}


#endif

