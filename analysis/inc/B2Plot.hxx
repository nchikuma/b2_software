#ifndef ____B2_Plot___
#define ____B2_Plot___

#include "FluxTuning.h"
#include "TwoDimRecon.hxx"
#include "ThreeDimRecon.hxx"
#include <TArrayF.h>

const int num_reweight = 161;

const float hittime_off[3] = {490.,470.,480.}; //ING, PM, WM

const int   IronPeneCut0     = 2; //pln ING
const int   IronPeneCut1     = 1; //pln PM
const int   IronPeneCut2     = 2; //pln WM
const float BeamTimingCut    = 100.; //ns
const float BeamTimingCut_MC = 300.; //ns

int   VetoCutPln0 = 1; //ING
int   VetoCutPln1 = 1; //PM
int   VetoCutPln2 = 4; //WM  //TopView (view:1,pln:1)

float FVcutX0     = 350.; //mm half  //ING
float FVcutY0     = 350.; //mm half
int   FVcutPlnZ0  = 8   ; //SideView
float FVcutX1     = 350.; //mm half  //PM
float FVcutY1     = 350.; //mm half
int   FVcutPlnZ1  = 15  ; //SideView
float FVcutX2     = 350.; //mm half  //WM
float FVcutY2     = 350.; //mm half
int   FVcutPlnZ2  = 15  ; // SideView (view:0,pln:5)
float ACCEPT_OFFSETZ = 0.; //mm
float ACCEPT_OFFSETY = 0.; //mm
float ACCEPT_OFFSETX = 0.; //mm


const int   INGStopPlnZ =    9;  //pln
const float INGStopXY   = 500.; //mm

double PION_MOM_TH  [3] = {0.2,0.1,0.2}; //GeV/c
double PION_ANG_TH  [3] = {70.,70.,91.}; //deg
double PROTON_MOM_TH[3] = {0.6,0.5,0.6}; //GeV/c
double PROTON_ANG_TH[3] = {70.,70.,91.}; //deg

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


const int NumHist  = 455;
const int NumHist2 = 115;
const int NumHist3 =  41;
const int NumSelec2 = 8;
const int NumSelec3 = 11;


void SetHist(int i,
    int &nbin, double &min, double &max,
    string &title, string &xtitle, string &ytitle,
    string &name
    )
{
  if(i<NumHist){
    name = Form("h1_%02d",i);
    if(i<=2){
      nbin   = 600; 
      min    = -300.;  //ns
      max    =  300.;  //ns
      if(i==0) title = "Hit Timing on 2D track: ING ";
      if(i==1) title = "Hit Timing on 2D track: PM ";
      if(i==2) title = "Hit Timing on 2D track: WM ";
      xtitle = "Hit Timing [ns]";
      ytitle = "";
    }
    else if(i<=5){
      nbin   =  12; 
      min    =  -1;  
      max    =  11;  
      if(i==3) title = "1st Track Stopping Position : ING ";
      if(i==4) title = "1st Track Stopping Position : PM ";
      if(i==5) title = "1st Track Stopping Position : WM ";
      xtitle = "Penetrated Iron Plane";
      ytitle = "";
    }
    else if(i<=8){
      nbin   =  80; 
      min    = -40.;  //deg
      max    =  40.;  //deg
      if(i==6) title = "2D Track Matching : WM-ING ";
      if(i==7) title = "2D Track Matching : PM-WM ";
      if(i==8) title = "2D Track Matching : PM-ING ";
      xtitle = "Angle Difference [deg]";
      ytitle = "";
    }
    else if(i<=11){
      nbin   =  400; 
      min    = -200.;  //mm
      max    =  200.;  //mm
      if(i== 9) title = "2D Track Matching : WM-ING ";
      if(i==10) title = "2D Track Matching : PM-WM ";
      if(i==11) title = "2D Track Matching : PM-ING ";
      xtitle = "Track Distance at the Middle [mm]";
      ytitle = "";
    }
    else if(i<=14){
      nbin   =  400; 
      min    = -200.;  //ns
      max    =  200.;  //ns
      if(i==12) title = "2D Track Matching : WM-ING ";
      if(i==13) title = "2D Track Matching : PM-WM ";
      if(i==14) title = "2D Track Matching : PM-ING ";
      xtitle = "Timing Difference [ns]";
      ytitle = "";
    }
    else if(i<=17){
      nbin   =  50 ; 
      min    = -25.;  //ns
      max    =  25.;  //ns
      if(i==15) title = "3D Track Matching : ING ";
      if(i==16) title = "3D Track Matching : PM ";
      if(i==17) title = "3D Track Matching : WM ";
      xtitle = "Differnce in Start Z Position [Plane]";
      ytitle = "";
    }
    else if(i<=20){
      nbin   = 600; 
      min    = -300.;  //ns
      max    =  300.;  //ns
      if(i==18) title = "Hit Timing: ING ";
      if(i==19) title = "Hit Timing: PM ";
      if(i==20) title = "Hit Timing: WM ";
      xtitle = "Hit Timing [ns]";
      ytitle = "";
    }
    else if(i<=23){
      nbin   = 600; 
      min    = -300.;  //ns
      max    =  300.;  //ns
      if(i==21) title = "Hit Timing (nhit>2 +/- 50ns): ING ";
      if(i==22) title = "Hit Timing (nhit>2 +/- 50ns): PM ";
      if(i==23) title = "Hit Timing (nhit>2 +/- 50ns): WM ";
      xtitle = "Hit Timing [ns]";
      ytitle = "";
    }
    else if(i<=26){
      nbin   = 600; 
      min    = -300.;  //ns
      max    =  300.;  //ns
      if(i==24) title = "Hit Timing (After Time Clustering): ING ";
      if(i==25) title = "Hit Timing (After Time Clustering): PM ";
      if(i==26) title = "Hit Timing (After Time Clustering): WM ";
      xtitle = "Hit Timing [ns]";
      ytitle = "";
    }
    else if(i<=29){
      nbin   = 600; 
      min    = -300.;  //ns
      max    =  300.;  //ns
      if(i==27) title = "Hit Timing Deference (After Time Clustering): ING ";
      if(i==28) title = "Hit Timing Deference (After Time Clustering): PM ";
      if(i==29) title = "Hit Timing Deference (After Time Clustering): WM ";
      xtitle = "Hit Timing Difference from its Reference [ns]";
      ytitle = "";
    }
    else if(i<=32){
      nbin   =  15; 
      min    =   0.;  //pln
      max    =  15.;  //pln
      if(i==30) title = "Vertexing : ING ";
      if(i==31) title = "Vertexing : PM ";
      if(i==32) title = "Vertexing : WM ";
      xtitle = "Difference in Z Position [Plane]";
      ytitle = "";
    }
    else if(i<=35){
      nbin   =  300; 
      min    =    0.;  //mm
      max    =  300.;  //mm
      if(i==33) title = "Vertexing : ING ";
      if(i==34) title = "Vertexing : PM ";
      if(i==35) title = "Vertexing : WM ";
      xtitle = "Difference in XY Position [mm]";
      ytitle = "";
    }
    else if(i<=38){
      nbin   = 5500; 
      min    = 2500.;  //ns
      max    = 8000.;  //ns
      if(i==36) title = "Timing of Hits on Track (After Vertexing) : ING ";
      if(i==37) title = "Timing of Hits on Track (After Vertexing) : PM ";
      if(i==38) title = "Timing of Hits on Track (After Vertexing) : WM ";
      xtitle = "Hit Timing [ns]";
      ytitle = "";
    }
    else if(i<=40){
      nbin   =  600; 
      min    = -300.;  //mm
      max    =  300.;  //mm
      if(i==39) title = "2D Track Matching with ING Hits : WM-ING ";
      if(i==40) title = "2D Track Matching with PM Hits : WM-PM ";
      xtitle = "Distance to Hit Position from Track [mm]";
      ytitle = "";
    }
    else if(i<=43){
      nbin   = 25; 
      min    = 0;  //pln
      max    = 25;  //pln
      if(i==41) title = "Reconstructed Vertex (After Vertexing) : ING ";
      if(i==42) title = "Reconstructed Vertex (After Vertexing) : PM ";
      if(i==43) title = "Reconstructed Vertex (After Vertexing) : WM ";
      xtitle = "Vertex Position Z [plane]";
      ytitle = "";
    }
    else if(i<=46){
      nbin   = 1500; 
      min    = -750.;  //mm
      max    =  750.;  //mm
      if(i==44) title = "Reconstructed Vertex (After Vertexing) : ING ";
      if(i==45) title = "Reconstructed Vertex (After Vertexing) : PM ";
      if(i==46) title = "Reconstructed Vertex (After Vertexing) : WM ";
      xtitle = "Vertex Position X [mm]";
      ytitle = "";
    }
    else if(i<=49){
      nbin   = 1500; 
      min    = -750.;  //mm
      max    =  750.;  //mm
      if(i==47) title = "Reconstructed Vertex (After Vertexing) : ING ";
      if(i==48) title = "Reconstructed Vertex (After Vertexing) : PM ";
      if(i==49) title = "Reconstructed Vertex (After Vertexing) : WM ";
      xtitle = "Vertex Position Y [mm]";
      ytitle = "";
    }
    else if(i<=52){
      nbin   =  90; 
      min    =   0.;  //deg
      max    =  90.;  //deg
      if(i==50) title = "1st Track Angle (After Vertexing) : ING ";
      if(i==51) title = "1st Track Angle (After Vertexing) : PM ";
      if(i==52) title = "1st Track Angle (After Vertexing) : WM ";
      xtitle = "Reconstructed Angle [deg]";
      ytitle = "";
    }
    else if(i<=55){
      nbin   = 600; 
      min    = -300.;  //ns
      max    =  300.;  //ns
      if(i==53) title = "Event Timing : ING ";
      if(i==54) title = "Event Timing : PM ";
      if(i==55) title = "Event Timing : WM ";
      xtitle = "Difference from Expected Beam Timing [ns]";
      ytitle = "";
    }
    else if(i<=58){
      nbin   = 25; 
      min    = 0;  //pln
      max    = 25;  //pln
      if(i==56) title = "Reconstructed Vertex (After Beam Timing Cut) : ING ";
      if(i==57) title = "Reconstructed Vertex (After Beam Timing Cut) : PM ";
      if(i==58) title = "Reconstructed Vertex (After Beam Timing Cut) : WM ";
      xtitle = "Vertex Position Z [plane]";
      ytitle = "";
    }
    else if(i<=61){
      nbin   = 1500; 
      min    = -750.;  //mm
      max    =  750.;  //mm
      if(i==59) title = "Reconstructed Vertex (After Beam Timing Cut) : ING ";
      if(i==60) title = "Reconstructed Vertex (After Beam Timing Cut) : PM ";
      if(i==61) title = "Reconstructed Vertex (After Beam Timing Cut) : WM ";
      xtitle = "Vertex Position X [mm]";
      ytitle = "";
    }
    else if(i<=64){
      nbin   = 1500; 
      min    = -750.;  //mm
      max    =  750.;  //mm
      if(i==62) title = "Reconstructed Vertex (After Beam Timing Cut) : ING ";
      if(i==63) title = "Reconstructed Vertex (After Beam Timing Cut) : PM ";
      if(i==64) title = "Reconstructed Vertex (After Beam Timing Cut) : WM ";
      xtitle = "Vertex Position Y [mm]";
      ytitle = "";
    }
    else if(i<=67){
      nbin   = 1500; 
      min    = -750.;  //mm
      max    =  750.;  //mm
      if(i==65) title = "Reconstructed Vertex (After Upst Veto Cut) : ING ";
      if(i==66) title = "Reconstructed Vertex (After Upst Veto Cut) : PM ";
      if(i==67) title = "Reconstructed Vertex (After Upst Veto Cut) : WM ";
      xtitle = "Vertex Position X [mm]";
      ytitle = "";
    }
    else if(i<=70){
      nbin   = 1500; 
      min    = -750.;  //mm
      max    =  750.;  //mm
      if(i==68) title = "Reconstructed Vertex (After Upst Veto Cut) : ING ";
      if(i==69) title = "Reconstructed Vertex (After Upst Veto Cut) : PM ";
      if(i==70) title = "Reconstructed Vertex (After Upst Veto Cut) : WM ";
      xtitle = "Vertex Position Y [mm]";
      ytitle = "";
    }
    else if(i<=73){
      nbin   =  90; 
      min    =   0.;  //deg
      max    =  90.;  //deg
      if(i==71) title = "1st Track Angle (After FV Cut) : ING ";
      if(i==72) title = "1st Track Angle (After FV Cut) : PM ";
      if(i==73) title = "1st Track Angle (After FV Cut) : WM ";
      xtitle = "Reconstructed Angle [deg]";
      ytitle = "";
    }
    else if(i<=76){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==74) title = "Expected Nu Flux by JNUBEAM : ING ";
      if(i==75) title = "Expected Nu Flux by JNUBEAM : PM ";
      if(i==76) title = "Expected Nu Flux by JNUBEAM : WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=79){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==77) title = "All Interaction in FV: ING ";
      if(i==78) title = "All Interaction in FV: PM ";
      if(i==79) title = "All Interaction in FV: WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=82){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==80) title = "CC-inclusive Interaction in FV : ING ";
      if(i==81) title = "CC-inclusive Interaction in FV : PM ";
      if(i==82) title = "CC-inclusive Interaction in FV : WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=85){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==83) title = "CC-inclusive Interaction in FV : ING ";
      if(i==84) title = "CC-inclusive Interaction in FV : PM ";
      if(i==85) title = "CC-inclusive Interaction in FV : WM ";
      xtitle = "Muon Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=88){
      nbin   =    180; 
      min    =     0.;  //GeV/c
      max    =   180.;  //GeV/c
      if(i==86) title = "CC-inclusive Interaction in FV : ING ";
      if(i==87) title = "CC-inclusive Interaction in FV : PM ";
      if(i==88) title = "CC-inclusive Interaction in FV : WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=91){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==89) title = "Selected CC-inc (After FV Cut): ING ";
      if(i==90) title = "Selected CC-inc (After FV Cut): PM ";
      if(i==91) title = "Selected CC-inc (After FV Cut): WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=94){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==92) title = "Selected CC-inc (After FV Cut): ING ";
      if(i==93) title = "Selected CC-inc (After FV Cut): PM ";
      if(i==94) title = "Selected CC-inc (After FV Cut): WM ";
      xtitle = "Muon Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=97){
      nbin   =    180; 
      min    =     0.;  //GeV/c
      max    =   180.;  //GeV/c
      if(i==95) title = "Selected CC-inc (After FV Cut): ING ";
      if(i==96) title = "Selected CC-inc (After FV Cut): PM ";
      if(i==97) title = "Selected CC-inc (After FV Cut): WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=100){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i== 98) title = "Selected CC-inc (After Track Angle Cut): ING ";
      if(i== 99) title = "Selected CC-inc (After Track Angle Cut): PM ";
      if(i==100) title = "Selected CC-inc (After Track Angle Cut): WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=103){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==101) title = "Selected CC-inc (After Track Angle Cut): ING ";
      if(i==102) title = "Selected CC-inc (After Track Angle Cut): PM ";
      if(i==103) title = "Selected CC-inc (After Track Angle Cut): WM ";
      xtitle = "Muon Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=106){
      nbin   =    180; 
      min    =     0.;  //GeV/c
      max    =   180.;  //GeV/c
      if(i==104) title = "Selected CC-inc (After Track Angle Cut): ING ";
      if(i==105) title = "Selected CC-inc (After Track Angle Cut): PM ";
      if(i==106) title = "Selected CC-inc (After Track Angle Cut): WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=109){
      nbin   =  1500; 
      min    =  -750;  //mm
      max    =   750;  //mm
      if(i==107) title = "Track Stop Position in INGRID (After Vertexing): ING ";
      if(i==108) title = "Track Stop Position in INGRID (After Vertexing): PM ";
      if(i==109) title = "Track Stop Position in INGRID (After Vertexing): WM ";
      xtitle = "Stop Position X [mm]";
      ytitle = "";
    }
    else if(i<=112){
      nbin   =  1500; 
      min    =  -750;  //mm
      max    =   750;  //mm
      if(i==110) title = "Track Stop Position in INGRID (After Vertexing): ING ";
      if(i==111) title = "Track Stop Position in INGRID (After Vertexing): PM ";
      if(i==112) title = "Track Stop Position in INGRID (After Vertexing): WM ";
      xtitle = "Stop Position Y [mm]";
      ytitle = "";
    }
    else if(i<=115){
      nbin   =  12; 
      min    =  -1;  //pln
      max    =  11;  //pln
      if(i==113) title = "Track Stop Position in INGRID (After Angle Cut): ING ";
      if(i==114) title = "Track Stop Position in INGRID (After Angle Cut): PM ";
      if(i==115) title = "Track Stop Position in INGRID (After Angle Cut): WM ";
      xtitle = "Stop Position Z [pln]";
      ytitle = "";
    }
    else if(i<=118){
      nbin   =  1500; 
      min    =  -750;  //mm
      max    =   750;  //mm
      if(i==116) title = "Track Stop Position in INGRID (After Angle Cut): ING ";
      if(i==117) title = "Track Stop Position in INGRID (After Angle Cut): PM ";
      if(i==118) title = "Track Stop Position in INGRID (After Angle Cut): WM ";
      xtitle = "Stop Position X [mm]";
      ytitle = "";
    }
    else if(i<=121){
      nbin   =  1500; 
      min    =  -750;  //mm
      max    =   750;  //mm
      if(i==119) title = "Track Stop Position in INGRID (After Angle Cut): ING ";
      if(i==120) title = "Track Stop Position in INGRID (After Angle Cut): PM ";
      if(i==121) title = "Track Stop Position in INGRID (After Angle Cut): WM ";
      xtitle = "Stop Position Y [mm]";
      ytitle = "";
    }
    else if(i<=124){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==122) title = "Selected CC-inc (After Side Escaping Cut): ING ";
      if(i==123) title = "Selected CC-inc (After Side Escaping Cut): PM ";
      if(i==124) title = "Selected CC-inc (After Side Escaping Cut): WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=127){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==125) title = "Selected CC-inc (After Side Escaping Cut): ING ";
      if(i==126) title = "Selected CC-inc (After Side Escaping Cut): PM ";
      if(i==127) title = "Selected CC-inc (After Side Escaping Cut): WM ";
      xtitle = "Muon Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=130){
      nbin   =    180; 
      min    =     0.;  //GeV/c
      max    =   180.;  //GeV/c
      if(i==128) title = "Selected CC-inc (After Side Escaping Cut): ING ";
      if(i==129) title = "Selected CC-inc (After Side Escaping Cut): PM ";
      if(i==130) title = "Selected CC-inc (After Side Escaping Cut): WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=133){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==131) title = "Selected CC-inc (After INGRID Stop Cut): ING ";
      if(i==132) title = "Selected CC-inc (After INGRID Stop Cut): PM ";
      if(i==133) title = "Selected CC-inc (After INGRID Stop Cut): WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=136){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==134) title = "Selected CC-inc (After INGRID Stop Cut): ING ";
      if(i==135) title = "Selected CC-inc (After INGRID Stop Cut): PM ";
      if(i==136) title = "Selected CC-inc (After INGRID Stop Cut): WM ";
      xtitle = "Muon Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=139){
      nbin   =    180; 
      min    =     0.;  //GeV/c
      max    =   180.;  //GeV/c
      if(i==137) title = "Selected CC-inc (After INGRID Stop Cut): ING ";
      if(i==138) title = "Selected CC-inc (After INGRID Stop Cut): PM ";
      if(i==139) title = "Selected CC-inc (After INGRID Stop Cut): WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=142){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==140) title = "Non-Selected CC-inc (After Side Escaping Cut): ING ";
      if(i==141) title = "Non-Selected CC-inc (After Side Escaping Cut): PM ";
      if(i==142) title = "Non-Selected CC-inc (After Side Escaping Cut): WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=145){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==143) title = "Non-Selected CC-inc (After Side Escaping Cut): ING ";
      if(i==144) title = "Non-Selected CC-inc (After Side Escaping Cut): PM ";
      if(i==145) title = "Non-Selected CC-inc (After Side Escaping Cut): WM ";
      xtitle = "Muon Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=148){
      nbin   =    180; 
      min    =     0.;  //GeV/c
      max    =   180.;  //GeV/c
      if(i==146) title = "Non-Selected CC-inc (After Side Escaping Cut): ING ";
      if(i==147) title = "Non-Selected CC-inc (After Side Escaping Cut): PM ";
      if(i==148) title = "Non-Selected CC-inc (After Side Escaping Cut): WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=151){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==149) title = "Non-Selected CC-inc (After INGRID Stop Cut): ING ";
      if(i==150) title = "Non-Selected CC-inc (After INGRID Stop Cut): PM ";
      if(i==151) title = "Non-Selected CC-inc (After INGRID Stop Cut): WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=154){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==152) title = "Non-Selected CC-inc (After INGRID Stop Cut): ING ";
      if(i==153) title = "Non-Selected CC-inc (After INGRID Stop Cut): PM ";
      if(i==154) title = "Non-Selected CC-inc (After INGRID Stop Cut): WM ";
      xtitle = "Muon Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=157){
      nbin   =    180; 
      min    =     0.;  //GeV/c
      max    =   180.;  //GeV/c
      if(i==155) title = "Non-Selected CC-inc (After INGRID Stop Cut): ING ";
      if(i==156) title = "Non-Selected CC-inc (After INGRID Stop Cut): PM ";
      if(i==157) title = "Non-Selected CC-inc (After INGRID Stop Cut): WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=160){
      nbin   =  90; 
      min    =   0.;  //deg
      max    =  90.;  //deg
      if(i==158) title = "1st Track Angle (After FV Cut) : ING ";
      if(i==159) title = "1st Track Angle (After FV Cut) : PM ";
      if(i==160) title = "1st Track Angle (After FV Cut) : WM ";
      xtitle = "Reconstructed Angle [deg]";
      ytitle = "";
    }
    else if(i<=163){
      nbin   =  90; 
      min    =   0.;  //deg
      max    =  90.;  //deg
      if(i==161) title = "1st Track Angle (After Side Escaping Cut) : ING ";
      if(i==162) title = "1st Track Angle (After Side Escaping Cut) : PM ";
      if(i==163) title = "1st Track Angle (After Side Escaping Cut) : WM ";
      xtitle = "Reconstructed Angle [deg]";
      ytitle = "";
    }
    else if(i<=166){
      nbin   =  90; 
      min    =   0.;  //deg
      max    =  90.;  //deg
      if(i==164) title = "1st Track Angle (After INGRID Stopping Cut) : ING ";
      if(i==165) title = "1st Track Angle (After INGRID Stopping Cut) : PM ";
      if(i==166) title = "1st Track Angle (After INGRID Stopping Cut) : WM ";
      xtitle = "Reconstructed Angle [deg]";
      ytitle = "";
    }
    else if(i<=169){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==167) title = "CC-inc (True, in FV && INGRID Stop): ING ";
      if(i==168) title = "CC-inc (True, in FV && INGRID Stop): PM ";
      if(i==169) title = "CC-inc (True, in FV && INGRID Stop): WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=172){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==170) title = "CC-inc (True, in FV && INGRID Stop): ING ";
      if(i==171) title = "CC-inc (True, in FV && INGRID Stop): PM ";
      if(i==172) title = "CC-inc (True, in FV && INGRID Stop): WM ";
      xtitle = "Muon Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=175){
      nbin   =    180; 
      min    =     0.;  //GeV/c
      max    =   180.;  //GeV/c
      if(i==173) title = "CC-inc (True, in FV && INGRID Stop): ING ";
      if(i==174) title = "CC-inc (True, in FV && INGRID Stop): PM ";
      if(i==175) title = "CC-inc (True, in FV && INGRID Stop): WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=178){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==176) title = "CC-inc (True, in FV && 3 INGRID Hits): ING ";
      if(i==177) title = "CC-inc (True, in FV && 3 INGRID Hits): PM ";
      if(i==178) title = "CC-inc (True, in FV && 3 INGRID Hits): WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=181){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==179) title = "CC-inc (True, in FV && 3 INGRID Hits): ING ";
      if(i==180) title = "CC-inc (True, in FV && 3 INGRID Hits): PM ";
      if(i==181) title = "CC-inc (True, in FV && 3 INGRID Hits): WM ";
      xtitle = "Muon Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=184){
      nbin   =    180; 
      min    =     0.;  //GeV/c
      max    =   180.;  //GeV/c
      if(i==182) title = "CC-inc (True, in FV && 3 INGRID Hits): ING ";
      if(i==183) title = "CC-inc (True, in FV && 3 INGRID Hits): PM ";
      if(i==184) title = "CC-inc (True, in FV && 3 INGRID Hits): WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=187){
      nbin   =  12; 
      min    =  -1;  //pln
      max    =  11;  //pln
      if(i==185) title = "Track Stop Position in INGRID (After Side Escaping Cut): ING ";
      if(i==186) title = "Track Stop Position in INGRID (After Side Escaping Cut): PM ";
      if(i==187) title = "Track Stop Position in INGRID (After Side Escaping Cut): WM ";
      xtitle = "Stop Position Z [pln]";
      ytitle = "";
    }
    else if(i<=190){
      nbin   = 25; 
      min    = 0;  //pln
      max    = 25;  //pln
      if(i==188) title = "Reconstructed Vertex (After FV Cut) : ING ";
      if(i==189) title = "Reconstructed Vertex (After FV Cut) : PM ";
      if(i==190) title = "Reconstructed Vertex (After FV Cut) : WM ";
      xtitle = "Vertex Position Z [plane]";
      ytitle = "";
    }
    else if(i<=193){
      nbin   = 1500; 
      min    = -750.;  //mm
      max    =  750.;  //mm
      if(i==191) title = "Reconstructed Vertex (After FV Cut) : ING ";
      if(i==192) title = "Reconstructed Vertex (After FV Cut) : PM ";
      if(i==193) title = "Reconstructed Vertex (After FV Cut) : WM ";
      xtitle = "Vertex Position X [mm]";
      ytitle = "";
    }
    else if(i<=196){
      nbin   = 1500; 
      min    = -750.;  //mm
      max    =  750.;  //mm
      if(i==194) title = "Reconstructed Vertex (After FV Cut) : ING ";
      if(i==195) title = "Reconstructed Vertex (After FV Cut) : PM ";
      if(i==196) title = "Reconstructed Vertex (After FV Cut) : WM ";
      xtitle = "Vertex Position Y [mm]";
      ytitle = "";
    }
    else if(i<=199){
      nbin   =  10; 
      min    =   0.;  //num
      max    =  10.;  //num 
      if(i==197) title = "Number of Tracks (After FV Cut) : ING ";
      if(i==198) title = "Number of Tracks (After FV Cut) : PM ";
      if(i==199) title = "Number of Tracks (After FV Cut) : WM ";
      xtitle = "Number of Tracks";
      ytitle = "";
    }
    else if(i<=202){
      nbin   =  10; 
      min    =   0.;  //num
      max    =  10.;  //num 
      if(i==200) title = "Number of Tracks (After Track Angle Cut) : ING ";
      if(i==201) title = "Number of Tracks (After Track Angle Cut) : PM ";
      if(i==202) title = "Number of Tracks (After Track Angle Cut) : WM ";
      xtitle = "Number of Tracks";
      ytitle = "";
    }
    else if(i<=205){
      nbin   =  10; 
      min    =   0.;  //num
      max    =  10.;  //num 
      if(i==203) title = "Number of Tracks (After ING Stop Cut) : ING ";
      if(i==204) title = "Number of Tracks (After ING Stop Cut) : PM ";
      if(i==205) title = "Number of Tracks (After ING Stop Cut) : WM ";
      xtitle = "Number of Tracks";
      ytitle = "";
    }
    else if(i<=208){
      nbin   =  90; 
      min    =   0.;  //deg
      max    =  90.;  //deg
      if(i==206) title = "1st Track Angle (After Track Angle & One Track Cut) : ING ";
      if(i==207) title = "1st Track Angle (After Track Angle & One Track Cut) : PM ";
      if(i==208) title = "1st Track Angle (After Track Angle & One Track Cut) : WM ";
      xtitle = "Reconstructed Angle [deg]";
      ytitle = "";
    }
    else if(i<=211){
      nbin   =   5; 
      min    =   0.;
      max    =   5.; //0:muon, 1:e/gamma, 2:p/n, 3:pi, 4:others
      if(i==209) title = "1st Track Hit PID (After FV Cut) : ING ";
      if(i==210) title = "1st Track Hit PID (After FV Cut) : PM ";
      if(i==211) title = "1st Track Hit PID (After FV Cut) : WM ";
      xtitle = "Particle ID";
      ytitle = "";
    }
    else if(i<=214){
      nbin   =   5; 
      min    =   0.;
      max    =   5.; //0:muon, 1:e/gamma, 2:p/n, 3:pi, 4:others
      if(i==212) title = "1st Track PID (After FV Cut) : ING ";
      if(i==213) title = "1st Track PID (After FV Cut) : PM ";
      if(i==214) title = "1st Track PID (After FV Cut) : WM ";
      xtitle = "Particle ID";
      ytitle = "";
    }
    else if(i<=217){
      nbin   =   5; 
      min    =   0.;
      max    =   5.; //0:muon, 1:e/gamma, 2:p/n, 3:pi, 4:others
      if(i==215) title = "1st Track PID (After Track Angle & One Track Cut) : ING ";
      if(i==216) title = "1st Track PID (After Track Angle & One Track Cut) : PM ";
      if(i==217) title = "1st Track PID (After Track Angle & One Track Cut) : WM ";
      xtitle = "Particle ID";
      ytitle = "";
    }
    else if(i<=220){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==218) title = "Selected CC-inc (After Track Angle & One Track Cut): ING ";
      if(i==219) title = "Selected CC-inc (After Track Angle & One Track Cut): PM ";
      if(i==220) title = "Selected CC-inc (After Track Angle & One Track Cut): WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=223){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==221) title = "Selected CC-inc (After Track Angle & One Track Cut): ING ";
      if(i==222) title = "Selected CC-inc (After Track Angle & One Track Cut): PM ";
      if(i==223) title = "Selected CC-inc (After Track Angle & One Track Cut): WM ";
      xtitle = "Muon Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=226){
      nbin   =    180; 
      min    =     0.;  //GeV/c
      max    =   180.;  //GeV/c
      if(i==224) title = "Selected CC-inc (After Track Angle & One Track Cut): ING ";
      if(i==225) title = "Selected CC-inc (After Track Angle & One Track Cut): PM ";
      if(i==226) title = "Selected CC-inc (After Track Angle & One Track Cut): WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=229){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==227) title = "CC-inclusive Interaction in FV, >400MeV/c: ING ";
      if(i==228) title = "CC-inclusive Interaction in FV, >400MeV/c: PM ";
      if(i==229) title = "CC-inclusive Interaction in FV, >400MeV/c: WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=232){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==230) title = "CC-inclusive Interaction in FV, >400MeV/c : ING ";
      if(i==231) title = "CC-inclusive Interaction in FV, >400MeV/c : PM ";
      if(i==232) title = "CC-inclusive Interaction in FV, >400MeV/c : WM ";
      xtitle = "Muon Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=235){
      nbin   =    180; 
      min    =     0.;  //GeV/c
      max    =   180.;  //GeV/c
      if(i==233) title = "CC-inclusive Interaction in FV, >400MeV/c : ING ";
      if(i==234) title = "CC-inclusive Interaction in FV, >400MeV/c : PM ";
      if(i==235) title = "CC-inclusive Interaction in FV, >400MeV/c : WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=238){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==236) title = "CC-inclusive Interaction in FV, <45deg: ING ";
      if(i==237) title = "CC-inclusive Interaction in FV, <45deg: PM ";
      if(i==238) title = "CC-inclusive Interaction in FV, <45deg: WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=241){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==239) title = "CC-inclusive Interaction in FV, <45deg : ING ";
      if(i==240) title = "CC-inclusive Interaction in FV, <45deg : PM ";
      if(i==241) title = "CC-inclusive Interaction in FV, <45deg : WM ";
      xtitle = "Muon Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=244){
      nbin   =    180; 
      min    =     0.;  //GeV/c
      max    =   180.;  //GeV/c
      if(i==242) title = "CC-inclusive Interaction in FV, <45deg : ING ";
      if(i==243) title = "CC-inclusive Interaction in FV, <45deg : PM ";
      if(i==244) title = "CC-inclusive Interaction in FV, <45deg : WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=247){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==245) title = "CC-inclusive Interaction in FV, <30deg: ING ";
      if(i==246) title = "CC-inclusive Interaction in FV, <30deg: PM ";
      if(i==247) title = "CC-inclusive Interaction in FV, <30deg: WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=250){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==248) title = "CC-inclusive Interaction in FV, <30deg : ING ";
      if(i==249) title = "CC-inclusive Interaction in FV, <30deg : PM ";
      if(i==250) title = "CC-inclusive Interaction in FV, <30deg : WM ";
      xtitle = "Muon Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=253){
      nbin   =    180; 
      min    =     0.;  //GeV/c
      max    =   180.;  //GeV/c
      if(i==251) title = "CC-inclusive Interaction in FV, <30deg : ING ";
      if(i==252) title = "CC-inclusive Interaction in FV, <30deg : PM ";
      if(i==253) title = "CC-inclusive Interaction in FV, <30deg : WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=256){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==254) title = "Selected CC-inc (After Acceptance Cut) : ING ";
      if(i==255) title = "Selected CC-inc (After Acceptance Cut) : PM ";
      if(i==256) title = "Selected CC-inc (After Acceptance Cut) : WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=259){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==257) title = "Selected CC-inc (After Acceptance Cut) : ING ";
      if(i==258) title = "Selected CC-inc (After Acceptance Cut) : PM ";
      if(i==259) title = "Selected CC-inc (After Acceptance Cut) : WM ";
      xtitle = "Muon Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=262){
      nbin   =    180; 
      min    =     0.;  //GeV/c
      max    =   180.;  //GeV/c
      if(i==260) title = "Selected CC-inc (After Acceptance Cut) : ING ";
      if(i==261) title = "Selected CC-inc (After Acceptance Cut) : PM ";
      if(i==262) title = "Selected CC-inc (After Acceptance Cut) : WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=265){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==263) title = "Selected CC-inc (After One Track Cut) : ING ";
      if(i==264) title = "Selected CC-inc (After One Track Cut) : PM ";
      if(i==265) title = "Selected CC-inc (After One Track Cut) : WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=268){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==266) title = "Selected CC-inc (After One Track Cut) : ING ";
      if(i==267) title = "Selected CC-inc (After One Track Cut) : PM ";
      if(i==268) title = "Selected CC-inc (After One Track Cut) : WM ";
      xtitle = "Muon Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=271){
      nbin   =    180; 
      min    =     0.;  //GeV/c
      max    =   180.;  //GeV/c
      if(i==269) title = "Selected CC-inc (After One Track Cut) : ING ";
      if(i==270) title = "Selected CC-inc (After One Track Cut) : PM ";
      if(i==271) title = "Selected CC-inc (After One Track Cut) : WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=274){
      nbin   =  90; 
      min    =   0.;  //deg
      max    =  90.;  //deg
      if(i==272) title = "1st Track Angle (After Acceptance Cut) : ING ";
      if(i==273) title = "1st Track Angle (After Acceptance Cut) : PM ";
      if(i==274) title = "1st Track Angle (After Acceptance Cut) : WM ";
      xtitle = "Reconstructed Angle [deg]";
      ytitle = "";
    }
    else if(i<=277){
      nbin   =  90; 
      min    =   0.;  //deg
      max    =  90.;  //deg
      if(i==275) title = "1st Track Angle (After One Track & Acceptance Cut) : ING ";
      if(i==276) title = "1st Track Angle (After One Track & Acceptance Cut) : PM ";
      if(i==277) title = "1st Track Angle (After One Track & Acceptance Cut) : WM ";
      xtitle = "Reconstructed Angle [deg]";
      ytitle = "";
    }
    else if(i<=280){
      nbin   =  90; 
      min    =   0.;  //deg
      max    =  90.;  //deg
      if(i==278) title = "1st Track Angle (After Acceptance Cut & OutsidePhaseSpace) : ING ";
      if(i==279) title = "1st Track Angle (After Acceptance Cut & OutsidePhaseSpace) : PM ";
      if(i==280) title = "1st Track Angle (After Acceptance Cut & OutsidePhaseSpace) : WM ";
      xtitle = "Reconstructed Angle [deg]";
      ytitle = "";
    }
    else if(i<=283){
      nbin   =  90; 
      min    =   0.;  //deg
      max    =  90.;  //deg
      if(i==281) title = "1st Track Angle (After One Track Cut & OutsidePhaseSpace) : ING ";
      if(i==282) title = "1st Track Angle (After One Track Cut & OutsidePhaseSpace) : PM ";
      if(i==283) title = "1st Track Angle (After One Track Cut & OutsidePhaseSpace) : WM ";
      xtitle = "Reconstructed Angle [deg]";
      ytitle = "";
    }
    else if(i<=286){
      nbin   =  10; 
      min    =   0.;  //num
      max    =  10.;  //num 
      if(i==284) title = "Number of Tracks (After Acceptance Cut) : ING ";
      if(i==285) title = "Number of Tracks (After Acceptance Cut) : PM ";
      if(i==286) title = "Number of Tracks (After Acceptance Cut) : WM ";
      xtitle = "Number of Tracks";
      ytitle = "";
    }
    else if(i<=289){
      nbin   =    180; 
      min    =     0.;  //GeV/c
      max    =   180.;  //GeV/c
      if(i==287) title = "CC-events in Phase Space : ING ";
      if(i==288) title = "CC-events in Phase Space : PM ";
      if(i==289) title = "CC-events in Phase Space : WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=292){
      nbin   =   5; 
      min    =   0.;
      max    =   5.; //0:muon, 1:e/gamma, 2:p/n, 3:pi, 4:others
      if(i==290) title = "1st Track PID (After Acceptance Cut) : ING ";
      if(i==291) title = "1st Track PID (After Acceptance Cut) : PM ";
      if(i==292) title = "1st Track PID (After Acceptance Cut) : WM ";
      xtitle = "Particle ID";
      ytitle = "";
    }
    else if(i<=295){
      nbin   =   5; 
      min    =   0.;
      max    =   5.; //0:muon, 1:e/gamma, 2:p/n, 3:pi, 4:others
      if(i==293) title = "1st Track PID (After One Track Cut) : ING ";
      if(i==294) title = "1st Track PID (After One Track Cut) : PM ";
      if(i==295) title = "1st Track PID (After One Track Cut) : WM ";
      xtitle = "Particle ID";
      ytitle = "";
    }
    else if(i<=298){
      nbin   =  2000; 
      min    =    0.;
      max    =   20.;
      if(i==296) title = "Cross Section for All CC Events : ING ";
      if(i==297) title = "Cross Section for All CC Events : PM ";
      if(i==298) title = "Cross Section for All CC Events : WM ";
      xtitle = "Total Cross Section [x1e-38 cm^2]";
      ytitle = "";
    }
    else if(i<=301){
      nbin   =  2000; 
      min    =    0.;
      max    =   20.;
      if(i==299) title = "Cross Section in Phase Space : ING ";
      if(i==300) title = "Cross Section in Phase Space : PM ";
      if(i==301) title = "Cross Section in Phase Space : WM ";
      xtitle = "Total Cross Section [x1e-38 cm^2]";
      ytitle = "";
    }
    else if(i<=304){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==302) title = "Nu Flux in FV : ING ";
      if(i==303) title = "Nu Flux in FV : PM ";
      if(i==304) title = "Nu Flux in FV : WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=307){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==305) title = "Nu Flux in Phase Space : ING ";
      if(i==306) title = "Nu Flux in Phase Space : PM ";
      if(i==307) title = "Nu Flux in Phase Space : WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=310){
      nbin   =    180; 
      min    =     0.;  //deg
      max    =   180.;  //deg
      if(i==308) title = "CC-events in Phase Space & On Target Material: ING ";
      if(i==309) title = "CC-events in Phase Space & On Target Material: PM ";
      if(i==310) title = "CC-events in Phase Space & On Target Material: WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=313){
      nbin   =    180; 
      min    =     0.;  //deg
      max    =   180.;  //deg
      if(i==311) title = "CC-events in Phase Space & Not On Target Material: ING ";
      if(i==312) title = "CC-events in Phase Space & Not On Target Material: PM ";
      if(i==313) title = "CC-events in Phase Space & Not On Target Material: WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=316){
      nbin   =   180; 
      min    =     0.;  //deg
      max    =   180.;  //deg
      if(i==314) title = "CC-events & On Target Material: ING ";
      if(i==315) title = "CC-events & On Target Material: PM ";
      if(i==316) title = "CC-events & On Target Material: WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=319){
      nbin   =  90; 
      min    =   0.;  //deg
      max    =  90.;  //deg
      if(i==317) title = "1st Track Angle (After Acceptance Cut & OutsideFV) : ING ";
      if(i==318) title = "1st Track Angle (After Acceptance Cut & OutsideFV) : PM ";
      if(i==319) title = "1st Track Angle (After Acceptance Cut & OutsideFV) : WM ";
      xtitle = "Reconstructed Angle [deg]";
      ytitle = "";
    }
    else if(i<=322){
      nbin   =  90; 
      min    =   0.;  //deg
      max    =  90.;  //deg
      if(i==320) title = "1st Track Angle (After OneTrk Cut & OutsideFV) : ING ";
      if(i==321) title = "1st Track Angle (After OneTrk Cut & OutsideFV) : PM ";
      if(i==322) title = "1st Track Angle (After OneTrk Cut & OutsideFV) : WM ";
      xtitle = "Reconstructed Angle [deg]";
      ytitle = "";
    }
    else if(i<=325){
      nbin   =  1000; 
      min    =    0.;  //pe
      max    =  100.;  //pe
      if(i==323) title = "Light Yield of Hits on 1st Track (After OneTrk Cut) : ING ";
      if(i==324) title = "Light Yield of Hits on 1st Track (After OneTrk Cut) : PM ";
      if(i==325) title = "Light Yield of Hits on 1st Track (After OneTrk Cut) : WM ";
      xtitle = "Light Yield [pe]";
      ytitle = "";
    }
    else if(i<=328){
      nbin   =  1000; 
      min    =    0.;  //pe
      max    =  100.;  //pe
      if(i==326) title = "Light Yield of Muon Hits : ING ";
      if(i==327) title = "Light Yield of Muon Hits : PM ";
      if(i==328) title = "Light Yield of Muon Hits : WM ";
      xtitle = "Light Yield [pe]";
      ytitle = "";
    }
    else if(i<=331){
      nbin   =  1000; 
      min    =    0.;  //pe
      max    =  100.;  //pe
      if(i==329) title = "Light Yield of Other Particle Hits : ING ";
      if(i==330) title = "Light Yield of Other Particle Hits : PM ";
      if(i==331) title = "Light Yield of Other Particle Hits : WM ";
      xtitle = "Light Yield [pe]";
      ytitle = "";
    }
    else if(i<=334){
      nbin   =   180; 
      min    =     0.;  //deg
      max    =   180.;  //deg
      if(i==332) title = "All Interaction in FV: ING ";
      if(i==333) title = "All Interaction in FV: PM ";
      if(i==334) title = "All Interaction in FV: WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=337){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==335) title = "Selected CC-inc (After Acceptance Cut from Phase Space) : ING ";
      if(i==336) title = "Selected CC-inc (After Acceptance Cut from Phase Space) : PM ";
      if(i==337) title = "Selected CC-inc (After Acceptance Cut from Phase Space) : WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=340){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==338) title = "Selected CC-inc (After Acceptance Cut from Phase Space) : ING ";
      if(i==339) title = "Selected CC-inc (After Acceptance Cut from Phase Space) : PM ";
      if(i==340) title = "Selected CC-inc (After Acceptance Cut from Phase Space) : WM ";
      xtitle = "Muon Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=343){
      nbin   =    180; 
      min    =     0.;  //GeV/c
      max    =   180.;  //GeV/c
      if(i==341) title = "Selected CC-inc (After Acceptance Cut from Phase Space) : ING ";
      if(i==342) title = "Selected CC-inc (After Acceptance Cut from Phase Space) : PM ";
      if(i==343) title = "Selected CC-inc (After Acceptance Cut from Phase Space) : WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=346){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==344) title = "Selected CC-inc (After OneTrk Cut from Phase Space) : ING ";
      if(i==345) title = "Selected CC-inc (After OneTrk Cut from Phase Space) : PM ";
      if(i==346) title = "Selected CC-inc (After OneTrk Cut from Phase Space) : WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=349){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==347) title = "Selected CC-inc (After OneTrk Cut from Phase Space) : ING ";
      if(i==348) title = "Selected CC-inc (After OneTrk Cut from Phase Space) : PM ";
      if(i==349) title = "Selected CC-inc (After OneTrk Cut from Phase Space) : WM ";
      xtitle = "Muon Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=352){
      nbin   =    180; 
      min    =     0.;  //deg
      max    =   180.;  //deg
      if(i==350) title = "Selected CC-inc (After OneTrk Cut from Phase Space) : ING ";
      if(i==351) title = "Selected CC-inc (After OneTrk Cut from Phase Space) : PM ";
      if(i==352) title = "Selected CC-inc (After OneTrk Cut from Phase Space) : WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=355){
      nbin   =    90; 
      min    =     0.;  //deg
      max    =    90.;  //deg
      if(i==353) title = "Selected Non-CC (After Acceptance Cut) : ING ";
      if(i==354) title = "Selected Non-CC (After Acceptance Cut) : PM ";
      if(i==355) title = "Selected Non-CC (After Acceptance Cut) : WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=358){
      nbin   =    90; 
      min    =     0.;  //deg
      max    =    90.;  //deg
      if(i==356) title = "Selected Non-CC (After OneTrk Cut) : ING ";
      if(i==357) title = "Selected Non-CC (After OneTrk Cut) : PM ";
      if(i==358) title = "Selected Non-CC (After OneTrk Cut) : WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=361){
      nbin   =    180; 
      min    =     0.;  //deg
      max    =   180.;  //deg
      if(i==359) title = "CC-events in FV from PS on CH (After Acceptance Cut): ING ";
      if(i==360) title = "CC-events in FV from PS on CH (After Acceptance Cut): PM ";
      if(i==361) title = "CC-events in FV from PS on CH (After Acceptance Cut): WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=364){
      nbin   =    180; 
      min    =     0.;  //deg
      max    =   180.;  //deg
      if(i==362) title = "CC-events in FV from PS on CH (After OneTrk Cut): ING ";
      if(i==363) title = "CC-events in FV from PS on CH (After OneTrk Cut): PM ";
      if(i==364) title = "CC-events in FV from PS on CH (After OneTrk Cut): WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=367){
      nbin   =    7; 
      min    =     0.;  //Interaction Type
      max    =     7.;  //Interaction Type
      if(i==365) title = "Interaction Type : CC-events in FV from PS : ING ";
      if(i==366) title = "Interaction Type : CC-events in FV from PS : PM ";
      if(i==367) title = "Interaction Type : CC-events in FV from PS : WM ";
      xtitle = "Interaction Type";
      ytitle = "";
    }
    else if(i<=370){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==368) title = "Pion Distribution : CC-events in FV from PS : ING ";
      if(i==369) title = "Pion Distribution : CC-events in FV from PS : PM ";
      if(i==370) title = "Pion Distribution : CC-events in FV from PS : WM ";
      xtitle = "Pion Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=373){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==371) title = "Proton Distribution (No Pion) : CC-events in FV from PS : ING ";
      if(i==372) title = "Proton Distribution (No Pion) : CC-events in FV from PS : PM ";
      if(i==373) title = "Proton Distribution (No Pion) : CC-events in FV from PS : WM ";
      xtitle = "Proton Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=376){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==374) title = "Proton Distribution : CC-events in FV from PS : ING ";
      if(i==375) title = "Proton Distribution : CC-events in FV from PS : PM ";
      if(i==376) title = "Proton Distribution : CC-events in FV from PS : WM ";
      xtitle = "Proton Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=379){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==377) title = "Pion Distribution : Selected as Multi-track events : ING ";
      if(i==378) title = "Pion Distribution : Selected as Multi-track events : PM ";
      if(i==379) title = "Pion Distribution : Selected as Multi-track events : WM ";
      xtitle = "Pion Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=382){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==380) title = "Proton Distribution (No Pion) : Selected as Multi-track events : ING ";
      if(i==381) title = "Proton Distribution (No Pion) : Selected as Multi-track events : PM ";
      if(i==382) title = "Proton Distribution (No Pion) : Selected as Multi-track events : WM ";
      xtitle = "Proton Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=385){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==383) title = "Proton Distribution : Selected as Multi-track events : ING ";
      if(i==384) title = "Proton Distribution : Selected as Multi-track events : PM ";
      if(i==385) title = "Proton Distribution : Selected as Multi-track events : WM ";
      xtitle = "Proton Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=388){
      nbin   =   8; 
      min    =   0.;  //bin
      max    =   8.;  //bin
      if(i==386) title = "All Interaction from Phase Space: ING ";
      if(i==387) title = "All Interaction from Phase Space: PM ";
      if(i==388) title = "All Interaction from Phase Space: WM ";
      xtitle = "Phase Space Bins";
      ytitle = "";
    }
    else if(i<=391){
      nbin   =   8; 
      min    =   0.;  //bin
      max    =   8.;  //bin
      if(i==389) title = "Selected Events from Phase Space: Acceptance Cut : ING ";
      if(i==390) title = "Selected Events from Phase Space: Acceptance Cut : PM ";
      if(i==391) title = "Selected Events from Phase Space: Acceptance Cut : WM ";
      xtitle = "Phase Space Bins";
      ytitle = "";
    }
    else if(i<=394){
      nbin   =   8; 
      min    =   0.;  //bin
      max    =   8.;  //bin
      if(i==392) title = "Selected Events from Phase Space: OneTrk Cut : ING ";
      if(i==393) title = "Selected Events from Phase Space: OneTrk Cut : PM ";
      if(i==394) title = "Selected Events from Phase Space: OneTrk Cut : WM ";
      xtitle = "Phase Space Bins";
      ytitle = "";
    }
    else if(i<=397){
      nbin   =   8; 
      min    =   0.;  //recon angle bin
      max    =   8.;  //recon angle bin
      if(i==395) title = "Number of Selected Events : Acceptance Cut : ING ";
      if(i==396) title = "Number of Selected Events : Acceptance Cut : PM ";
      if(i==397) title = "Number of Selected Events : Acceptance Cut : WM ";
      xtitle = "Reconstructed Angle Bins";
      ytitle = "";
    }
    else if(i<=400){
      nbin   =   8; 
      min    =   0.;  //recon angle bin
      max    =   8.;  //recon angle bin
      if(i==398) title = "Number of Selected Events : OneTrk Cut : ING ";
      if(i==399) title = "Number of Selected Events : OneTrk Cut : PM ";
      if(i==400) title = "Number of Selected Events : OneTrk Cut : WM ";
      xtitle = "Reconstructed Angle Bins";
      ytitle = "";
    }
    else if(i<=403){
      nbin   =   8; 
      min    =   0.;  //recon angle bin
      max    =   8.;  //recon angle bin
      if(i==401) title = "BG CC out of FV : Acceptance Cut : ING ";
      if(i==402) title = "BG CC out of FV : Acceptance Cut : PM ";
      if(i==403) title = "BG CC out of FV : Acceptance Cut : WM ";
      xtitle = "Reconstructed Angle Bins";
      ytitle = "";
    }
    else if(i<=406){
      nbin   =   8; 
      min    =   0.;  //recon angle bin
      max    =   8.;  //recon angle bin
      if(i==404) title = "BG CC out of FV : OneTrk Cut : ING ";
      if(i==405) title = "BG CC out of FV : OneTrk Cut : PM ";
      if(i==406) title = "BG CC out of FV : OneTrk Cut : WM ";
      xtitle = "Reconstructed Angle Bins";
      ytitle = "";
    }
    else if(i<=409){
      nbin   =   8; 
      min    =   0.;  //recon angle bin
      max    =   8.;  //recon angle bin
      if(i==407) title = "BG NC : Acceptance Cut : ING ";
      if(i==408) title = "BG NC : Acceptance Cut : PM ";
      if(i==409) title = "BG NC : Acceptance Cut : WM ";
      xtitle = "Reconstructed Angle Bins";
      ytitle = "";
    }
    else if(i<=412){
      nbin   =   8; 
      min    =   0.;  //recon angle bin
      max    =   8.;  //recon angle bin
      if(i==410) title = "BG NC : OneTrk Cut : ING ";
      if(i==411) title = "BG NC : OneTrk Cut : PM ";
      if(i==412) title = "BG NC : OneTrk Cut : WM ";
      xtitle = "Reconstructed Angle Bins";
      ytitle = "";
    }
    else if(i<=415){
      nbin   =   8; 
      min    =   0.;  //bin
      max    =   8.;  //bin
      if(i==413) title = "All Interaction on CH : ING ";
      if(i==414) title = "All Interaction on CH : PM ";
      if(i==415) title = "All Interaction on CH : WM ";
      xtitle = "Phase Space Bins";
      ytitle = "";
    }
    else if(i<=418){
      nbin   =   8; 
      min    =   0.;  //bin
      max    =   8.;  //bin
      if(i==416) title = "Selected Events on CH : Acceptance Cut : ING ";
      if(i==417) title = "Selected Events on CH : Acceptance Cut : PM ";
      if(i==418) title = "Selected Events on CH : Acceptance Cut : WM ";
      xtitle = "Phase Space Bins";
      ytitle = "";
    }
    else if(i<=421){
      nbin   =   8; 
      min    =   0.;  //bin
      max    =   8.;  //bin
      if(i==419) title = "Selected Events on CH : OneTrk Cut : ING ";
      if(i==420) title = "Selected Events on CH : OneTrk Cut : PM ";
      if(i==421) title = "Selected Events on CH : OneTrk Cut : WM ";
      xtitle = "Phase Space Bins";
      ytitle = "";
    }
    else if(i<=427){
      nbin   =    180; 
      min    =     0.;  //deg
      max    =   180.;  //deg
      if(i==422) title = "Reconstructed Angle (SideView): FV Cut : ING ";
      if(i==423) title = "Reconstructed Angle (SideView): FV Cut : PM ";
      if(i==424) title = "Reconstructed Angle (SideView): FV Cut : WM ";
      if(i==425) title = "Reconstructed Angle (TopView): FV Cut : ING ";
      if(i==426) title = "Reconstructed Angle (TopView): FV Cut : PM ";
      if(i==427) title = "Reconstructed Angle (TopView): FV Cut : WM ";
      if(i<=424) xtitle = "Reconstructed Angle X [deg]";
      if(i>=425) xtitle = "Reconstructed Angle Y [deg]";
      ytitle = "";
    }
    else if(i<=430){
      nbin   =  10; 
      min    =   0.;  //num
      max    =  10.;  //num 
      if(i==428) title = "Number of Tracks (After Pre-Selection) : ING ";
      if(i==429) title = "Number of Tracks (After Pre-Selection) : PM ";
      if(i==430) title = "Number of Tracks (After Pre-Selection) : WM ";
      xtitle = "Number of Tracks";
      ytitle = "";
    }
    else if(i<=433){
      nbin   =    180; 
      min    =     0.;  //deg
      max    =   180.;  //deg
      if(i==431) title = "Reconstructed Angle (One-Track): ING ";
      if(i==432) title = "Reconstructed Angle (One-Track): PM ";
      if(i==433) title = "Reconstructed Angle (One-Track): WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=436){
      nbin   = 25; 
      min    = 0;  //pln
      max    = 25;  //pln
      if(i==434) title = "Reconstructed Vertex Z (One-Track: 0-30deg) : ING ";
      if(i==435) title = "Reconstructed Vertex Z (One-Track: 0-30deg) : PM ";
      if(i==436) title = "Reconstructed Vertex Z (One-Track: 0-30deg) : WM ";
      xtitle = "Vertex Position Z [plane]";
      ytitle = "";
    }
    else if(i<=439){
      nbin   = 1500; 
      min    = -750.;  //mm
      max    =  750.;  //mm
      if(i==437) title = "Reconstructed Vertex X (One-Track: 0-30deg) : ING ";
      if(i==438) title = "Reconstructed Vertex X (One-Track: 0-30deg) : PM ";
      if(i==439) title = "Reconstructed Vertex X (One-Track: 0-30deg) : WM ";
      xtitle = "Vertex Position X [mm]";
      ytitle = "";
    }
    else if(i<=442){
      nbin   = 1500; 
      min    = -750.;  //mm
      max    =  750.;  //mm
      if(i==440) title = "Reconstructed Vertex Y (One-Track: 0-30deg) : ING ";
      if(i==441) title = "Reconstructed Vertex Y (One-Track: 0-30deg) : PM ";
      if(i==442) title = "Reconstructed Vertex Y (One-Track: 0-30deg) : WM ";
      xtitle = "Vertex Position Y [mm]";
      ytitle = "";
    }
    else if(i<=445){
      nbin   = 1000; 
      min    =   0.;  //GeV
      max    =  10.;  //GeV
      if(i==443) title = "Neutrino Energy (One-Track: 0-30deg) : ING ";
      if(i==444) title = "Neutrino Energy (One-Track: 0-30deg) : PM ";
      if(i==445) title = "Neutrino Energy (One-Track: 0-30deg) : WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "";
    }
    else if(i<=448){
      nbin   = 1000; 
      min    =   0.;  //GeV/c
      max    =  10.;  //GeV/c
      if(i==446) title = "Muon Momentum (One-Track: 0-30deg) : ING ";
      if(i==447) title = "Muon Momentum (One-Track: 0-30deg) : PM ";
      if(i==448) title = "Muon Momentum (One-Track: 0-30deg) : WM ";
      xtitle = "Muon Momentum [GeV/c]";
      ytitle = "";
    }
    else if(i<=451){
      nbin   =    180; 
      min    =     0.;  //deg
      max    =   180.;  //deg
      if(i==449) title = "Muon Angle (One-Track: 0-30deg) : ING ";
      if(i==450) title = "Muon Angle (One-Track: 0-30deg) : PM ";
      if(i==451) title = "Muon Angle (One-Track: 0-30deg) : WM ";
      xtitle = "Muon Angle [deg]";
      ytitle = "";
    }
    else if(i<=454){
      nbin   =  12; 
      min    =  -1;  
      max    =  11;  
      if(i==452) title = "Penetrated Iron Plate (One-Track: 0-30deg) : ING ";
      if(i==453) title = "Penetrated Iron Plate (One-Track: 0-30deg) : PM ";
      if(i==454) title = "Penetrated Iron Plate (One-Track: 0-30deg) : WM ";
      xtitle = "Penetrated Iron Plane";
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
    if(i<=5){
      nbin1   = 1500; 
      min1    = -750.; //mm
      max1    =  750.; //mm //along Z axis
      nbin2   = 1500; 
      min2    = -750.; //mm
      max2    =  750.; //mm //along XY axis
      if(i==0) title = "Vertex : ING : Side";
      if(i==1) title = "Vertex : ING : Top";
      if(i==2) title = "Vertex : PM : Side";
      if(i==3) title = "Vertex : PM : Top";
      if(i==4) title = "Vertex : WM : Side";
      if(i==5) title = "Vertex : WM : Top";
      xtitle = "Vertex Z [mm]";
      if(i%2==0) ytitle = "Vertex Y [mm]";
      if(i%2==1) ytitle = "Vertex X [mm]";
    }
    else if(i<=29){
      if(i>=6+NumSelec2*3){
        cout << "ID exceeds the selection number" << endl;
        cout << " i=" << i << " NumSelec2=" << NumSelec2 << endl;
        exit(1);
      }
      nbin1   =  500; 
      min1    =   0.; // GeV/c
      max1    =   5.; // GeV/c
      nbin2   =  180; 
      min2    =   0.; // deg
      max2    = 180.; // deg 
      title   = "Selected Events";
      if(i== 6) title += " : ING : All CC-inc in FV";
      if(i== 7) title += " : ING : After FV Cut";
      if(i== 8) title += " : ING : After Track Angle Cut";
      if(i== 9) title += " : ING : After Side Escaping Cut";
      if(i==10) title += " : ING : After INGRID Stop Cut";
      if(i==11) title += " : ING : 3 INGRID Hits (True)";
      if(i==12) title += " : ING : Stop in INGRID (True)";
      if(i==13) title += " : ING : After TrkAng and OneTrk Cut";
      if(i==14) title += " : PM : All CC-inc in FV";
      if(i==15) title += " : PM : After FV Cut";
      if(i==16) title += " : PM : After Track Angle Cut";
      if(i==17) title += " : PM : After Side Escaping Cut";
      if(i==18) title += " : PM : After INGRID Stop Cut";
      if(i==19) title += " : PM : 3 INGRID Hits (True)";
      if(i==20) title += " : PM : Stop in INGRID (True)";
      if(i==21) title += " : PM : After TrkAng and OneTrk Cut";
      if(i==22) title += " : WM : All CC-inc in FV";
      if(i==23) title += " : WM : After FV Cut";
      if(i==24) title += " : WM : After Track Angle Cut";
      if(i==25) title += " : WM : After Side Escaping Cut";
      if(i==26) title += " : WM : After INGRID Stop Cut";
      if(i==27) title += " : WM : 3 INGRID Hits (True)";
      if(i==28) title += " : WM : Stop in INGRID (True)";
      if(i==29) title += " : WM : After TrkAng and OneTrk Cut";
      xtitle = "True Muon Momentum [GeV/c]";
      ytitle = "True Muon Angle [deg]";
    }
    else if(i<=41){
      nbin1   =  2000; 
      min1    = -1500.; //cm
      max1    =   500.; //cm //along Z axis
      nbin2   =  2000; 
      min2    = -1000.; //cm
      max2    =  1000.; //cm //along XY axis
      if(i==30) title = "True Vertex (After Vertexing): ING : Side";
      if(i==31) title = "True Vertex (After Vertexing): ING : Top";
      if(i==32) title = "True Vertex (After Vertexing): PM : Side";
      if(i==33) title = "True Vertex (After Vertexing): PM : Top";
      if(i==34) title = "True Vertex (After Vertexing): WM : Side";
      if(i==35) title = "True Vertex (After Vertexing): WM : Top";
      if(i==36) title = "True Vertex (After FV Cut): ING : Side";
      if(i==37) title = "True Vertex (After FV Cut): ING : Top";
      if(i==38) title = "True Vertex (After FV Cut): PM : Side";
      if(i==39) title = "True Vertex (After FV Cut): PM : Top";
      if(i==40) title = "True Vertex (After FV Cut): WM : Side";
      if(i==41) title = "True Vertex (After FV Cut): WM : Top";
      xtitle = "Vertex Z [cm]";
      if(i%2==0) ytitle = "Vertex Y [cm]";
      if(i%2==1) ytitle = "Vertex X [cm]";
    }
    else if(i<=47){
      nbin1   =   45;
      min1    =    0.; //deg
      max1    =   45.; //deg
      nbin2   =   45; 
      min2    =    0.; //deg
      max2    =   45.; //deg
      if(i==42) title = "Resolution [Angle] (After FV Cut): ING";
      if(i==43) title = "Resolution [Angle] (After FV Cut): PM";
      if(i==44) title = "Resolution [Angle] (After FV Cut): WM";
      if(i==45) title = "Resolution [Angle] (After FV & OneTrk Cut): ING";
      if(i==46) title = "Resolution [Angle] (After FV & OneTrk Cut): PM";
      if(i==47) title = "Resolution [Angle] (After FV & OneTrk Cut): WM";
      xtitle = "True track angle [deg]";
      ytitle = "Reconstructed track angle [deg]";
    }
    else if(i<=53){
      nbin1   =  1500;
      min1    =  -750.; //mm
      max1    =   750.; //mm
      nbin2   =  1500; 
      min2    =  -750.; //mm
      max2    =   750.; //mm
      if(i==48) title = "Resolution [Vertex] (After UpstVeto Cut): X : ING";
      if(i==49) title = "Resolution [Vertex] (After UpstVeto Cut): X : PM" ;
      if(i==50) title = "Resolution [Vertex] (After UpstVeto Cut): X : WM" ;
      if(i==51) title = "Resolution [Vertex] (After UpstVeto Cut): Y : ING";
      if(i==52) title = "Resolution [Vertex] (After UpstVeto Cut): Y : PM" ;
      if(i==53) title = "Resolution [Vertex] (After UpstVeto Cut): Y : WM" ;
      xtitle = "True track vertex [mm]";
      ytitle = "Reconstructed Vertex [mm]";
    }
    else if(i<=54){
      nbin1   =   53; //18+24+11
      min1    =   0.; //pln
      max1    =  53.; //pln
      nbin2   =   53; 
      min2    =   0.; //pln
      max2    =  53.; //pln
      if(i==54) title = "Resolution [Vertex] (After Vertexing): Z";
      xtitle = "True track vertex [pln]";
      ytitle = "Reconstructed Vertex [pln]";
    }
    else if(i<=60){
      nbin1   =   45;
      min1    =    0.; //deg
      max1    =   45.; //deg
      nbin2   =   45; 
      min2    =    0.; //deg
      max2    =   45.; //deg
      if(i==55) title = "U Matrix [Angle] (After Acceptance Cut): ING";
      if(i==56) title = "U Matrix [Angle] (After Acceptance Cut): PM";
      if(i==57) title = "U Matrix [Angle] (After Acceptance Cut): WM";
      if(i==58) title = "U Matrix [Angle] (After One Track Cut): ING";
      if(i==59) title = "U Matrix [Angle] (After One Track Cut): PM";
      if(i==60) title = "U Matrix [Angle] (After One Track Cut): WM";
      xtitle = "True track angle [deg]";
      ytitle = "Reconstructed track angle [deg]";
    }
    else if(i<=66){
      nbin1   =  500; 
      min1    =   0.; // GeV/c
      max1    =   5.; // GeV/c
      nbin2   =  180; 
      min2    =   0.; // deg
      max2    = 180.; // deg 
      title   = "Selected Events";
      if(i==61) title += " : ING : After Acceptance Cut"; 
      if(i==62) title += " : PM  : After Acceptance Cut"; 
      if(i==63) title += " : WM  : After Acceptance Cut"; 
      if(i==64) title += " : ING : After One Track Cut"; 
      if(i==65) title += " : PM  : After One Track Cut"; 
      if(i==66) title += " : WM  : After One Track Cut"; 
      xtitle = "True Muon Momentum [GeV/c]";
      ytitle = "True Muon Angle [deg]";
    }
    else if(i<=69){
      nbin1   =  1000; 
      min1    =    0.;
      max1    =   10.;
      nbin2   =  2000; 
      min2    =    0.;
      max2    =   20.;
      title = "Cross Section vs Neturino Energy";
      if(i==67) title += " CC Interaction : ING ";
      if(i==68) title += " CC Interaction : PM ";
      if(i==69) title += " CC Interaction : WM ";
      xtitle = "Neutrino Energy [GeV]";
      ytitle = "Total Cross Section [x1e-38 cm^2]";
    }
    else if(i<=81){
      nbin1   =   45;
      min1    =    0.; //deg
      max1    =   45.; //deg
      nbin2   =   45; 
      min2    =    0.; //deg
      max2    =   45.; //deg
      if(i==70) title = "P Matrix [Angle] (After Acceptance Cut): ING";
      if(i==71) title = "P Matrix [Angle] (After Acceptance Cut): PM";
      if(i==72) title = "P Matrix [Angle] (After Acceptance Cut): WM";
      if(i==73) title = "P Matrix [Angle] (After One Track Cut): ING";
      if(i==74) title = "P Matrix [Angle] (After One Track Cut): PM";
      if(i==75) title = "P Matrix [Angle] (After One Track Cut): WM";
      if(i==76) title = "P Matrix on CH [Angle] (After Acceptance Cut): ING";
      if(i==77) title = "P Matrix on CH [Angle] (After Acceptance Cut): PM";
      if(i==78) title = "P Matrix on CH [Angle] (After Acceptance Cut): WM";
      if(i==79) title = "P Matrix on CH [Angle] (After One Track Cut): ING";
      if(i==80) title = "P Matrix on CH [Angle] (After One Track Cut): PM";
      if(i==81) title = "P Matrix on CH [Angle] (After One Track Cut): WM";
      xtitle = "Reconstructed track angle [deg]";
      ytitle = "True track angle [deg]";
    }
    else if(i<=84){
      nbin1   =    7;
      min1    =    0.; //deg
      max1    =    7.; //deg
      nbin2   =    4; 
      min2    =    0.; //deg
      max2    =    4.; //deg
      if(i==82) title = "Interaction Type (CC events in FV from PS): ING";
      if(i==83) title = "Interaction Type (CC events in FV from PS): PM";
      if(i==84) title = "Interaction Type (CC events in FV from PS): WM";
      xtitle = "Interaction Type";
      ytitle = "Interaction Type After FSI";
    }
    else if(i<=99){
      nbin1   =  500; 
      min1    =   0.; // GeV/c
      max1    =   5.; // GeV/c
      nbin2   =  180; 
      min2    =   0.; // deg
      max2    = 180.; // deg
      if     (i<=90){title   = "CC events in FV from PS";}
      else if(i<=99){title   = "Selected Multi-track Events";}
      if(i==85) title += " : Pion Distribution : ING "; 
      if(i==86) title += " : Pion Distribution : PM  "; 
      if(i==87) title += " : Pion Distribution : WM  "; 
      if(i==88) title += " : Proton Distribution (No Pion) : ING "; 
      if(i==89) title += " : Proton Distribution (No Pion) : PM  "; 
      if(i==90) title += " : Proton Distribution (No Pion) : WM  "; 
      if(i==94) title += " : Pion Distribution : ING "; 
      if(i==95) title += " : Pion Distribution : PM  "; 
      if(i==96) title += " : Pion Distribution : WM  "; 
      if(i==97) title += " : Proton Distribution (No Pion) : ING "; 
      if(i==98) title += " : Proton Distribution (No Pion) : PM  "; 
      if(i==99) title += " : Proton Distribution (No Pion) : WM  "; 
      if(i<=87){
        xtitle = "True Pion Momentum [GeV/c]";
        ytitle = "True Pion Angle [deg]";
      }else if(i<=90){
        xtitle = "True Proton Momentum [GeV/c]";
        ytitle = "True Proton Angle [deg]";
      }else if(i<=96){
        xtitle = "True Pion Momentum [GeV/c]";
        ytitle = "True Pion Angle [deg]";
      }else if(i<=99){
        xtitle = "True Proton Momentum [GeV/c]";
        ytitle = "True Proton Angle [deg]";
      }
    }
    else if(i<=105){
      nbin1   =   8; 
      min1    =   0.; // phase space bin 
      max1    =   8.; // phase space bin 
      nbin2   =   8; 
      min2    =   0.; // reconstructed angle bin 
      max2    =   8.; // reconstructed angle bin 
      title   = "U Matrix (8 bins x 8 bins)";
      if(i==100) title += " : Acceptance Cut : ING "; 
      if(i==101) title += " : Acceptance Cut : PM  "; 
      if(i==102) title += " : Acceptance Cut : WM  "; 
      if(i==103) title += " : One Track Cut : ING "; 
      if(i==104) title += " : One Track Cut : PM  "; 
      if(i==105) title += " : One Track Cut : WM  "; 
      xtitle = "True Phase Space Bins";
      ytitle = "Reconstructed Track Angle Bins";
    }
    else if(i<=111){
      nbin1   =   8; 
      min1    =   0.; // reconstructed angle bin 
      max1    =   8.; // reconstructed angle bin 
      nbin2   =   8; 
      min2    =   0.; // phase space bin 
      max2    =   8.; // phase space bin 
      title   = "P Matrix (8 bins x 8 bins)";
      if(i==106) title += " : Acceptance Cut : ING "; 
      if(i==107) title += " : Acceptance Cut : PM  "; 
      if(i==108) title += " : Acceptance Cut : WM  "; 
      if(i==109) title += " : One Track Cut : ING "; 
      if(i==110) title += " : One Track Cut : PM  "; 
      if(i==111) title += " : One Track Cut : WM  "; 
      xtitle = "Reconstructed Track Angle Bins";
      ytitle = "True Phase Space Bins";
    }
    else if(i<=114){
      nbin1   =   25; 
      min1    =    0.; // reconstructed vertex pln
      max1    =   25.; // reconstructed vertex pln  
      nbin2   =  180; 
      min2    =    0.; // reconstructed angle 
      max2    =  180.; // reconstructed angle 
      title   = "Recon. Angle vs VertexZ";
      if(i==112) title += " : FV Cut : ING "; 
      if(i==113) title += " : FV Cut : PM  "; 
      if(i==114) title += " : FV Cut : WM  "; 
      xtitle = "Reconstructed Vertex Z [pln]";
      ytitle = "Reconstricted Angle [deg]";
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


int GetTrueStartPln(int mod,double posz)
{
  if(mod<0){return -1;}
  int vtx_pln = -1;
  int maxpln = plnmax(mod);
  for(int pln=0;pln<maxpln;pln++){
    double sciwid = sciwidth(mod,1,pln);
    if(posz<zposi(mod,1,pln)+sciwid){
      vtx_pln = pln;
      break;
    }
  }
  return vtx_pln;
}


void fHitTimeWindow(vector<Hit> &hit,double cTdcRsr=50.,int hit_th=5)
{

  Int_t nhit = hit.size();
  int hitmod = hit[0].mod;

  for(int i=0; i<nhit; i++){
    int ncount = 0;
    for(int j=0;j<nhit; j++){
      if(i==j){continue;}
      if(fabs(hit[i].time-hit[j].time) < cTdcRsr){
        ncount++;
      }
    }
    if(ncount>hit_th){ hit[i].used = 1; }
  }
}


bool fTimeClustering(vector<Hit> &hit,double cTdcRsr=50)
{

  Int_t nhit = hit.size();
  int maxhit = 5; //Maximum for clustered hits
  int hitmod = hit[0].mod;

  if(nhit<=maxhit) return false;

  vector<Hit> tmp(nhit);
  for(int i=0;i<nhit;i++){ 
    tmp[i] = hit[i];
    tmp[i].id = i;
  }

  nhit = tmp.size();
  for(int i=0; i<nhit-maxhit; i++){
    if(fabs( tmp[i].time - tmp[i+maxhit].time) < cTdcRsr){
      long  basetime = 0.;    
      float highpe   = 0.;
      float sumpe    = 0.;

      //Extract largest pe event to be used for basetime
      for(int j=0; j<maxhit+1; j++){
        if(tmp[i+j].pe > highpe){
          basetime = tmp[i+j].time;
          highpe   = tmp[i+j].pe;
          sumpe    = sumpe + tmp[i+j].pe;
        }
      }

      if(sumpe<cPeCut) continue;

      vector<Hit>::iterator it;
      Int_t ncount=0;
      for(it=tmp.begin(); it!=tmp.end(); it++){
        if(fabs( basetime - it->time) < cTdcRsr){
          ncount++;
        }
      }
      if(ncount<=maxhit) continue;

      //Clstering
      double reftime;
      for(it=tmp.begin(); it!=tmp.end(); it++){
        if(fabs( basetime - it->time) < cTdcRsr){
          if(highpe<it->pe){
            highpe  = it->pe;
            reftime = it->time;
          }
        }
      }
      for(it=tmp.begin(); it!=tmp.end(); it++){
        if(fabs( basetime - it->time) < cTdcRsr){
          hit[it->id].id  = 1; //clustered 
          hit[it->id].tdc = reftime; //used as a reference time
          it = tmp.erase(it);
          it--;
        }
      }
    }
  }
}



void SetHist3(int i,
    int &nbin, double &min, double &max,
    string &title, string &xtitle, string &ytitle,
    string &name
    )
{
  name = Form("h3_%02d",i);
  if(i<=NumHist3){
    nbin   =  1; 
    min    =  0; 
    max    =  1;
    if(i== 0) title = "POT";
    if(i== 1) title = "Number of Spills";
    if(i== 2) title = "Number of Event: ING : 1. Vertexing";
    if(i== 3) title = "Number of Event: ING : 2. Iron Penetration";
    if(i== 4) title = "Number of Event: ING : 3. Beam Timing Cut";
    if(i== 5) title = "Number of Event: ING : 4. Upstream Veto Cut";
    if(i== 6) title = "Number of Event: ING : 5. Fiducial Volume Cut";
    if(i== 7) title = "Number of Event: ING : 6. Acceptance Cut";
    if(i== 8) title = "";
    if(i== 9) title = "";
    if(i==10) title = "";
    if(i==11) title = "Number of Event: ING : 7. ";
    if(i==12) title = "Number of Event: ING : 8. Track Angle Cut";
    if(i==13) title = "Number of Event: PM : 1. Vertexing";
    if(i==14) title = "Number of Event: PM : 2. Iron Penetration";
    if(i==15) title = "Number of Event: PM : 3. Beam Timing Cut";
    if(i==16) title = "Number of Event: PM : 4. Upstream Veto Cut";
    if(i==17) title = "Number of Event: PM : 5. Fiducial Volume Cut";
    if(i==18) title = "Number of Event: PM : 6. Acceptance Cut";
    if(i==19) title = "";
    if(i==20) title = "";
    if(i==21) title = "";
    if(i==22) title = "Number of Event: PM : 7. ";
    if(i==23) title = "Number of Event: PM : 8. Track Angle Cut";
    if(i==24) title = "Number of Event: WM : 1. Vertexing";
    if(i==25) title = "Number of Event: WM : 2. Iron Penetration";
    if(i==26) title = "Number of Event: WM : 3. Beam Timing Cut";
    if(i==27) title = "Number of Event: WM : 4. Upstream Veto Cut";
    if(i==28) title = "Number of Event: WM : 5. Fiducial Volume Cut";
    if(i==29) title = "Number of Event: WM : 6. Acceptance Cut";
    if(i==30) title = "";
    if(i==31) title = "";
    if(i==32) title = "";
    if(i==33) title = "Number of Event: WM : 7. ";
    if(i==34) title = "Number of Event: WM : 8. Track Angle Cut";
    if(i==35) title = "Neutrino Flux : ING";
    if(i==36) title = "Neutrino Flux : PM";
    if(i==37) title = "Neutrino Flux : WM";
    if(i==38) title = "Number of Event in Phase Space : ING";
    if(i==39) title = "Number of Event in Phase Space : PM";
    if(i==40) title = "Number of Event in Phase Space : WM";
    xtitle = "";
    ytitle = title;
  }
  else{
    cout << "Unexpected ID for 1D (data) histogram." << endl;
    cout << " i=" << i << endl;
    exit(1);
  }
}

#endif

