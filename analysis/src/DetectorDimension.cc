#ifndef _DetectorDemension_C
#define _DetectorDemension_C

//#define DESIGND_VALUE

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include "DetectorDimension.hh"
#include "Const.hh"

using namespace std;

//Apply measurement results of scintillator dimension
DetectorDimension::DetectorDimension(){

  //===============================
  //===   INGRID Water Module   ===
  //===============================
  f = TFile::Open(Form("%s/dimension/iwm/Scintillator_position.root",MY_PATH),"READ");
  t = (TTree*) f->Get("tree");
  double xy1, xy3, xy0;
  int view, pln, ch;
  t->SetBranchAddress("xy1",  &xy1);
  t->SetBranchAddress("xy3",  &xy3);
  t->SetBranchAddress("view", &view);
  t->SetBranchAddress("pln",  &pln);
  t->SetBranchAddress("ch",   &ch);

  for(int i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    if(ch==0){
        //xy0 = (xy1 + xy3) / 2. + 0.1; //Offset
        xy0 = 0.3; //Offset
    }
    if(ch<40){
      position_xy[view][pln][ch] = (xy1 + xy3) / 2.; //mm
      position_xy[view][pln][ch] = position_xy[view][pln][ch] - xy0;
      //std::cout << "position_xy:  " << position_xy[view][pln][ch] << std::endl;
    }
    else{
      position_xy[view][pln][ch] = 0;
    }
  }
  f->Close();

  f = TFile::Open(Form("%s/dimension/iwm/z_position.root",MY_PATH),"READ");
  t = (TTree*) f->Get("tree");
  double z;
  double z0;
  t->SetBranchAddress("z",    &z);
  t->SetBranchAddress("view", &view);
  t->SetBranchAddress("pln",  &pln);

  for(int i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    if(i==0){
        z0 = z; //Offset
    }
    for(int ch_temp=0; ch_temp<80; ch_temp++){
    	position_z[view][pln][ch_temp] = z; //mm
    	position_z[view][pln][ch_temp] = position_z [view][pln][ch_temp] - z0;
    }
  }
  f->Close();

  std::cout << "DetectorDimension:: APPLY MEASURED SCINTILLATOR DIMENSION" << std::endl;

  //===================
  //===   WAGASCI   ===
  //===================
  f = TFile::Open(Form("%s/dimension/wagasci/Scintillator_position.root",MY_PATH),"READ");
  t = (TTree*) f->Get("tree");
  double wxy1, wxy3, wxy0;
  int wview, wpln, wch;
  t->SetBranchAddress("xy1",  &wxy1);
  t->SetBranchAddress("xy3",  &wxy3);
  t->SetBranchAddress("view", &wview);
  t->SetBranchAddress("pln",  &wpln);
  t->SetBranchAddress("ch",   &wch);

  for(int i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    if(wch==0){
        //wxy0 = (wxy1 + wxy3) / 2. + 0.1; //Offset
        wxy0 = 0.3; //Offset
    }
    if(wch<40){
      wposition_xy[wview][wpln][wch] = (wxy1 + wxy3) / 2.; //mm
      wposition_xy[wview][wpln][wch] = wposition_xy[wview][wpln][wch] - wxy0;
      //std::cout << "wposition_xy: " << wposition_xy[wview][wpln][wch] << std::endl;
    }
    else{
      wposition_xy[wview][wpln][wch] = 0;
    }
  }
  f->Close();

  f = TFile::Open(Form("%s/dimension/wagasci/z_position.root",MY_PATH),"READ");
  t = (TTree*) f->Get("tree");
  double wz;
  double wz0;
  t->SetBranchAddress("z",    &wz);
  t->SetBranchAddress("view", &wview);
  t->SetBranchAddress("pln",  &wpln);

  for(int i=0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    if(i==0){
        wz0 = wz; //Offset
    }
    for(int ch_temp=0; ch_temp<80; ch_temp++){
    	wposition_z[wview][wpln][ch_temp] = wz; //mm
    	wposition_z[wview][wpln][ch_temp] = wposition_z [wview][wpln][ch_temp] - wz0;
    }
  }
  f->Close();
  std::cout << "APPLY MEASURED SCINTILLATOR DIMENSION" << std::endl;


}

void DetectorDimension::GetPos(int mod, int pln, int view, int ch, double *x, double *y, double *z){

  double Offset_x, Offset_y, Offset_z;
  double x_tmp, y_tmp, z_tmp;

  if(mod<0 || pln<0 || ch<0){
    std::cout << "Bad mod#. mod>=0 && pln>=0 && ch>=0.";
    exit(1);
  }

  this->GetPosInMod(mod,pln,view,ch,&x_tmp,&y_tmp,&z_tmp);

  // INGRID
  if(mod<NUMINGMOD){
    if(mod<=MOD_INGRID_H){
      Offset_x = C_INGHMotherPosX + C_INGStart + C_INGSpace*mod;
      Offset_y = C_INGHMotherPosY;
      Offset_z = C_INGHMotherPosZ;
    }else{
      Offset_x = C_INGVMotherPosX;
      Offset_y = C_INGVMotherPosY + C_INGStart + C_INGSpace*(mod-7);
      Offset_z = C_INGVMotherPosZ;
    }
  }
  // ProtonModule
  else if(mod==MOD_PM){
    Offset_x = C_PMMotherPosX;
    Offset_y = C_PMMotherPosY;
    Offset_z = C_PMMotherPosZ;
  }
  // OnAxis WaterModule
  else if(mod==MOD_ONAXIS_WM){
    Offset_x = C_PMMotherPosX; //same as Proton Module
    Offset_y = C_PMMotherPosY;
    Offset_z = C_PMMotherPosZ;
  }
  //WaterModule at B2
  else if(mod==MOD_B2_WM){
    Offset_x = C_B2MotherPosX + C_B2WMPosX;
    Offset_y = C_B2MotherPosY + C_B2WMPosY;
    Offset_z = C_B2MotherPosZ + C_B2WMPosZ;
  }
  //ProtonModule(CHModule) at B2
  else if(mod==MOD_B2_CH){
    Offset_x = C_B2MotherPosX + C_B2CHPosX;
    Offset_y = C_B2MotherPosY + C_B2CHPosY;
    Offset_z = C_B2MotherPosZ + C_B2CHPosZ;
  }
  //Downstream INGRID at B2
  else if(mod==MOD_B2_INGRID){
    Offset_x = C_B2MotherPosX + C_B2INGPosX;
    Offset_y = C_B2MotherPosY + C_B2INGPosY;
    Offset_z = C_B2MotherPosZ + C_B2INGPosZ;
  }
  else{
    Offset_x = 0.;
    Offset_y = 0.;
    Offset_z = 0.;
  }

  *x = x_tmp + Offset_x;
  *y = y_tmp + Offset_y;
  *z = z_tmp + Offset_z;

#ifdef DEBUG_DETDIM
  std::cout << "###DEBUG### mod:" << mod << " pln:" << pln <<  " view:" << view << " ch:" << ch
	    << " x:" << *x << " y:" << *y << " z:" << *z << std::endl;
#endif
}


void DetectorDimension::GetPosInMod(int mod, int pln, int view, int ch, double *x, double *y, double *z){

  double x_tmp,y_tmp,z_tmp;

  if(mod<0 || pln<0 || ch<0){
    std::cout << "Bad mod#. mod>=0 && pln>=0 && ch>=0.";
    exit(1);
  }
  // INGRID
  if(mod<NUMINGMOD){
    this->GetPosING(mod,pln,view,ch,&x_tmp,&y_tmp,&z_tmp);
  }
  // ProtonModule
  else if(mod==MOD_PM){
    this->GetPosPM(pln,view,ch,&x_tmp,&y_tmp,&z_tmp);
  }
  // OnAxis WaterModule & WaterModule at B2
  else if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
    int grid;
    if(ch<C_WMNumXYLayerCh)                     grid = 0;
    else if (ch<C_WMNumXYLayerCh+C_WMNumGridCh) grid = 1;
    else                                        grid = 2;
    this->GetPosWM(mod,pln,view,ch,grid,&x_tmp,&y_tmp,&z_tmp);
  }
  else if(mod==MOD_B2_CH){
    this->GetPosPM(pln,view,ch,&x_tmp,&y_tmp,&z_tmp);
  }
  // Downstream INGRID at B2
  else if(mod==MOD_B2_INGRID){
    this->GetPosING(8,pln,view,ch,&x_tmp,&y_tmp,&z_tmp);
  }

  *x = x_tmp;
  *y = y_tmp;
  *z = z_tmp;

#ifdef DEBUG_DETDIM
  std::cout << "###DEBUG### mod:" << mod << " pln:" << pln <<  " view:" << view << " ch:" << ch
	    << " x:" << *x << " y:" << *y << " z:" << *z << std::endl;
#endif
}


//INGRID
void DetectorDimension::GetPosING(int mod, int pln, int view, int ch, double *x, double *y, double *z){

  //DEBUG
  //cout << "GetPosING(mod:"<<mod<<",pln:"<<pln<<",view:"<<view<<",ch:"<<ch<<");" << endl;

  //INGRID horizontal module
  if(mod<=MOD_INGRID_H){

    //if(mod==0) INGNumVetoPln = 4;
    //else       INGNumVetoPln = 3;
    INGNumVetoPln = 4;

    if(pln<C_INGNumPln){
      if(ch>=C_INGNumCh){
	std::cout << "Bad channel arrangement for INGRID. ch is 0-23.";
	exit(1);
      }
      if(mod==3) *z = C_INGPlnStart + C_INGScintiThick*view + (C_INGPlnDist-0.2)*pln + 1.;
      else       *z = C_INGPlnStart + C_INGScintiThick*view + C_INGPlnDist*pln;
      if(view==0){ *x = 0.; *y = C_INGChStart + C_INGScintiWidth*ch;}
      else       { *x = C_INGChStart + C_INGScintiWidth*ch; *y = 0.;}
    }
    else if(pln < C_INGNumPln + INGNumVetoPln){
      if(ch>=C_INGNumVetoCh){
	std::cout << "Bad channel arrangement for INGRID veto. ch is 0-21.";
	exit(1);
      }
      *z = C_INGVetoStartZ + C_INGScintiWidth*ch;
      if      (pln==11&&view==1){ *x = C_INGVetoPos4X; *y = C_INGVetoPos4Y;} //rihgt (only for mod==0)
      else if (pln==12&&view==1){ *x = C_INGVetoPos3X; *y = C_INGVetoPos3Y;} //left
      else if (pln==13&&view==0){ *x = C_INGVetoPos2X; *y = C_INGVetoPos2Y;} //bottom
      else if (pln==14&&view==0){ *x = C_INGVetoPos1X; *y = C_INGVetoPos1Y;} //top
    }
    else{
      std::cout << "mod: " << mod <<  " pln:"<< pln << " view:"<<view<<" ch:"<<ch << std::endl;
      std::cout << "Bad pln/channel arrangement for INGRID. pln is 0-10. ch is 0-23.";
      exit(1);
    }
  }

  //INGRID vertical module
  else if(mod<NUMINGMOD){ 

    //if(mod==7) INGNumVetoPln = 4;
    //else       INGNumVetoPln = 3;
    INGNumVetoPln = 4;

    if(pln<C_INGNumPln){
      if(ch>=C_INGNumCh){
	std::cout << "Bad channel arrangement for INGRID. ch is 0-23.";
	exit(1);
      }
	    
      *z = C_INGPlnStart + C_INGScintiThick*view + C_INGPlnDist*pln;
      if(view==0){ *x = 0.; *y = C_INGChStart + C_INGScintiWidth*ch;}
      else       { *x = C_INGChStart + C_INGScintiWidth*ch; *y = 0.;}
    }
    else if(pln < C_INGNumPln + INGNumVetoPln){
      if(ch>=C_INGNumVetoCh){
	std::cout << "Bad channel arrangement for INGRID veto. ch is 0-21.";
	exit(1);
      }
      *z = C_INGVetoStartZ + C_INGScintiWidth*ch;
      if      (pln==11&&view==1){ *x = C_INGVetoPos4X; *y = C_INGVetoPos4Y;} //right
      else if (pln==12&&view==1){ *x = C_INGVetoPos3X; *y = C_INGVetoPos3Y;} //left
      else if (pln==13&&view==0){ *x = C_INGVetoPos2X; *y = C_INGVetoPos2Y;} //bottom
      else if (pln==14&&view==0){ *x = C_INGVetoPos1X; *y = C_INGVetoPos1Y;} //top
    }
    else{
      std::cout << "mod: " << mod <<  " pln:"<< pln << " view:"<<view<<" ch:"<<ch << std::endl;
      std::cout << "Bad pln/channel arrangement for INGRID. pln is 0-10. ch is 0-23.";
      exit(1);
    }

  }

}

//Proton Module
void DetectorDimension::GetPosPM(int pln, int view, int ch, double *x, double *y, double *z){

  //DEBUG
  //cout << "GetPosPM(pln:"<<pln<<",view:"<<view<<",ch:"<<ch<<");" << endl;

  if(pln==0){
    if(ch>=C_PMNumCh1){
      std::cout << "Bad channel arrangement for PM (pln=0). ch is 0-23.";
      exit(1);
    }
    if(view==0){ *x = 0; *y = C_INGChStart + C_INGScintiWidth*ch;}
    else       { *x = C_INGChStart + C_INGScintiWidth*ch; *y = 0;}
    *z = C_PMPlnStart + C_PMPlnDist*view; //view:0 horizontal, 1 vertical
  }
  else if(pln<C_PMNumPln){
    *z = C_PMPlnStart + C_PMPlnDist*(2*pln+view-1) + C_PMPlnDist_First;
    if(ch<C_PMNumCh_side){
      if(view==0){ *x = 0.; *y = C_PMChStart + C_INGScintiWidth*ch;}
      else       { *x = C_PMChStart + C_INGScintiWidth*ch; *y = 0.;}
    }else if(ch < C_PMNumCh_side + C_PMNumCh_mid){
      if(view==0){ *x = 0.; *y = C_PMChStart + C_INGScintiWidth*(C_PMNumCh_side-0.5) + C_PMScintiWidth*(0.5+ch-C_PMNumCh_side);}
      else       { *x = C_PMChStart + C_INGScintiWidth*(C_PMNumCh_side-0.5) + C_PMScintiWidth*(0.5+ch-C_PMNumCh_side); *y = 0.;}
    }else if(ch<C_PMNumCh){
      if(view==0){ *x = 0; *y = C_PMChStart + C_INGScintiWidth*(C_PMNumCh_side-0.5) 
			     + C_PMScintiWidth*C_PMNumCh_mid + C_INGScintiWidth*(0.5+ch-C_PMNumCh_mid-C_PMNumCh_side);}
      else       { *y = 0; *x = C_PMChStart + C_INGScintiWidth*(C_PMNumCh_side-0.5) 
			     + C_PMScintiWidth*C_PMNumCh_mid + C_INGScintiWidth*(0.5+ch-C_PMNumCh_mid-C_PMNumCh_side);}
    }
    else{
      std::cout << " pln:"<< pln << " view:"<<view<<" ch:"<<ch << std::endl;
      std::cout << "Bad channel arrangment for PM (pln>0). ch is 0-31.";
      exit(1);
    }
  }
  else if(pln < C_PMNumPln + C_PMNumVetoPln){
    if(ch>=C_PMNumVetoCh){
      std::cout << "Bad channel arrangement for PM veto. ch is 0-16.";
      exit(1);
    }
    *z = C_PMVetoStartZ + C_INGScintiWidth*ch;
    if     (pln==18&&view==0){ *x = C_PMVetoPos2X; *y = C_PMVetoPos2Y; } //bottom
    else if(pln==19&&view==1){ *x = C_PMVetoPos3X; *y = C_PMVetoPos3Y; } //left
    else if(pln==20&&view==0){ *x = C_PMVetoPos1X; *y = C_PMVetoPos1Y; } //top
    else if(pln==21&&view==1){ *x = C_PMVetoPos4X; *y = C_PMVetoPos4Y; } //right
  }
  else{
    std::cout <<  " pln:"<< pln << " view:"<<view<<" ch:"<<ch << std::endl;
    std::cout << "Bad pln arrangment for PM. pln is 0-21." << std::endl;
    exit(1);
  }

}

void DetectorDimension::GetPosWM(int mod, int pln, int view, int ch, int grid, double *x, double *y, double *z){

  //DEBUG
  //cout << "GetPosWM(mod:"<<mod<< ",pln:"<<pln<<",view:"<<view<<",ch:"<<ch<<",grid:"<<grid<<");" << endl;

  //Measured
  //===============================
  //===   INGRID Water Module   ===
  //===============================
  double posxy,posz;
  if(view==0){
    *x = 0.;
    if(grid==0 && ch<C_WMNumXYLayerCh){
      if(mod==MOD_ONAXIS_WM){
        posxy = position_xy[view][pln][ch]; 
        posz  = position_z [view][pln][ch];
      }
      else if(mod==MOD_B2_WM){
        posxy = wposition_xy[view][pln][ch]; 
        posz  = wposition_z [view][pln][ch];
      }
      *y = C_WMChStart  + posxy;
      *z = C_WMPlnStart + posz ;
    }
    else if(grid==1 && ch>=C_WMNumXYLayerCh && ch<C_WMNumXYLayerCh+C_WMNumGridCh){
      if(mod==MOD_ONAXIS_WM){
        posz  = position_z [0][pln][ch];
      }
      else if(mod==MOD_B2_WM){
        posz  = wposition_z [0][pln][ch];
      }
      *y = C_WMGridChStart + C_WMScintiSlitStep*(ch-C_WMNumXYLayerCh);
      *z = C_WMPlnStart + C_WMScintiThick/2. + C_WMScintiWidth/2. + posz;
    }
    else if(grid==2 && ch>=C_WMNumXYLayerCh+C_WMNumGridCh && ch<C_WMNumCh){
      if(mod==MOD_ONAXIS_WM){
        posz  = position_z [1][pln][ch];
      }
      else if(mod==MOD_B2_WM){
        posz  = wposition_z [1][pln][ch];
      }
      *y = C_WMGridChStart + C_WMScintiSlitStep*(ch-C_WMNumXYLayerCh-C_WMNumGridCh);
      *z = C_WMPlnStart + C_WMScintiThick/2. + C_WMScintiWidth/2. + posz;
    }
    else{
      std::cout << " pln:"<< pln << " view:"<<view<<" ch:"<<ch << " grid:" << grid << std::endl;
      std::cout << "Bad arrangement for INGRID Water Module. (view0)" << std::endl;
      exit(1);
    }
  }
  else if(view==1){
    *y = 0.;
    if(grid==0 && ch<C_WMNumXYLayerCh){
      if(mod==MOD_ONAXIS_WM){
        posxy = position_xy[view][pln][ch]; 
        posz  = position_z [view][pln][ch];
      }
      else if(mod==MOD_B2_WM){
        posxy = wposition_xy[view][pln][ch]; 
        posz  = wposition_z [view][pln][ch];
      }
      *x = C_WMChStart  + posxy;
      *z = C_WMPlnStart + posz;
    }
    else if(grid==1 && ch>=C_WMNumXYLayerCh && ch<C_WMNumXYLayerCh+C_WMNumGridCh){
      if(mod==MOD_ONAXIS_WM){
        posz  = position_z [0][pln][ch];
      }
      else if(mod==MOD_B2_WM){
        posz  = wposition_z [0][pln][ch];
      }
      *x = C_WMGridChStart + C_WMScintiSlitStep*(ch-C_WMNumXYLayerCh);
      *z = C_WMPlnStart + C_WMScintiThick/2. + C_WMScintiWidth/2. + posz;
    }
    else if(grid==2 && ch>=C_WMNumXYLayerCh+C_WMNumGridCh && ch<C_WMNumCh){
      if(mod==MOD_ONAXIS_WM){
        posz  = position_z [1][pln][ch];
      }
      else if(mod==MOD_B2_WM){
        posz  = wposition_z [1][pln][ch];
      }
      *x = C_WMGridChStart + C_WMScintiSlitStep*(ch-C_WMNumXYLayerCh-C_WMNumGridCh);
      *z = C_WMPlnStart + C_WMScintiThick/2. + C_WMScintiWidth/2. + posz;
    }
    else{
      std::cout << " pln:"<< pln << " view:"<<view<<" ch:"<<ch <<" grid:"<<grid << std::endl;
      std::cout << "Bad arrangement for INGRID Water Module. (view1)" << std::endl;
      exit(1);
    }
  }
  else{
    std::cout << " pln:"<< pln << " view:"<<view<<" ch:"<<ch <<" grid:"<<grid<<std::endl;
    std::cout << "Bad arrangment for INGRID Water Module" << std::endl;
    exit(1);
  }

#ifdef DESIGND_VALUE
  //Designed value (Default)
  if(view==0){
    *x = 0.;
    if(grid==0 && ch<C_WMNumXYLayerCh){
      *y = C_WMChStart + C_WMScintiWidth*ch;
      *z = C_WMPlnStart + C_WMPlnDist*pln;
    }
    else if(grid==1 && ch>=C_WMNumXYLayerCh && ch<C_WMNumXYLayerCh+C_WMNumGridCh){
      *y = C_WMGridChStart + C_WMScintiSlitStep*(ch-C_WMNumXYLayerCh);
      *z = C_WMPlnStart + C_WMPlnDist*pln + C_WMScintiThick/2. + C_WMScintiWidth/2.+0.1;
    }
    else if(grid==2 && ch>=C_WMNumXYLayerCh+C_WMNumGridCh && ch<C_WMNumCh){
      *y = C_WMGridChStart + C_WMScintiSlitStep*(ch-C_WMNumXYLayerCh-C_WMNumGridCh);
      *z = C_WMPlnStart + C_WMPlnDist*pln + C_WMLayerDist + C_WMScintiThick/2. + C_WMScintiWidth/2.+0.1;
    }
    else{
      std::cout << " pln:"<< pln << " view:"<<view<<" ch:"<<ch << " grid:" << grid << std::endl;
      std::cout << "Bad arrangement for INGRID Water Module. (view0)" << std::endl;
      exit(1);
    }
  }
  else if(view==1){
    *y = 0.;
    if(grid==0 && ch<C_WMNumXYLayerCh){
      *x = C_WMChStart + C_WMScintiWidth*ch;
      *z = C_WMPlnStart + C_WMPlnDist*pln + C_WMLayerDist;
    }
    else if(grid==1 && ch>=C_WMNumXYLayerCh && ch<C_WMNumXYLayerCh+C_WMNumGridCh){
      *x = C_WMGridChStart + C_WMScintiSlitStep*(ch-C_WMNumXYLayerCh);
      *z = C_WMPlnStart + C_WMPlnDist*pln + C_WMScintiThick/2. + C_WMScintiWidth/2.+0.1;
    }
    else if(grid==2 && ch>=C_WMNumXYLayerCh+C_WMNumGridCh && ch<C_WMNumCh){
      *x = C_WMGridChStart + C_WMScintiSlitStep*(ch-C_WMNumXYLayerCh-C_WMNumGridCh);
      *z = C_WMPlnStart + C_WMPlnDist*pln + C_WMLayerDist + C_WMScintiThick/2. + C_WMScintiWidth/2.+0.1;
    }
    else{
      std::cout << " pln:"<< pln << " view:"<<view<<" ch:"<<ch <<" grid:"<<grid << std::endl;
      std::cout << "Bad arrangement for INGRID Water Module. (view1)" << std::endl;
      exit(1);
    }
  }
  else{
    std::cout << " pln:"<< pln << " view:"<<view<<" ch:"<<ch <<" grid:"<<grid<<std::endl;
    std::cout << "Bad arrangment for INGRID Water Module" << std::endl;
    exit(1);
  }
#endif

#ifdef DEBUG_DETDIM
  std::cout << " pln:"<< pln << " view:"<<view<<" ch:"<<ch <<" grid:"<< grid << std::endl;
  std::cout << "   => X:" << *x << " Y:" << *y <<" Z:" << *z << std::endl;
#endif

}


bool DetectorDimension::GetWMGridCh(int pln, int view, int ch, int *grid, int *gridch){
  if(ch>=0 && ch<40){
    *grid=0;
    *gridch=ch;
  }
  else if(ch>=40 && ch < 60){
    *grid=1;
    *gridch=ch-40;
  }
  else if(ch>=60 && ch < 80){
    *grid=2;
    *gridch=ch-60;
  }
  else{ 
    std::cout << "DetectorDimension  chnum exceed 80 or less than 0" << "\n";
    return false;
  }

  return true;
}


bool DetectorDimension::GetWMGridCellID(int mod, int pln, int view, int ch, double posx, double posy, double posz,
					int* gridcell_id_x1, int* gridcell_id_x2, int* gridcell_id_y1, int* gridcell_id_y2){

#ifdef DEBUG_DETDIM
  std::cout << " ===== GetWMGridCellID ======== " << std::endl;
#endif
  int grid, gridch;
  double posx_ch,posy_ch,posz_ch;
  this->GetWMGridCh(pln,view,ch,&grid,&gridch);
  if(view==0){
    if(grid==0){
      if     (gridch== 0){ *gridcell_id_x1 = -1; *gridcell_id_x2 =  0; }
      else if(gridch==39){ *gridcell_id_x1 = 19; *gridcell_id_x2 = -1; }
      else{
        *gridcell_id_x1 = (gridch-1)/2;
        *gridcell_id_x2 = (gridch-1)/2+1;
      }
      for(int i=0;i<20;i++){
	if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
	  this->GetPosWM(mod, pln, 1, i+40, 1, &posx_ch, &posy_ch, &posz_ch);
	}
	if(posx_ch > posx){
	  *gridcell_id_y1 = i-1;
	  *gridcell_id_y2 = i;
	  break;
	}
	if(i==19){
          *gridcell_id_y1 = 19;
          *gridcell_id_y2 = -1;
        }
      }
    }
    else if(grid==1||grid==2){
      *gridcell_id_x1= gridch;
      *gridcell_id_x2= -1;
      for(int i=0;i<20;i++){
	if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
	  this->GetPosWM(mod, pln, 1, i+40, 1, &posx_ch, &posy_ch, &posz_ch);
	}
	if(posx_ch > posx){
	  *gridcell_id_y1 = i-1;
	  *gridcell_id_y2 = i;
	  break;
	}
	if(i==19){
          *gridcell_id_y1 = 19;
          *gridcell_id_y2 = -1;
        }
      }
    }
  }
  else if(view==1){
    if(grid==0){
      if     (gridch== 0){ *gridcell_id_y1 = -1; *gridcell_id_y2 =  0; }
      else if(gridch==39){ *gridcell_id_y1 = 19; *gridcell_id_y2 = -1; }
      else{
        *gridcell_id_y1 = (gridch-1)/2;
        *gridcell_id_y2 = (gridch-1)/2+1;
      }
      for(int i=0;i<20;i++){
	if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
	  this->GetPosWM(mod, pln, 0, i+40, 1, &posx_ch, &posy_ch, &posz_ch);
	}
	if(posy_ch > posy){
	  *gridcell_id_x1 = i-1;
	  *gridcell_id_x2 = i;
	  break;
	}
	if(i==19){
          *gridcell_id_x1 = 19;
          *gridcell_id_x2 = -1;
        }
      }
    }
    else if(grid==1||grid==2){
      *gridcell_id_y1= gridch;
      *gridcell_id_y2= -1;
      for(int i=0;i<20;i++){
	if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
	  this->GetPosWM(mod,pln, 0, i+40, 1, &posx_ch, &posy_ch, &posz_ch);
	}
	if(posy_ch > posy){
	  *gridcell_id_x1 = i-1;
	  *gridcell_id_x2 = i;
	  break;
	}
	if(i==19){
          *gridcell_id_x1 = 19;
          *gridcell_id_x2 = -1;
        }
      }
    }
  }

  return true;

}


void DetectorDimension::GetScintiID(int mod, int view, int pln, int gridcell_id, int* cross_n, int* ch){
  if(view==0){
    if(gridcell_id==0){
      *cross_n=2;
      ch[0] = 0;
      ch[1] = 40;
    }
    else if(gridcell_id>=1 && gridcell_id<=19){
      *cross_n=4;
      ch[0] = 2*gridcell_id -1;
      ch[1] = 2*gridcell_id;
      ch[2] = 40+gridcell_id-1;
      ch[3] = 40+gridcell_id;
    }
    else if(gridcell_id==20){
      *cross_n=2;
      ch[0] = 39;
      ch[1] = 59;
    }
    else if(gridcell_id==21){
      *cross_n=1;
      ch[0] = 60;
    }
    else if(gridcell_id>=22 && gridcell_id<=40){
      *cross_n=2;
      ch[0] = 40+gridcell_id-2;
      ch[1] = 40+gridcell_id-1;
    }
    else if(gridcell_id==41){
      *cross_n=1;
      ch[0] = 79;
    }
  }
  else if(view==1){
    if(gridcell_id==0){
      *cross_n=1;
      ch[0] = 40;
    }
    else if(gridcell_id>=1 && gridcell_id<=19){
      *cross_n=2;
      ch[0] = 40+gridcell_id-1;
      ch[1] = 40+gridcell_id;
    }
    else if(gridcell_id==20){
      *cross_n=1;
      ch[0] = 59;
    }
    else if(gridcell_id==21){
      *cross_n=2;
      ch[0] = 0;
      ch[1] = 60;
    }
    else if(gridcell_id>=22 && gridcell_id<=40){
      *cross_n=4;
      ch[0] = 2*(gridcell_id-21) -1;
      ch[1] = 2*(gridcell_id-21);
      ch[2] = 60+(gridcell_id-21)-1;
      ch[3] = 60+(gridcell_id-21);
    }
    else if(gridcell_id==41){
      *cross_n=2;
      ch[0] = 39;
      ch[1] = 79;
    }
  }
}

bool DetectorDimension::GetReconPlnCh(int mod, int view, int pln, int ch, int axis, int* reconpln, int* reconch){
  if(mod<NUMINGMOD||mod==MOD_B2_INGRID){ //INGRID
    *reconpln= pln;
    *reconch = ch;
  }
  else if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){ //watermodule
    int grid,gridch;
    this->GetWMGridCh(pln, view, ch, &grid, &gridch);
    if(axis==0){
      if(view==0){
	if(grid==0){
	  *reconpln=pln*3;
	  *reconch =gridch;
	}
	else if(grid==1){
	  *reconpln=pln*3+1;
	  *reconch =gridch;
	}
	else if(grid==2){
	  *reconpln=pln*3+2;
	  *reconch =gridch;
	}
	else return false;
      }
      else if(view==1){
	if(grid==0){
	  *reconpln=pln*3+1;
	  *reconch =gridch;
	}
	else if(grid==1){
	  *reconpln=pln*3;
	  *reconch =gridch;
	}
	else if(grid==2){
	  *reconpln=pln*3+2;
	  *reconch =gridch;
	}
	else return false;
      }
      else return false;
    }
    else if(axis==1){
      if(view==0){
	if(grid==0){
	  if(gridch%2==0)*reconpln=(int) (gridch*1.5);
	  else if(gridch%2==1)*reconpln=(int) (gridch*1.5+0.5);
	  else return false;
	  *reconch =pln;
	}
	else if(grid==1){
	  *reconpln=gridch*3+1;
	  *reconch =pln*2;
	}
	else if(grid==2){
	  *reconpln=gridch*3+1;
	  *reconch =pln*2+1;
	}
	else return false;
      }
      else if(view==1){
	if(grid==0){
	  if(gridch%2==0)*reconpln=(int) (gridch*1.5);
	  else if(gridch%2==1)*reconpln=(int) (gridch*1.5+0.5);
	  *reconch =pln;
	}
	else if(grid==1){
	  *reconpln=gridch*3+1;
	  *reconch =pln*2;
	}
	else if(grid==2){
	  *reconpln=gridch*3+1;
	  *reconch =pln*2+1;
	}
	else return false;
      }
      else return false;
    }
    else return false;
  }
  else if(mod==MOD_PM||mod==MOD_B2_CH){ //ProtonModule
    *reconpln= pln;
    *reconch = ch;
  }
  else return false;
	
  return true;
}


bool DetectorDimension::GetRawPlnChGrid(int mod, int view, int reconpln, int reconch, int axis, int* pln, int* gridch, int* grid){
  if(mod<NUMINGMOD||mod==MOD_B2_INGRID){ //INGRID
    *pln    = reconpln;
    *gridch = -1;
    *grid   = -1;
  }
  else if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){ //watermodule

    if(axis==0){
      *pln=reconpln/3;
      *gridch=reconch;
      if(view==0){
	if(reconpln%3==0){
	  *grid=0;
	}
	else if(reconpln%3==1){
	  *grid=1;
	}
	else if(reconpln%3==2){
	  *grid=2;
	}
      }
      else{ //view==1
	if(reconpln%3==0){
	  *grid=1;
	}
	else if(reconpln%3==1){
	  *grid=0;
	}
	else if(reconpln%3==2){
	  *grid=2;
	}
      }
    }
    else if(axis==1){
      if(reconpln%3==1){
	*pln=reconch/2;
	*gridch=reconpln/3;
	if(reconch%2==0)     *grid=1;
	else if(reconch%2==1)*grid=2;
      }
      else{
	*pln=reconch;
	*gridch = (int) (reconpln/1.5);
	*grid=0;
      }
    }
    else return false;
  }
  else if(mod==MOD_PM||mod==MOD_B2_CH){ //ProtonModule
    *pln    = reconpln;
    *gridch = -1;
    *grid   = -1;
  }
  else return false;
	
  return true;
}

int DetectorDimension::GetChMax(int mod, int view, int pln, int axis){
  if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){ //WaterModule
    if(axis==0){
      if(view==0){
	if(pln%3==0)return C_WMNumXYLayerCh;
	else        return C_WMNumGridCh;
      }
      else if(view==1){
	if(pln%3==1)return C_WMNumXYLayerCh;
	else        return C_WMNumGridCh;
      }
      else return 0;
    }
    else if(axis==1){
      if(view==0){
	if(pln%3==1)return C_WMNumPln*2;
	else        return C_WMNumPln;
      }
      else if(view==1){
	if(pln%3==1)return C_WMNumPln*2;
	else        return C_WMNumPln;
      }
      else return 0;
    }
    else return 0;
  }
  else if(mod==MOD_PM||mod==MOD_B2_CH){ //ProtonModule
    if(pln==0) return C_PMNumCh1;
    else       return C_PMNumCh;
  }
  else if(mod<NUMINGMOD||mod==MOD_B2_INGRID){ //INGRID
    return C_INGNumCh;
  }
  else return 0;
}

int DetectorDimension::GetPlnMax(int mod, int view, int pln, int axis){
  if(mod<NUMINGMOD||mod==MOD_B2_INGRID){ //INGRID
    return C_INGNumPln;
  }
  else if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){ //WaterModule
    if(axis==0){
      return C_WMNumPln*3;
    }
    else if(axis==1){
      return C_WMNumXYLayerCh + C_WMNumGridCh;
    }
    else return 0;
  }
  else if(mod==MOD_PM||mod==MOD_B2_CH){ //ProtonModule
    return C_PMNumPln;
  }
  else return 0;
}


bool DetectorDimension::GetPos_TwoDimRecon(int mod, int view, int reconpln, int reconch, int axis, double* posxy, double* posz){
  int gridch,grid,pln;
  double x_tmp, y_tmp, z_tmp;

  //DEBUG
  //std::cout << "GetPos_TwoDimRecon("<<mod<<","<<view<<","<<reconpln<<","<<reconch<<","<<axis<<")"<<std::endl;

  if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
    this->GetRawPlnChGrid(mod, view, reconpln, reconch, axis, &pln, &gridch, &grid);

    //DEBUG
    //std::cout << "GetRawPlnChGrid("<<mod<<","<<view<<","<<reconpln<<","<<reconch<<","<<axis
    //	  <<","<<pln<<","<<gridch<<","<<grid<<");" << std::endl;

    int ch;
    if(grid==0) ch = gridch;
    else ch = gridch + 20*(grid+1);
    if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
      this->GetPosWM(mod,pln,view,ch,grid,&x_tmp,&y_tmp,&z_tmp);
    }	
    *posz = z_tmp;
    if(view==0) *posxy = y_tmp; else *posxy = x_tmp;
    return true;
  }
  else return false;
}

float DetectorDimension::GetScintiWidth(int mod, int view, int reconpln, int reconch, int axis){
  int gridch,grid,pln;
  if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
    this->GetRawPlnChGrid(mod, view, reconpln, reconch, axis, &pln, &gridch, &grid);
    if(axis==0){
      if(grid==0) return C_WMScintiWidth/2.;
      else        return C_WMScintiThick/2.;
    }
    else if(axis==1){
      if(grid==0) return C_WMScintiThick/2.;
      else        return C_WMScintiWidth/2.;
    }
  }
  else if(mod<NUMINGMOD ||mod==MOD_B2_INGRID) return C_INGScintiWidth/2.;
  else if(mod==MOD_PM||mod==MOD_B2_CH){
    if(pln==0||reconch<8||reconch>=24) return C_INGScintiWidth/2.;
    else                               return C_PMScintiWidth/2.;
  }
  else return 0.0;
}

float DetectorDimension::GetScintiThick(int mod, int view, int reconpln, int reconch, int axis){
  int gridch,grid,pln;
  if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
    this->GetRawPlnChGrid(mod, view, reconpln, reconch, axis, &pln, &gridch, &grid);
    if(axis==0){
      if(grid==0) return C_WMScintiThick/2.;
      else        return C_WMScintiWidth/2.;
    }
    else if(axis==1){
      if(grid==0) return C_WMScintiWidth/2.;
      else        return C_WMScintiThick/2.;
    }
  }
  else if(mod<NUMINGMOD ||mod==MOD_B2_INGRID) return C_INGScintiThick/2.;
  else if(mod==MOD_PM||mod==MOD_B2_CH){
    if(pln==0||reconch<8||reconch>=24) return C_INGScintiThick/2.;
    else                               return C_PMScintiThick/2.;
  }

  else return 0.0;
}

#endif
