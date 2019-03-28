#ifndef _SETUP_H
#define _SETUP_H


#include<iostream>
#include<sstream>
#include "INGRID_Ch_config.cxx"
#include "B2ING_Ch_config.cxx"
#include "PM_Ch_config.cxx"
#include "WM_Ch_config.cxx"
#include "INGRID_BadCh_mapping.cxx"
#include "INGRID_Dimension.cxx"

#define NumMod     17
#define NumTPL     11
#define NumVETO    4
#define NumTFB     24
#define NumCh      64
#define NumTPLCh   48 
#define NumVETOCh  22
#define UseNumCh   48
#define NumCyc     23
#define UseNumTFB  11
#define VetoNumTFB 2
#define UseNumCh   48
#define LayerNumCh 80
#define NumLayer   2
#define GATE       480 //nsec

//for Proton Module
#define NumTPL_PM     18
#define NumVETO_PM    4
#define NumTFB_PM     22
#define NumCh_PM      64
#define NumTPLCh_PM   64 
#define NumVETOCh_PM  17
#define UseNumCh_PM   64
#define UseNumTFB_PM  18
#define VetoNumTFB_PM 3
#define UseNumCh_PM   64
#define LayerNumCh_PM 32
#define NumLayer_PM   2

#define MaxNum_BadCh 100

// Data directory
string ingrid_dir
= "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data";
string data_midas_dir  = Form("%s/data_midas" ,ingrid_dir.c_str());
string data_dst_dir    = Form("%s/data_dst"   ,ingrid_dir.c_str());
string data_cosmic_dir = Form("%s/data_cosmic",ingrid_dir.c_str());
string data_calib_dir  = Form("%s/data_calib" ,ingrid_dir.c_str());

Int_t BadCh_Map[MaxNum_BadCh][4];
Int_t Num_BadCh;
bool  BadCh_Ana = false;

INGRID_Ch_config * fINGRID_Ch_config;
PM_Ch_config     * fPM_Ch_config    ;
WM_Ch_config     * fWM_Ch_config    ;
B2ING_Ch_config  * fB2ING_Ch_config ;
INGRID_Dimension * fINGRID_Dimension;


bool isPM(int rmm, int tfb){
  if( rmm==4 && tfb>=26)
    return true;
  else
    return false;
}

bool isWM(int rmm, int tfb){
  if( rmm==4 && tfb>-1 && tfb<20) return true; //hosomi 160501
  else return false;
}

bool isB2ING(int rmm, int tfb){
  if( 
      (rmm==4 && tfb>=20 && tfb<=25)||
      (rmm==3 && tfb>=43 && tfb<=47)||
      (rmm==1 && tfb>=44 && tfb<=45) )
  {
    return true;
  }
  else return false;
}

bool isPM(int mod){
  if( mod==16 )
    return true;
  else
    return false;
}

bool isWM(int mod){
  if( mod==15 )
    return true;
  else
    return false;
}

bool isB2ING(int mod){
  if( mod==14 )
    return true;
  else
    return false;
}

// _________________________________________________________
void set_BadCh(bool opt=true)
{
  if(opt){
    string cardfilename 
      = Form("%s/ingrid_info/card.txt",ingrid_dir.c_str());
    INGRID_BadCh_mapping* fINGRID_BadCh_mapping = new INGRID_BadCh_mapping();
    fINGRID_BadCh_mapping->readbad(cardfilename,&Num_BadCh,BadCh_Map);
    delete fINGRID_BadCh_mapping;
  }
  else{
    cout << "Bad channel masks are opened for this analysis." << endl;
    Num_BadCh = 0;
  }
  cout << "Number of bad channels : " << Num_BadCh << endl;
  if(Num_BadCh>0){
    cout << "-----------------" << endl;
    cout << " mod view pln ch " << endl;
    for(int ich=0;ich<MaxNum_BadCh&&ich<Num_BadCh;ich++){
      cout << " " << setw(3) << BadCh_Map[ich][0];
      cout << " " << setw(4) << BadCh_Map[ich][1];
      cout << " " << setw(3) << BadCh_Map[ich][2];
      cout << " " << setw(2) << BadCh_Map[ich][3];
      cout << endl;
    }
    cout << "-----------------" << endl;
  }
}

// _________________________________________________________
bool is_BadCh(int mod,int view,int pln,int ch)
{
  for(int ich=0;ich<MaxNum_BadCh&&ich<Num_BadCh;ich++){
    if(
        (BadCh_Map[ich][0]==mod )&&
        (BadCh_Map[ich][1]==view)&&
        (BadCh_Map[ich][2]==pln )&&
        (BadCh_Map[ich][3]==ch  ))
    {
      return true;
    }
  }
  return false;
}

//__________________________________________________________

bool isWAGASCIconfig(int *rmm,int *tfb,int *trip,int *trip_ch,int *mod,int *view,int *plane,int *ch)
{
  bool read = false;

  int irmm     = *rmm;
  int itfb     = *tfb;
  int itrip    = *trip;
  int itrip_ch = *trip_ch;
  int imod, iview, ipln, ich;

  if( isPM(irmm, itfb) ){
    if(fPM_Ch_config
        ->channel_configuration(&irmm,&itfb,&itrip,&itrip_ch,&imod,&iview,&ipln,&ich))
    { read = true; }
  }
  else if( isWM(irmm, itfb) ){
    if(fWM_Ch_config
        ->channel_configuration(&irmm,&itfb,&itrip,&itrip_ch,&imod,&iview,&ipln,&ich))
    { read = true; }
  }
  else if( isB2ING(irmm, itfb) ){
    if(fB2ING_Ch_config
        ->channel_configuration(&irmm,&itfb,&itrip,&itrip_ch,&imod,&iview,&ipln,&ich))
    { read = true; }
  }
  else{
    if(fINGRID_Ch_config
        ->channel_configuration(&irmm,&itfb,&itrip,&itrip_ch,&imod,&iview,&ipln,&ich))
    { if(imod==3){read = true;} }
  }

  *mod   = imod;
  *view  = iview;
  *plane = ipln;
  *ch    = ich;

  return read;
}
//__________________________________________________________

#endif
