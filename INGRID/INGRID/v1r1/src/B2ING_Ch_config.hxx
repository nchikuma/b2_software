#ifndef _B2ING_Ch_config_H
#define _B2ING_Ch_config_H

#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;

#define nB2INGRMM 3
#define nB2INGRMMch 48
#define nB2INGTFB 13

// RMM, RMMch
const int B2ING_TFB[nB2INGTFB][2] = 
{
  {4,20},{4,21},{4,22},{4,23},{4,24},{4,25},
  {3,43},{3,44},{3,45},{3,46},{3,47},
  {1,44},{1,45}
};

class B2ING_Ch_config{
private:
  int ReadoutId[nB2INGTFB][2];
  //ReadoutId[][0]:RMM 0~4
  //ReadoutId[][1]:TFB 0~47

  int PhysId[nB2INGTFB][2];
  bool tpl_flag[nB2INGTFB];
  //PhysId[][0]:Module = 14 for B2 INGRID
  //PhysId[][1]:TPL 0~10 VETO 11~12
  //tpl[]=true:tracking plane
  //tpl[]=false:VETO plane

  bool horizontal;
  bool vertical;
  bool offdiag;
  bool full_veto;

  int index;
public:
  B2ING_Ch_config(){
    for(int iTFB=0;iTFB<nB2INGTFB;iTFB++){
      ReadoutId[iTFB][0]=B2ING_TFB[iTFB][0];
      ReadoutId[iTFB][1]=B2ING_TFB[iTFB][1];
      PhysId   [iTFB][0]=14;
      PhysId   [iTFB][1]=iTFB;
      if(iTFB==11||iTFB==12) tpl_flag[iTFB] = false;
      else                   tpl_flag[iTFB] = true;
    }
  };
  ~B2ING_Ch_config(){};


  bool channel_configuration(int *rmm,int *tfb,int *trip,int *trip_ch,int *mod,int *plane,int *ch, bool *tpl, bool *veto);
  //rmm=0~4, tfb=0~47, trip=0~3, trip_ch=0~15
  //mod=0~13, plane=0~10(tpl), 0~4(veto), ch=0~47,0~21
  bool channel_configuration(int *rmm,int *tfb,int *trip,int *trip_ch,int *mod,int *view,int *plane,int *ch);
  //new channel_configuration for new data structure
  //rmm=0~4, tfb=0~47, trip=0~3, trip_ch=0~15
  //mod=0~13, plane=0~10(tpl), 0~4(veto), ch=0~47,0~21
  bool active_veto(int *mod,int *plane);
  bool active_pln(int mod, int view, int plane);
  bool get_global_pln(int *mod,int *plane, bool *tpl, bool *veto, int *global_pln);
  bool get_global_ch(int *mod, int *plane,int *ch, bool *tpl, bool *veto, int *global_ch);
  bool get_global_ch(int *mod, int *view, int *plane,int *ch, int *global_ch);
  //for new data structure
  bool get_global_ch(int mod, int view, int plane,int ch, int *global_ch);

};
#endif

