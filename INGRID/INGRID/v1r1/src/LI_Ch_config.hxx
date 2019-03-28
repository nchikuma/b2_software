#ifndef _LI_Ch_config_H
#define _LI_Ch_config_H

#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;

class LI_Ch_config{
private:
  int ReadoutId[300][2];
  //ReadoutId[][0]:RMM 0~4
  //ReadoutId[][1]:TFB 0~47

  int PhysId[300][2];
  bool tpl_flag[300];
  //PhysId[][0]:Module 0~14
  //PhysId[][1]:TPL 0~10 VETO 0~3
  //tpl[]=true:tracking plane
  //tpl[]=false:VETO plane

  bool horizontal;
  bool vertical;
  bool full_veto;

  int index;
public:
  LI_Ch_config(){
    //make a list of RMM and TFB
    index=0;
    for(int numrmm=1;numrmm<2;numrmm++){//k
      //for(int numtfb=43;numtfb<44;numtfb++){//k
      //for(int numtfb=1; numtfb<2 ;numtfb++){//added for mppc test 2016/1/11
      for(int numtfb=14;numtfb<15;numtfb++){//added for mppc test 2016/2/24
	ReadoutId[index][0]=numrmm;
	ReadoutId[index][1]=numtfb;
	index++;
      }//numtfb
    }//numrmm
    
    //make a list of nummod and numplane
    index=0;
  for(int nummod=17;nummod<18;nummod++){//k
    for(int numplane=0;numplane<1;numplane++){//k
	  PhysId[index][0]=nummod;
	PhysId[index][1]=numplane;
	tpl_flag[index]=true;
	index++;
      }//numplane
    }//nummod
  };
  ~LI_Ch_config(){};


  bool channel_configuration(int *rmm,int *tfb,int *trip,int *trip_ch,int *mod,int *plane,int *ch, bool *tpl, bool *veto);
  //rmm=0~4, tfb=0~47, trip=0~3, trip_ch=0~15
  //mod=0~13, plane=0~10(tpl), 0~4(veto), ch=0~47,0~21
  bool channel_configuration(int *rmm,int *tfb,int *trip,int *trip_ch,int *mod,int *view,int *plane,int *ch);
  //new channel_configuration for new data structure
  //rmm=0~4, tfb=0~47, trip=0~3, trip_ch=0~15
  //mod=0~13, plane=0~10(tpl), 0~4(veto), ch=0~47,0~21
  bool active_veto(int *mod,int *plane);
  bool get_global_pln(int *mod,int *plane, bool *tpl, bool *veto, int *global_pln);
  bool get_global_ch(int *mod, int *plane,int *ch, bool *tpl, bool *veto, int *global_ch);
  bool get_global_ch(int *mod, int *view, int *plane,int *ch, int *global_ch);
  //for new data structure
  bool get_global_ch(int mod, int view, int plane,int ch, int *global_ch);

};
#endif
