#ifndef _WM_Ch_config_H
#define _WM_Ch_config_H

#include<iostream>
#include<sstream>
#include<fstream>
using namespace std;

class WM_Ch_config{
private:
  int ReadoutId[300][2];
  //ReadoutId[][0]:RMM 0~4
  //ReadoutId[][1]:TFB 0~47

  int PhysId[1280][2]; //corresponds 1280 ch of WModule, hosomi 160430
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
  WM_Ch_config(){
    const char* mapping_file__pln_ch =
      "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/INGRID/INGRID/v1r1/src/WM_Ch_mapping__pln_ch.txt";
    std::fstream ifs(mapping_file__pln_ch);

    index=0;
    if(!ifs){
      cout << "There is no such a file: " << mapping_file__pln_ch << endl;
    }
    else{
      int numpln;
      int numch;
      while(!ifs.eof()){
        ifs >> numpln >> numch;
        if(index>=1280){break;}
        PhysId[index][0]=numpln;
        PhysId[index][1]=numch;
        index++;
      }
    }

    while(index<1280){
      PhysId[index][0] = -1;
      PhysId[index][1] = -1;
      index++;
    }

  };
  ~WM_Ch_config(){};


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
