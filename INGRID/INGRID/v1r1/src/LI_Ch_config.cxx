#ifndef _LI_Ch_config_C
#define _LI_Ch_config_C

#include<iostream>
#include<sstream>
#include<fstream>
#include"LI_Ch_config.hxx"
using namespace std;

bool LI_Ch_config::channel_configuration(int *rmm,int *tfb,int *trip,int *trip_ch,int *mod,int *plane,int *ch,bool *tpl,bool *veto){

  //for(int i=0;i<300;i++){
  for(int i=0;i<1;i++){
    if(ReadoutId[i][0]==*rmm&&ReadoutId[i][1]==*tfb){
        *mod=PhysId[i][0];
	*plane=PhysId[i][1];
	*tpl=true;*veto=false;
	if(*trip>1){
		*ch=(*trip-2)*16+(*trip_ch);
		//if(*ch>7)return false;
		//if(*ch<3 || *ch>18)return false;
		//if(*ch<0 || *ch>15)return false;
	}
	else{
		*ch=*trip*16+(*trip_ch);
		//if(*ch>15)return false;
		//if(*ch<4 || *ch>11)return false;
	}

	return true;

    }
  }
  return false;  
}

bool LI_Ch_config::active_veto(int *mod,int *plane){
  return true;
}

bool LI_Ch_config::get_global_pln(int *mod,int *plane, bool *tpl, bool *veto, int *global_pln){
  *global_pln=0;
  int tempmod=*mod;
  int temppln=*plane;
  if(*veto&&!this->active_veto(&tempmod, &temppln))return false;
  for(int nummod=17;nummod<18;nummod++){
    for(int pln=0;pln<1;pln++){
      if((*mod)==nummod&&(*plane)==pln&&*tpl)return true;
      (*global_pln)++; 
    }//tpl
  }//nummod
  return false;
}

bool LI_Ch_config::get_global_ch(int *mod, int *plane, int *ch, bool *tpl, bool *veto, int *global_ch){
  *global_ch=0;
  int tempmod=*mod;
  int temppln=*plane;
  if(plane!=0)return false;

  for(int pln=0;pln<1;pln++){
    for(int numch=0;numch<16;numch++){
      if((*plane)==pln&&(*ch)==numch&&*tpl)return true;
      (*global_ch)++;
    }//nuch
    //for(int numch=16*2;numch<16*2+8;numch++){
    //  if((*plane)==pln&&(*ch)==numch&&*tpl)return true;
    //  (*global_ch)++;
    //}//nuch

  }//pln

  return false;
}


bool LI_Ch_config::channel_configuration(int *rmm,int *tfb,int *trip,int *trip_ch,int *mod, int *view, int *plane,int *ch){
  bool tpl_flag;
  bool veto_flag;
  if(*trip>1)*view=1;
  else *view=0;
  if(!this->channel_configuration(rmm,tfb,trip,trip_ch,mod,plane,ch,&tpl_flag,&veto_flag))return false;
  else if(*plane==0)return true;
  else return false;
}

bool LI_Ch_config::get_global_ch(int *mod, int *view, int *plane,int *ch, int *global_ch){
  *global_ch=0;

  if(*plane>0)return false;

  for(int pln=0;pln<1;pln++){
    for(int numch=0;numch<32;numch++){
    //for(int numch=4;numch<12;numch++){
      if((*plane)==pln&&(*ch)==numch&&*view==0)return true;
      (*global_ch)++;
    }
    for(int numch=0;numch<32;numch++){
    //for(int numch=3;numch<15;numch++){
      if((*plane)==pln&&(*ch)==numch&&*view==1)return true;
      (*global_ch)++;
    }//nuch
 }//pln

  return false;
}

bool LI_Ch_config::get_global_ch(int mod, int view, int plane,int ch, int *global_ch){
  *global_ch=0;

  if(plane>0)return false;

  for(int pln=0;pln<1;pln++){
    for(int numch=0;numch<32;numch++){
      if( plane==pln && ch==numch&&view==0)return true;
      (*global_ch)++;
    }
    for(int numch=0;numch<32;numch++){
      if(plane==pln && ch==numch &&view==1)return true;
      (*global_ch)++;
    }//nuch
  }//pln

  return false;

}



#endif
