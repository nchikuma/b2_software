#ifndef _B2ING_Ch_config_C
#define _B2ING_Ch_config_C

#include<iostream>
#include<sstream>
#include<fstream>
#include"B2ING_Ch_config.hxx"
using namespace std;

bool B2ING_Ch_config::channel_configuration(int *rmm,int *tfb,int *trip,int *trip_ch,int *mod,int *plane,int *ch,bool *tpl,bool *veto){

  for(int i=0;i<nB2INGTFB;i++){
    if(ReadoutId[i][0]==*rmm&&ReadoutId[i][1]==*tfb){
      horizontal=false;
      vertical  =false;
      offdiag   =true;
      full_veto =false;
      *mod=PhysId[i][0];

      //read physical channel
      if(tpl_flag[i]){//tracking plane
	*tpl=true;*veto=false;
	*plane=PhysId[i][1];
	*ch=(*trip)*16+(*trip_ch);
	if(*ch>47){return false;}
	if(*ch<24){ *ch=23-*ch; }
	return true;
      }//tracking plane end

      if(!tpl_flag[i]){//VETO plane
	*tpl=false;*veto=true;
	if(PhysId[i][1]==11){//*TFB==11
	  if(*trip==0){
	    *ch=*trip_ch;
	    *plane=11+3;
	    *ch=21-*ch;
	    return true;
	  }
	  if(*trip==1){
	    if(*trip_ch<6){
	      *ch=(*trip_ch)+16;
	      *plane=11+3;
	      *ch=21-*ch;
	      return true;
	    }
	    if(*trip_ch==6||*trip_ch==7)return false;
	    if(*trip_ch>=8){
	      *ch=*trip_ch-8;
	      *plane=11+0;
	      if(horizontal&&(!(full_veto))){return false;}
	      return true;
	    }
	  }
	  if(*trip==2){
	    *ch=*trip_ch+8;
	    *plane=11+0;
	    if(horizontal&&(!(full_veto))){return false;}
	    if(*ch>21)return false;

	    return true;
	  }
	  if(*trip==3){
	    return false;
	  }
	}//*TFB==11
	if(PhysId[i][1]==12){//*TFB==12
	  if(*trip==0){
	    *ch=*trip_ch;
	    *plane=11+2;
	    //replace physical channel(add 090928 by Otani)
	    *ch=21-*ch;
	    if(vertical&&(!(full_veto))){return false;}
	    return true;
	  }
	  if(*trip==1){
	    if(*trip_ch<6){
	      *ch=*trip_ch+16;
	      *plane=11+2;
	      //replace physical channel(add 090928 by Otani)
	      *ch=21-*ch;
	      if(vertical&&(!(full_veto))){return false;}
	      return true;
	    }
	    if(*trip_ch==6||*trip_ch==7)return false;
	    if(*trip_ch>=8){
	      *ch=*trip_ch-8;
	      *plane=11+1;
	      return true;
	    }
	  }
	  if(*trip==2){
	    *ch=*trip_ch+8;
	    *plane=11+1;
	    if(*ch>21)return false;
	    return true;
	  }
	  if(*trip==3){
	    return false;
	  }
	}//*TFB==12
      }//end of veto plane

    }
  }
  return false;  
}

bool B2ING_Ch_config::active_pln(int mod, int view ,int plane){
  if(plane<11) return true;
  else{
    if     ( view == 0 ){
      if( plane == 13 || plane == 14 )
       return true;
      else
       return false;
    }
    else if( view == 1 ){
      if( plane == 11 || plane == 12 )
        return true;
      else
        return false;
    }
  }
  return false;
}

bool B2ING_Ch_config::active_veto(int *mod,int *plane){
  if(*plane>=11&&*plane<=14) return true;
  else return false;
}

bool B2ING_Ch_config::get_global_pln(int *mod,int *plane, bool *tpl, bool *veto, int *global_pln){
  *global_pln=0;
  int tempmod=*mod;
  int temppln=*plane;
  if(*veto&&!this->active_veto(&tempmod, &temppln))return false;
  int nummod =14;
  for(int pln=0;pln<11;pln++){
    if((*mod)==nummod&&(*plane)==pln&&*tpl)return true;
    (*global_pln)++; 
  }//tpl
  for(int pln=11;pln<15;pln++){
    if((*mod)==nummod&&(*plane)==pln&&*veto)return true;
    if(this->active_veto(&nummod, &pln))(*global_pln)++; 
  }//veto
  return false;
}

bool B2ING_Ch_config::get_global_ch(int *mod, int *plane, int *ch, bool *tpl, bool *veto, int *global_ch){
  *global_ch=0;
  int tempmod=*mod;
  int temppln=*plane;
  if(*veto&&!this->active_veto(&tempmod, &temppln))return false;

  for(int pln=0;pln<11;pln++){
    for(int numch=0;numch<48;numch++){
      if((*plane)==pln&&(*ch)==numch&&*tpl)return true;
      (*global_ch)++;
    }//nuch
  }//pln
  for(int pln=11;pln<15;pln++){
    if(this->active_veto(&tempmod, &pln)){
      for(int numch=0;numch<22;numch++){
	if((*plane)==pln&&(*ch)==numch&&*veto)return true;
	(*global_ch)++;
      }//numch
    }//active
  }//pln

  return false;
}


bool B2ING_Ch_config::channel_configuration(int *rmm,int *tfb,int *trip,int *trip_ch,int *mod, int *view, int *plane,int *ch){
  bool tpl_flag;
  bool veto_flag;
  if(!this->channel_configuration(rmm,tfb,trip,trip_ch,mod,plane,ch,&tpl_flag,&veto_flag))return false;
  if(veto_flag){
    *plane=(*plane);
    if(*plane==11||*plane==12){
      *view=1;
      return true;
    }
    else if(*plane==13||*plane==14){
      *view=0;
      return true;
    }
    return false;
  }
  else if(tpl_flag){
    if(*ch>=24){
      *view=1;
      *ch=47-(*ch);
      return true;
    }
    else{
      *view=0;
      *ch=23-(*ch);
      return true;
    }
  }
  return true;
}

bool B2ING_Ch_config::get_global_ch(int *mod, int *view, int *plane,int *ch, int *global_ch){
  *global_ch=0;
  int temppln=*plane;
  int tempmod=*mod;
  if((*plane>10)&&!this->active_veto(&tempmod, &temppln))return false;
  for(int pln=0;pln<11;pln++){
    for(int numch=0;numch<24;numch++){
      if((*plane)==pln&&(*ch)==numch&&(*view)==0)return true;
      (*global_ch)++;
    }//numchX
    for(int numch=0;numch<24;numch++){
      if((*plane)==pln&&(*ch)==numch&&(*view)==1)return true;
      (*global_ch)++;
    }//numchY
  }//pln
  for(int pln=11;pln<15;pln++){
    int temp=pln;
    int v;
    if(pln==11||pln==12)v=1;
    else v=0;
    if(this->active_veto(&tempmod, &temp)){
      for(int numch=0;numch<22;numch++){
	if(*view==v&&(*plane)==pln&&(*ch)==numch)return true;
	(*global_ch)++;
      }
    }
  }
  return false;

}

bool B2ING_Ch_config::get_global_ch(int mod, int view, int plane,int ch, int *global_ch){
  *global_ch=0;
  int temppln=plane;
  int tempmod=mod;

  if((plane>10)&&!this->active_veto(&tempmod, &temppln))return false;
  for(int pln=0;pln<11;pln++){
    for(int numch=0;numch<24;numch++){
      if((plane)==pln&&(ch)==numch&&(view)==0)return true;
      (*global_ch)++;
    }//numchX
    for(int numch=0;numch<24;numch++){
      if((plane)==pln&&(ch)==numch&&(view)==1)return true;
      (*global_ch)++;
    }//numchY
  }//pln
  for(int pln=11;pln<15;pln++){
    int temp=pln;
    int v;
    if(pln>=11&&pln<=14)v=1;
    else v=0;
    if(this->active_veto(&tempmod, &temp)){
      for(int numch=0;numch<22;numch++){
	if(view==v&&(plane)==pln&&(ch)==numch)return true;
	(*global_ch)++;
      }
    }
  }
  return false;

}



#endif
