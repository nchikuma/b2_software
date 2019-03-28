#ifndef __B2_PE_HXX__
#define __B2_PE_HXX__

#include "TwoDimRecon.hxx"

const int START_UNIXT = 1507993200; //Sun Oct 15 00:00:00 JST 2017
const int STOP_UNIXT  = 1527778800; //Sun Jun  1 00:00:00 JST 2018
const int RANGE_UNIXT = STOP_UNIXT - START_UNIXT;
const int NBIN_UNIXT  = (int)(RANGE_UNIXT/86400);

/*
double calc_pathlength(double slopex,double slopey,int mod,int view,int pln,int ch)
{
  const double ANGLE_LIMIT = acos(0.001)*180/PI; // 0.1%
  const double ANGLE_ABORT = -300.;
  double angleX, angleY, angleZ;
  double angle1=0., angle2=0.;
  double path = 0.;
  angleX = 180./PI*acos(slopey/sqrt(1.+pow(slopex,2)+pow(slopey,2))); //Angle from Y-axis (0,1,0)
  if(angleX>90) angleX = 180. - angleX;
  angleY = 180./PI*acos(slopex/sqrt(1.+pow(slopex,2)+pow(slopey,2))); //Angle from X-axis (1,0,0)
  if(angleY>90) angleY = 180. - angleY;
  angleZ = 180./PI*atan(sqrt(pow(slopex,2)+pow(slopey,2))); //Angle from Z-axis (0,0,1)
  
  if(!(fabs(angleX)<ANGLE_LIMIT && fabs(angleY)<ANGLE_LIMIT && fabs(angleZ)<ANGLE_LIMIT)){
    angleX = ANGLE_ABORT; angleY = ANGLE_ABORT; angleZ = ANGLE_ABORT;
  }

  //Classify: Grid or Plane
  double scinti_width, scinti_thick;
  if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
    if(view==0&&ch< 40){angle1=angleZ;angle2=fabs(180./PI*angleX);}
    if(view==1&&ch< 40){angle1=angleZ;angle2=fabs(180./PI*angleY);}
    if(view==0&&ch>=40){angle1=angleY;angle2=fabs(180./PI*angleX);}
    if(view==1&&ch>=40){angle1=angleX;angle2=fabs(180./PI*angleY);}
    scinti_width = C_WMScintiWidth;
    scinti_thick = C_WMScintiThick;
  }
  else{
    if(view==0){angle1=angleZ;angle2=fabs(180./PI*angleX);}
    if(view==1){angle1=angleZ;angle2=fabs(180./PI*angleY);}
    if( (mod==MOD_PM||mod==MOD_B2_CH)&&
        (pln==0||ch<8||ch>=24))
    {
      scinti_width = C_PMScintiWidth;
      scinti_thick = C_PMScintiThick;
    }
    else{
      scinti_width = C_INGScintiWidth;
      scinti_thick = C_INGScintiThick;
    }
  }
  //angle1 = PI/180.*(angle1+0.1); angle2 = PI/180.*(angle2+0.1);
  angle1 = PI/180.*angle1; angle2 = PI/180.*angle2;

  //Pathlength
  if(tan(angle2)<scinti_width/scinti_thick){
    path = (
        scinti_thick*tan(angle2)*tan(angle2)/sin(angle1)
        - scinti_thick*tan(angle2)/cos(angle1)
        + scinti_width/cos(angle1)
        )
      /( scinti_width+scinti_thick*tan(angle2) );
  }
  else{
    path = (scinti_width*tan(angle2)/sin(angle1)) 
      / (scinti_width+scinti_thick*tan(angle2));
  }
  if(path>500.) path = 500.;
  if(path<0.  ) path = 0.;

  return path;
}

double calc_pathlength_wg(double slopex,double slopey,int mod,int view,int pln,int ch)
{
  double path=-1.;
  double x,y,z;

  double norm = sqrt(1.+slopex*slopex+slopey*slopey);
  //double cos_zen = slopey/norm;
  //double cos_azi = 1./norm;
  if(norm<=0.){return -1.;}
  double cos_z = 1./norm;
  double cos_y = slopey/norm;
  double cos_x = slopex/norm;

  if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
    if(pln>=0&&pln<8){
      if     (view==0&&ch< 40){z=C_WMScintiThick; y=C_WMScintiWidth; x=C_WMScintiLength;}
      else if(view==0&&ch>=40){y=C_WMScintiThick; z=C_WMScintiWidth; x=C_WMScintiLength;}
      else if(view==1&&ch< 40){z=C_WMScintiThick; x=C_WMScintiWidth; y=C_WMScintiLength;}
      else if(view==1&&ch>=40){x=C_WMScintiThick; z=C_WMScintiWidth; y=C_WMScintiLength;}
      else{return -1.;}
    }
    else{return -1.;}
  }
  else if(mod>=0&&mod<NUMINGMOD){
    if(pln>=0&&pln<11){
      if     (view==0&&ch<24){z=C_INGScintiThick; y=C_INGScintiWidth; x=C_INGScintiLength;}
      else if(view==1&&ch<24){z=C_INGScintiThick; x=C_INGScintiWidth; y=C_INGScintiLength;}
      else{return -1.;}
    }
    else{return -1.;}
  }
  else if(mod==MOD_PM||mod==MOD_B2_CH){
    if(pln==0){
      if     (view==0&&ch<24){z=C_INGScintiThick; y=C_INGScintiWidth; x=C_INGScintiLength;}
      else if(view==1&&ch<24){z=C_INGScintiThick; x=C_INGScintiWidth; y=C_INGScintiLength;}
      else{return -1.;}
    }
    else if(pln>=1&&pln<18){
      if((ch>=0&&ch<8)||(ch>=24&&ch<32))
      {
        if     (view==0){z=C_INGScintiThick; y=C_INGScintiWidth; x=C_INGScintiLength;}
        else if(view==1){z=C_INGScintiThick; x=C_INGScintiWidth; y=C_INGScintiLength;}
        else{return -1.;}
      }
      else if(ch>=8&&ch<24){
        if     (view==0){z=C_PMScintiThick; y=C_PMScintiWidth; x=C_PMScintiLength;}
        else if(view==1){z=C_PMScintiThick; x=C_PMScintiWidth; y=C_PMScintiLength;}
        else{return -1.;}
      }
      else{return -1.;}
    }
  }
  else{return -1.;}
 

  double volume = x*y*z;
  //double sin_zen,sin_azi;
  //sin_zen = sqrt(1.0-cos_zen*cos_zen);
  //sin_azi = sqrt(1.0-cos_azi*cos_azi);
  //double crosssection=fabs(x*y*cos_zen)+fabs(x*z*sin_zen*cos_azi)+fabs(y*z*sin_zen*sin_azi);
  double crosssection=fabs(x*y*cos_z)+fabs(x*z*cos_y)+fabs(y*z*cos_x);

  if(crosssection>0.){ path = volume/crosssection; }else{ return -1.; }

  if(path>500.0) path =500.0;

  return path; 
}
*/

#endif
