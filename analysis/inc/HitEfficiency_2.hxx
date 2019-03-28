#ifndef __HIT_EFF_2_HXX__
#define __HIT_EFF_2_HXX__

#include "TwoDimRecon.hxx"

/*
bool get_hitpos_in_scinti(
    double intcptx,double slopex,double intcpty,double slopey,
    int hitmod,int hitview,int hitpln,int hitch,
    double &posx,double &posy)
{
  double chx,chy,chz;
  detdim->GetPosInMod(hitmod,hitpln,hitview,hitch,&chx,&chy,&chz);
  double dist;
  if(hitview==0){dist = fabs(chy - intcpty - slopey*chz)/sqrt(1.+slopey*slopey);}
  else          {dist = fabs(chx - intcptx - slopex*chz)/sqrt(1.+slopex*slopex);}

  if(hitmod<=NUMINGMOD||hitmod==MOD_PM){
    if(hitview==0){posx = intcptx+slopex*chz; posy = intcpty+slopey*chz-chy;}
    else          {posx = intcpty+slopey*chz; posy = intcptx+slopex*chz-chx;}
  }
  else if(hitmod==MOD_B2_WM||hitmod==MOD_ONAXIS_WM){
    if(hitview==0){
      if(hitch<40){ posx = intcptx+slopex*chz; posy = intcpty+slopey*chz  -chy; }
      else        { posx = intcptx+slopex*chz; posy = (chy-intcpty)/slopey-chz ;}
    }
    else{
      if(hitch<40){ posx = intcpty+slopey*chz; posy = intcptx+slopex*chz-chx; }
      else        { posx = intcpty+slopey*chz; posy = chz-(chx-intcptx)/slopex;}
    }
  }

  //cout 
  //  << " mod:"  << hitmod
  //  << " view:" << hitview
  //  << " pln:"  << hitpln
  //  << " ch:"   << hitch
  //  << " x:"    << chx
  //  << " y:"    << chy
  //  << " z:"    << chz
  //  << " dist:" << dist
  //  << " posx:" << posx
  //  << " posy:" << posy
  //  << endl;


  return true;

}
*/

#endif
