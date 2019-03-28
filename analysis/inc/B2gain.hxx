#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TSpectrum.h"
#include "TTree.h"
#include "TArc.h"
#include "TBox.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TText.h"
#include "TGaxis.h"
#include "TChain.h"
#include "TColor.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <stdio.h>


#include "Const.hh"
#include "INGRID_BadCh_mapping.hh"
#include "EVENTSUMMARY.h"


void plot_gain_history(string filename,bool BadCh_Ana);


bool getchannel_id(int mod,int id,int& view, int& pln, int& ch)
{
  int all   = 0;
  int npln  = 0;
  int nch   = 0;
  int nch2  = 0;
  int ntrkch = 0;
  int nvpln = 0;
  int nvch  = 0;
  int tmp;
  if(mod==3){
     all   = 594;
     npln  =  11; 
     nch   =  24; 
     ntrkch = npln*nch*2;
     nvpln =   3; 
     nvch  =  22;
     if(id<all){
       if(id<ntrkch){
         view = id/(ntrkch/2); 
         pln  = id%(ntrkch/2)/nch;
         ch   = id%(ntrkch/2)%nch;
       }
       else{
         pln  = (id-ntrkch)/nvch+npln;
         ch   = (id-ntrkch)%nvch;
         if(pln==11||pln==12) view=1; else view=0;
       }
     }
     else{return false;}
  }
  else if(mod==14){
     all   = 616;
     npln  =  11; 
     nch   =  24; 
     ntrkch = npln*nch*2;
     nvpln =   4; 
     nvch  =  22;
     if(id<all){
       if(id<ntrkch){
         view = id/(ntrkch/2); 
         pln  = id%(ntrkch/2)/nch;
         ch   = id%(ntrkch/2)%nch;
       }
       else{
         pln  = (id-ntrkch)/nvch+npln;
         ch   = (id-ntrkch)%nvch;
         if(pln==11||pln==12) view=1; else view=0;
       }
     }
     else{return false;}
  }
  else if(mod==15){
     all   = 1280;
     npln  =    8; 
     nch   =   80; 
     ntrkch = npln*nch*2;
     nvpln =   0; 
     nvch  =   0;
     if(id<all){
       view = id/(ntrkch/2); 
       pln  = id%(ntrkch/2)/nch;
       ch   = id%(ntrkch/2)%nch;
     }
     else{return false;}
  }
  else if(mod==16){
     all   = 1204;
     npln  =   18; 
     nch   =   32; 
     ntrkch = npln*nch*2-8*2;
     nvpln =    4; 
     nvch  =   17;
     if(id<all){
       if(id<ntrkch){
         view = id/(ntrkch/2); 
         tmp  = id%(ntrkch/2);
         if(tmp<24){
           pln = 0;
           ch  = tmp;
         }
         else{
           pln = (tmp-24)/nch+1;
           ch  = (tmp-24)%nch;
         }
       }
       else{
         pln  = (id-ntrkch)/nvch+npln;
         ch   = (id-ntrkch)%nvch;
         if(pln==18||pln==20) view=0; else view=1;
       }
     }
     else{return false;}
  }
  else{
    return false;
  }
  return true;
}



bool getid_channel(int mod,int view, int pln, int ch,int& id)
{
  int all   = 0;
  int npln  = 0;
  int nch   = 0;
  int nch2  = 0;
  int ntrkch = 0;
  int nvpln = 0;
  int nvch  = 0;
  int tmp;
  if(mod==3){
     all   = 594;
     npln  =  11; 
     nch   =  24; 
     ntrkch = npln*nch*2;
     nvpln =   3; 
     nvch  =  22;
     if(pln<npln){
       id = view*npln*nch + pln*nch + ch;
     }
     else{
       id = (pln-1-npln)*nvch + ch + ntrkch;
     }
     if(id<all){return true ;}
     else      {return false;}
  }
  else if(mod==14){
     all   = 616;
     npln  =  11; 
     nch   =  24; 
     ntrkch = npln*nch*2;
     nvpln =   4; 
     nvch  =  22;
     if(pln<npln){
       id = view*npln*nch + pln*nch + ch;
     }
     else{
       id = (pln-npln)*nvch + ch + ntrkch;
     }

     if(id<all){return true ;}
     else      {return false;}
  }
  else if(mod==15){
     all   = 1280;
     npln  =    8; 
     nch   =   80; 
     ntrkch = npln*nch*2;
     nvpln =   0; 
     nvch  =   0;
     id = view*npln*nch + pln*nch + ch;

     if(id<all){return true ;}
     else      {return false;}
  }
  else if(mod==16){
     all   = 1204;
     npln  =   18; 
     nch   =   32; 
     ntrkch = npln*nch*2-8*2;
     nvpln =    4; 
     nvch  =   17;
     if(pln<1){
       id = view*ntrkch/2 + ch;
     }
     else if(pln<npln){
       id = view*ntrkch/2 + (pln-1)*nch + ch + 24;
     }
     else{
       id = (pln-npln)*nvch + ch + ntrkch;
     }

     if(id<all){return true ;}
     else      {return false;}
  }
  else{
    return false;
  }
  return true;

}
