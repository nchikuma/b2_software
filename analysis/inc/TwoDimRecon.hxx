#ifndef __TWODIM_RECON_HXX__
#define __TWODIM_RECON_HXX__

//#define PMEDGE
//#define DEBUG_TWODIMRECON

//C++ libraly
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <deque>
#include <vector>
#include <sys/stat.h>
#include <unistd.h>
#include <complex>

//ROOT
#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TEventList.h>
#include <TBranch.h>
#include <TH1.h>
#include <TH1F.h>
#include <TApplication.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TH2F.h>
#include <TString.h>
#include <TSpectrum.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TMarker.h>
#include <TText.h>
#include <TMath.h>
#include <TBox.h>
#include <TLatex.h>
#include <TSystem.h>

//
#include "Const.hh"
#include "ReconConst.hh"
#include "DetectorDimension.hh"
#include "INGRID_BadCh_mapping.hh"

//INGRID library
#include "EVENTSUMMARY.h"

INGRID_BadCh_mapping *badch;
DetectorDimension *detdim;

using namespace std;

const int Cview = 2;
const int Cpln  = 200;
const int Cch   = 200;
const int Ccls  = 200;
const int Ccell = 5000;
const int Chit  = 5000;
const int Ctra  = 5000;
//const int Cdist = 4;
const int Cdist = 6;
const int Cmod  = 30;
const int Cdir  = 3;
const int MAX_CELLVALUE = 80;

const int    nactpln_threshold    = 3;
const int    nhit_threshold       = 5;
const int    nhit_threshold_WM    = 3;
const double layerpe_threshold    = 6.5;
double hitpe_threshold_ING  = 3.5; //3.5;
double hitpe_threshold_WM   = 2.5;
double hitpe_threshold_PM   = 3.5; //3.5;

const double wg_hittime_corr = 12550;


const double joint_th            = 25; //35.;
const int    celldist_max        = 6;
const double cellcenter_diffmax  = 120.;
const double cellcenter_diff1max = 150.;
const double neibor_th           = 3.;  //WM
const double neibor_th1          = 3.;  //PM/ING for DIST=DIST2=0
const double neibor_th2          = 1.5; //PM/ING
const double neibor_diffmax      = 60.;
const double upstpe_th           = 4.5;
const double upstpe_th1          = 1.5;
const double axis_border         = 60.; //deg

int   cTdcRsr = 50;  //nsec
float cPeCut  = 1.5; //18

//__________________________________________________________


class Hit{
  public:
    Int_t   id; //Between Hit class and IngridHitSummary
    Int_t   mod;
    Int_t   pln;
    Int_t   view;
    Int_t   ch;
    Int_t   used;
    Float_t pe; //p.e.
    Float_t pe_cross; //p.e.
    Float_t lope; //p.e.
    Long_t  tdc;  
    Long_t  time; //nsec
    double  posxy;
    double  posz;

    void clear(){
      id       =  -1;
      mod      =  -1;
      pln      =  -1;
      view     =  -1;
      ch       =  -1;
      used     =  -1;
      pe       =  -1e-5;
      pe_cross =  -1e-5;
      lope     =  -1e-5;
      tdc      =  -1;
      time     =  -1;
      posxy    =  -1.e-5;
      posz     =  -1.e-5;
    }
};

vector<Hit> allhit;
vector<Hit> allhit_for_disp;
vector<Hit> hitcls;
vector<Hit> hitcls_for_joint;

class Track{
  public:
    Int_t   mod;
    Int_t   view;
    Int_t   ipln;
    Int_t   fpln;
    Float_t ixy;
    Float_t fxy;
    Float_t iz;
    Float_t fz;
    Float_t slope;
    Float_t intcpt;
    Float_t ang;
    Float_t clstime;
    Bool_t  veto;
    Bool_t  edge;
    Bool_t  stop;
    Float_t vetodist;
    vector<Int_t>  hitid;
    vector<Bool_t> isohit;

    void clear(){
      mod      =  -1;
      view     =  -1;
      ipln     =  -1;
      fpln     =  -1;
      ixy      =  -1.e-5;
      fxy      =  -1.e-5;
      iz       =  -1.e-5;
      fz       =  -1.e-5;
      slope    =  -1.e-5;
      intcpt   =  -1.e-5;
      ang      =  -1.e-5;
      clstime  =  -1.e-5;
      veto     = false;
      edge     = false;
      stop     = false;
      vetodist =  -1.e-5;
      hitid.clear();
      isohit.clear();
    }
};

vector<Track> alltrack;


int chmax(int mod, int view, int pln, int axis=0){
  int maxch = detdim->GetChMax(mod, view, pln, axis);
  return maxch;
}


int plnmax(int mod, int view=0, int pln=0, int axis=0){
  int maxpln = detdim->GetPlnMax( mod, view, pln, axis);
  return maxpln;
}


float zposi(int mod,int view,int pln, int ch=0, int axis=0){
  if(mod<0){return -1e+5;}

  double posiz, posixy, x_tmp, y_tmp;

  if(mod==MOD_ONAXIS_WM || mod==MOD_B2_WM){
    detdim->GetPos_TwoDimRecon(mod,view,pln,ch,axis,&posixy,&posiz);
    if(axis==0) posiz = posiz;
    else        posiz = posixy;
  }
  else{
    detdim->GetPosInMod(mod,pln,view,ch,&x_tmp,&y_tmp,&posiz);
    posiz = posiz;
  }

  return posiz;
}

float xyposi(int mod, int view, int pln,int ch, int axis=0){
  if(mod<0){return -1e+5;}

  double posixy, posiz, x_tmp, y_tmp;

  if(mod==MOD_ONAXIS_WM || mod==MOD_B2_WM){
    detdim->GetPos_TwoDimRecon(mod,view,pln,ch,axis,&posixy,&posiz);
    if(axis==0) posixy = posixy;
    else        posixy = posiz;
  }
  else{
    detdim->GetPosInMod(mod,pln,view,ch,&x_tmp,&y_tmp,&posiz);
    if(view==SideView) posixy = y_tmp;
    else               posixy = x_tmp;
  }

  return posixy;
}

//Veto
double vposiz(int mod,int ch){
  double posiz, x_tmp, y_tmp;

  if(mod==MOD_PM || mod==MOD_B2_CH){
    detdim->GetPosPM(19,0,ch,&x_tmp,&y_tmp,&posiz);
  }
  else if(mod==MOD_ONAXIS_WM || mod==MOD_B2_WM){
    //Temporarily the same as PM
    detdim->GetPosPM(19,0,ch,&x_tmp,&y_tmp,&posiz);
  }
  else if(mod<NUMINGMOD || mod==MOD_B2_INGRID){
    detdim->GetPosING(0,12,0,ch,&x_tmp,&y_tmp,&posiz);
  }
  else{
    posiz = 0.;
  }
  return posiz;
}

//Veto
double vposixy(int mod,int pln){
  double posixy, x_tmp, y_tmp, z_tmp;

  int view;
  if(mod<NUMINGMOD || mod==MOD_B2_INGRID){
    if(pln==11 || pln==12) view=1;
    else                   view=0;

    detdim->GetPosING(mod,pln,view,0,&x_tmp,&y_tmp,&z_tmp);

    if(view==0) posixy = y_tmp;
    else        posixy = x_tmp;
  }
  else if(mod==MOD_PM || mod==MOD_B2_CH){
    if(pln==18||pln==20) view=0;
    else                 view=1;
    detdim->GetPosPM(pln,view,0,&x_tmp,&y_tmp,&z_tmp);

    if(view==0) posixy = y_tmp;
    else        posixy = x_tmp;
  }
  else if(mod==MOD_ONAXIS_WM || mod==MOD_B2_WM){
    //Temporarily the same as PM
    if(pln==18 || pln==20) view=0;
    else                   view=1;
    detdim->GetPosPM(pln,view,0,&x_tmp,&y_tmp,&z_tmp);

    if(view==0) posixy = y_tmp;
    else        posixy = x_tmp;
  }
  else posixy = 0.;

  return posixy;
}

float sciwidth(int mod, int view, int pln, int ch=0, int axis=0, int flag=0){
  float sciw = detdim->GetScintiWidth(mod,view,pln,ch,axis);

  if(flag==1) sciw = 0.5;

  return sciw;
};


float scithick(int mod, int view, int pln, float xy, int axis=0, int flag=0){
  float scith;
  int ch;

  if(mod==MOD_PM && pln!=0){
    ch = (int) xy/C_INGScintiWidth;
    if(ch<8)       ch = ch;
    else if(ch<16) ch = ch*2-8;
    else           ch = ch+8;
  }
  else{
    ch = 0;
  }

  scith = detdim->GetScintiThick(mod,view,pln,ch,axis);

  if(flag==1) scith = 0.5;

  return scith;
};


float Yerr(int mod,int view, int pln,int ch, int axis, int edge){
  float yerr;

  if(edge==0) yerr = -sciwidth(mod,view,pln,ch,axis);
  else        yerr =  sciwidth(mod,view,pln,ch,axis);

  // For grid scintillator in axis==0
  // and plane scinti in axis==1
  double grid_err = 12.5;
  if(mod==MOD_B2_WM || mod==MOD_ONAXIS_WM){
    if(axis==0){
      if((view==0&&pln%3!=0)||(view==1&&pln%3!=1)){
        if(edge==0) yerr = -grid_err;
        else        yerr =  grid_err;
      }
    }
    else{
      if(pln%3!=1){
        if(edge==0) yerr = -grid_err;
        else        yerr =  grid_err;
      }
    }
  }


  return yerr;
};


bool actpln[Cview][Cpln];
int  nactpln;
int fNactpln(int mod){
  nactpln = 0;
  memset(actpln,false,sizeof(actpln));
  for(int i=0; i<(int)hitcls.size(); i++){
    if(hitcls[i].pln<plnmax(mod, 0, 0, 0)){
      actpln[hitcls[i].view][hitcls[i].pln] = true;
    }
  }

  for(int i=0; i<plnmax(mod, 0, 0, 0); i++){
    if(actpln[0][i] && actpln[1][i]){
      nactpln++;
    }
    else{
      actpln[0][i] = false;
      actpln[1][i] = false;
    }
  }
  return nactpln;
};

float fLayerpe(int mod){
  float layerpe = 0;
  for(int i=0; i<(int)hitcls.size(); i++){
    if(hitcls[i].pln<plnmax(mod, 0, 0, 0)){
      if(actpln[hitcls[i].view][hitcls[i].pln]){
        layerpe += hitcls[i].pe;
      }
    }
  }

  if(nactpln==0) layerpe = 0;
  else           layerpe = layerpe/nactpln;

  return layerpe;
};


bool withtime(const Hit& left, const Hit& right){
  return left.time < right.time;
};

void fSortTime(vector<Hit> &a){
  std::stable_sort(a.begin(), a.end(), withtime);
};


bool fFindTimeClster(vector<Hit> &hit, vector<Hit> &hitclster, 
    Long_t &ctime, bool cosmic=false){

  Int_t nhit = hit.size();
  int maxhit = 5; //Maximum for clustered hits
  int hitmod = hit[0].mod;
  if(hitmod==MOD_B2_WM){ cTdcRsr = 580; }
  else                 { cTdcRsr =  50; }

  if(nhit<=maxhit) return false;

  //For neutrino beam measurement
  if(!cosmic){
    for(Int_t i=0; i<nhit-maxhit; i++){
      if(fabs( hit[i].time - hit[i+maxhit].time) < cTdcRsr){
        long  basetime = 0.;    
        float highpe   = 0.;
        float sumpe    = 0.;

        //Extract largest pe event to be used for basetime
        for(int j=0; j<maxhit+1; j++){
          if(hit[i+j].pe > highpe){
            basetime = hit[i+j].time;
            highpe   = hit[i+j].pe;
            sumpe    = sumpe + hit[i+j].pe;
          }
        }
        if(hit[0].mod!=MOD_ONAXIS_WM && hit[0].mod!=MOD_B2_WM){
          if(sumpe<cPeCut) continue;
        }

        vector<Hit>::iterator it;
        Int_t ncount=0;
        for(it=hit.begin(); it!=hit.end(); it++){
          if(fabs( basetime - it->time) < cTdcRsr){
            ncount++;
          }
        }
        if(ncount<=maxhit) continue;

        //Clstering
        hitcls.clear(); //Reset hit clster
        for(it=hit.begin(); it!=hit.end(); it++){
          if(fabs( basetime - it->time) < cTdcRsr){
            hitclster.push_back(*it);
            it = hit.erase(it);
            it--;
          }
        }

        //Caluculate time of clster with maximum p.e.
        highpe = 0.;
        for(int j=0; j<(int)hitclster.size(); j++){
          if(highpe<hitclster[j].pe){
            highpe = hitclster[j].pe;
            ctime  = hitclster[j].time;
          }
        }

        return true;
      }
    }
    return false;
  }
  //For cosmic-rays
  else{
    //Clstering
    vector<Hit>::iterator it;
    hitcls.clear(); //Reset hit clster
    for(it=hit.begin(); it!=hit.end(); it++){
      hitclster.push_back(*it);
      it = hit.erase(it);
      it--;
    }

    //Caluculate time of clster with maximum p.e.
    float highpe = 0.;
    for(int j=0; j<(int)hitclster.size(); j++){
      if(highpe < hitclster[j].pe){
        highpe = hitclster[j].pe;
        ctime  = hitclster[j].time;
      }
    }
    return true;
  }
}


int   view,id[Cview][Chit], pln[Cview][Chit], ch[Cview][Chit];
int   hit[Cview], hitnum[Cview][Cpln][Cch], used[Cview][Chit];
float pe[Cview][Chit];
int   PLN, PLN2, CH;
int   CL, CL2, CLHIT, CLHIT2, CLHIT3;
int   CELL, CELL2, NEI;
int   TRA, TRA2, TRACELL, TRACL, TRACL2;
int   HIT, DIST, DIST2, DIST3, dummy, TMP;
bool  hitcl[Cview][Cpln][Ccls];
int   clchi[Cview][Cpln][Ccls], clchf[Cview][Cpln][Ccls], ncl[Cview][Cpln], numcl[Cview][Cpln][Ccls], clhit[Cview][Cpln][Ccls][Cch];
float clcenter[Cview][Cpln][Ccls], clpe[Cview][Cpln][Ccls];
int   clused[Cview][Cpln][Ccls];
int   cellu[Cview][Cpln][Ccell][Cdist], celld[Cview][Cpln][Ccell][Cdist], ncell[Cview][Cpln][Cdist];
int   value[Cview][Cpln][Ccell][Cdist], nvalue[Cview][Cpln][Ccell][Cdist];
bool  neibor[Cview][Cpln][Ccell][Cdist][Cdist];
int   neiu[Cview][Cpln][Ccell][Cdist][Cdist], neid[Cview][Cpln][Ccell][Cdist][Cdist], nnei[Cview][Cpln][Cdist][Cdist];
int   track_cell[Cview][Ccell][Cpln], track_pln[Cview][Ccell][Cpln], track_dist[Cview][Ccell][Cpln];
int   ntracell[Cview][Ccell], ntracl[Cview][Ccell], ntracl2[Cview][Ccell], ntrack[Cview], ntrack2[Cview], ntrack3[Cview];
int   trank[Cview][Cpln][Cpln], ncltra[Cview][Cpln], rank[Cview][Ccell];
bool  ttrack[Cview][Ccell], ttrack2[Cview][Ccell];
int   plane[Cview][Ctra][Cpln], clus[Cview][Ctra][Cpln];
bool  vetowtracking[Cview][Ctra], edgewtracking[Cview][Ctra], track_stop[Cview][Ctra];
float vetodist[Cview][Ctra];
bool  ovcl[Cview][Cpln][Ccls];
float dis;
float XX[Cview][Ctra][Cpln], YY[Cview][Ctra][Cpln], Xe[Cview][Ctra][Cpln], Yel[Cview][Ctra][Cpln], Yeh[Cview][Ctra][Cpln];
float x[Cdir], y[Cdir], xe[Cdir], yel[Cdir], yeh[Cdir];
float chi2, tmp;
float par[Cview][Ctra][2];
float trape[Cview][Ctra];
float ulimit;
bool  upst, dwst;
int   nhit;
int   ue, shita;
int   hit_n[Cview][Chit], hit_id[Cview][Chit][Chit];
bool  isohit[Cview][Ctra][Chit];



bool fTracking(int mod, int axis=0, int N_long=1){

#ifdef DEBUG_TWODIMRECON
  if(mod!=14){return false;}

  cout
    << endl
    << " =============================================="
    << endl
    << " fTracking :"
    << " mod:" << mod
    << " axis:" << axis
    << " N_long:" << N_long
    << endl;

#endif

  memset(hit   , 0,sizeof(hit   ));
  memset(hitnum,-1,sizeof(hitnum));
  if(mod!=MOD_ONAXIS_WM && mod!=MOD_B2_WM){
    alltrack.clear();
  }
  int hitclssize = hitcls.size();
  for(int i=0; i<hitclssize; i++){
    view = hitcls[i].view;
    id  [view][hit[view]]   = hitcls[i].id;
    pln [view][hit[view]]   = hitcls[i].pln;
    ch  [view][hit[view]]   = hitcls[i].ch;
    pe  [view][hit[view]]   = hitcls[i].pe;
    used[view][hit[view]]   = hitcls[i].used;
    //pdg[view][hit[view]]  = hitcls[i].pdg;
    if(pln[view][hit[view]]<plnmax(mod,view,pln[view][hit[view]],axis)){
      hitnum[view][pln[view][hit[view]]][ch[view][hit[view]]] = hit[view];
    }
    hit[view]++;
  }

#ifdef DEBUG_TWODIMRECON
  for(int VIEW=0; VIEW<2; VIEW++){
#else
  for(int VIEW=0; VIEW<2; VIEW++){
#endif

    //===============
    //Define clusters
#ifdef DEBUG_TWODIMRECON
    cout << "---- Define clusters ----------- " <<endl;
#endif 

    memset(ncl  [VIEW],0,sizeof(ncl  [VIEW]));
    memset(numcl[VIEW],0,sizeof(numcl[VIEW]));
    memset(clhit[VIEW],0,sizeof(clhit[VIEW]));
    for(PLN=0; PLN<plnmax(mod,VIEW,PLN,axis); PLN++){
      CH = 0;

      while(1){
        if(hitnum[VIEW][PLN][CH]>=0){
          clchi[VIEW][PLN][ncl[VIEW][PLN]] = CH;
          clhit[VIEW][PLN][ncl[VIEW][PLN]][numcl[VIEW][PLN][ncl[VIEW][PLN]]] = hitnum[VIEW][PLN][CH];
          while(1){
            numcl[VIEW][PLN][ncl[VIEW][PLN]]++;

            if(((mod==MOD_ONAXIS_WM || mod==MOD_B2_WM) && !((axis==0&&view==0&&PLN%3==0) || (axis==0&&view==1&&PLN%3==1))))
            { break; }
            if(CH==chmax(mod,VIEW,PLN,axis)-1) break;
            if(hitnum[VIEW][PLN][CH+1]<0)      break;
            CH++;
            clhit[VIEW][PLN][ncl[VIEW][PLN]][numcl[VIEW][PLN][ncl[VIEW][PLN]]] = hitnum[VIEW][PLN][CH];
          }
          clchf[VIEW][PLN][ncl[VIEW][PLN]] = CH;
          ncl[VIEW][PLN]++;
        }
        if(CH==chmax(mod,VIEW,PLN,axis)-1) break;
        CH++;
      }
    }

    //==============================
    //Fill pe and center of clusters
#ifdef DEBUG_TWODIMRECON
    cout << "---- Fill pe and center of clusters ----------- " <<endl;
#endif 
    memset(clpe[VIEW]    ,0,sizeof(clpe[VIEW]));
    memset(clcenter[VIEW],0,sizeof(clcenter[VIEW]));
    for(PLN=0; PLN<plnmax(mod,VIEW,PLN,axis); PLN++){
      for(CL=0; CL<ncl[VIEW][PLN]; CL++){
        for(CLHIT=0; CLHIT<numcl[VIEW][PLN][CL]; CLHIT++){
          clpe[VIEW][PLN][CL] += pe[VIEW][clhit[VIEW][PLN][CL][CLHIT]];
          clcenter[VIEW][PLN][CL] += pe[VIEW][clhit[VIEW][PLN][CL][CLHIT]]*xyposi(mod,VIEW,PLN,clchi[VIEW][PLN][CL]+CLHIT,axis);
        }
        clcenter[VIEW][PLN][CL] = clcenter[VIEW][PLN][CL]/clpe[VIEW][PLN][CL];
      }
    }

    //==================
    //Used clusters
#ifdef DEBUG_TWODIMRECON
    cout << "---- Used clusters ----------- " <<endl;
#endif 
    memset(clused[VIEW],0,sizeof(clused[VIEW]));
    for(PLN=0; PLN<plnmax(mod,VIEW,PLN,axis); PLN++){
      for(CL=0; CL<ncl[VIEW][PLN]; CL++){
        for(CLHIT=0; CLHIT<numcl[VIEW][PLN][CL]; CLHIT++){
          if(used[VIEW][clhit[VIEW][PLN][CL][CLHIT]]==1){
            clused[VIEW][PLN][CL] = 1;
          }
        }
      }
    }

    //==================
    //Define cells
#ifdef DEBUG_TWODIMRECON
    cout << "---- Define cells ----------- " <<endl;
#endif 

    memset(ncell[VIEW],0,sizeof(ncell[VIEW]));
    for(PLN=0; PLN<plnmax(mod,VIEW,PLN,axis)-1; PLN++){
      for(CL=0; CL<ncl[VIEW][PLN]; CL++){
        for(DIST=0; DIST<celldist_max; DIST++){
          if(mod!=MOD_ONAXIS_WM && mod!=MOD_B2_WM){
           if(DIST>1){ continue;}
          }
          else{
            if     (DIST==4&&axis==0&&VIEW==1&&PLN%3==0){continue;}
            else if(DIST==4&&axis==0&&VIEW==0&&PLN%3==2){continue;}
            else if(DIST==4&&axis==1&&         PLN%3==0){continue;}
            else if(DIST==5&&axis==0&&VIEW==1&&PLN%3!=1){continue;}
            else if(DIST==5&&axis==0&&VIEW==0&&PLN%3!=0){continue;}
            else if(DIST==5&&axis==1&&         PLN%3!=1){continue;}
          }
          if(PLN+DIST>plnmax(mod,VIEW,PLN,axis)-2){ continue; }

          for(CL2=0; CL2<ncl[VIEW][PLN+DIST+1]; CL2++){
            if((mod==MOD_ONAXIS_WM || mod==MOD_B2_WM) && 
                fabs( clcenter[VIEW][PLN][CL]-clcenter[VIEW][PLN+DIST+1][CL2])>cellcenter_diffmax)
            { 
              continue;
            }
            else if(DIST==1 &&
                fabs( clcenter[VIEW][PLN][CL] - clcenter[VIEW][PLN+2][CL2] )>cellcenter_diff1max)
            {
              continue;
            }
            if(clused[VIEW][PLN][CL]==1 || clused[VIEW][PLN+DIST+1][CL2]==1){ continue; }
            cellu[VIEW][PLN][ncell[VIEW][PLN][DIST]][DIST] = CL;
            celld[VIEW][PLN][ncell[VIEW][PLN][DIST]][DIST] = CL2;
            ncell[VIEW][PLN][DIST]++;
#ifdef DEBUG_TWODIMRECON
            cout << "view=" << VIEW << " pln=" << PLN << " dist=" << DIST << endl;
#endif
          }
        }
      }
    }

    //===================
    //Define neiborhoods
#ifdef DEBUG_TWODIMRECON
    cout << "---- Define neiborhoods ----------- " <<endl;
#endif 
    memset(nnei[VIEW],0,sizeof(nnei[VIEW]));
    for(PLN=0; PLN<plnmax(mod,VIEW,PLN,axis)-2; PLN++){
      for(DIST=0; DIST<celldist_max; DIST++){
        if(PLN-DIST<0){ continue; }

        for(CELL=0; CELL<ncell[VIEW][PLN-DIST][DIST]; CELL++){
          for(DIST2=0; DIST2<celldist_max; DIST2++){
            if(PLN+DIST2>plnmax(mod,VIEW,PLN,axis)-3){ continue;}

            if((mod!=MOD_ONAXIS_WM && mod!=MOD_B2_WM) && DIST==1&&DIST2==1){ continue; }
            if((mod!=MOD_ONAXIS_WM && mod!=MOD_B2_WM) && DIST>1           ){ continue; }
            if((mod!=MOD_ONAXIS_WM && mod!=MOD_B2_WM) && DIST2>1          ){ continue; }

            for(CELL2=0; CELL2<ncell[VIEW][PLN+1][DIST2]; CELL2++){
              if(celld[VIEW][PLN-DIST][CELL][DIST]==cellu[VIEW][PLN+1][CELL2][DIST2]){
                x[0] = zposi(mod,VIEW,PLN-DIST,0,axis);
                x[1] = zposi(mod,VIEW,PLN+1,0,axis);
                x[2] = zposi(mod,VIEW,PLN+2+DIST2,0,axis);

                y[0] = clcenter[VIEW][PLN-DIST][cellu[VIEW][PLN-DIST][CELL][DIST]];
                y[1] = clcenter[VIEW][PLN+1][celld[VIEW][PLN-DIST][CELL][DIST]];
                y[2] = clcenter[VIEW][PLN+2+DIST2][celld[VIEW][PLN+1][CELL2][DIST2]];
                if(DIST>2 && DIST2>2 && (fabs(y[2]-y[0])>neibor_diffmax)){ continue; }

                xe[0] = scithick(mod,VIEW,PLN-DIST,y[0],axis,0);
                xe[1] = scithick(mod,VIEW,PLN+1,y[1],axis,0);
                xe[2] = scithick(mod,VIEW,PLN+2+DIST2,y[2],axis,0);

                yel[0] = -Yerr(mod,VIEW,PLN-DIST   ,clchi[VIEW][PLN-DIST   ][cellu[VIEW][PLN-DIST][CELL ][DIST ]],axis,0);
                yel[1] = -Yerr(mod,VIEW,PLN+1      ,clchi[VIEW][PLN+1      ][celld[VIEW][PLN-DIST][CELL ][DIST ]],axis,0);
                yel[2] = -Yerr(mod,VIEW,PLN+2+DIST2,clchi[VIEW][PLN+2+DIST2][celld[VIEW][PLN+1   ][CELL2][DIST2]],axis,0);
                yeh[0] =  Yerr(mod,VIEW,PLN-DIST   ,clchf[VIEW][PLN-DIST   ][cellu[VIEW][PLN-DIST][CELL ][DIST ]],axis,1);
                yeh[1] =  Yerr(mod,VIEW,PLN+1      ,clchf[VIEW][PLN+1      ][celld[VIEW][PLN-DIST][CELL ][DIST ]],axis,1);
                yeh[2] =  Yerr(mod,VIEW,PLN+2+DIST2,clchf[VIEW][PLN+2+DIST2][celld[VIEW][PLN+1   ][CELL2][DIST2]],axis,1);

                TGraphAsymmErrors *graph = new TGraphAsymmErrors(3,x,y,xe,xe,yel,yeh);
                TF1 *f = new TF1("f","[0]+[1]*x");
                f->SetParameters(y[0]-x[0]*(y[2]-y[0])/(x[2]-x[0]),(y[2]-y[0])/(x[2]-x[0]));

                graph->Fit("f","Q");
                chi2 = f->GetChisquare();
#ifdef DEBUG_TWODIMRECON
                double fitpara0 = f->GetParameter(0);
                double fitpara1 = f->GetParameter(1);
                double rms = 0.;
                double xx1 = x[1]-x[0];
                double xx2 = x[2]-x[1];
                double yy1 = y[1]-y[0];
                double yy2 = y[2]-y[1];
                double cos = (yy1*yy2+xx1*xx2)/sqrt(xx1*xx1+yy1*yy1)/sqrt(xx2*xx2+yy2*yy2);
                for(int i=0;i<3;i++){rms += pow(y[i]-fitpara0-x[i]*fitpara1,2);}
                rms = sqrt(rms);
                cout
                  << " view=" << VIEW << " axis=" << axis
                  << " pln={" << PLN-DIST<<","<< PLN+1<<","<<PLN+2+DIST2<<"}"
                  << " x={"   << x[0]<<","<<x[1]<<","<<x[2]<<"}"
                  << " y={"   << y[0]<<","<<y[1]<<","<<y[2]<<"}"
                  << " chi2=" << chi2
                  << " rms="  << rms
                  << " cos="  << cos;
#endif
                graph->Delete();
                f->Delete();

                if(mod==MOD_ONAXIS_WM || mod==MOD_B2_WM){
                  ulimit = neibor_th;
                }
                else{
                  if(DIST==0&&DIST2==0) ulimit = neibor_th1;
                  else                  ulimit = neibor_th2;
                }

                if(chi2 < ulimit){
                //if(rms < 100.){
                //if(cos > 0.60){
                  neiu[VIEW][PLN][nnei[VIEW][PLN][DIST][DIST2]][DIST][DIST2] = CELL;
                  neid[VIEW][PLN][nnei[VIEW][PLN][DIST][DIST2]][DIST][DIST2] = CELL2;
                  nnei[VIEW][PLN][DIST][DIST2]++;
#ifdef DEBUG_TWODIMRECON
                  cout <<  " OK ";
#endif
                }
#ifdef DEBUG_TWODIMRECON
                  cout <<  endl;
#endif

              }
            }
          }
        }
      }
    }

    //===========================
    //Define value of cells
#ifdef DEBUG_TWODIMRECON
    cout << "----------- Define value of cells -----------------" << endl;
#endif
    memset(value [VIEW],0,sizeof(value [VIEW]));
    memset(nvalue[VIEW],0,sizeof(nvalue[VIEW]));
    for(int jndex=0; jndex<MAX_CELLVALUE; jndex++){
      for(PLN=0; PLN<plnmax(mod,VIEW,PLN,axis)-2; PLN++){
        for(DIST=0; DIST<celldist_max; DIST++){
          for(DIST2=0; DIST2<celldist_max; DIST2++){
            if(PLN-DIST <0                          ){ continue; }
            if(PLN+DIST2>plnmax(mod,VIEW,PLN,axis)-3){ continue; }

            for(NEI=0; NEI<nnei[VIEW][PLN][DIST][DIST2]; NEI++){
              if(value   [VIEW][PLN-DIST][neiu[VIEW][PLN][NEI][DIST][DIST2]][DIST ]
                  ==value[VIEW][PLN+1   ][neid[VIEW][PLN][NEI][DIST][DIST2]][DIST2])
              {
                nvalue   [VIEW][PLN+1][neid[VIEW][PLN][NEI][DIST][DIST2]][DIST2]
                  = value[VIEW][PLN+1][neid[VIEW][PLN][NEI][DIST][DIST2]][DIST2]+1;
              }
            }
          }
        }
      }
      for(PLN=0; PLN<plnmax(mod,VIEW,PLN,axis)-1; PLN++){
        for(DIST=0; DIST<celldist_max; DIST++){
          if(PLN+DIST>plnmax(mod,VIEW,PLN,axis)-2){ continue;}
          for(CELL=0; CELL<ncell[VIEW][PLN][DIST]; CELL++){
            value[VIEW][PLN][CELL][DIST] = nvalue[VIEW][PLN][CELL][DIST];
          }
        }
      }
    }
#ifdef DEBUG_TWODIMRECON
    cout 
      << " view=" << VIEW
      << " axis=" << axis
      << endl;
    for(PLN=0; PLN<plnmax(mod,VIEW,PLN,axis)-1; PLN++){
      for(DIST=0; DIST<celldist_max; DIST++){
        if(PLN+DIST>plnmax(mod,VIEW,PLN,axis)-2){ continue;}
        for(CELL=0; CELL<ncell[VIEW][PLN][DIST]; CELL++){
          cout 
            << " {" << PLN
            << ", " << DIST
            << ", " << value[VIEW][PLN][CELL][DIST]
            << "}";
        }
      }
    }
    cout << endl;
#endif


    //==============
    //Define tracks
#ifdef DEBUG_TWODIMRECON
    cout << "---- Define tracks -----------" <<endl;
#endif 
    ntrack[VIEW] = 0;
    memset(ntracell[VIEW],   0,sizeof(ntracell[VIEW]));
    memset(trape   [VIEW],   0,sizeof(trape   [VIEW]));
    memset(neibor  [VIEW],true,sizeof(neibor  [VIEW]));
    //for(PLN=1; PLN<plnmax(mod,VIEW,PLN,axis)-1; PLN++){ 
    for(PLN=plnmax(mod,VIEW,PLN,axis)-2; PLN>=1; PLN--){
      for(DIST=0; DIST<celldist_max; DIST++){
        if(PLN+DIST>plnmax(mod,VIEW,PLN,axis)-2){ continue; }

        for(CELL=0; CELL<ncell[VIEW][PLN][DIST]; CELL++){
          if(value[VIEW][PLN][CELL][DIST]>0){
            if(PLN+DIST==plnmax(mod,VIEW,PLN,axis)-2){
              upst = true;
            }
            else{
              upst = true;
              for(DIST2=0; DIST2<celldist_max; DIST2++){
                if((mod!=MOD_ONAXIS_WM && mod!=MOD_B2_WM)&&DIST==1&&DIST2==1){ continue; }
                if((mod!=MOD_ONAXIS_WM && mod!=MOD_B2_WM)&&DIST >1          ){ continue; }
                if((mod!=MOD_ONAXIS_WM && mod!=MOD_B2_WM)&&DIST2>1          ){ continue; }

                for(NEI=0; NEI<nnei[VIEW][PLN+DIST][DIST][DIST2]; NEI++){
                  if(CELL==neiu[VIEW][PLN+DIST][NEI][DIST][DIST2]){ upst = false; }
                  if(CELL==neiu[VIEW][PLN+DIST][NEI][DIST][DIST2] 
                      && neibor[VIEW][PLN+DIST][NEI][DIST][DIST2])
                  {
                    upst = true;
                    break;
                  }
                }
                if(!upst){break;}
              }
            }

            if(upst){
              track_pln [VIEW][ntrack[VIEW]][0] = PLN;
              track_cell[VIEW][ntrack[VIEW]][0] = CELL;
              track_dist[VIEW][ntrack[VIEW]][0] = DIST;

              //1st cluster info
              plane[VIEW][ntrack[VIEW]][0] = track_pln[VIEW][ntrack[VIEW]][0]+track_dist[VIEW][ntrack[VIEW]][0]+1;
              clus [VIEW][ntrack[VIEW]][0]  
                = celld[VIEW]
                       [track_pln [VIEW][ntrack[VIEW]][0]]
                       [track_cell[VIEW][ntrack[VIEW]][0]]
                       [track_dist[VIEW][ntrack[VIEW]][0]];
              XX   [VIEW][ntrack[VIEW]][0] = zposi(mod,VIEW,plane[VIEW][ntrack[VIEW]][0],0,axis);
              YY   [VIEW][ntrack[VIEW]][0] = clcenter[VIEW][plane[VIEW][ntrack[VIEW]][0]][clus[VIEW][ntrack[VIEW]][0]];
              Yel  [VIEW][ntrack[VIEW]][0] = -Yerr(mod,VIEW,plane[VIEW][ntrack[VIEW]][0],
                                                   clchi[VIEW]
                                                        [plane[VIEW][ntrack[VIEW]][0]]
                                                        [clus [VIEW][ntrack[VIEW]][0]],
                                                   axis,0);
              Yeh  [VIEW][ntrack[VIEW]][0] = +Yerr(mod,VIEW,plane[VIEW][ntrack[VIEW]][0],
                                                   clchf[VIEW]
                                                        [plane[VIEW][ntrack[VIEW]][0]]
                                                        [clus [VIEW][ntrack[VIEW]][0]],
                                                   axis,1);
              Xe   [VIEW][ntrack[VIEW]][0] = scithick(mod,VIEW,plane[VIEW][ntrack[VIEW]][0],YY[VIEW][ntrack[VIEW]][0],axis,0);
              trape[VIEW][ntrack[VIEW]]   += clpe[VIEW][plane[VIEW][ntrack[VIEW]][0]][clus[VIEW][ntrack[VIEW]][0]];

              //2nd cluster info
              plane[VIEW][ntrack[VIEW]][1] =  track_pln[VIEW][ntrack[VIEW]][0];
              clus [VIEW][ntrack[VIEW]][1] 
                =  cellu[VIEW]
                        [track_pln [VIEW][ntrack[VIEW]][0]]
                        [track_cell[VIEW][ntrack[VIEW]][0]]
                        [track_dist[VIEW][ntrack[VIEW]][0]];
              XX   [VIEW][ntrack[VIEW]][1] = zposi(mod,VIEW,plane[VIEW][ntrack[VIEW]][1],0,axis);
              YY   [VIEW][ntrack[VIEW]][1] = clcenter[VIEW][plane[VIEW][ntrack[VIEW]][1]][clus[VIEW][ntrack[VIEW]][1]];
              Yel  [VIEW][ntrack[VIEW]][1] = -Yerr(mod,VIEW,plane[VIEW][ntrack[VIEW]][1],
                                                        clchi[VIEW]
                                                             [plane[VIEW][ntrack[VIEW]][1]]
                                                             [clus [VIEW][ntrack[VIEW]][1]],
                                                        axis,0);
              Yeh  [VIEW][ntrack[VIEW]][1] = +Yerr(mod,VIEW,plane[VIEW][ntrack[VIEW]][1],
                                                 clchf[VIEW]
                                                      [plane[VIEW][ntrack[VIEW]][1]]
                                                      [clus [VIEW][ntrack[VIEW]][1]],
                                                 axis,1);
              Xe   [VIEW ][ntrack[VIEW]][1] = scithick(mod,VIEW,plane[VIEW][ntrack[VIEW]][1],YY[VIEW][ntrack[VIEW]][1],axis,0);
              trape[VIEW][ntrack[VIEW]]    += clpe[VIEW][plane[VIEW][ntrack[VIEW]][1]][clus[VIEW][ntrack[VIEW]][1]];

              ntracell[VIEW][ntrack[VIEW]] = 1;
              PLN2 = PLN-1;
              double chi2old = 3.;
              if(mod==MOD_ONAXIS_WM || mod==MOD_B2_WM){
                chi2old=neibor_th;
              }

              while(PLN2>=0){
#ifdef DEBUG_TWODIMRECON
                cout << "   PLN2=" << PLN2 << endl;
#endif
                dwst = true;
                DIST3 = track_dist[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]-1];
                for(DIST2=0; DIST2<celldist_max; DIST2++){
                  if(PLN2-DIST2<0){ continue; }

                  for(TMP=0; TMP<2; TMP++){
                    for(NEI=0; NEI<nnei[VIEW][PLN2][DIST2][DIST3]; NEI++){
                      if( neid[VIEW][PLN2][NEI][DIST2][DIST3]
                          ==track_cell[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]-1])
                      {
                        if(   value[VIEW][PLN2-DIST2][neiu[VIEW][PLN2][NEI][DIST2][DIST3]][DIST2]+1
                            ==value[VIEW][PLN2+1    ][neid[VIEW][PLN2][NEI][DIST2][DIST3]][DIST3] 
                            || TMP!=0)
                        {
                          if(!dwst)
                          {
                            float Xetmp  = Xe [VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1];
                            float YYtmp  = YY [VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1];
                            float Yeltmp = Yel[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1];
                            float Yehtmp = Yeh[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1];

                            YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]
                              = clcenter[VIEW][PLN2-DIST2][cellu[VIEW][PLN2-DIST2][neiu[VIEW][PLN2][NEI][DIST2][DIST3]][DIST2]];
                            Yel[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]
                              = -Yerr(mod,VIEW,PLN2-DIST2,
                                  clchi[VIEW][PLN2-DIST2][cellu[VIEW][PLN2-DIST2][neiu[VIEW][PLN2][NEI][DIST2][DIST3]][DIST2]],
                                  axis,0);
                            Yeh[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]
                              = +Yerr(mod,VIEW,PLN2-DIST2,
                                  clchf[VIEW][PLN2-DIST2][cellu[VIEW][PLN2-DIST2][neiu[VIEW][PLN2][NEI][DIST2][DIST3]][DIST2]],
                                  axis,1);
                            Xe[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]
                              = scithick(mod,VIEW,PLN2-DIST2,YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1],axis,0);

                            TGraphAsymmErrors *graph2 = new TGraphAsymmErrors((ntracell[VIEW][ntrack[VIEW]]+2),
                                XX [VIEW][ntrack[VIEW]],
                                YY [VIEW][ntrack[VIEW]],
                                Xe [VIEW][ntrack[VIEW]],
                                Xe [VIEW][ntrack[VIEW]],
                                Yel[VIEW][ntrack[VIEW]],
                                Yeh[VIEW][ntrack[VIEW]]);
                            TF1 *f2 = new TF1("f2","[0]+[1]*x");
                            f2->SetParameters(YY[VIEW][ntrack[VIEW]][0] 
                                - XX[VIEW][ntrack[VIEW]][0]
                                *(YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1] - YY[VIEW][ntrack[VIEW]][0])
                                /(XX[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1] - XX[VIEW][ntrack[VIEW]][0]),
                                 (YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1] - YY[VIEW][ntrack[VIEW]][0])
                                /(XX[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1] - XX[VIEW][ntrack[VIEW]][0])
                                );
                            graph2->Fit("f2","Q");
                            tmp = f2->GetChisquare();
#ifdef DEBUG_TWODIMRECON
                            cout << "---- Define tracks (!dwst)" <<endl;
                            double fitpara0 = f2->GetParameter(0);
                            double fitpara1 = f2->GetParameter(1);
                            double rms = 0.;
                            cout << "{";
                            for(int i=0;i<ntracell[VIEW][ntrack[VIEW]];i++){
                              rms += pow( YY[VIEW][ntrack[VIEW]][i]-fitpara0
                                         -XX[VIEW][ntrack[VIEW]][i]*fitpara1,2);
                              cout << " " << plane[VIEW][ntrack[VIEW]][i];
                            }
                            cout << "}" << endl;
                            rms = sqrt(rms);
                            cout
                              << " view="   << VIEW << " axis=" << axis
                              << " track="  << ntrack[VIEW]
                              << " intcpt=" << fitpara0
                              << " slope="  << fitpara1
                              << " chi2="   << chi2
                              << " rms="    << rms
                              << endl;
#endif
                            graph2->Delete();
                            f2->Delete();

                            float upstpe1 
                              = clpe[VIEW]
                                    [PLN2-DIST2]
                                    [cellu[VIEW]
                                          [PLN2-DIST2]
                                          [track_cell[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]]]
                                          [DIST2]];
                            float upstpe2 
                              = clpe[VIEW]
                                    [PLN2-DIST2]
                                    [cellu[VIEW][PLN2-DIST2][neiu[VIEW][PLN2][NEI][DIST2][DIST3]][DIST2]];

                            if( ((mod <NUMINGMOD||mod==MOD_B2_INGRID) && (tmp>chi2||upstpe2<upstpe_th) && upstpe1>=upstpe_th) ||
                                ((mod==MOD_PM   ||mod==MOD_B2_CH    ) && (tmp>chi2||upstpe2<upstpe_th) && upstpe1>=upstpe_th) ||
                                ((mod==MOD_ONAXIS_WM||mod==MOD_B2_WM) && (tmp>chi2)                    && upstpe1>=upstpe_th1))
                            {
                              Xe [VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=Xetmp;
                              YY [VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=YYtmp;
                              Yel[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=Yeltmp;
                              Yeh[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]=Yehtmp;
                              continue;
                            }

                          }//if(dwst)

                          track_cell[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]] = neiu[VIEW][PLN2][NEI][DIST2][DIST3];
                          track_pln [VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]] = PLN2-DIST2;
                          track_dist[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]] = DIST2;

                          plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1] 
                            = track_pln[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]];

                          clus[VIEW ][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1] 
                            = cellu[VIEW]
                                   [track_pln [VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]]]
                                   [track_cell[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]]]
                                   [track_dist[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]]];

                          XX[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]    
                            = zposi(mod,VIEW,plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1],0,axis);			  
                          YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]    
                            = clcenter[VIEW]
                                      [plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]]
                                      [clus [VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]];

                          Yel[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1] 
                            = -Yerr(mod,
                                VIEW,
                                plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1],
                                clchi[VIEW]
                                     [plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]]
                                     [clus [VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]],
                                axis,
                                0);
                          Yeh[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]
                            = +Yerr(mod,
                                VIEW,
                                plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1],
                                clchf[VIEW]
                                     [plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]]
                                     [clus [VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]],
                                axis,
                                1);
                          Xe[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]
                            = scithick(mod,VIEW,plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1],
                                YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1],axis,0);

                          TGraphAsymmErrors *graph1 = new TGraphAsymmErrors((ntracell[VIEW][ntrack[VIEW]]+2),
                              XX[VIEW][ntrack[VIEW]],
                              YY[VIEW][ntrack[VIEW]],
                              Xe[VIEW][ntrack[VIEW]],
                              Xe[VIEW][ntrack[VIEW]],
                              Yel[VIEW][ntrack[VIEW]],
                              Yeh[VIEW][ntrack[VIEW]]);
                          TF1 *f1 = new TF1("f1","[0]+[1]*x");
                          f1->SetParameters(YY[VIEW][ntrack[VIEW]][0]
                              -XX[VIEW][ntrack[VIEW]][0]
                              *(YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]-YY[VIEW][ntrack[VIEW]][0])
                              /(XX[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]-XX[VIEW][ntrack[VIEW]][0]),
                               (YY[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]-YY[VIEW][ntrack[VIEW]][0])
                              /(XX[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]-XX[VIEW][ntrack[VIEW]][0])
                              );

                          graph1->Fit("f1","Q");
                          chi2 = f1->GetChisquare();
#ifdef DEBUG_TWODIMRECON
                          double fitpara0 = f1->GetParameter(0);
                          double fitpara1 = f1->GetParameter(1);
                          double ndf      = f1->GetNDF();
                          double rms = 0.;
                          cout << "{";
                          for(int i=0;i<ntracell[VIEW][ntrack[VIEW]]+2;i++){
                            rms += pow( YY[VIEW][ntrack[VIEW]][i]-fitpara0
                                       -XX[VIEW][ntrack[VIEW]][i]*fitpara1,2);
                            cout << " " << plane[VIEW][ntrack[VIEW]][i];
                          }
                          cout << "}" << endl;
                          rms = sqrt(rms);
                          cout
                            << " view="     << VIEW << " axis=" << axis
                            << " track="    << ntrack[VIEW]
                            << " ntracell=" << ntracell[VIEW][ntrack[VIEW]]
                            << " intcpt="   << fitpara0
                            << " slope="    << fitpara1
                            << " chi2="     << chi2
                            << " chi2/ntracel="<<  chi2/ntracell[VIEW][ntrack[VIEW]]
                            << " chi2old="  << chi2old
                            << " rms="      << rms
                            << " ndf="      << ndf
                            << endl;
#endif
                          graph1->Delete();
                          f1->Delete();


                          if(((mod<NUMINGMOD || mod==MOD_B2_INGRID) &&
                                ( (chi2/ntracell[VIEW][ntrack[VIEW]]>3 && ntracell[VIEW][ntrack[VIEW]]==1)|| 
                                  (chi2/ntracell[VIEW][ntrack[VIEW]]>2 && ntracell[VIEW][ntrack[VIEW]]>1 )|| 
                                 ((chi2/ntracell[VIEW][ntrack[VIEW]]-chi2old)>1.5 -(ntracell[VIEW][ntrack[VIEW]]-1)*0.02)))
                              ||
                              ((mod==MOD_PM || mod==MOD_B2_CH) &&
                               ( (chi2/ntracell[VIEW][ntrack[VIEW]]>3 && ntracell[VIEW][ntrack[VIEW]]==1)||
                                 (chi2/ntracell[VIEW][ntrack[VIEW]]>2 && ntracell[VIEW][ntrack[VIEW]]>1)|| 
                                ((chi2/ntracell[VIEW][ntrack[VIEW]]-chi2old)>1.5 -(ntracell[VIEW][ntrack[VIEW]]-1)*0.015)))
                              ||
                              ((mod==MOD_ONAXIS_WM || mod==MOD_B2_WM) &&
                               ((chi2/ntracell[VIEW][ntrack[VIEW]]>neibor_th&&ntracell[VIEW][ntrack[VIEW]]==1)|| 
                                ((chi2/ntracell[VIEW][ntrack[VIEW]]-chi2old)> 1.6+(ntracell[VIEW][ntrack[VIEW]]-1)*0.3))) 
                            )
                          {
                            if(value[VIEW][PLN2-DIST2][neiu[VIEW][PLN2][NEI][DIST2][DIST3]][DIST2]>2 && TMP==0)
                            {
                              neibor[VIEW][PLN2][NEI][DIST2][DIST3] = false;
                            }
#ifdef DEBUG_TWODIMRECON
                            cout << "    ::::: neibor is changed to false." 
                              << " pln:"     << PLN2
                              << " dist_u:"  << DIST2
                              << " dist_d:"  << DIST3
                              << " value:"   << value[VIEW][PLN2-DIST2][neiu[VIEW][PLN2][NEI][DIST2][DIST3]][DIST2]
                              << endl;
#endif
                            continue;
                          }
                          dwst = false;
                        }//if(value0 = value1 + 1)
                      }//if(neid = trackcell)
                    }//for(NEI)
                    if(!dwst) break;
                  }//for(TMP)
                  if(!dwst) break;
                }
                if(dwst) break;

                PLN2    = PLN2 - track_dist[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]]-1;
                chi2old = chi2 / ntracell[VIEW][ntrack[VIEW]];

                trape[VIEW][ntrack[VIEW]] 
                  += clpe[VIEW]
                         [plane[VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]]
                         [clus [VIEW][ntrack[VIEW]][ntracell[VIEW][ntrack[VIEW]]+1]];

                ntracell[VIEW][ntrack[VIEW]]++;
              }//while(PLN2>=0)
              if(ntracell[VIEW][ntrack[VIEW]]<N_long) continue;

              ntracl[VIEW][ntrack[VIEW]] = ntracell[VIEW][ntrack[VIEW]]+1;
              ntrack[VIEW]++;
#ifdef DEBUG_TWODIMRECON
              cout << "      >>> Track : " << ntrack[VIEW]-1 << endl;
#endif
            } //if(upst)
          }//if(value>0)
        }//for(CELL)
      }//for(DIST)
    }//for(PLN)

    //=============
    //Track rank
    memset(trank [VIEW],0,sizeof(trank [VIEW]));
    memset(rank  [VIEW],0,sizeof(rank  [VIEW]));
    memset(ncltra[VIEW],0,sizeof(ncltra[VIEW]));
    ntrack2[VIEW] = 0;

    for(TRACL=MAX_CELLVALUE; TRACL>0; TRACL--){
      for(TRA=0; TRA<ntrack[VIEW]; TRA++){
        if(ntracl[VIEW][TRA]==TRACL){
          trank[VIEW][TRACL][ncltra[VIEW][TRACL]] = TRA;
          ncltra[VIEW][TRACL]++;
        }
      }

      for(CL=0; CL<ncltra[VIEW][TRACL]; CL++){
        for(CL2=CL+1; CL2<ncltra[VIEW][TRACL]; CL2++){
          if(trape[VIEW][trank[VIEW][TRACL][CL]]<trape[VIEW][trank[VIEW][TRACL][CL2]]){
            dummy                   = trank[VIEW][TRACL][CL];
            trank[VIEW][TRACL][CL]  = trank[VIEW][TRACL][CL2];
            trank[VIEW][TRACL][CL2] = dummy;
          }
        }
        rank[VIEW][ntrack2[VIEW]] = trank[VIEW][TRACL][CL];
        ntrack2[VIEW]++;
      }
    }



    //======================
    //True track selection
    memset(ttrack [VIEW],true ,sizeof(ttrack [VIEW]));
    memset(hitcl  [VIEW],false,sizeof(hitcl  [VIEW]));
    memset(ntracl2[VIEW],   0 ,sizeof(ntracl2[VIEW]));
    memset(ovcl   [VIEW],false,sizeof(ovcl   [VIEW]));
    ntrack2[VIEW] = 0;
#ifdef DEBUG_TWODIMRECON
      cout << "---------- True track selection ----" << endl;
#endif
    for(TRA=0; TRA<ntrack[VIEW]; TRA++){
      TRA2 = rank[VIEW][TRA];

      for(TRACL=0; TRACL<ntracl[VIEW][TRA2]; TRACL++){
        if(  ! hitcl [VIEW][plane[VIEW][TRA2][TRACL]][clus[VIEW][TRA2][TRACL]]
            && clused[VIEW][plane[VIEW][TRA2][TRACL]][clus[VIEW][TRA2][TRACL]]==0)
        {
          ntracl2[VIEW][TRA2]++;
        }
        else if((mod==MOD_ONAXIS_WM || mod==MOD_B2_WM) && TRACL!=ntracl[VIEW][TRA2]-1)
        {
          ntracl2[VIEW][TRA2] = ntracl2[VIEW][TRA2]-1;
        }
      }

      if(ntracl2[VIEW][TRA2]>N_long)
      {
        for(TRACL=0; TRACL<ntracl[VIEW][TRA2]; TRACL++){
          if(hitcl[VIEW][plane[VIEW][TRA2][TRACL]][clus[VIEW][TRA2][TRACL]]){
            ovcl[VIEW][plane[VIEW][TRA2][TRACL]][clus[VIEW][TRA2][TRACL]] = true;
          }
          hitcl[VIEW][plane[VIEW][TRA2][TRACL]][clus[VIEW][TRA2][TRACL]] = true;
        }
        ntrack2[VIEW]++;
      }
      else{
        ttrack[VIEW][TRA2] = false;
      }
#ifdef DEBUG_TWODIMRECON
      cout << " Track " << TRA2  
        << " : " << ttrack[VIEW][TRA2] 
        << " : pe=" << trape[VIEW][TRA2]
        << endl;
#endif
    }
    
    //===================
    //Fit tracks
    for(TRA=0; TRA<ntrack[VIEW]; TRA++){
      if(!ttrack[VIEW][TRA]) continue;
      TGraphAsymmErrors *g = new TGraphAsymmErrors(ntracell[VIEW][TRA]+1,
          XX [VIEW][TRA],
          YY [VIEW][TRA],
          Xe [VIEW][TRA],
          Xe [VIEW][TRA],
          Yel[VIEW][TRA],
          Yeh[VIEW][TRA]);
      TF1 *f = new TF1("f","[0]+[1]*x");
      f->SetParameters(YY[VIEW][TRA][0] 
          - XX[VIEW][TRA][0]
          *(YY[VIEW][TRA][ntracell[VIEW][TRA]]-YY[VIEW][TRA][0])
          /(XX[VIEW][TRA][ntracell[VIEW][TRA]]-XX[VIEW][TRA][0]),
           (YY[VIEW][TRA][ntracell[VIEW][TRA]]-YY[VIEW][TRA][0])
          /(XX[VIEW][TRA][ntracell[VIEW][TRA]]-XX[VIEW][TRA][0]));
      g->Fit("f","Q");
      TF1 *func = g->GetFunction("f");
      par[VIEW][TRA][0] = func->GetParameter(0);
      par[VIEW][TRA][1] = func->GetParameter(1);
#ifdef DEBUG_TWODIMRECON
      cout << "--- Fit Track : " << TRA << endl;
      cout 
        << " view="   << VIEW
        << " intcpt=" << par[VIEW][TRA][0]
        << " slope="  << par[VIEW][TRA][1]
        << " pln={";
      for(int i=0;i<ntracell[VIEW][TRA]+1;i++){
        cout << " " << plane[VIEW][TRA][i];
      }
      cout <<"}"<<endl;
      cout << endl;
      cout << "XX={";
      for(int i=0;i<ntracell[VIEW][TRA]+1;i++){cout<<" "<<XX[VIEW][TRA][i];}
      cout << "}"<<endl;
      cout << "YY={";
      for(int i=0;i<ntracell[VIEW][TRA]+1;i++){cout<<" "<<YY[VIEW][TRA][i];}
      cout << "}"<<endl;
      cout << "Xe={";
      for(int i=0;i<ntracell[VIEW][TRA]+1;i++){cout<<" "<<Xe[VIEW][TRA][i];}
      cout << "}"<<endl;
      cout << "Yel={";
      for(int i=0;i<ntracell[VIEW][TRA]+1;i++){cout<<" "<<Yel[VIEW][TRA][i];}
      cout << "}"<<endl;
      cout << "Yeh={";
      for(int i=0;i<ntracell[VIEW][TRA]+1;i++){cout<<" "<<Yeh[VIEW][TRA][i];}
      cout << "}"<<endl;

#endif

      g->Delete();
      f->Delete();
    }
    if(axis==0&&fabs(atan(par[VIEW][TRA][1]))>axis_border    ){continue;}
    if(axis==1&&fabs(atan(par[VIEW][TRA][1]))>90.-axis_border){continue;}

    //=============================================
    // Remove duplicate track
#ifdef DEBUG_TWODIMRECON
    cout << "------------ Remove duplicate track ----- " <<endl;
#endif
    for(TRA=0; TRA<ntrack[VIEW]-1; TRA++){
      if(ttrack[VIEW][TRA]){
        double startX = XX[VIEW][TRA][0];
        double stopX  = XX[VIEW][TRA][ntracell[VIEW][TRA]];
        for(TRA2=TRA+1; TRA2<ntrack[VIEW]; TRA2++){
          if(ttrack[VIEW][TRA2]){
            double diff_angle  = fabs(atan(par[VIEW][TRA][1])-atan(par[VIEW][TRA2][1]))*180./PI;
            if(diff_angle<10.){
              int TRA3,TRA4;
              if(rank[VIEW][TRA]<rank[VIEW][TRA2]){TRA3=TRA2;TRA4=TRA ;}
              else                                {TRA3=TRA ;TRA4=TRA2;}
              double hit_dist = 0.;
              int num_closehit = 0;
              int num_awayhit  = 0;
              for(int i=0;i<ntracell[VIEW][TRA3]+2;i++){
                double hitx_tmp = XX[VIEW][TRA3][i];
                double hity_tmp = YY[VIEW][TRA3][i]; 
                hit_dist = fabs(hity_tmp-par[VIEW][TRA4][0]
                    -hitx_tmp*par[VIEW][TRA4][1])/sqrt(1.+par[VIEW][TRA4][1]*par[VIEW][TRA4][1]);
                if(hitx_tmp>startX&&hitx_tmp<stopX&&hit_dist<joint_th){
                  num_closehit++; 
                }else{ 
                  num_awayhit++ ; 
                }
#ifdef DEBUG_TWODIMRECON
              cout << "       hit_dist=" << hit_dist << endl;
#endif
              }
              if(num_closehit>num_awayhit){
                ttrack[VIEW][TRA3]=false;
#ifdef DEBUG_TWODIMRECON
              cout << "  >>> TRA:" << TRA3 << " is removed."  << endl;
#endif

              }
            }
          }
        }
      }
    }

    //==================================
    //Reset hitcl for not selected track
#ifdef DEBUG_TWODIMRECON
    cout << "------------ Reset hitcl for not selected track ------ " << endl;
#endif

    for(TRA=0; TRA<ntrack[VIEW]; TRA++){
      if(!ttrack[VIEW][TRA]){
        for(TRACL=0; TRACL<ntracl[VIEW][TRA]; TRACL++){
          hitcl[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]] = false;
        }
      }
    }
    for(TRA=0; TRA<ntrack[VIEW]; TRA++){
      if(ttrack[VIEW][TRA]){
        for(TRACL=0; TRACL<ntracl[VIEW][TRA]; TRACL++){
          hitcl[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]] = true;
        }
      }
    }

    //==========================================
    //Upstream veto & edge channel cut
#ifdef DEBUG_TWODIMRECON
    cout << "------------ Upstream veto & edge channel cut ------ " << endl;
#endif

    memset(vetowtracking[VIEW],false,sizeof(vetowtracking[VIEW]));
    memset(edgewtracking[VIEW],false,sizeof(edgewtracking[VIEW]));
    memset(track_stop   [VIEW],true ,sizeof(track_stop   [VIEW]));
    for(TRA=0; TRA<ntrack[VIEW]; TRA++){
      if(!ttrack[VIEW][TRA]) continue;

      if(   ((mod<NUMINGMOD || mod==MOD_B2_INGRID ) && plane[VIEW][TRA][ntracl[VIEW][TRA]-1]==0) 
          ||((mod==MOD_ONAXIS_WM || mod==MOD_B2_WM) && plane[VIEW][TRA][ntracl[VIEW][TRA]-1]<=2) 
          ||((mod==MOD_PM || mod==MOD_B2_CH       ) && plane[VIEW][TRA][ntracl[VIEW][TRA]-1]<=1)
          )
      {
        vetowtracking[VIEW][TRA] = true;
      }
      if(plane[VIEW][TRA][0]==plnmax(mod,VIEW,plane[VIEW][TRA][0],axis)-1){
        track_stop[VIEW][TRA] = false;
      }

#ifdef PMEDGE
      for(TRACL=ntracl[VIEW][TRA]-2; TRACL<ntracl[VIEW][TRA]; TRACL++){
        for(CLHIT=0; CLHIT<numcl[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]]; CLHIT++){
          if(ch[VIEW][clhit[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]][CLHIT]]==0){
            edgewtracking[VIEW][TRA] = true;
          }
          else if(ch[VIEW][clhit[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]][CLHIT]]
              == chmax(mod,VIEW,plane[VIEW][TRA][TRACL],axis)-1)
          {
            edgewtracking[VIEW][TRA] = true;
          }
        }
      }

      for(TRACL=0; TRACL<2; TRACL++){
        for(CLHIT=0; CLHIT<numcl[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]]; CLHIT++){
          if(ch[VIEW][clhit[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]][CLHIT]]==0){
            track_stop[VIEW][TRA] = false;
          }
          else if(ch[VIEW][clhit[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]][CLHIT]]
              ==chmax(mod,VIEW,plane[VIEW][TRA][TRACL],axis)-1)
          {
            track_stop[VIEW][TRA] = false;
          }
        }
      }
#else
      if(par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]],0,axis))+par[VIEW][TRA][0]<200){
        edgewtracking[VIEW][TRA] = true;
      }
      if(par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]],0,axis))+par[VIEW][TRA][0]>1000){
        edgewtracking[VIEW][TRA] = true;
      }

      if(par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][0],0,axis))+par[VIEW][TRA][0]<100){
        track_stop[VIEW][TRA] = false;
      }
      if(par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][0],0,axis))+par[VIEW][TRA][0]>1100){
        track_stop[VIEW][TRA] = false;
      }
#endif
    }

    //=======================
    //Side veto
#ifdef DEBUG_TWODIMRECON
    cout << "------------ Side Veto ------ " << endl;
#endif

    for(TRA=0; TRA<ntrack[VIEW]; TRA++){
      vetodist[VIEW][TRA] = 1e5;
      if(!ttrack[VIEW][TRA]) continue;

      for(HIT=0; HIT<hit[VIEW]; HIT++){
        if(pln[VIEW][HIT]>=plnmax(mod,VIEW,pln[VIEW][HIT],axis)){

          if(par[VIEW][TRA][1]>0){
            if(vposixy(mod,pln[VIEW][HIT])>0) continue;
            if((mod==MOD_PM||mod==MOD_B2_CH) 
                && XX[VIEW][TRA][ntracl[VIEW][TRA]-1]-(0-par[VIEW][TRA][0])/par[VIEW][TRA][1]>46*2
                )
            { continue;}
            if((mod <NUMINGMOD||mod==MOD_B2_INGRID) 
                && XX[VIEW][TRA][ntracl[VIEW][TRA]-1]-(0-par[VIEW][TRA][0])/par[VIEW][TRA][1]>107*2
                )
            { continue;}
          }
          else if(par[VIEW][TRA][1]<0){
            if(vposixy(mod,pln[VIEW][HIT])<0) continue;
            if((mod==MOD_PM||mod==MOD_B2_CH) 
                && XX[VIEW][TRA][ntracl[VIEW][TRA]-1]-(1200-par[VIEW][TRA][0])/par[VIEW][TRA][1]>46*2
                )
            { continue;}
            if((mod <NUMINGMOD||mod==MOD_B2_INGRID) 
                && XX[VIEW][TRA][ntracl[VIEW][TRA]-1]-(1200-par[VIEW][TRA][0])/par[VIEW][TRA][1]>107*2
                )
            { continue; }
          }
          else continue;

          if(mod==MOD_ONAXIS_WM || mod==MOD_B2_WM){ continue; }

          dis = fabs(par[VIEW][TRA][1]*vposiz(mod,ch[VIEW][HIT])-vposixy(mod,pln[VIEW][HIT])+par[VIEW][TRA][0])
            /sqrt((par[VIEW][TRA][1])*(par[VIEW][TRA][1])+1);

          if(dis<80)                  vetowtracking[VIEW][TRA] = true;
          if(dis<vetodist[VIEW][TRA]) vetodist[VIEW][TRA] = dis;
        }
      }
    }

    //================
    //FC with veto
#ifdef DEBUG_TWODIMRECON
    cout << "------------ FC with Veto ------ " << endl;
#endif

    for(TRA=0; TRA<ntrack[VIEW]; TRA++){
      if(!ttrack[VIEW][TRA]) continue;

      for(HIT=0; HIT<hit[VIEW]; HIT++){
        if(pln[VIEW][HIT]>=plnmax(mod,VIEW,pln[VIEW][HIT],axis)){

          if(par[VIEW][TRA][1]>0){
            if(vposixy(mod,pln[VIEW][HIT])<0) continue;
            if((mod==MOD_PM||mod==MOD_B2_CH) && XX[VIEW][TRA][0]-(1200-par[VIEW][TRA][0])/par[VIEW][TRA][1]>46*2)
            { continue;}
            if((mod <NUMINGMOD||mod==MOD_B2_INGRID) && XX[VIEW][TRA][0]-(1200-par[VIEW][TRA][0])/par[VIEW][TRA][1]>107*2)
            { continue;}
          }
          else if(par[VIEW][TRA][1]<0){
            if(vposixy(mod,pln[VIEW][HIT])>0) continue;
            if((mod==MOD_PM||mod==MOD_B2_CH) && XX[VIEW][TRA][0]-(0-par[VIEW][TRA][0])/par[VIEW][TRA][1]>46*2)
            { continue; }
            if((mod <NUMINGMOD||mod==MOD_B2_INGRID) && XX[VIEW][TRA][0]-(0-par[VIEW][TRA][0])/par[VIEW][TRA][1]>107*2)
            { continue; }
          }
          else continue;

          if(mod==MOD_ONAXIS_WM || mod==MOD_B2_WM) continue;
          dis = fabs(par[VIEW][TRA][1]*vposiz(mod,ch[VIEW][HIT])-vposixy(mod,pln[VIEW][HIT])+par[VIEW][TRA][0])
            /sqrt((par[VIEW][TRA][1])*(par[VIEW][TRA][1])+1);

          if(dis<80) track_stop[VIEW][TRA] = false;
        }
      }
    }

    //==========
    //Track hits
#ifdef DEBUG_TWODIMRECON
    cout << "------------ Track Hits ------ " << endl;
#endif

    for(TRA=0; TRA<ntrack[VIEW]; TRA++){
      if(ttrack[VIEW][TRA]){
        TMP = 0;
        for(TRACL=0; TRACL<ntracl[VIEW][TRA]; TRACL++){
          for(CLHIT=0; CLHIT<numcl[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]]; CLHIT++){
            hit_id[VIEW][TRA][TMP] = id[VIEW][clhit[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]][CLHIT]];

            if(ovcl[VIEW][plane[VIEW][TRA][TRACL]][clus[VIEW][TRA][TRACL]]){
              isohit[VIEW][TRA][TMP] = false;
            }
            else{	    
              isohit[VIEW][TRA][TMP] = true;
            }
            TMP++;
          }
        }
        hit_n[VIEW][TRA] = TMP;
      }
    }


    for(int iFit=0;iFit<1;iFit++){
      if(mod!=MOD_ONAXIS_WM && mod!=MOD_B2_WM && iFit!=0){continue;}
      //=======================
      //Merge hit around track
#ifdef DEBUG_TWODIMRECON
      cout << "------------ Merge hit around track ----- " <<endl;
#endif
      int num_mergedhit = 0;
      if(mod==MOD_ONAXIS_WM || mod==MOD_B2_WM){
        for(TRA=0; TRA<ntrack[VIEW]; TRA++){
          if(ttrack[VIEW][TRA]){
            TMP = hit_n[VIEW][TRA];
            for(PLN=plane[VIEW][TRA][ntracell[VIEW][TRA]]-4; PLN<=plane[VIEW][TRA][0]+4; PLN++){
              if(PLN<0||PLN>=plnmax(mod,VIEW,PLN,axis)){continue;}
              for(CL=0; CL<ncl[VIEW][PLN]; CL++){
                if(!hitcl[VIEW][PLN][CL]){
                  double dist_hit = fabs(
                      par[VIEW][TRA][0] + par[VIEW][TRA][1]*zposi(mod,VIEW,PLN,0,axis)
                      -clcenter[VIEW][PLN][CL])
                    / sqrt(par[VIEW][TRA][1]*par[VIEW][TRA][1]+1.);
#ifdef DEBUG_TWODIMRECON
                  cout
                    << " pln=" << PLN
                    << " dist=" << dist_hit
                    << endl;
#endif
                  if(dist_hit<joint_th)
                  {
                    for(CLHIT=0; CLHIT<numcl[VIEW][PLN][CL]; CLHIT++){
                      hit_id[VIEW][TRA][TMP] = id[VIEW][clhit[VIEW][PLN][CL][CLHIT]];
                      isohit[VIEW][TRA][TMP] = true;
                      hitcl [VIEW][PLN][CL]  = true;
                      TMP++;
                      num_mergedhit++;
                    }
                  }
                }
              }
              hit_n[VIEW][TRA] = TMP;
            }
          }
        }
      }


      //=============================================
      //Arrange reconstructed hits order in  a track
      for(TRA=0; TRA<ntrack[VIEW]; TRA++){
        if(ttrack[VIEW][TRA]){
          int id_temp;
          //for(int TMP=0; TMP<hit_n[VIEW][TRA]; TMP++){
          //  id_temp[TMP] = hit_id[VIEW][TRA][hit_n[VIEW][TRA]-1-TMP];
          //}
          //for(int TMP=0; TMP<hit_n[VIEW][TRA]; TMP++){
          //  hit_id[VIEW][TRA][TMP] = id_temp[TMP];
          //}
          for(int TMP=0; TMP<hit_n[VIEW][TRA]-1; TMP++){
            for(int TMP1=TMP+1; TMP1<hit_n[VIEW][TRA]; TMP1++){
              int pln1 = -1;
              for(int i=0; i<hitclssize; i++){
                if(hitcls[i].id==hit_id[VIEW][TRA][TMP]){
                  pln1 = hitcls[i].pln;
                  break;
                }
              }
              int pln2 = -1;
              for(int i=0; i<hitclssize; i++){
                if(hitcls[i].id==hit_id[VIEW][TRA][TMP1]){
                  pln2 = hitcls[i].pln;
                  break;
                }
              }
              if(pln1>pln2){
                id_temp = hit_id[VIEW][TRA][TMP];
                hit_id[VIEW][TRA][TMP] = hit_id[VIEW][TRA][TMP1];
                hit_id[VIEW][TRA][TMP1] = id_temp;
              }
            }
          }
        }
      }

      // ===================================================
      // Reset start/stop pln & Fit the all hits in a track
#ifdef DEBUG_TWODIMRECON
    cout << "--- Reset start/stop pln and Fit all hits  " << endl;
#endif
      if((mod==MOD_ONAXIS_WM || mod==MOD_B2_WM) && num_mergedhit>0){
        for(TRA=0; TRA<ntrack[VIEW]; TRA++){
          if(ttrack[VIEW][TRA]){
            int plnhitid=0;
            int plnlast =-1;
            int npln = 0;
            double yytmp1=0.,yytmp2=0.;
            for(int TMP=0; TMP<hit_n[VIEW][TRA]; TMP++){
              int pln1=-1,ch1=-1;
              for(int i=0; i<hitclssize; i++){
                if(hitcls[i].id==hit_id[VIEW][TRA][TMP]){
                  pln1 = hitcls[i].pln;
                  ch1  = hitcls[i].ch ;
                  break;
                }
              }
              if(pln1==plnlast){plnhitid++;}
              else             {plnhitid=0;npln++;}
              if(pln1<plane[VIEW][TRA][ntracell[VIEW][TRA]]){plane[VIEW][TRA][ntracell[VIEW][TRA]]=pln1;}
              if(pln1>plane[VIEW][TRA][0]                  ){plane[VIEW][TRA][0]                  =pln1;}
              if(plnhitid==0){
                XX [VIEW][TRA][npln-1] = zposi   (mod,VIEW,pln1,  0,axis);
                YY [VIEW][TRA][npln-1] = xyposi  (mod,VIEW,pln1,ch1,axis);
                Xe [VIEW][TRA][npln-1] = scithick(mod,VIEW,pln1,YY[VIEW][TRA][pln1],axis,0);
                Yel[VIEW][TRA][npln-1] = sciwidth(mod,VIEW,pln1,ch1,axis);
                Yeh[VIEW][TRA][npln-1] = sciwidth(mod,VIEW,pln1,ch1,axis);
                yytmp1 = YY[VIEW][TRA][npln-1];
                yytmp2 = YY[VIEW][TRA][npln-1];
              }
              else{
                if(yytmp1>xyposi(mod,VIEW,pln1,ch1,axis)){yytmp1=xyposi(mod,VIEW,pln1,ch1,axis);}
                if(yytmp2<xyposi(mod,VIEW,pln1,ch1,axis)){yytmp2=xyposi(mod,VIEW,pln1,ch1,axis);}
                YY [VIEW][TRA][npln-1] *= plnhitid;
                YY [VIEW][TRA][npln-1] += xyposi  (mod,VIEW,pln1,ch1,axis);
                YY [VIEW][TRA][npln-1] /= (plnhitid+1);
                Yel[VIEW][TRA][npln-1] = sciwidth(mod,VIEW,pln1,ch1,axis)+(yytmp2-yytmp1)/2;
                Yeh[VIEW][TRA][npln-1] = sciwidth(mod,VIEW,pln1,ch1,axis)+(yytmp2-yytmp1)/2;
              }
              plnlast = pln1;
#ifdef DEBUG_TWODIMRECON
              cout << " " << pln1;
#endif
            }
            
#ifdef DEBUG_TWODIMRECON
            cout << endl;
            cout << "XX={";
            for(int i=0;i<npln;i++){cout<<" "<<XX[VIEW][TRA][i];}
            cout << "}"<<endl;
            cout << "YY={";
            for(int i=0;i<npln;i++){cout<<" "<<YY[VIEW][TRA][i];}
            cout << "}"<<endl;
            cout << "Xe={";
            for(int i=0;i<npln;i++){cout<<" "<<Xe[VIEW][TRA][i];}
            cout << "}"<<endl;
            cout << "Yel={";
            for(int i=0;i<npln;i++){cout<<" "<<Yel[VIEW][TRA][i];}
            cout << "}"<<endl;
            cout << "Yeh={";
            for(int i=0;i<npln;i++){cout<<" "<<Yeh[VIEW][TRA][i];}
            cout << "}"<<endl;
#endif
            TGraphAsymmErrors *g = new TGraphAsymmErrors(npln,
                XX [VIEW][TRA],
                YY [VIEW][TRA],
                Xe [VIEW][TRA],
                Xe [VIEW][TRA],
                Yel[VIEW][TRA],
                Yeh[VIEW][TRA]);
            TF1 *f = new TF1("f","[0]+[1]*x");
            f->SetParameters(YY[VIEW][TRA][0] 
                - XX[VIEW][TRA][0]
                *(YY[VIEW][TRA][npln]-YY[VIEW][TRA][0])
                /(XX[VIEW][TRA][npln]-XX[VIEW][TRA][0]),
                 (YY[VIEW][TRA][npln]-YY[VIEW][TRA][0])
                /(XX[VIEW][TRA][npln]-XX[VIEW][TRA][0]));
            g->Fit("f","Q");
            TF1 *func = g->GetFunction("f");
            par[VIEW][TRA][0] = func->GetParameter(0);
            par[VIEW][TRA][1] = func->GetParameter(1);
            delete g;
            delete f;
          }
        }
      }
    }

#ifdef DEBUG_TWODIMRECON
    for(TRA=0; TRA<ntrack[VIEW]; TRA++){
      if(ttrack[VIEW][TRA]){
        cout 
          << " TRA :"  << TRA
          << " view="   << VIEW
          << " intcpt=" << par[VIEW][TRA][0]
          << " slope="  << par[VIEW][TRA][1]
          << endl;
      }
    }
#endif


    //===================================
    //Reconstructed track information
#ifdef DEBUG_TWODIMRECON
    cout<< "------- Reconstructed track information---- " << endl;
#endif

    Track track;
    for(TRA=0; TRA<ntrack[VIEW]; TRA++){
      if(ttrack[VIEW][TRA]){
#ifdef DEBUG_TWODIMRECON
        cout << " Track : " << TRA << endl;
#endif
        track.view     = VIEW;
        track.fpln     = plane[VIEW][TRA][0];
        track.ipln     = plane[VIEW][TRA][ntracell[VIEW][TRA]];
        track.fxy      = par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][0],0,axis))+par[VIEW][TRA][0];
        track.ixy      = par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]],0,axis))+par[VIEW][TRA][0];
        track.fz       = zposi(mod,VIEW,plane[VIEW][TRA][0],0,axis);
        track.iz       = zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]],0,axis);
        track.intcpt   = par[VIEW][TRA][0];
        track.slope    = par[VIEW][TRA][1];
        track.ang      = atan(par[VIEW][TRA][1])*180./PI;
        track.veto     = vetowtracking[VIEW][TRA];
        track.edge     = edgewtracking[VIEW][TRA];
        track.stop     = track_stop[VIEW][TRA];
        track.vetodist = vetodist[VIEW][TRA];
        track.hitid.clear();
        track.isohit.clear();

        for(TMP=0; TMP<hit_n[VIEW][TRA]; TMP++){
          track.hitid .push_back(hit_id[VIEW][TRA][TMP]);
          track.isohit.push_back(isohit[VIEW][TRA][TMP]);
        }

        if(axis==1){
          if(track.slope>=0){
            track.fz  = par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][0],0,axis))+par[VIEW][TRA][0];
            track.iz  = par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]],0,axis))+par[VIEW][TRA][0];
            track.fxy = zposi(mod,VIEW,plane[VIEW][TRA][0],0,axis);
            track.ixy = zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]],0,axis);

            int pln, grid, gridch;
            int fch  = ch[VIEW][clhit[VIEW][plane[VIEW][TRA][0]][clus[VIEW][TRA][0]][0]];
            int fpln = plane[VIEW][TRA][0];
            detdim->GetRawPlnChGrid(mod, VIEW, fpln, fch, axis, &pln, &gridch, &grid);
            if(VIEW==0){
              track.fpln = pln*3 + grid;
            }
            else if(VIEW==1){
              if(grid==0) track.fpln = pln*3+1;
              if(grid==1) track.fpln = pln*3+0;
              if(grid==2) track.fpln = pln*3+2;
            }

            int ich  = ch[VIEW][clhit[VIEW][plane[VIEW][TRA][ntracell[VIEW][TRA]]][clus[VIEW][TRA][ntracell[VIEW][TRA]]][0]];
            int ipln = plane[VIEW][TRA][ntracell[VIEW][TRA]];
            detdim->GetRawPlnChGrid(mod, VIEW, ipln, ich, axis, &pln, &gridch, &grid);
            if(VIEW==0){
              track.ipln = pln*3 + grid;
            }
            else if(VIEW==1){
              if(grid==0) track.ipln = pln*3+1;
              if(grid==1) track.ipln = pln*3+0;
              if(grid==2) track.ipln = pln*3+2;
            }
          }
          else{
            track.iz  = par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][0],0,axis))+par[VIEW][TRA][0];
            track.fz  = par[VIEW][TRA][1]*(zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]],0,axis))
              + par[VIEW][TRA][0];
            track.ixy = zposi(mod,VIEW,plane[VIEW][TRA][0],0,axis); 
            track.fxy = zposi(mod,VIEW,plane[VIEW][TRA][ntracell[VIEW][TRA]],0,axis);

            int pln, grid, gridch;
            int ich = ch[VIEW][clhit[VIEW][plane[VIEW][TRA][0]][clus[VIEW][TRA][0]][0]];
            int ipln = plane[VIEW][TRA][0];
            detdim->GetRawPlnChGrid(mod, VIEW, ipln, ich, axis, &pln, &gridch, &grid);
            if(VIEW==0){
              track.ipln = pln*3 + grid;
            }
            else if(VIEW==1){
              if(grid==0) track.ipln = pln*3+1;
              if(grid==1) track.ipln = pln*3+0;
              if(grid==2) track.ipln = pln*3+2;
            }

            int fch  = ch[VIEW][clhit[VIEW][plane[VIEW][TRA][ntracell[VIEW][TRA]]][clus[VIEW][TRA][ntracell[VIEW][TRA]]][0]];
            int fpln = plane[VIEW][TRA][ntracell[VIEW][TRA]];
            detdim->GetRawPlnChGrid(mod, VIEW, fpln, fch, axis, &pln, &gridch, &grid);
            if(VIEW==0){
              track.fpln = pln*3 + grid;
            }
            else if(VIEW==1){
              if(grid==0) track.fpln = pln*3+1;
              if(grid==1) track.fpln = pln*3+0;
              if(grid==2) track.fpln = pln*3+2;
            }
          }

          if(fabs(track.slope)<1e-8) track.slope = 1e+8;
          else                       track.slope = 1./track.slope; //y=ax+b, a'=-1/a

          track.intcpt = -track.intcpt*track.slope; //y=ax+b, b'=-b/a
          track.ang    = atan(track.slope)*180./PI;

        }

#ifdef DEBUG_TWODIMRECON
        cout
          << " mod="   << mod
          << " view="  << VIEW
          << " axis="  << axis
          << " ipos={" << track.iz << "," << track.ixy << "}"
          << " fpos={" << track.fz << "," << track.fz  << "}"
          << " intcpt="<< track.intcpt
          << " slope=" << track.slope
          << " ang="   << track.ang
          <<endl;
#endif

        //Angle cut
        if(mod==MOD_ONAXIS_WM || mod==MOD_B2_WM){
          if(axis==0 && (track.ang>axis_border || track.ang<-axis_border)){ continue; }
        }
        alltrack.push_back(track);

#ifdef DEBUG_TWODIMRECON
        cout << " ----> Filled." << endl;
#endif

      }
    }
  }

  if(alltrack.size()>0) return true;
  else                  return false;
}; //fTracking


bool fTracking_WM(int mod)
{

#ifdef DEBUG_TWODIMRECON
  cout << "fTracking_WM::Start" << endl;
#endif

  alltrack.clear();
  int N_long = 5;
  vector<Hit> hit_temp;

  //Reconstructed long track along z axis
  hit_temp = hitcls;
  int hitclssize= (int)hitcls.size();
  for(int i=0; i<hitclssize; i++){
    int reconpln,reconch;
    detdim->GetReconPlnCh(mod, hitcls[i].view, hitcls[i].pln, hitcls[i].ch, 0, &reconpln, &reconch );
    hitcls[i].pln  = reconpln;
    hitcls[i].ch   = reconch;
    hitcls[i].used = 0;
  }
  fTracking(mod,0,N_long);

  //Reconstructed long track along xy axis
  for(int i=0; i<(int)alltrack.size(); i++){
    for(int j=1; j<(int)alltrack[i].hitid.size(); j++){
      vector<Hit>::iterator it;
      for(it=hit_temp.begin(); it!=hit_temp.end(); it++){
        if(it->id==alltrack[i].hitid[j]){
          it->used = 1;
        }
      }
    }
  }

  hitcls = hit_temp;
  for(int i=0; i<hitclssize; i++){
    int reconpln,reconch;
    detdim->GetReconPlnCh(mod, hitcls[i].view, hitcls[i].pln, hitcls[i].ch, 1, &reconpln, &reconch );
    hitcls[i].pln = reconpln;
    hitcls[i].ch  = reconch;
  }
  fTracking(mod,1,N_long);

  //Reconstructed track along z axis
  for(int i=0; i<(int)alltrack.size(); i++){
    for(int j=1; j<(int)alltrack[i].hitid.size(); j++){
      vector<Hit>::iterator it;

      for(it=hit_temp.begin(); it!=hit_temp.end(); it++){
        if(it->id==alltrack[i].hitid[j]){
          it->used = 1;
        }
      }
    }
  }

  hitcls = hit_temp;
  for(int i=0; i<hitclssize; i++){
    int reconpln,reconch;
    detdim->GetReconPlnCh(mod, hitcls[i].view, hitcls[i].pln, hitcls[i].ch, 0, &reconpln, &reconch );
    hitcls[i].pln = reconpln;
    hitcls[i].ch  = reconch;
  }
  fTracking(mod,0,1);

  //Reconstructed track along xy axis
  for(int i=0; i<(int)alltrack.size(); i++){
    for(int j=1; j<(int)alltrack[i].hitid.size(); j++){
      vector<Hit>::iterator it;
      for(it=hit_temp.begin(); it!=hit_temp.end(); it++){
        if(it->id==alltrack[i].hitid[j]){
          it->used = 1;
        }
      }
    }
  }

  hitcls = hit_temp;
  for(int i=0; i<hitclssize; i++){
    int reconpln, reconch;
    detdim->GetReconPlnCh(mod, hitcls[i].view, hitcls[i].pln, hitcls[i].ch, 1, &reconpln, &reconch );
    hitcls[i].pln = reconpln;
    hitcls[i].ch  = reconch;
  }
  fTracking(mod,1,1);

  hit_temp.clear();

  return true;
};

bool fFindNearestTracks(int* a, int* b, float angle_th, float dist_th){

  /* Function to find two nearest tracks to be connected for water module.
   * Conditions:
   * (1) Angle between two tracks should be less than "angle_th".
   * (2) One track's position should be lower than the other track.
   * (3) Two tracks should be near close each other.
   *     (extending two tracks, the distance between two tracks at the middle point should be less than "dist_th".)
   * If more than one candidate is found, calculate joilik=sqrt((angle/angle_th)^2+(dist/dist_th)^2) for each candidate and select minimum one.
   * If no candidate is found, return false.
   */

  const int Ntrack = (int)alltrack.size();
  float joilik = 1e-6;
  bool candidate = false;
  for(int i=0; i<Ntrack; i++){
    for(int j=0; j<i; j++){
      if(alltrack[i].view!=alltrack[j].view) continue;

      float angle;
      if(fabs(alltrack[i].ang-alltrack[j].ang)<90.0){
        angle = fabs(alltrack[i].ang-alltrack[j].ang);
      }
      else{
        angle = 180.0-fabs(alltrack[i].ang-alltrack[j].ang);
      }

#ifdef DEBUG_TWODIMRECON
      cout << "---- fFindNearestTracks " << endl;
      cout << "angle=" << angle << endl;
#endif

      if(angle>angle_th) continue; //Condidtion (1)

      vector<vector<vector<float> > > pos(2, vector<vector<float> >(2, vector<float>(2, 0.0)));;
      pos[0][0][0] = alltrack[i].ixy;
      pos[0][0][1] = alltrack[i].iz;
      pos[0][1][0] = alltrack[i].fxy;
      pos[0][1][1] = alltrack[i].fz;
      pos[1][0][0] = alltrack[j].ixy;
      pos[1][0][1] = alltrack[j].iz;
      pos[1][1][0] = alltrack[j].fxy;
      pos[1][1][1] = alltrack[j].fz;

      /*  Let theta to be mean of two tracks,
       *  if -45 < theta < 45 then fitaxis=0 (fit along z-axis),
       *  else fitaxis=1 (fit along xy-axis).
       */

      int fitaxis=(fabs(alltrack[i].ang-alltrack[j].ang)>90.0 || fabs(alltrack[i].ang+alltrack[j].ang)>90.0);
      if(fitaxis){
        swap(pos[0][0][0], pos[0][0][1]);
        swap(pos[0][1][0], pos[0][1][1]);
        swap(pos[1][0][0], pos[1][0][1]);
        swap(pos[1][1][0], pos[1][1][1]);
      }

      //Swap positions in order
      if(pos[0][0][1]>pos[0][1][1]) swap(pos[0][0], pos[0][1]);
      if(pos[1][0][1]>pos[1][1][1]) swap(pos[1][0], pos[1][1]);
      if(pos[0][0][1]>pos[1][0][1]) swap(pos[0]   , pos[1]   );
      if(pos[0][1][1]>pos[1][0][1]) continue; //Condidion (2)

      //Note: dist is distance between two tracks at the middle points
      float dist = fabs(((pos[0][0][0]-pos[0][1][0])/(pos[0][0][1]-pos[0][1][1])
            +(pos[1][0][0]-pos[1][1][0])/(pos[1][0][1]-pos[1][1][1]))
          *(pos[0][1][1]-pos[1][0][1])/2.0-pos[0][1][0]+pos[1][0][0]);

#ifdef DEBUG_TWODIMRECON
      cout << "dist=" << dist << endl;
#endif


      if(dist>dist_th) continue; //Condition (3)

      //If more than one candidate is found, ...
      if(!candidate || joilik>sqrt(angle*angle/angle_th/angle_th+dist*dist/dist_th/dist_th)){
        *a        = i;
        *b        = j;
        joilik    = sqrt(angle*angle/angle_th/angle_th+dist*dist/dist_th/dist_th);
        candidate = true;
      }
    }
  }

  return candidate;
};

void fGetHitFromTrack(vector<Hit> &hit, const vector<Track> &track){
  vector<Hit>::iterator it;
  for(int i=0; i<(int)track.size(); i++){
    for(int j=0; j<(int)track[i].hitid.size(); j++){
      for(it=hitcls_for_joint.begin(); it!=hitcls_for_joint.end(); it++){
        if(it->id==track[i].hitid[j]){
          hit.push_back(*it);
          break;
        }
      }
    }
  }
};

void fFitAllHit(
    const vector<Hit> &hit, vector<float> &par, int axis, 
    bool disp=false, int mod=MOD_ONAXIS_WM)
{

#ifdef DEBUG_TWODIMRECON
  cout <<"fFitAllHit::Start" << endl;
#endif

  vector<float> x, y, xe, ye;
  const int Nhits = (int)hit.size();
  for(int i=0; i<Nhits; i++){
    int reconpln, reconch;
    detdim->GetReconPlnCh(mod, hit[i].view, hit[i].pln, hit[i].ch, axis, &reconpln, &reconch);
    x.push_back (zposi   (mod, hit[i].view, reconpln, reconch, axis));
    y.push_back (xyposi  (mod, hit[i].view, reconpln, reconch, axis));
    xe.push_back(scithick(mod, hit[i].view, reconpln, y.at(i), axis)/sqrt(3));
    ye.push_back(sciwidth(mod, hit[i].view, reconpln, reconch, axis)/sqrt(3));

#ifdef DEBUG_TWODIMRECON
    cout
      << " mod="  << mod
      << " view=" << hit[i].view
      << " pln="  << hit[i].pln
      << " ch="   << hit[i].ch
      << endl;
#endif

  }

  TGraphErrors* g = new TGraphErrors(Nhits, &x[0], &y[0], &xe[0], &ye[0]);
  TF1 *f = new TF1("f", "[0]+[1]*x");
  const float deltamin = 0.1;
  const float deltax   = (fabs(x[Nhits-1]-x[0])<deltamin ? deltamin : (x[Nhits-1]-x[0])); //Avoid initial parameter of 0
  const float deltay   = (fabs(y[Nhits-1]-y[0])<deltamin ? deltamin : (y[Nhits-1]-y[0]));

  f->SetParameters(y[0]-x[0]*deltay/deltax, deltay/deltax);
  g->Fit("f", "WQ"); //Ignore errors (Any better solutions?)

  if(disp){
    TCanvas* c_allhit = new TCanvas("c_allhit", "c_allhit");
    g->Draw("AP");
    c_allhit->Update();
    getchar();
    delete c_allhit;
  }
  par.push_back(f->GetParameter(0));
  par.push_back(f->GetParameter(1));

#ifdef DEBUG_TWODIMRECON
  cout
    << " par0=" << f->GetParameter(0)
    << " par1=" << f->GetParameter(1)
    << " ang="  << f->GetParameter(1)*180./PI
    << endl;
#endif 

  f->Delete();
  g->Delete();
};

Track fGetConnectedTrack(
    const vector<Hit> &hit, const vector<Track> &track, 
    const vector<float> &par, int axis, int mod)
{
  Track tr;
  tr.clear();
  tr.mod  = mod;
  tr.view = hit[0].view;

  for(int i=0; i<(int)track.size(); i++){
    if(track[i].veto) tr.veto = true;
    if(track[i].edge) tr.edge = true;
    if(track[i].stop) tr.stop = true;
  }

  if(axis==0){
    tr.intcpt = par[0];
    tr.slope  = par[1];
    tr.ang    = atan(tr.slope)*180./PI;

    //Initialization
    int tmppln, tmpch;
    detdim->GetReconPlnCh(mod, hit[0].view, hit[0].pln, hit[0].ch, 0, &tmppln, &tmpch);
    tr.ipln = tmppln; //hit[0].pln;
    tr.fpln = tmppln; //hit[0].pln;
    tr.iz   = zposi(mod, hit[0].view, tmppln, tmpch);
    tr.fz   = zposi(mod, hit[0].view, tmppln, tmpch);
    tr.ixy  = par[0]+par[1]*tr.iz;
    tr.fxy  = par[0]+par[1]*tr.fz;

    for(int i=0; i<(int)hit.size(); i++){
      int reconpln, reconch;
      detdim->GetReconPlnCh(mod, hit[i].view, hit[i].pln, hit[i].ch, 0, &reconpln, &reconch);

      if(tr.iz>zposi(mod, hit[i].view, reconpln, reconch)){
        tr.ipln = reconpln; //hit[i].pln;
        tr.iz   = zposi(mod, hit[i].view, reconpln, reconch);
        tr.ixy  = par[0]+par[1]*tr.iz;
      }
      else if(tr.fz<zposi(mod, hit[i].view, reconpln, reconch)){
        tr.fpln = reconpln; //hit[i].pln;
        tr.fz   = zposi(mod, hit[i].view, reconpln, reconch);
        tr.fxy  = par[0]+par[1]*tr.fz;
      }
    }
  }
  else{
    if(par[1]==0){
      tr.intcpt = -par[0]*1e+8;
      tr.slope  = 1e+8;
      tr.ang    = 90.0;
    }
    else{
      tr.intcpt = -par[0]/par[1];
      tr.slope  = 1.0/par[1];
      tr.ang    = atan(tr.slope)*180.0/PI;
    }

    //Initialization
    int tmppln, tmpch;
    detdim->GetReconPlnCh(mod, hit[0].view, hit[0].pln, hit[0].ch, 0, &tmppln, &tmpch);
    tr.ipln = tmppln; //hit[0].pln;
    tr.fpln = tmppln; //hit[0].pln;
    tr.ixy  = xyposi(mod, hit[0].view, tmppln, tmpch);
    tr.fxy  = xyposi(mod, hit[0].view, tmppln, tmpch);
    tr.iz   = par[0]+par[1]*tr.ixy;
    tr.fz   = par[0]+par[1]*tr.fxy;

    for(int i=0; i<(int)hit.size(); i++){
      int reconpln, reconch;
      detdim->GetReconPlnCh(mod, hit[i].view, hit[i].pln, hit[i].ch, 0, &reconpln, &reconch);

      if( (par[1]>0) != (tr.ixy<xyposi(mod, hit[i].view, reconpln, reconch)) ){
        tr.ipln = reconpln; //hit[i].pln;
        tr.ixy  = xyposi(mod, hit[i].view, reconpln, reconch);
        tr.iz   = par[0]+par[1]*tr.ixy;
      }
      else if( (par[1]>0) != (tr.fxy>xyposi(mod, hit[i].view, reconpln, reconch)) ){
        tr.fpln = reconpln; //hit[i].pln;
        tr.fxy  = xyposi(mod, hit[i].view, reconpln, reconch);
        tr.fz   = par[0]+par[1]*tr.fxy;
      }
    }
  }

  for(int i=0; i<(int)hit.size(); i++){
    tr.hitid.push_back(hit[i].id);
    bool isohit = true;

    //Here alltrack does not include tracks to be connected.
    for(int j=0; j<(int)alltrack.size(); j++){ 
      for(int k=0; k<(int)alltrack[j].hitid.size(); k++){
        if(hit[i].id==alltrack[j].hitid[k]){
          isohit = false;
          break;
        }
      }
      if(!isohit) break;
    }
    tr.isohit.push_back(isohit);
  }

#ifdef DEBUG_TWODIMRECON
  cout
    << " New trk:"
    << " ipos={" << tr.iz << "," << tr.ixy << "}"
    << " fpos={" << tr.fz << "," << tr.fxy << "}"
    << " ipln="  << tr.ipln
    << " fpln="  << tr.fpln
    << " slope=" << tr.slope
    << " ang="   << tr.ang
    << endl;
#endif

  return tr;
};

bool fConnectTracks(int mod, float angle_th = 20., float dist_th = 50.){
  /* Connect break two tracks due to gap between layers in water module.
   * Procedure:
   * - Find nearest tracks
   * - Get hit info from two tracks to be connected
   * - Fit all hits to be connected
   * - Get connected trak info
   * If no candidate is found, return false
   */

  int a, b;
  if(!fFindNearestTracks(&a, &b, angle_th, dist_th)){
    return false;
  }

  vector<Hit> hit_connect;
  vector<Track> track_connect;
  track_connect.push_back(alltrack[a]);
  track_connect.push_back(alltrack[b]);
  fGetHitFromTrack(hit_connect, track_connect);

  if((int)hit_connect.size()<2){
    return false;
  }

  /*  Let theta to be mean of two tracks,
   *  if -45 < theta < 45 then fitaxis=0 (fit along z-axis),
   *  else fitaxis=1 (fit along xy-axis).
   */

  int fitaxis = (fabs(alltrack[a].ang-alltrack[b].ang)>90.0 || fabs(alltrack[a].ang+alltrack[b].ang)>90.0);
  vector<float> par;
  fFitAllHit(hit_connect, par, fitaxis, false,  mod);

  alltrack.erase(alltrack.begin()+a);
  alltrack.erase(alltrack.begin()+b);
  alltrack.push_back(fGetConnectedTrack(hit_connect, track_connect, par, fitaxis, mod));

  return true;
};

void fGetTrackInfo(char* buf, int ievt, int Nconnect){
  int ntracks = (int)alltrack.size();
  int nview = 0;

  for(int i=0; i<ntracks; i++){
    nview += alltrack[i].view;
  }

  sprintf(buf, "%d\t%d\t%d\t%d", ievt, ntracks, nview, Nconnect);
}



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
  double x=0.,y=0.,z=0.;

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
  else{
    x=0.; y=0.; z=0.;
    return -1.;
  }
 

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


void Add_RandomNoise(int id,int targetmod,int cyc,EventSummary* evt,double mean){
  //cout << " ====== Add_RandomNoise :  mod=" << targetmod << " cyc=" << cyc << " ======== " << endl;
  int mod;
  int num_channel;
  int norm_view;
  int norm_pln;
  int offset_ch;
  int num_hits = (int)allhit.size();
  if(id==0){mod=21;num_channel=1280;norm_view=640;norm_pln=80;offset_ch=0;}
  else     {mod=16;num_channel=1136;norm_view=568;norm_pln=32;offset_ch=8;}

  if(mod!=targetmod){ return; }


  int num_randhit;
  if     (id==0){num_randhit =  1;}
  else if(id==1){num_randhit =  5;}
  else if(id==2){num_randhit = 10;}
  else if(id==3){num_randhit = 15;}
  else if(id==4){num_randhit = 20;}
  else if(id==5){num_randhit = 25;}
  else{ return; }


  int cnt=0;
  while(cnt<num_randhit){
    int ch_id = rand()%num_channel;
    int view = ch_id/norm_view;
    int pln  = (ch_id%norm_view+offset_ch)/norm_pln;
    int ch   = (ch_id%norm_view+offset_ch)%norm_pln;
    if(mod==16&&pln==0) ch -= offset_ch;

    bool is_in_cls = false;
    for(int j=0;j<num_hits;j++){
      if(
          (mod  == allhit[j].mod ) &&
          (view == allhit[j].view) &&
          (pln  == allhit[j].pln ) &&
          (ch   == allhit[j].ch  ))
      {
        is_in_cls = true;
        break;
      }
    }
    //cout << "DEBUG "
    //  << " cnt:"   << cnt 
    //  << " ch_id:" << ch_id 
    //  << " num_hits:" << num_hits
    //  << endl;
    if(!is_in_cls){
      Hit randhit;
      randhit.clear();
      randhit.mod  = mod;
      randhit.id   = (int)evt->NModHits(mod,cyc);
      randhit.view = view;
      randhit.pln  = pln;
      randhit.ch   = ch;
      randhit.used = 0;
      randhit.pe   = mean;
      randhit.time = 0.;
      allhit.push_back(randhit);

      double posx,posy,posz,posxy;
      detdim->GetPosInMod(mod,pln,view,ch,&posx,&posy,&posz);
      if(view==0){posxy=posy;}else{posxy=posx;}

      HitSummary* hitsum = new HitSummary();
      hitsum->Clear("C");
      hitsum->mod  = randhit.mod;
      hitsum->view = randhit.view;
      hitsum->pln  = randhit.pln;
      hitsum->ch   = randhit.ch;
      hitsum->xy   = posxy;
      hitsum->z    = posz;
      hitsum->time = randhit.time;
      hitsum->pe   = randhit.pe;
      hitsum->lope = randhit.pe;
      evt->AddModHit(hitsum,mod,4);
      delete hitsum;
      
      //cout
      //  << " mod:"  << randhit.mod
      //  << " view:" << randhit.view
      //  << " pln:"  << randhit.pln
      //  << " ch:"   << randhit.ch
      //  << " time:" << randhit.time
      //  << endl;
    }
    cnt++;
  }
}


void Add_Crosstalk(double crosstalk,int targetmod,int cyc,EventSummary* evt){
  if(targetmod!=MOD_B2_WM){ return; }

  //cout << " ====== Add_Crosstalk :  mod=" << targetmod << " cyc=" << cyc << " ======== " << endl;
  int num_hits = (int)allhit.size();
  for(int j=0;j<num_hits;j++){
    if((allhit[j].mod==MOD_B2_WM)&&(allhit[j].ch>=40))
    {
      int    id   = allhit[j].id;
      int    mod  = allhit[j].mod;
      int    view = allhit[j].view;
      int    pln  = allhit[j].pln;
      int    ch   = allhit[j].ch;
      double pe   = allhit[j].pe;
      double time = allhit[j].time;

      HitSummary* hitsum = evt->GetModHit(id,mod,cyc);
      int gridx1 = hitsum->gridcell_id_x1;
      int gridx2 = hitsum->gridcell_id_x2;
      int gridy1 = hitsum->gridcell_id_y1;
      int gridy2 = hitsum->gridcell_id_y2;

      int grid_ch[2];
      if     (view==0){grid_ch[0]=gridy1;grid_ch[1]=gridy2;}
      else if(view==1){grid_ch[0]=gridx1;grid_ch[1]=gridx2;}

      double pe_crosstalk = pe*crosstalk;
      if(pe_crosstalk>hitpe_threshold_WM){

        for(int igrid=0;igrid<2;igrid++){
          if(grid_ch[igrid]<0||grid_ch[igrid]>=20){continue;}

          Hit hit_crosstalk;
          hit_crosstalk.clear();
          hit_crosstalk.mod  = mod;
          hit_crosstalk.id   = (int)evt->NModHits(mod,cyc);
          hit_crosstalk.view = 1-view;
          hit_crosstalk.pln  = pln;
          hit_crosstalk.ch   = grid_ch[igrid]+40+20*igrid;
          hit_crosstalk.used = 0;
          hit_crosstalk.pe   = pe_crosstalk;
          hit_crosstalk.time = time;
      
          bool is_in_cls = false;
          for(int j=0;j<num_hits;j++){
            if(
                (hit_crosstalk.mod  == allhit[j].mod ) &&
                (hit_crosstalk.view == allhit[j].view) &&
                (hit_crosstalk.pln  == allhit[j].pln ) &&
                (hit_crosstalk.ch   == allhit[j].ch  ))
            {
              is_in_cls = true;
              break;
            }
          }
          if(is_in_cls){continue;}

          allhit.push_back(hit_crosstalk);

          double posx,posy,posz,posxy;
          detdim->GetPosInMod(
              hit_crosstalk.mod,
              hit_crosstalk.pln,
              hit_crosstalk.view,
              hit_crosstalk.ch,
              &posx,&posy,&posz);
          if(view==0){posxy=posy;}else{posxy=posx;}

          HitSummary* hitsum_crosstalk = new HitSummary();
          hitsum_crosstalk->Clear("C");
          hitsum_crosstalk->mod  = hit_crosstalk.mod;
          hitsum_crosstalk->view = hit_crosstalk.view;
          hitsum_crosstalk->pln  = hit_crosstalk.pln;
          hitsum_crosstalk->ch   = hit_crosstalk.ch;
          hitsum_crosstalk->xy   = posxy;
          hitsum_crosstalk->z    = posz;
          hitsum_crosstalk->time = hit_crosstalk.time;
          hitsum_crosstalk->pe   = hit_crosstalk.pe;
          hitsum_crosstalk->lope = hit_crosstalk.pe;
          evt->AddModHit(hitsum_crosstalk,mod,4);
          delete hitsum_crosstalk;

        }
      }
    }
  }
}


#endif
