#ifndef __THREEDIMRECOFUNC_HXX__
#define __THREEDIMRECOFUNC_HXX__

#include "ThreeDimRecon_Class.hxx"
#include "TwoDimRecon.hxx"


int max_dif_cyc = 0;

float PE(float pe,float lope, int mod, int view, int pln, int ch){
  float Pe;
  if(1.422*lope+pe<100) Pe = pe;
  else                  Pe = 1.422*lope;
  if((mod==MOD_PM || mod==MOD_B2_CH) && pln>0 && ch>=8 && ch<24){
    Pe = Pe/2;
  }

  return Pe;
};

int pdg2num(int pdg){
  int num;
  if     (abs(pdg)==13)   num = 0;
  else if(abs(pdg)==211)  num = 1;
  else if(abs(pdg)==321)  num = 2;
  else if(abs(pdg)==2212) num = 3;
  else if(abs(pdg)==2112) num = 4;
  else if(abs(pdg)==11)   num = 5;
  else                    num = 6;

  return num;
};

int num2pdg(int *trkpdg){
  int num=0;
  int pdg;

  for(int i=1; i<6; i++){
    if(trkpdg[num]<=trkpdg[i]) num = i;
  }

  if(trkpdg[num]==0) num = 6;

  if     (num==0) pdg = 13;
  else if(num==1) pdg = 211;
  else if(num==2) pdg = 321;
  else if(num==3) pdg = 2212;
  else if(num==4) pdg = 2112;
  else if(num==5) pdg = 11;
  else if(num==6) pdg = 0;

  return pdg;
};

bool withendPM(const TrackPM& left, const TrackPM& right){
  if(left.fpln != right.fpln)      return left.fpln > right.fpln;
  else if(left.ipln != right.ipln) return left.ipln < right.ipln;
  else                             return fabs(left.ang) < fabs(right.ang);
};

bool withendWM(const TrackWM& left, const TrackWM& right){
  if(left.fpln != right.fpln)      return left.fpln > right.fpln;
  else if(left.ipln != right.ipln) return left.ipln < right.ipln;
  else                             return fabs(left.ang) < fabs(right.ang);
};

bool withcenter(const TrackIng& left, const TrackIng& right){
  if(abs(left.mod-INGMODNUM_mid) != abs(right.mod-INGMODNUM_mid)){
    return abs(left.mod-INGMODNUM_mid) < abs(right.mod-INGMODNUM_mid);
  }
  else if(left.ipln != right.ipln){
    return left.ipln < right.ipln;
  }
  else if(left.fpln != right.fpln){
    return left.fpln > right.fpln;
  }
  else{
    return fabs(left.ang) < fabs(right.ang);
  }
};

void fWMSortTrack(vector<TrackWM> &a){
  std::stable_sort(a.begin(), a.end(), withendWM);
};

void fPMSortTrack(vector<TrackPM> &a){
  std::stable_sort(a.begin(), a.end(), withendPM);
};

void fIngSortTrack(vector<TrackIng> &a){
  std::stable_sort(a.begin(), a.end(), withcenter);
};

double CLi(float pe, bool scibar){
  return 1.;
};

double TTCL; //For mucl
int    nCL;  //For mucl

void initMuCL(){
  TTCL = 1;
  nCL  = 0;
};

void addMuCL(vector<Hits> &allhit, float ang,int ipln, float intcpt, float slope){
  bool  hitpln[Cpln];
  int   stype[2][Cpln];
  float plnpe[Cpln];
  float corr = 1;

  memset(hitpln,false,sizeof(hitpln));
  memset(plnpe,0,sizeof(plnpe));
  memset(stype,0,sizeof(stype));

  for(int k=0; k<(int)allhit.size(); k++){
    if(!allhit[k].isohit)     continue; //Test
    if(allhit[k].pln==ipln+1) continue; //Test
    if(allhit[k].pln==ipln)   continue;

    if(allhit[k].view==0){
      corr=exp(-((intcpt+slope*zposi(allhit[k].mod,
				     allhit[k].view,
				     allhit[k].pln,
				     allhit[k].ch)))/2417);
    }
    else{
      corr=exp(-(1200-(intcpt+slope*zposi(allhit[k].mod,
					  allhit[k].view,
					  allhit[k].pln,
					  allhit[k].ch)))/2417);
    }
    
    hitpln[allhit[k].pln] = true;
    plnpe[allhit[k].pln] += PE(allhit[k].pe,
			       allhit[k].lope,
			       allhit[k].mod,
			       allhit[k].view,
			       allhit[k].pln,
			       allhit[k].ch)
      *cos(ang*PI/180.)/corr;

    if(allhit[k].pln>0 && allhit[k].ch>=8 && allhit[k].ch<24){
      stype[0][allhit[k].pln]++;
    }
    else{
      stype[1][allhit[k].pln]++;
    }
  }

  for(int i=0; i<Cpln; i++){
    if(hitpln[i]){
      nCL++;
      if(stype[1][i]<=stype[0][i]){
        TTCL *= CLi(plnpe[i],true);
      }
      else{
        TTCL *= CLi(plnpe[i],false);
      }
    }
  }
};

void addMuCLRe(vector<Hits> &allhit, float ang,int ipln, float intcpt, float slope, int upln){
  bool  hitpln[Cpln];
  int   stype[2][Cpln];
  float plnpe[Cpln];
  float corr = 1;

  memset(hitpln,false,sizeof(hitpln));
  memset(plnpe,0,sizeof(plnpe));
  memset(stype,0,sizeof(stype));

  for(int k=0; k<(int)allhit.size(); k++){
    if(!allhit[k].isohit)     continue; //Test
    if(allhit[k].pln==ipln+1) continue; //Test
    if(allhit[k].pln==ipln)   continue;
    if(allhit[k].pln<=upln)   continue; //Remove overlap hits

    if(allhit[k].view==0){
      corr = exp(-((intcpt+slope*zposi(allhit[k].mod,
				       allhit[k].view,
				       allhit[k].pln,
				       allhit[k].ch)))/2417);
    }
    else{
      corr = exp(-(1200-(intcpt+slope*zposi(allhit[k].mod,
					    allhit[k].view,
					    allhit[k].pln,
					    allhit[k].ch)))/2417);
    }

    hitpln[allhit[k].pln] = true;
    plnpe[allhit[k].pln] += PE(allhit[k].pe,
			       allhit[k].lope,
			       allhit[k].mod,
			       allhit[k].view,
			       allhit[k].pln,
			       allhit[k].ch)
      *cos(ang*PI/180.)/corr;

    if(allhit[k].pln>0 && allhit[k].ch>=8 && allhit[k].ch<24){
      stype[0][allhit[k].pln]++;
    }
    else{
      stype[1][allhit[k].pln]++;
    }
  }

  for(int i=0; i<Cpln; i++){
    if(hitpln[i]){
      nCL++;
      if(stype[1][i]<=stype[0][i]){
        TTCL *= CLi(plnpe[i],true);
      }
      else{
        TTCL *= CLi(plnpe[i],false);
      }
    }
  }
};


void addIngMuCL(vector<Hits> &allhit, float ang,float intcpt, float slope){
  bool  hitpln[Cpln];
  float plnpe[Cpln];
  float corr = 1;

  memset(hitpln,false,sizeof(hitpln));
  memset(plnpe,0,sizeof(plnpe));

  for(int k=0; k<(int)allhit.size(); k++){
    if(!allhit[k].isohit) continue; //Test

    if(allhit[k].view==0){
      corr = exp(-((intcpt+slope*zposi(allhit[k].mod,
				       allhit[k].view,
				       allhit[k].pln,
				       allhit[k].ch)))/2417);
    }
    else{
      corr = exp(-(1200-(intcpt+slope*zposi(allhit[k].mod,
					    allhit[k].view,
					    allhit[k].pln,
					    allhit[k].ch)))/2417);
    }

    hitpln[allhit[k].pln] = true;
    plnpe[allhit[k].pln] += PE(allhit[k].pe,
			       allhit[k].lope,
			       allhit[k].mod,
			       allhit[k].view,
			       allhit[k].pln,
			       allhit[k].ch)
      *cos(ang*PI/180.)/corr;
  }

  for(int i=0; i<Cpln; i++){
    if(hitpln[i]){
      nCL++;
      TTCL *= CLi(plnpe[i],false);
    }
  }
};

double calcMuCL(){
  double lncli = -log(TTCL);
  double mucl = 0;
  double kaijo = 1;

  for(int m=0; m<nCL; m++){
    if(m!=0) kaijo *= m;
    mucl += pow(lncli,m)/kaijo;
  };

  mucl = TTCL*mucl;

  return mucl;
};

int Ingmod(float hintcpt, float hslope, float vintcpt, float vslope){
  int   incmod = -1;
  float ypos = zposi(3,0,0,0)*hslope+hintcpt;
  float xpos = zposi(3,1,0,0)*vslope+vintcpt;
  if(ypos>-750 && ypos<750){
    if(xpos>-5250 && xpos<5250){
      incmod = (xpos+5250.)/1500.;
    }
  }

  return incmod;
};

void addIngHitAng(int incmod,Trk &trk, TrackWM &htrk, TrackWM &vtrk, int ivi, int evi){
  return;

  /*
  if(incmod<0 || incmod>=7) return;
  bool hitpln[2];
  int hitch[2];
  float hitpe[2];
  float slope, intcpt;
  float xing, ying, pldist;
  float xip[40], yip[40], xeip[40], yeip[40];
  int ntp = 0;

  for(int view=ivi; view<=evi; view++){
    memset(hitpln,false,sizeof(hitpln));
    memset(hitch,0,sizeof(hitch));
    memset(hitpe,0,sizeof(hitpe));

    if(view==0){
      slope = trk.slopex;
      intcpt = trk.intcptx;
    }
    else{
      slope=trk.slopey;
      intcpt=trk.intcpty-(incmod-3)*1500;
    }
    for(int pln=0; pln<=1; pln++){
      for(int ch=0; ch<24; ch++){
        xing = zposi(incmod,view,pln);
        ying = xyposi(incmod,pln,ch);
        pldist = fabs(slope*xing-ying+intcpt)/sqrt(1+slope*slope);
        if(
          nonrechits[incmod][view][pln][ch]>2.5 && 
          nonrechits[incmod][view][pln][ch]>hitpe[pln] && pldist<75)
        {
          hitch[pln]  = ch;
          hitpe[pln]  = nonrechits[incmod][view][pln][ch];
          hitpln[pln] = true;
        }
      }
    }
    if(hitpln[0]){
      ntp=0;

      int pst,ped;
      if(view==0){
        pst = htrk.ipln;
        ped = htrk.fpln;
      }
      else{
        pst = vtrk.ipln;
        ped = vtrk.fpln;
      }

      for(int p=pst; p<=ped; p++){
        xip[ntp]=zposi(15,view,p);
        yip[ntp]=intcpt+slope*xip[ntp];
        if(yip[ntp]>400&&yip[ntp]<800){
          xeip[ntp] = 6.5;
          yeip[ntp] = 12.5;
        }
        else{
          xeip[ntp] = 5.;
          yeip[ntp] = 25.;
        }


        xeip[ntp] = 1.5;
        yeip[ntp] = 12.5;

        ntp++;
      }

      for(int p=0; p<=1; p++){
        if(hitpln[p]){
          xip[ntp]  = zposi(incmod,view,p);
          yip[ntp]  = xyposi(incmod,p,hitch[p]);
          xeip[ntp] = 5.;
          yeip[ntp] = 25.;
          ntp++;
        }
      }

      TGraphErrors *gip = new TGraphErrors(ntp,xip,yip,xeip,yeip);
      TF1 *fip = new TF1("fip","[0]+[1]*x");
      fip->SetParameters(
          yip[0]-xip[0]*(yip[ntp-1]-yip[0])/(xip[ntp-1]-xip[0]),
          (yip[ntp-1]-yip[0])/(xip[ntp-1]-xip[0]));
      gip->Fit("fip","Q");

      TF1 *funcip    = gip->GetFunction("fip");
      float intcptip = funcip->GetParameter(0);
      float slopeip  = funcip->GetParameter(1);

      gip->Delete();
      fip->Delete();


      if(view==0){
        trk.slopex=slopeip;
        trk.intcptx=intcptip;
        trk.thetax=atan(slopeip)/PI*180.;
        trk.angle=180./PI*atan(
            sqrt(pow(tan(trk.thetax*PI/180.),2)+pow(tan(trk.thetay*PI/180.),2)));
      }
      else{
        trk.slopey=slopeip;
        trk.intcpty=intcptip+(incmod-3)*1500;
        trk.thetay=atan(slopeip)/PI*180.;
        trk.angle=180./PI*atan(
            sqrt(pow(tan(trk.thetax*PI/180.),2)+pow(tan(trk.thetay*PI/180.),2)));
      }	
    }
  }
  */
};

void addIngHitMuCL(
    int incmod, int view, float ang,float intcpt, 
    float slope,float intcpt2, float slope2)
{
  return;

  /*
  bool  hitpln[2];
  float plnpe[2];
  float corr = 1;
  float xing, ying, pldist;
  memset(hitpln,false,sizeof(hitpln));
  memset(plnpe,0,sizeof(plnpe));

  if(incmod<0 || incmod>=7) return;
  for(int pln=0; pln<=1; pln++){
    if(view==0) corr = exp(-((intcpt2+slope2*zposi(incmod,view,pln)))/2417);
    else        corr = exp(-(1200-(intcpt2+slope2*zposi(incmod,view,pln)))/2417);

    for(int ch=0; ch<24; ch++){
      xing   = zposi(incmod,view,pln);
      ying   = xyposi(incmod,pln,ch);
      pldist = fabs(slope*xing-ying+intcpt)/sqrt(1+slope*slope);
      if(nonrechits[incmod][view][pln][ch]>2.5 && pldist<75){
	hitpln[pln]  = true;
	plnpe[pln]  += nonrechits[incmod][view][pln][ch]*cos(ang*PI/180.)/corr;
      }
    }
  }

  for(int i=0; i<2; i++){
    if((i==0 && hitpln[0]) || (i==1 && hitpln[0] && hitpln[1])){
      nCL++;
      TTCL *= CLi(plnpe[i],false);
    }
  }
  */
};


float veract(
    HitSummary* inghitsum, int wmmod, int xpln, int ypln, 
    float xch, float ych, 
    int diffpln=0., float diffxy=0., int usepe=0)
{

  int   view = inghitsum->view;
  int   pln  = inghitsum->pln;
  int   ch   = inghitsum->ch;
  float pe   = inghitsum->pe;

  if(pe < hitpe_threshold_WM)             pe = 0;
  if(usepe==0 && pe > hitpe_threshold_WM) pe = 1;

  int reconpln, reconch;
  int axis = 0;
  detdim -> GetReconPlnCh( wmmod, view, pln, ch, axis, &reconpln, &reconch );
  pln = reconpln;
  ch  = reconch;

  if(view==0){
    if(abs(pln-xpln)>diffpln||fabs(xyposi(wmmod,view,pln,ch,axis)-xch)>diffxy){
      pe=0;
    }
  }
  else{
    if(abs(pln-ypln)>diffpln||fabs(xyposi(wmmod,view,pln,ch,axis)-ych)>diffxy){
      pe=0;
    }
  }

  return pe;
};


float veract(
    int wmmod, int xpln, int ypln, float xch, float ych, 
    int diffpln=0, float diffxy=0, int usepe=0)
{

  float pe;
  float totpe=0;
  int   reconpln, reconch;
  for(int view=0; view<Cview; view++){
    for(int pln=0; pln<C_WMNumPln; pln++){
      for(int ch=0; ch<C_WMNumCh; ch++){
        DetectorDimension *fdim = new DetectorDimension();
        fdim->GetReconPlnCh(wmmod,view,pln,ch,0,&reconpln,&reconch);
        delete fdim;

        int axis = 0;
        pe = nonrechits[wmmod][view][reconpln][reconch];
        if(pe<hitpe_threshold_WM)             pe = 0;
        if(usepe==0 && pe>hitpe_threshold_WM) pe = 1;
        if(view==0){
          if(abs(pln-xpln)>diffpln||fabs(xyposi(wmmod,view,reconpln,reconch,axis)-xch)>diffxy){
            pe = 0;
	  }
        } 
        else{
          if(abs(pln-ypln)>diffpln||fabs(xyposi(wmmod,view,reconpln,reconch,axis)-ych)>diffxy){
            pe = 0;
	  }
        }
        totpe += pe;
      }
    }
  }

  return totpe;
};


void fAddPE(
    vector<Hits> &allhit, float ang, float &totalpe, int &totalhit, 
    int *trkpdg,int ipln, float intcpt, float slope)
{
  bool  hitpln[Cpln];
  float corr = 1;
  memset(hitpln,false,sizeof(hitpln));
  for(int k=0; k<(int)allhit.size(); k++){
    if(!allhit[k].isohit)     continue; //Test
    if(allhit[k].pln==ipln+1) continue; //Test
    if(allhit[k].pln==ipln)   continue;    

    if(allhit[k].view==0){
      corr = exp(-((intcpt+slope*zposi(allhit[k].mod,
				       allhit[k].view,
				       allhit[k].pln,
				       allhit[k].ch)))/2417);
    }
    else{
      corr = exp(-(1200-(intcpt+slope*zposi(allhit[k].mod,
					    allhit[k].view,
					    allhit[k].pln,
					    allhit[k].ch)))/2417);
    }

    totalpe += PE(allhit[k].pe,
		  allhit[k].lope,
		  allhit[k].mod,
		  allhit[k].view,
		  allhit[k].pln,
		  allhit[k].ch)
      *cos(ang*PI/180.)/corr;

    trkpdg[pdg2num(allhit[k].pdg)]++;
    hitpln[allhit[k].pln] = true;
  }

  for(int i=0; i<Cpln; i++){
    if(hitpln[i]) totalhit++;
  }
};

void GetNonRecHits(int wmmod,int ing_start, int ing_end,int pmmod, EventSummary* evt, int cyc)
{
  HitSummary* hitsum;
  int         imod, view, pln, ch;
  float       pe, lope;
  int         nwmhit,ninghit,npmhit,pdg;
  memset(nonrechits     ,0,sizeof(nonrechits     ));
  memset(nonrechits_lope,0,sizeof(nonrechits_lope));
  memset(nonrechits_pdg ,0,sizeof(nonrechits_pdg ));
  memset(nonrechits_id  ,0,sizeof(nonrechits_id  ));

  //------- Fill all hits -> Remove used hits -------
  //WaterModule
  if(wmmod==MOD_B2_WM){max_dif_cyc=1;}else{max_dif_cyc=0;}
  for(int icyc=0;icyc<=max_dif_cyc;icyc++){
    nwmhit = evt -> NModHits(wmmod, cyc+icyc);
    for(int ihit=0; ihit<nwmhit; ihit++){
      hitsum = (HitSummary*) (evt -> GetModHit(ihit,wmmod,cyc+icyc));
      view = hitsum->view;
      ch   = hitsum->ch;
      pln  = hitsum->pln;
      pe   = hitsum->pe;
      lope = hitsum->lope;

      if(badch->is_BadCh(hitsum->mod,view,pln,ch)) continue;
      if((hitsum -> NSimHits()) > 0){
        pdg = hitsum -> GetSimHit(0)->pdg;
      }
      else{
	pdg = 0;
      }

      int reconpln, reconch;
      int axis = 0;
      detdim -> GetReconPlnCh( wmmod, view, pln, ch, axis, &reconpln, &reconch );
      pln = reconpln;
      ch  = reconch;

      nonrechits     [wmmod][view][pln][ch] = pe;
      nonrechits_lope[wmmod][view][pln][ch] = lope;
      nonrechits_pdg [wmmod][view][pln][ch] = pdg;
      nonrechits_id  [wmmod][view][pln][ch] = ihit;

    }
  }

  //INGRID
  for(int mod=ing_start; mod<=ing_end; mod++){
    ninghit = evt -> NModHits(mod, cyc);
    for(int ihit=0; ihit<ninghit; ihit++){
      hitsum = (HitSummary*)(evt -> GetModHit(ihit,mod,cyc));

      view = hitsum->view;
      ch   = hitsum->ch;
      pln  = hitsum->pln;
      pe   = hitsum->pe;
      lope = hitsum->lope;

      if(badch->is_BadCh(hitsum->mod,view,pln,ch)) continue;

      if((hitsum -> NSimHits()) > 0){
        pdg = hitsum -> GetSimHit(0)->pdg;
      }
      else{
	pdg = 0;
      }
      nonrechits     [mod][view][pln][ch] = pe;
      nonrechits_lope[mod][view][pln][ch] = lope;
      nonrechits_pdg [mod][view][pln][ch] = pdg;
      nonrechits_id  [mod][view][pln][ch] = ihit;
    }
  }

  // Proton Module
  if(pmmod>0){
    npmhit = evt -> NModHits(pmmod, cyc);
    for(int ihit=0; ihit<npmhit; ihit++){
      hitsum = (HitSummary*)(evt -> GetModHit(ihit,pmmod,cyc));
      view = hitsum->view;
      ch   = hitsum->ch;
      pln  = hitsum->pln;
      pe   = hitsum->pe;
      lope = hitsum->lope;
      if(badch->is_BadCh(hitsum->mod,view,pln,ch)){ continue;}
      if((hitsum -> NSimHits()) > 0){
        pdg = hitsum -> GetSimHit(0)->pdg;
      }
      else{
        pdg = 0;
      }
      nonrechits     [pmmod][view][pln][ch] = pe;
      nonrechits_lope[pmmod][view][pln][ch] = lope;
      nonrechits_pdg [pmmod][view][pln][ch] = pdg;
      nonrechits_id  [pmmod][view][pln][ch] = ihit;
    }
  }

  //------- Remove the hits already used for tracking -------
  //WaterModule
  int hn,vn;
  hn = hwmtrack.size();
  vn = vwmtrack.size();
  for(int i=0; i<hn; i++){
    int nhit = hwmtrack[i].hit.size();
    for(int j=0; j<nhit; j++){
      view = hwmtrack[i].hit[j].view;
      pln  = hwmtrack[i].hit[j].pln;
      ch   = hwmtrack[i].hit[j].ch;
      nonrechits     [wmmod][view][pln][ch] = 0;
      nonrechits_lope[wmmod][view][pln][ch] = 0;
      nonrechits_pdg [wmmod][view][pln][ch] = 0;
      nonrechits_id  [wmmod][view][pln][ch] = 0;
    }
  }

  for(int i=0; i<vn; i++){
    int nhit = vwmtrack[i].hit.size();
    for(int j=0; j<nhit; j++){
      view = vwmtrack[i].hit[j].view;
      pln  = vwmtrack[i].hit[j].pln;
      ch   = vwmtrack[i].hit[j].ch;
      nonrechits     [wmmod][view][pln][ch] = 0;
      nonrechits_lope[wmmod][view][pln][ch] = 0;
      nonrechits_pdg [wmmod][view][pln][ch] = 0;
      nonrechits_id  [wmmod][view][pln][ch] = 0;
    }
  }

  //INGRID
  hn = hingtrack.size();
  vn = vingtrack.size();
  for(int i=0; i<hn; i++){
    int nhit = hingtrack[i].hit.size();
    for(int j=0; j<nhit; j++){
      imod = hingtrack[i].hit[j].mod;
      view = hingtrack[i].hit[j].view;
      pln  = hingtrack[i].hit[j].pln;
      ch   = hingtrack[i].hit[j].ch;
      nonrechits     [imod][view][pln][ch] = 0;
      nonrechits_lope[imod][view][pln][ch] = 0;
      nonrechits_pdg [imod][view][pln][ch] = 0;
      nonrechits_id  [imod][view][pln][ch] = 0;
    }
  }

  for(int i=0; i<vn; i++){
    int nhit = vingtrack[i].hit.size();
    for(int j=0; j<nhit; j++){
      imod = vingtrack[i].hit[j].mod;
      view = vingtrack[i].hit[j].view;
      pln  = vingtrack[i].hit[j].pln;
      ch   = vingtrack[i].hit[j].ch;
      nonrechits     [imod][view][pln][ch] = 0;
      nonrechits_lope[imod][view][pln][ch] = 0;
      nonrechits_pdg [imod][view][pln][ch] = 0;
      nonrechits_id  [imod][view][pln][ch] = 0;
    }
  }

  //Proton Module
  hn = hpmtrack.size();
  vn = vpmtrack.size();
  for(int i=0; i<hn; i++){
    int nhit = hpmtrack[i].hit.size();
    for(int j=0; j<nhit; j++){
      imod = hpmtrack[i].hit[j].mod;
      view = hpmtrack[i].hit[j].view;
      pln  = hpmtrack[i].hit[j].pln;
      ch   = hpmtrack[i].hit[j].ch;
      nonrechits     [imod][view][pln][ch] = 0;
      nonrechits_lope[imod][view][pln][ch] = 0;
      nonrechits_pdg [imod][view][pln][ch] = 0;
      nonrechits_id  [imod][view][pln][ch] = 0;
    }
  }

  for(int i=0; i<vn; i++){
    int nhit = vpmtrack[i].hit.size();
    for(int j=0; j<nhit; j++){
      imod = vpmtrack[i].hit[j].mod;
      view = vpmtrack[i].hit[j].view;
      pln  = vpmtrack[i].hit[j].pln;
      ch   = vpmtrack[i].hit[j].ch;
      nonrechits     [imod][view][pln][ch] = 0;
      nonrechits_lope[imod][view][pln][ch] = 0;
      nonrechits_pdg [imod][view][pln][ch] = 0;
      nonrechits_id  [imod][view][pln][ch] = 0;
    }
  } 
};


void SortArray(int nsize,float *array_key,float *array1,float *array2,float *array3){
  float tmp;
  for(int i=0; i<nsize; i++){
    for(int j=i+1; j<nsize; j++){
      if(array_key[i]>array_key[j]){
        tmp          = array_key[i];
        array_key[i] = array_key[j];
        array_key[j] = tmp;
        tmp          = array1[i];
        array1[i]    = array1[j];
        array1[j]    = tmp;
        tmp          = array2[i];
        array2[i]    = array2[j];
        array2[j]    = tmp;
        tmp          = array3[i];
        array3[i]    = array3[j];
        array3[j]    = tmp;
      }
    }
  }
};


/*
double calc_pathlength(double anglex,double angley,int mod,int view,int pln,int ch)
{
  //double anglex, angley;
  double angleX, angleY, angleZ;
  double angle1=0., angle2=0.;
  //double x1, x2, y1, y2;
  //double z1, z2;
  double path;
  //Thetax, Thetay
  //anglex = analyzed_trk[i].trk[t].thetax;
  //angley = analyzed_trk[i].trk[t].thetay;
  //Angle from each axis
  angleX = 180./PI*acos(tan(angley) / sqrt(1 + pow(tan(anglex),2) + pow(tan(angley),2)));
  if(angleX>90) angleX = 180. - angleX;
  angleY = 180./PI*acos(tan(anglex) / sqrt(1 + pow(tan(anglex),2) + pow(tan(angley),2)));
  if(angleY>90) angleY = 180. - angleY;
  angleZ = 180./PI*atan(sqrt(pow(tan(anglex),2) + pow(tan(angley),2)));
  if(!(angleX<85 && angleX>-85 && angleY<85 && angleY>-85)){
    angleX = -300; angleY = -300; angleZ = -300;
  }

  //Classify: Grid or Plane
  double scinti_width, scinti_thick;
  if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
    if(view==0&&ch< 40){angle1=angleZ;angle2=fabs(180./PI*anglex);}
    if(view==1&&ch< 40){angle1=angleZ;angle2=fabs(180./PI*angley);}
    if(view==0&&ch>=40){angle1=angleY;angle2=fabs(180./PI*anglex);}
    if(view==1&&ch>=40){angle1=angleX;angle2=fabs(180./PI*angley);}
    scinti_width = C_WMScintiWidth;
    scinti_thick = C_WMScintiThick;
  }
  else{
    if(view==0){angle1=angleZ;angle2=fabs(180./PI*anglex);}
    if(view==1){angle1=angleZ;angle2=fabs(180./PI*angley);}
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
  angle1 = PI/180.*(angle1+0.1); angle2 = PI/180.*(angle2+0.1);

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
*/

#endif
