#ifndef __THREEDIMRECO_HXX__
#define __THREEDIMRECO_HXX__

//#define USEINGTRK

//#define DEBUG_THREEDIMRECON

//INGRID
#include "TwoDimRecon.hxx"
#include "ThreeDimRecon_Functions.hxx"
#include "ThreeDimRecon_Class.hxx"

bool NOTJOINT=false;
bool USE_NON_3DTRK=false;

int cyc = -1;

//Vertexing : |diffZ_x| + |diffZ_y| < pln_th
int   pln_th0  = 2; //Vertexing ING
int   pln_th1  = 3; //Vertexing PM
int   pln_th2  = 9; //Vertexing WM
float ch_th0   = 150.; //Vertexing ING
float ch_th1   = 150.; //Vertexing PM
float ch_th2   = 150.; //Vertexing WM

//3D track match : diff_z <diff_th
int   diff_th0 = 3; //3D track match ING
int   diff_th1 = 3; //3D track match PM
int   diff_th2 = 9; //3D track match WM

//2D track match b/w modules
float ang_th0   =  35.; //2D track match WM-ING
float ang_th1   =  35.; //2D track match PM-WM
float ang_th2   =  35.; //2D track match PM-ING
float pos_th0   = 150.; //2D track match WM-ING
float pos_th1   = 150.; //2D track match PM-WM
float pos_th2   = 150.; //2D track match PM-ING
float hitpos_th2  =  75.; //2D track match WM-INGhit
//float pos_th_hit3  =  75.; //2D track match WM-PMhit

const float PEth_extrahits = 3.5;

bool fIngWMJoint(vector<TrackIng> &itrk, vector<TrackWM> &wmtrk, bool vertical, int mod){

#ifdef DEBUG_THREEDIMRECON
  cout << " -------------- fIngWMJoint " << vertical << " --------------- " << endl;
#endif

  float diff_ang  = -1e-5;
  float diff_pos  = -1e-5;
  float joilik    = -1e-5;
  int   joitra    = -1;
  bool  jointed   = false;
  bool  hasingtrk = false;

  float ang_th = ang_th0;
  float pos_th = pos_th0;

  for(int j=0; j<(int)itrk.size(); j++){
    if(itrk[j].mod!=INGMODNUM_mid){continue;}
    if(itrk[j].ipln>2) continue;

    jointed=false;
    for(int i=0; i<(int)wmtrk.size(); i++){

      double zpos_at_joint = -1.;
      if(MODE_DET==0){
        zpos_at_joint = ( C_INGHMotherPosZ + zposi(INGMODNUM_mid,0,0,0)
			  + C_PMMotherPosZ + zposi(WMMODNUM,1,C_WMNumPln*3-1,0) ) /2.;
        diff_pos = ( (wmtrk[i].intcpt + wmtrk[i].slope*(zpos_at_joint-C_PMMotherPosZ))
		     -(itrk[j].intcpt + itrk[j].slope*(zpos_at_joint-C_INGHMotherPosZ)) );
        if(vertical) diff_pos += C_PMMotherPosX - C_INGHMotherPosX;
        else         diff_pos += C_PMMotherPosY - C_INGHMotherPosY;
      }
      else if(MODE_DET==1){	      
        zpos_at_joint = (  C_B2INGPosZ + zposi(INGMODNUM_mid,0,0,0)
			 + C_B2WMPosZ  + zposi(WMMODNUM,1,C_WMNumPln*3-1,0) ) /2.;
        diff_pos = ( (wmtrk[i].intcpt + wmtrk[i].slope*(zpos_at_joint-C_B2WMPosZ))
		     -(itrk[j].intcpt + itrk[j].slope*(zpos_at_joint-C_B2INGPosZ)) );
        if(vertical) diff_pos += C_B2WMPosX - C_B2INGPosX;
        else         diff_pos += C_B2WMPosY - C_B2INGPosY;
      }

      diff_ang = wmtrk[i].ang-itrk[j].ang;

#ifdef DEBUG_THREEDIMRECON
      cout 
        << " wmtrk" << i
        << " itrk"  << j
        << " diff_ang=" << diff_ang
        << " diff_pos=" << diff_pos
        << " zpos_at_joint=" << zpos_at_joint
        << endl;
#endif


      if(fabs(diff_ang)<ang_th && fabs(diff_pos)<pos_th){
        if(jointed){
          if(joilik > sqrt( fabs(diff_ang)*fabs(diff_ang)/ang_th/ang_th
			   +fabs(diff_pos)*fabs(diff_pos)/pos_th/pos_th)){
            wmtrk[joitra].ing_trk=false;
	  }
          else continue;
        }

        if(wmtrk[i].ing_trk == false){
          wmtrk[i].ing_imod  = itrk[j].mod;
          wmtrk[i].ing_fmod  = itrk[j].mod;
          wmtrk[i].ing_ipln  = itrk[j].ipln;
          wmtrk[i].ing_fpln  = itrk[j].fpln;
          wmtrk[i].ing_trk   = true;
          wmtrk[i].ing_stop  = itrk[j].stop;	  
          wmtrk[i].ing_num   = j;
          wmtrk[i].diff_pos  = diff_pos;
          wmtrk[i].diff_ang  = diff_ang;
          wmtrk[i].diff_time = wmtrk[i].clstime - itrk[j].clstime;

          if(itrk[j].fpln==C_INGNumPln-1){
            wmtrk[i].iron_pene  = itrk[j].fpln - itrk[j].ipln-1;
	  }
          else{
            wmtrk[i].iron_pene  = itrk[j].fpln - itrk[j].ipln;
	  }
        }
	//If wmtrk[i] already has jointed ingtrk
        else{
          if(abs(itrk[j].mod-INGMODNUM_mid) <= abs(wmtrk[i].ing_fmod-INGMODNUM_mid)) continue;
          if(itrk[j].ipln < wmtrk[i].ing_fpln) continue;

          wmtrk[i].ing_fmod   = itrk[j].mod;
          wmtrk[i].ing_fpln   = itrk[j].fpln;
          wmtrk[i].ing_stop   = itrk[j].stop;
          wmtrk[i].diff_pos   = diff_pos;
          wmtrk[i].diff_ang   = diff_ang;
          wmtrk[i].diff_time  = wmtrk[i].clstime - itrk[j].clstime;

          if(itrk[j].fpln==C_INGNumPln-1){
            wmtrk[i].iron_pene  += itrk[j].fpln - itrk[j].ipln-1;
	  }
          else{
            wmtrk[i].iron_pene  += itrk[j].fpln - itrk[j].ipln;
	  }
        }

        jointed   = true;
        joitra    = i;
        joilik    = sqrt( fabs(diff_ang)*fabs(diff_ang)/ang_th/ang_th
			 +fabs(diff_pos)*fabs(diff_pos)/pos_th/pos_th);
        hasingtrk = true;

#ifdef DEBUG_THREEDIMRECON
        cout << " Jointed" << endl;
#endif

      }
    }
  }

//#ifdef USEINGTRK 
//
//  for(int i=0; i<(int)wmtrk.size(); i++){
//    if(wmtrk[i].ing_trk){
//      int inum = wmtrk[i].ing_num;
//      float xip[Cpln], yip[Cpln], xeip[Cpln], yeip[Cpln];
//      int ntp = 0;
//      
//      for(int hitnum=0; hitnum<(int)wmtrk[i].hit.size(); hitnum++){
//        int view  = wmtrk[i].hit[hitnum].view;
//        int pln   = wmtrk[i].hit[hitnum].pln;
//        int ch    = wmtrk[i].hit[hitnum].ch;
//        xip[ntp]  = zposi   (mod,view,pln,ch,0);
//        yip[ntp]  = xyposi  (mod,view,pln,ch,0);
//        xeip[ntp] = scithick(mod,view,pln,yip[ntp],0);
//        yeip[ntp] = sciwidth(mod,view,pln,ch,0);
//
//        if(MODE_DET==0){
//          xip[ntp] += C_PMMotherPosZ;
//          if(vertical) yip[ntp] += C_PMMotherPosX;
//          else         yip[ntp] += C_PMMotherPosY;
//        }
//        else if(MODE_DET==1){
//          xip[ntp] += C_B2WMPosZ;
//          if(vertical) yip[ntp] += C_B2WMPosX;
//          else         yip[ntp] += C_B2WMPosY;
//        }
//        ntp++;
//      }
//
//      for(int hitnum=0; hitnum<(int)itrk[inum].hit.size(); hitnum++){
//        int ingmod  = itrk[inum].mod ;
//        int view    = itrk[inum].hit[hitnum].view;
//        int pln     = itrk[inum].hit[hitnum].pln ;
//        int ch      = itrk[inum].hit[hitnum].ch  ;
//        xip[ntp] = zposi (ingmod,view,pln,ch,0);
//        yip[ntp] = xyposi(ingmod,view,pln,ch,0);
//        xeip[ntp] = scithick(ingmod,view,pln,yip[ntp],0);
//        yeip[ntp] = sciwidth(ingmod,view,pln,ch,0);
//
//        if(MODE_DET==0){
//          xip[ntp] += C_INGHMotherPosZ;
//          if(vertical) yip[ntp] += C_INGHMotherPosX;
//          else         yip[ntp] += C_INGHMotherPosY;
//        }
//        else if(MODE_DET==1){
//          xip[ntp] += C_B2INGPosZ;
//          if(vertical) yip[ntp] += C_B2INGPosX;
//          else         yip[ntp] += C_B2INGPosY;
//        }
//
//        if(vertical){
//          yip[ntp] += (ingmod-INGMODNUM_mid)*C_INGSpace;
//        }
//
//        ntp++;
//      }
//
//      SortArray(ntp,xip,yip,xeip,yeip); 
//      int cnt=0;
//      float last_x=-1e+5;
//      float max_y =-1e+5;
//      float min_y = 1e+5;
//      for(int itp=0;itp<ntp;itp++){
//        if(last_x==xip[itp]){
//          if(max_y<yip[itp]){max_y=yip[itp];}
//          if(min_y>yip[itp]){min_y=yip[itp];}
//          cnt++;
//          yip [itp-1] *= cnt;
//          yip [itp-1] += yip[itp];
//          yip [itp-1] /= cnt+1;
//          yeip[itp-1] = max_y-min_y+yeip[itp];
//          for(int jtp=itp;jtp<ntp-1;jtp++){
//            xip [jtp] = xip [jtp+1];
//            yip [jtp] = yip [jtp+1];
//            xeip[jtp] = xeip[jtp+1];
//            yeip[jtp] = yeip[jtp+1];
//          }
//          ntp--;
//          itp--;
//        }
//        else{
//          cnt=0;
//          last_x=-1e+5;
//          max_y =-1e+5;
//          min_y = 1e+5;
//          last_x = xip[itp];
//        }
//      }
//
//      TGraphErrors *gip = new TGraphErrors(ntp,xip,yip,xeip,yeip);
//      TF1 *fip = new TF1("fip","[0]+[1]*x");
//      fip->SetParameters( yip[0]-xip[0]*(yip[ntp-1]-yip[0])/(xip[ntp-1]-xip[0]),
//			 (yip[ntp-1]-yip[0])/(xip[ntp-1]-xip[0]));
//      gip->Fit("fip","Q");
//
//      TF1 *funcip    = gip->GetFunction("fip");
//      float intcptip = funcip->GetParameter(0);
//      float slopeip  = funcip->GetParameter(1);
//#ifdef DEBUG_THREEDIMRECON
//      cout << " Fit:" 
//        << " slope=" << slopeip
//        << " ang="   << atan(slopeip)*180./PI
//        << endl;
//#endif
//      gip->Delete();
//      fip->Delete();
//
//      wmtrk[i].ang    = atan(slopeip)*180./PI;
//      wmtrk[i].slope  = slopeip;
//      wmtrk[i].intcpt = wmtrk[i].ixy-wmtrk[i].slope*wmtrk[i].iz;
//      //wmtrk[i].intcpt = intcptip;
//      //if     (MODE_DET==0) wmtrk[i].intcpt += slopeip*C_PMMotherPosZ;
//      //else if(MODE_DET==1) wmtrk[i].intcpt += slopeip*C_B2WMPosZ;
//
//      itrk[inum].ang    = atan(slopeip)*180./PI;
//      itrk[inum].slope  = slopeip;
//      itrk[inum].intcpt = itrk[inum].ixy-itrk[inum].slope*itrk[inum].iz;
//      //itrk[inum].intcpt = intcptip;
//      //if     (MODE_DET==0) itrk[inum].intcpt += slopeip*C_INGHMotherPosZ;
//      //else if(MODE_DET==1) itrk[inum].intcpt += slopeip*C_B2INGPosZ;
//
//    }
//  }
//#endif

  return hasingtrk;
};


bool fIngPMJoint(vector<TrackIng> &itrk, vector<TrackWM> &wmtrk, vector<TrackPM> &pmtrk, bool vertical, int mod){

  if(MODE_DET!=1){return false;}

#ifdef DEBUG_THREEDIMRECON
  cout << " -------------- fIngPMJoint " << vertical << " --------------- " << endl;
#endif

  float diff_ang  = -1e-5;
  float diff_pos  = -1e-5;
  float joilik    = -1e-5;
  int   joitra    = -1;
  bool  jointed   = false;
  bool  hasingtrk = false;

  float ang_th = ang_th2;
  float pos_th = pos_th2;

  for(int j=0; j<(int)itrk.size(); j++){
    if(itrk[j].mod!=INGMODNUM_mid){continue;}
    if(itrk[j].ipln>2) continue;
    bool wmjointed = false;
    for(int k=0;k<(int)wmtrk.size();k++){
      if(wmtrk[k].ing_num==j){wmjointed=true;break;}
    }
    if(wmjointed){continue;}

    jointed=false;
    for(int i=0; i<(int)pmtrk.size(); i++){
      if(pmtrk[i].wm_trk ){continue;}
      if(pmtrk[i].ing_trk){continue;}

      double zpos_at_joint = -1.;
      zpos_at_joint = (  C_B2INGPosZ + zposi(INGMODNUM_mid,0,0,0)
      		 + C_B2CHPosZ  + zposi(PMMODNUM,1,C_PMNumPln-1,0) ) /2.;
      diff_pos = ( (pmtrk[i].intcpt + pmtrk[i].slope*(zpos_at_joint-C_B2CHPosZ))
      	     -(itrk[j].intcpt + itrk[j].slope*(zpos_at_joint-C_B2INGPosZ)) );
      if(vertical) diff_pos += C_B2CHPosX - C_B2INGPosX;
      else         diff_pos += C_B2CHPosY - C_B2INGPosY;

      diff_ang = pmtrk[i].ang-itrk[j].ang;

#ifdef DEBUG_THREEDIMRECON
      cout 
        << " pmtrk" << i
        << " itrk"  << j
        << " diff_ang=" << diff_ang
        << " diff_pos=" << diff_pos
        << " zpos_at_joint=" << zpos_at_joint
        << endl;
#endif

      if(fabs(diff_ang)<ang_th && fabs(diff_pos)<pos_th){
        if(jointed){
          if(joilik > sqrt( fabs(diff_ang)*fabs(diff_ang)/ang_th/ang_th
			   +fabs(diff_pos)*fabs(diff_pos)/pos_th/pos_th)){
            pmtrk[joitra].ing_trk=false;
	  }
          else continue;
        }

        if(pmtrk[i].ing_trk == false){
          pmtrk[i].ing_imod  = itrk[j].mod;
          pmtrk[i].ing_fmod  = itrk[j].mod;
          pmtrk[i].ing_ipln  = itrk[j].ipln;
          pmtrk[i].ing_fpln  = itrk[j].fpln;
          pmtrk[i].ing_trk   = true;
          pmtrk[i].ing_stop  = itrk[j].stop;	  
          pmtrk[i].ing_num   = j;
          pmtrk[i].diff_pos  = diff_pos;
          pmtrk[i].diff_ang  = diff_ang;
          pmtrk[i].diff_time = pmtrk[i].clstime - itrk[j].clstime;

          if(itrk[j].fpln==C_INGNumPln-1){
            pmtrk[i].iron_pene  = itrk[j].fpln - itrk[j].ipln-1;
	  }
          else{
            pmtrk[i].iron_pene  = itrk[j].fpln - itrk[j].ipln;
	  }
        }
	//If pmtrk[i] already has jointed ingtrk
        else{
          if(abs(itrk[j].mod-INGMODNUM_mid) <= abs(pmtrk[i].ing_fmod-INGMODNUM_mid)) continue;
          if(itrk[j].ipln < pmtrk[i].ing_fpln) continue;

          pmtrk[i].ing_fmod   = itrk[j].mod;
          pmtrk[i].ing_fpln   = itrk[j].fpln;
          pmtrk[i].ing_stop   = itrk[j].stop;
          pmtrk[i].diff_pos   = diff_pos;
          pmtrk[i].diff_ang   = diff_ang;
          pmtrk[i].diff_time  = pmtrk[i].clstime - itrk[j].clstime;

          if(itrk[j].fpln==C_INGNumPln-1){
            pmtrk[i].iron_pene  += itrk[j].fpln - itrk[j].ipln-1;
	  }
          else{
            pmtrk[i].iron_pene  += itrk[j].fpln - itrk[j].ipln;
	  }
        }

        jointed   = true;
        joitra    = i;
        joilik    = sqrt( fabs(diff_ang)*fabs(diff_ang)/ang_th/ang_th
			 +fabs(diff_pos)*fabs(diff_pos)/pos_th/pos_th);
        hasingtrk = true;

#ifdef DEBUG_THREEDIMRECON
        cout << " Jointed" << endl;
#endif

      }
    }
  }

//#ifdef USEINGTRK 
//
//  for(int i=0; i<(int)pmtrk.size(); i++){
//    if(pmtrk[i].ing_trk&&!pmtrk[i].wm_trk){
//      int inum = pmtrk[i].ing_num;
//      float xip[Cpln], yip[Cpln], xeip[Cpln], yeip[Cpln];
//      int ntp = 0;
//      
//      for(int hitnum=0; hitnum<(int)pmtrk[i].hit.size(); hitnum++){
//        int view  = pmtrk[i].hit[hitnum].view;
//        int pln   = pmtrk[i].hit[hitnum].pln;
//        int ch    = pmtrk[i].hit[hitnum].ch;
//        xip [ntp] = zposi   (mod,view,pln,ch,0);
//        yip [ntp] = xyposi  (mod,view,pln,ch,0);
//        xeip[ntp] = scithick(mod,view,pln,yip[ntp],0);
//        yeip[ntp] = sciwidth(mod,view,pln,ch,0);
//
//        xip[ntp] += C_B2CHPosZ;
//        if(vertical) yip[ntp] += C_B2CHPosX;
//        else         yip[ntp] += C_B2CHPosY;
//        ntp++;
//      }
//
//      for(int hitnum=0; hitnum<(int)itrk[inum].hit.size(); hitnum++){
//        int ingmod  = itrk[inum].mod ;
//        int view    = itrk[inum].hit[hitnum].view;
//        int pln     = itrk[inum].hit[hitnum].pln ;
//        int ch      = itrk[inum].hit[hitnum].ch  ;
//        xip [ntp] = zposi   (ingmod,view,pln,ch,0);
//        yip [ntp] = xyposi  (ingmod,view,pln,ch,0);
//        xeip[ntp] = scithick(ingmod,view,pln,yip[ntp],0);
//        yeip[ntp] = sciwidth(ingmod,view,pln,ch,0);
//
//        xip[ntp] += C_B2INGPosZ;
//        if(vertical) yip[ntp] += C_B2INGPosX;
//        else         yip[ntp] += C_B2INGPosY;
//
//        if(vertical){
//          yip[ntp] += (ingmod-INGMODNUM_mid)*C_INGSpace;
//        }
//
//        ntp++;
//      }
//
//      SortArray(ntp,xip,yip,xeip,yeip); 
//      int cnt=0;
//      float last_x=-1e+5;
//      float max_y =-1e+5;
//      float min_y = 1e+5;
//      for(int itp=0;itp<ntp;itp++){
//        if(last_x==xip[itp]){
//          if(max_y<yip[itp]){max_y=yip[itp];}
//          if(min_y>yip[itp]){min_y=yip[itp];}
//          cnt++;
//          yip [itp-1] *= cnt;
//          yip [itp-1] += yip[itp];
//          yip [itp-1] /= cnt+1;
//          yeip[itp-1] = max_y-min_y+yeip[itp];
//          for(int jtp=itp;jtp<ntp-1;jtp++){
//            xip [jtp] = xip [jtp+1];
//            yip [jtp] = yip [jtp+1];
//            xeip[jtp] = xeip[jtp+1];
//            yeip[jtp] = yeip[jtp+1];
//          }
//          ntp--;
//          itp--;
//        }
//        else{
//          cnt=0;
//          last_x=-1e+5;
//          max_y =-1e+5;
//          min_y = 1e+5;
//          last_x = xip[itp];
//        }
//      }
//
//      TGraphErrors *gip = new TGraphErrors(ntp,xip,yip,xeip,yeip);
//      TF1 *fip = new TF1("fip","[0]+[1]*x");
//      fip->SetParameters( yip[0]-xip[0]*(yip[ntp-1]-yip[0])/(xip[ntp-1]-xip[0]),
//			 (yip[ntp-1]-yip[0])/(xip[ntp-1]-xip[0]));
//      gip->Fit("fip","Q");
//
//      TF1 *funcip    = gip->GetFunction("fip");
//      float intcptip = funcip->GetParameter(0);
//      float slopeip  = funcip->GetParameter(1);
//#ifdef DEBUG_THREEDIMRECON
//      cout << " Fit:" 
//        << " slope=" << slopeip
//        << " ang="   << atan(slopeip)*180./PI
//        << endl;
//#endif
//      gip->Delete();
//      fip->Delete();
//
//      pmtrk[i]   .ang    = atan(slopeip)*180./PI;
//      pmtrk[i]   .slope  = slopeip;
//      pmtrk[i]   .intcpt = pmtrk[i].ixy-pmtrk[i].slope*pmtrk[i].iz;
//      itrk [inum].ang    = atan(slopeip)*180./PI;
//      itrk [inum].slope  = slopeip;
//      itrk [inum].intcpt = itrk[inum].ixy-itrk[inum].slope*itrk[inum].iz;
//    }
//  }
//#endif

  return hasingtrk;
};

bool fIngHitWMJoint( vector<TrackIng> &itrk, vector<TrackWM> &wmtrk, vector<TrackPM> &pmtrk, bool vertical, int mod){


#ifdef DEBUG_THREEDIMRECON
  cout << " -------------- fIngHitWMJoint " << vertical << " --------------- " << endl;
#endif

  float diff_pos;
  //float joilik = -1e-5;
  int   joitra = -1;
  bool  jointed;
  bool  hasingtrk = false;
  int view;
  if(vertical) view = 1;
  else         view = 0;

  for(int incmod=INGMODNUM_start; incmod<=INGMODNUM_end; incmod++){
    for(int i=0; i<(int)wmtrk.size(); i++){
      if( wmtrk[i].ing_trk ) continue;
      if(!wmtrk[i].pm_match) continue;
      jointed = false;
      
      for(int pln=0; pln<2; pln++){
        for(int ch=0; ch<C_INGNumCh; ch++){
          
          if(nonrechits[incmod][view][pln][ch]<PEth_extrahits) continue;
          if(incmod!=INGMODNUM_mid) continue;

          double zpos  = zposi (incmod,view,pln,ch);
          double xypos = xyposi(incmod,view,pln,ch);

          if     (MODE_DET==0) zpos += C_INGHMotherPosZ - C_PMMotherPosZ;
          else if(MODE_DET==1) zpos += C_B2INGPosZ      - C_B2WMPosZ;
          if(view==1){
            xypos += C_INGSpace*(incmod-INGMODNUM_mid);
          }
          if(view==0&&MODE_DET==1) xypos += C_B2INGPosY - C_B2WMPosY;
          if(view==1&&MODE_DET==1) xypos += C_B2INGPosX - C_B2WMPosX;

          diff_pos = (wmtrk[i].intcpt+wmtrk[i].slope*zpos) - xypos;

#ifdef DEBUG_THREEDIMRECON
      cout 
        << " wmtrk" << i
        << " ihit:"
        << " pln=" << pln
        << " ch="  << ch
        << " diff_pos=" << diff_pos
        << endl;
#endif

          if(fabs(diff_pos)<hitpos_th2){
            if(!jointed){

              wmtrk[i].ing_imod = incmod;
              wmtrk[i].ing_fmod = incmod;
              wmtrk[i].ing_ipln = 0;
              wmtrk[i].ing_fpln = pln;
              wmtrk[i].ing_trk  = true;
              if(ch<4 || ch>20) wmtrk[i].ing_stop = false;
              else              wmtrk[i].ing_stop = true;
              wmtrk[i].ing_num    = (int)itrk.size();
              wmtrk[i].iron_pene  = pln;

              wmtrk[i].diff_pos   = -1000.;
              wmtrk[i].diff_time  = -1000.;
              wmtrk[i].diff_ang   = -1000.;

              TrackIng  ingtrack;
              ingtrack.clear();
              ingtrack.mod    = incmod;
              ingtrack.view   = view;
              ingtrack.ipln   = pln;
              ingtrack.fpln   = pln;
              ingtrack.ixy    = xyposi (incmod,view,pln,ch);
              ingtrack.fxy    = xyposi (incmod,view,pln,ch);
              ingtrack.iz     = zposi  (incmod,view,pln,ch);
              ingtrack.fz     = zposi  (incmod,view,pln,ch);
              ingtrack.slope  = wmtrk[i].slope;
              ingtrack.intcpt = ingtrack.ixy-ingtrack.slope*ingtrack.iz;
              ingtrack.ang    = wmtrk[i].ang;
              ingtrack.clstime= wmtrk[i].clstime;
              ingtrack.veto   = true;
              if(ch==0||ch==C_INGNumCh-1) ingtrack.edge = true;
              else                        ingtrack.edge = false;
              if(ch<4 || ch>20) ingtrack.stop = false;
              else              ingtrack.stop = true;

              Hits      hits;
              hits.clear();
              hits.cyc      = cyc;
              hits.mod      = incmod;
              hits.view     = view;
              hits.pln      = pln;
              hits.ch       = ch;
              hits.pe       = nonrechits     [incmod][view][pln][ch];
              hits.lope     = nonrechits_lope[incmod][view][pln][ch];
              hits.pdg      = nonrechits_pdg [incmod][view][pln][ch];
              hits.hit_id   = nonrechits_id  [incmod][view][pln][ch];
              hits.isohit   = false;
              ingtrack.hit.push_back(hits);
              itrk.push_back(ingtrack);

              int psize = pmtrk.size();
              for(int ip=0;ip<psize;ip++){
                if(pmtrk[ip].wm_num==i){
                  pmtrk[ip].ing_trk   = wmtrk[i].ing_trk;
                  pmtrk[ip].ing_num   = wmtrk[i].ing_num;
                  pmtrk[ip].ing_imod  = wmtrk[i].ing_imod;
                  pmtrk[ip].ing_fmod  = wmtrk[i].ing_fmod;
                  pmtrk[ip].ing_ipln  = wmtrk[i].ing_ipln;
                  pmtrk[ip].ing_fpln  = wmtrk[i].ing_fpln;
                  pmtrk[ip].ing_stop  = wmtrk[i].ing_stop;
                  pmtrk[ip].iron_pene = wmtrk[i].iron_pene;
                }
              }
            }
            else{
              wmtrk[i].ing_fpln  = pln;
              wmtrk[i].iron_pene = pln;

              if(ch<4 || ch>20) wmtrk[i].ing_stop = false;
              else              wmtrk[i].ing_stop = true;

              if(itrk[wmtrk[i].ing_num].fpln<pln){
                itrk[wmtrk[i].ing_num].fpln =pln;
                itrk[wmtrk[i].ing_num].fxy  =xyposi(incmod,view,pln,ch);
              }
              if(itrk[wmtrk[i].ing_num].ipln>pln){
                itrk[wmtrk[i].ing_num].ipln=pln;
                itrk[wmtrk[i].ing_num].ixy =xyposi(incmod,view,pln,ch);
              }
              wmtrk[i].diff_pos   = -1000.;
              wmtrk[i].diff_time  = -1000.;
              wmtrk[i].diff_ang   = -1000.;

              Hits hits;
              hits.clear();
              hits.cyc    = cyc;
              hits.mod    = incmod;
              hits.view   = view;
              hits.pln    = pln;
              hits.ch     = ch;
              hits.pe     = nonrechits     [incmod][view][pln][ch];
              hits.lope   = nonrechits_lope[incmod][view][pln][ch];
              hits.pdg    = nonrechits_pdg [incmod][view][pln][ch];
              hits.hit_id = nonrechits_id  [incmod][view][pln][ch];
              hits.isohit = false;
              itrk[wmtrk[i].ing_num].hit.push_back(hits);
            }

            jointed   = true;
            joitra    = i;
            hasingtrk = true;
            nonrechits     [incmod][view][pln][ch] = 0;
            nonrechits_lope[incmod][view][pln][ch] = 0;
            nonrechits_pdg [incmod][view][pln][ch] = 0;
            nonrechits_id  [incmod][view][pln][ch] = 0;

#ifdef DEBUG_THREEDIMRECON
            cout << "  Jointed " << endl;
#endif

          }
        }
      }
    }
  }

//#ifdef USEINGTRK  
//  for(int i=0; i<(int)wmtrk.size(); i++){
//    if(wmtrk[i].ing_trk){
//      int   inum = wmtrk[i].ing_num;
//      float xip[Cpln], yip[Cpln], xeip[Cpln], yeip[Cpln];
//      int   ntp = 0;
//
//      for(int hitnum=0; hitnum<(int)wmtrk[i].hit.size(); hitnum++){
//        int view = wmtrk[i].hit[hitnum].view;
//        int pln  = wmtrk[i].hit[hitnum].pln ;
//        int ch   = wmtrk[i].hit[hitnum].ch  ;
//        xip[ntp] = zposi (mod,view,pln,ch,0);
//        yip[ntp] = xyposi(mod,view,pln,ch,0);
//
//        if(MODE_DET==0){
//          xip[ntp] += C_PMMotherPosZ;
//          if(vertical) yip[ntp] += C_PMMotherPosX;
//          else         yip[ntp] += C_PMMotherPosY;
//        }
//        else if(MODE_DET==1){
//          xip[ntp] += C_B2WMPosZ;
//          if(vertical) yip[ntp] += C_B2WMPosX;
//          else         yip[ntp] += C_B2WMPosY;
//        }
//        xeip[ntp]=scithick(mod,view,pln,yip[ntp],0);
//        yeip[ntp]=sciwidth(mod,view,pln,ch,0);
//        ntp++;
//      }
//
//      for(int hitnum=0; hitnum<(int)itrk[inum].hit.size(); hitnum++){
//        int ingmod = itrk[inum].mod ;
//        int view   = itrk[inum].hit[hitnum].view;
//        int pln    = itrk[inum].hit[hitnum].pln ;
//        int ch     = itrk[inum].hit[hitnum].ch  ;
//        xip[ntp]=zposi (ingmod,view,pln,ch,0);
//        yip[ntp]=xyposi(ingmod,view,pln,ch,0);
//
//        if(MODE_DET==0){
//          xip[ntp] += C_INGHMotherPosZ;
//          if(vertical) yip[ntp] += C_INGHMotherPosX;
//          else         yip[ntp] += C_INGHMotherPosY;
//        }
//        else if(MODE_DET==1){
//          xip[ntp] += C_B2INGPosZ;
//          if(vertical) yip[ntp] += C_B2INGPosX;
//          else         yip[ntp] += C_B2INGPosY;
//        }
//
//        if(view==1){
//          yip[ntp]+=(ingmod-INGMODNUM_mid)*C_INGSpace;
//        }
//        xeip[ntp]=scithick(ingmod,view,pln,yip[ntp],0);
//        yeip[ntp]=sciwidth(ingmod,view,pln,ch,0);
//
//        ntp++;
//      }
//
//      SortArray(ntp,xip,yip,xeip,yeip);
//      int cnt=0;
//      float last_x=-1e+5;
//      float max_y =-1e+5;
//      float min_y = 1e+5;
//      for(int itp=0;itp<ntp;itp++){
//        if(last_x==xip[itp]){
//          if(max_y<yip[itp]){max_y=yip[itp];}
//          if(min_y>yip[itp]){min_y=yip[itp];}
//          cnt++;
//          yip [itp-1] *= cnt;
//          yip [itp-1] += yip[itp];
//          yip [itp-1] /= cnt+1;
//          yeip[itp-1] = max_y-min_y+yeip[itp];
//          for(int jtp=itp;jtp<ntp-1;jtp++){
//            xip [jtp] = xip [jtp+1];
//            yip [jtp] = yip [jtp+1];
//            xeip[jtp] = xeip[jtp+1];
//            yeip[jtp] = yeip[jtp+1];
//          }
//          ntp--;
//          itp--;
//        }
//        else{
//          cnt=0;
//          last_x=-1e+5;
//          max_y =-1e+5;
//          min_y = 1e+5;
//          last_x = xip[itp];
//        }
//      }
//
//      TGraphErrors *gip = new TGraphErrors(ntp,xip,yip,xeip,yeip);
//      TF1 *fip = new TF1("fip","[0]+[1]*x");
//      fip->SetParameters( yip[0]-xip[0]*(yip[ntp-1]-yip[0])/(xip[ntp-1]-xip[0]),
//			 (yip[ntp-1]-yip[0])/(xip[ntp-1]-xip[0]));
//      gip->Fit("fip","Q");
//
//      TF1 *funcip    = gip->GetFunction("fip");
//      float intcptip = funcip->GetParameter(0);
//      float slopeip  = funcip->GetParameter(1);
//
//      gip->Delete();
//      fip->Delete();
//
//
//#ifdef DEBUG_THREEDIMRECON
//      cout << " Fit:" 
//        << " slope=" << slopeip
//        << " ang="   << atan(slopeip)*180./PI
//        << endl;
//#endif
//
//      wmtrk[i].ang    = atan(slopeip)*180./PI;
//      wmtrk[i].slope  = slopeip;
//      wmtrk[i].intcpt = wmtrk[i].ixy-wmtrk[i].slope*wmtrk[i].iz;
//      itrk[inum].ang    = atan(slopeip)*180./PI;
//      itrk[inum].slope  = slopeip;
//      itrk[inum].intcpt = itrk[inum].ixy-itrk[inum].slope*itrk[inum].iz;
//    }
//  }
//#endif

  return hasingtrk;
};

//bool fIngHitPMJoint( vector<TrackIng> &itrk, vector<TrackPM> &pmtrk, bool vertical, int mod){
//
//  if(MODE_DET!=1){return false;}
//
//#ifdef DEBUG_THREEDIMRECON
//  cout << " -------------- fIngHitPMJoint " << vertical << " --------------- " << endl;
//#endif
//
//  float diff_pos;
//  //float joilik = -1e-5;
//  int   joitra = -1;
//  bool  jointed;
//  bool  hasingtrk = false;
//  int view;
//  if(vertical) view = 1;
//  else         view = 0;
//
//  for(int incmod=INGMODNUM_start; incmod<=INGMODNUM_end; incmod++){
//    for(int i=0; i<(int)pmtrk.size(); i++){
//      if(pmtrk[i].ing_trk) continue;
//      if(pmtrk[i].wm_trk ) continue;
//      if(pmtrk[i].fpln<C_PMNumPln-2) continue;
//      jointed = false;
//      
//      for(int pln=0; pln<2; pln++){
//        for(int ch=0; ch<C_INGNumCh; ch++){
//          
//          if(nonrechits[incmod][view][pln][ch]<PEth_extrahits) continue;
//          if(incmod!=INGMODNUM_mid) continue;
//
//          double zpos  = zposi (incmod,view,pln,ch);
//          double xypos = xyposi(incmod,view,pln,ch);
//
//          zpos += C_B2INGPosZ - C_B2CHPosZ;
//          if(view==1){
//            xypos += C_INGSpace*(incmod-INGMODNUM_mid);
//          }
//          if(view==1) xypos += C_B2INGPosX - C_B2CHPosX;
//
//          diff_pos = (pmtrk[i].intcpt+pmtrk[i].slope*zpos) - xypos;
//
//#ifdef DEBUG_THREEDIMRECON
//      cout 
//        << " pmtrk" << i
//        << " ihit:"
//        << " pln=" << pln
//        << " ch="  << ch
//        << " diff_pos=" << diff_pos
//        << endl;
//#endif
//
//          if(fabs(diff_pos)<hitpos_th2){
//            if(!jointed){
//
//              pmtrk[i].ing_imod = incmod;
//              pmtrk[i].ing_fmod = incmod;
//              pmtrk[i].ing_ipln = pln;
//              pmtrk[i].ing_fpln = pln;
//              pmtrk[i].ing_trk  = true;
//              if(ch<4 || ch>20) pmtrk[i].ing_stop = false;
//              else              pmtrk[i].ing_stop = true;
//              int itrksize = (int)itrk.size();
//              pmtrk[i].ing_num   = itrksize;
//              pmtrk[i].iron_pene = pln;
//
//              pmtrk[i].diff_pos  = -1000.;
//              pmtrk[i].diff_time = -1000.;
//              pmtrk[i].diff_ang  = -1000.;
//
//              TrackIng  ingtrack;
//              ingtrack.clear();
//              ingtrack.mod    = incmod;
//              ingtrack.view   = view;
//              ingtrack.ipln   = pln;
//              ingtrack.fpln   = pln;
//              ingtrack.ixy    = xyposi (incmod,view,pln,ch);
//              ingtrack.fxy    = xyposi (incmod,view,pln,ch);
//              ingtrack.iz     = zposi  (incmod,view,pln,ch);
//              ingtrack.fz     = zposi  (incmod,view,pln,ch);
//              ingtrack.slope  = pmtrk[i].slope;
//              ingtrack.intcpt = ingtrack.ixy-ingtrack.slope*ingtrack.iz;
//              ingtrack.ang    = pmtrk[i].ang;
//              ingtrack.clstime= pmtrk[i].clstime;
//              ingtrack.veto   = true;
//              if(ch==0||ch==C_INGNumCh-1) ingtrack.edge = true;
//              else                        ingtrack.edge = false;
//              if(ch<4 || ch>20) ingtrack.stop = false;
//              else              ingtrack.stop = true;
//
//              Hits      hits;
//              hits.clear();
//              hits.cyc      = cyc;
//              hits.mod      = incmod;
//              hits.view     = view;
//              hits.pln      = pln;
//              hits.ch       = ch;
//              hits.pe       = nonrechits     [incmod][view][pln][ch];
//              hits.lope     = nonrechits_lope[incmod][view][pln][ch];
//              hits.pdg      = nonrechits_pdg [incmod][view][pln][ch];
//              hits.hit_id   = nonrechits_id  [incmod][view][pln][ch];
//              hits.isohit   = false;
//              ingtrack.hit.push_back(hits);
//
//              itrk.push_back(ingtrack);
//            }
//            else{
//              if(pmtrk[i].ing_ipln>pln) pmtrk[i].ing_ipln = pln;
//              if(pmtrk[i].ing_fpln<pln) pmtrk[i].ing_fpln = pln;
//              pmtrk[i].iron_pene = pln;
//
//              if(ch<4 || ch>20) pmtrk[i].ing_stop = false;
//              else              pmtrk[i].ing_stop = true;
//
//              if(itrk[pmtrk[i].ing_num].fpln<pln){
//                itrk[pmtrk[i].ing_num].fpln=pln;
//                itrk[pmtrk[i].ing_num].fxy =xyposi(incmod,view,pln,ch);
//              }else if(itrk[pmtrk[i].ing_num].fpln==pln){
//                itrk[pmtrk[i].ing_num].fxy += xyposi(incmod,view,pln,ch);
//                itrk[pmtrk[i].ing_num].fxy /= 2.;
//              }
//              if(itrk[pmtrk[i].ing_num].ipln>pln){
//                itrk[pmtrk[i].ing_num].ipln=pln;
//                itrk[pmtrk[i].ing_num].ixy =xyposi(incmod,view,pln,ch);
//              }else if(itrk[pmtrk[i].ing_num].ipln==pln){
//                itrk[pmtrk[i].ing_num].ixy += xyposi(incmod,view,pln,ch);
//                itrk[pmtrk[i].ing_num].ixy /= 2.;
//              }
//              pmtrk[i].diff_pos  = -1000.;
//              pmtrk[i].diff_time = -1000.;
//              pmtrk[i].diff_ang  = -1000.;
//
//              Hits hits;
//              hits.clear();
//              hits.cyc    = cyc;
//              hits.mod    = incmod;
//              hits.view   = view;
//              hits.pln    = pln;
//              hits.ch     = ch;
//              hits.pe     = nonrechits     [incmod][view][pln][ch];
//              hits.lope   = nonrechits_lope[incmod][view][pln][ch];
//              hits.pdg    = nonrechits_pdg [incmod][view][pln][ch];
//              hits.hit_id = nonrechits_id  [incmod][view][pln][ch];
//              hits.isohit = false;
//              itrk[pmtrk[i].ing_num].hit.push_back(hits);
//            }
//
//            jointed   = true;
//            joitra    = i;
//            hasingtrk = true;
//            nonrechits     [incmod][view][pln][ch] = 0;
//            nonrechits_lope[incmod][view][pln][ch] = 0;
//            nonrechits_pdg [incmod][view][pln][ch] = 0;
//            nonrechits_id  [incmod][view][pln][ch] = 0;
//
//#ifdef DEBUG_THREEDIMRECON
//      cout << "  Jointed " << endl;
//#endif
//
//
//          }
//        }
//      }
//    }
//  }
//
//#ifdef USEINGTRK  
//  for(int i=0; i<(int)pmtrk.size(); i++){
//    if(pmtrk[i].ing_trk&&!pmtrk[i].wm_trk){
//      int   inum = pmtrk[i].ing_num;
//      float xip[Cpln], yip[Cpln], xeip[Cpln], yeip[Cpln];
//      int   ntp = 0;
//
//      for(int hitnum=0; hitnum<(int)pmtrk[i].hit.size(); hitnum++){
//        int view = pmtrk[i].hit[hitnum].view;
//        int pln  = pmtrk[i].hit[hitnum].pln ;
//        int ch   = pmtrk[i].hit[hitnum].ch  ;
//        xip[ntp] = zposi (mod,view,pln,ch,0);
//        yip[ntp] = xyposi(mod,view,pln,ch,0);
//
//        xip[ntp] += C_B2CHPosZ;
//        if(vertical) yip[ntp] += C_B2CHPosX;
//        else         yip[ntp] += C_B2CHPosY;
//        xeip[ntp]=scithick(mod,view,pln,yip[ntp],0);
//        yeip[ntp]=sciwidth(mod,view,pln,ch,0);
//        ntp++;
//      }
//      for(int hitnum=0; hitnum<(int)itrk[inum].hit.size(); hitnum++){
//        int ingmod = itrk[inum].mod ;
//        int view   = itrk[inum].hit[hitnum].view;
//        int pln    = itrk[inum].hit[hitnum].pln ;
//        int ch     = itrk[inum].hit[hitnum].ch  ;
//        xip[ntp]=zposi (ingmod,view,pln,ch,0);
//        yip[ntp]=xyposi(ingmod,view,pln,ch,0);
//
//        xip[ntp] += C_B2INGPosZ;
//        if(vertical) yip[ntp] += C_B2INGPosX;
//        else         yip[ntp] += C_B2INGPosY;
//        xeip[ntp]=scithick(ingmod,view,pln,yip[ntp],0);
//        yeip[ntp]=sciwidth(ingmod,view,pln,ch,0);
//        ntp++;
//      }
//
//      SortArray(ntp,xip,yip,xeip,yeip);
//      int cnt=0;
//      float last_x=-1e+5;
//      float max_y =-1e+5;
//      float min_y = 1e+5;
//      for(int itp=0;itp<ntp;itp++){
//        if(last_x==xip[itp]){
//          if(max_y<yip[itp]){max_y=yip[itp];}
//          if(min_y>yip[itp]){min_y=yip[itp];}
//          cnt++;
//          yip [itp-1] *= cnt;
//          yip [itp-1] += yip[itp];
//          yip [itp-1] /= cnt+1;
//          yeip[itp-1] = max_y-min_y+yeip[itp];
//          for(int jtp=itp;jtp<ntp-1;jtp++){
//            xip [jtp] = xip [jtp+1];
//            yip [jtp] = yip [jtp+1];
//            xeip[jtp] = xeip[jtp+1];
//            yeip[jtp] = yeip[jtp+1];
//          }
//          ntp--;
//          itp--;
//        }
//        else{
//          cnt=0;
//          last_x=-1e+5;
//          max_y =-1e+5;
//          min_y = 1e+5;
//          last_x = xip[itp];
//        }
//      }
//
//      TGraphErrors *gip = new TGraphErrors(ntp,xip,yip,xeip,yeip);
//      TF1 *fip = new TF1("fip","[0]+[1]*x");
//      fip->SetParameters( yip[0]-xip[0]*(yip[ntp-1]-yip[0])/(xip[ntp-1]-xip[0]),
//			 (yip[ntp-1]-yip[0])/(xip[ntp-1]-xip[0]));
//      gip->Fit("fip","Q");
//
//      TF1 *funcip    = gip->GetFunction("fip");
//      float intcptip = funcip->GetParameter(0);
//      float slopeip  = funcip->GetParameter(1);
//
//      gip->Delete();
//      fip->Delete();
//
//
//#ifdef DEBUG_THREEDIMRECON
//      cout << " Fit:" 
//        << " slope=" << slopeip
//        << " ang="   << atan(slopeip)*180./PI
//        << endl;
//#endif
//
//      pmtrk[i].ang    = atan(slopeip)*180./PI;
//      pmtrk[i].slope  = slopeip;
//      pmtrk[i].intcpt = pmtrk[i].ixy-pmtrk[i].slope*pmtrk[i].iz;
//
//      itrk[inum].ang    = atan(slopeip)*180./PI;
//      itrk[inum].slope  = slopeip;
//      itrk[inum].intcpt = itrk[inum].ixy-itrk[inum].slope*itrk[inum].iz;
//
//    }
//  }
//#endif
//
//  return hasingtrk;
//};

//bool fWMHitIngJoint( vector<TrackIng> &itrk, vector<TrackWM> &wmtrk, bool vertical, int mod){
//
//#ifdef DEBUG_THREEDIMRECON
//  cout << " -------------- fWMHitWMJoint " << vertical << " --------------- " << endl;
//#endif
//
//  float diff_pos;
//  bool  jointed;
//  bool  hasingtrk = false;
//  int view;
//  if(vertical) view = 1;
//  else         view = 0;
//
//
//  int used_nonreco[2][2][C_PMNumCh];
//  for(int ii=0;ii<2;ii++){
//    for(int jj=0;jj<2;jj++){
//      for(int kk=0;kk<C_WMNumCh;kk++){
//        used_nonreco[ii][jj][kk] = 0;
//      }
//    }
//  }
//
//  for(int i=0; i<(int)itrk.size(); i++){
//#ifdef DEBUG_THREEDIMRECON
//      cout 
//        << " imtrk" << i
//        << " ipln=" << itrk[i].ipln
//        << endl;
//#endif
//
//    if( itrk[i].ipln>2  ){continue;}
//    bool wmjointed = false;
//    for(int ii=0;ii<(int)wmtrk.size();ii++){
//      if(wmtrk[ii].ing_num==i){wmjointed=true;break;}
//    }
//    if(wmjointed){continue;}
//
//    jointed = false;
//    for(int pln=C_WMNumPln*3-1; pln>=C_WMNumPln*3-7; pln--){
//      for(int ch=0; ch<40; ch++){
//
//        int rawpln,rawch,grid;
//        if(!detdim->GetRawPlnChGrid(mod,view,pln,ch,0,&rawpln,&rawch,&grid)){continue;}
//        if(grid==1) rawch+=40;
//        else if(grid==2) rawch+=60;
//
//
//        if(nonrechits[mod][view][pln][ch]<PEth_extrahits){ continue;}
//        if(used_nonreco[view][pln-C_WMNumPln+2][ch]>2){continue;}
//        double zpos  = zposi (mod,view,pln,ch);
//        double xypos = xyposi(mod,view,pln,ch);
//        zpos += C_B2WMPosZ -C_B2INGPosZ;
//        if(view==1){ xypos += C_B2WMPosX-C_B2INGPosX; }
//        else       { xypos += C_B2WMPosY-C_B2INGPosY; }
//
//        diff_pos = (itrk[i].intcpt+itrk[i].slope*zpos) - xypos;
//
//#ifdef DEBUG_THREEDIMRECON
//      cout 
//        << " wmhit:"
//        << " pln=" << pln
//        << " ch="  << ch
//        << " diff_pos=" << diff_pos
//        << " pe="  << nonrechits[mod][view][pln][ch]
//        << endl;
//#endif
//
//
//        if(fabs(diff_pos)<pos_th_hit3){
//          used_nonreco[view][pln-C_WMNumPln+2][ch]++;
//          if(jointed){
//
//            Hits hits;
//            hits.clear();
//            hits.cyc    = cyc;
//            hits.mod    = mod;
//            hits.view   = view;
//            hits.pln    = rawpln;
//            hits.ch     = rawch;
//            hits.pe     = nonrechits     [mod][view][pln][ch];
//            hits.lope   = nonrechits_lope[mod][view][pln][ch];
//            hits.pdg    = nonrechits_pdg [mod][view][pln][ch];
//            hits.hit_id = nonrechits_id  [mod][view][pln][ch];
//            hits.isohit = false;
//            int lptrk = (int)wmtrk.size()-1;
//            wmtrk[lptrk].hit.push_back(hits);
//            if(wmtrk[lptrk].fpln<pln){
//              wmtrk[lptrk].fpln = pln;
//              wmtrk[lptrk].fxy  = xyposi(mod,view,pln,ch);
//              wmtrk[lptrk].fz   = zposi (mod,view,pln,ch);
//              if(grid==0&&(ch==39||ch== 0)) wmtrk[lptrk].edge = true;
//              if(grid==1&&(ch==59||ch==40)) wmtrk[lptrk].edge = true;
//              if(grid==2&&(ch==79||ch==60)) wmtrk[lptrk].edge = true;
//            }else if(wmtrk[lptrk].fpln==pln){
//              wmtrk[lptrk].fxy += xyposi(mod,view,pln,ch);
//              wmtrk[lptrk].fxy /= 2.;
//            }
//            if(wmtrk[lptrk].ipln>pln){
//              wmtrk[lptrk].ipln = pln;
//              wmtrk[lptrk].ixy  = xyposi(mod,view,pln,ch);
//              wmtrk[lptrk].iz   = zposi (mod,view,pln,ch);
//            }else if(wmtrk[lptrk].ipln==pln){
//              wmtrk[lptrk].ixy  += xyposi(mod,view,pln,ch);
//              wmtrk[lptrk].ixy  /= 2;
//            }
//
//          }
//          else{
//            jointed   = true;
//
//            TrackWM wtrack;
//            wtrack.clear();
//            wtrack.view      = view;
//            wtrack.ipln      = pln;
//            wtrack.fpln      = pln;
//            wtrack.ixy       = xyposi(mod,view,pln,ch);
//            wtrack.fxy       = xyposi(mod,view,pln,ch);
//            wtrack.iz        = zposi (mod,view,pln,ch);
//            wtrack.fz        = zposi (mod,view,pln,ch);
//            wtrack.slope     = itrk[i].slope;
//            wtrack.intcpt    = wtrack.ixy - wtrack.slope*wtrack.iz;
//            wtrack.ang       = itrk[i].ang;
//            wtrack.clstime   = itrk[i].clstime;
//            wtrack.veto      = false;
//            if     (grid==0&&(ch==39||ch== 0)) wtrack.edge=true;
//            else if(grid==1&&(ch==59||ch==40)) wtrack.edge=true;
//            else if(grid==2&&(ch==79||ch==60)) wtrack.edge=true;
//            else                               wtrack.edge=false;
//            wtrack.stop      = false;
//            wtrack.ing_imod  = INGMODNUM_mid;
//            wtrack.ing_fmod  = INGMODNUM_mid;
//            wtrack.ing_ipln  = itrk[i].ipln;
//            wtrack.ing_fpln  = itrk[i].fpln;
//            wtrack.ing_trk   = true;
//            wtrack.ing_num   = i;
//            wtrack.pm_stop   = false;
//            wtrack.ing_stop  = itrk[i].stop;
//
//            Hits hits;
//            hits.clear();
//            hits.cyc      = cyc;
//            hits.mod      = mod;
//            hits.view     = view;
//            hits.pln      = rawpln;
//            hits.ch       = rawch;
//            hits.pe       = nonrechits     [mod][view][pln][ch];
//            hits.lope     = nonrechits_lope[mod][view][pln][ch];
//            hits.pdg      = nonrechits_pdg [mod][view][pln][ch];
//            hits.hit_id   = nonrechits_id  [mod][view][pln][ch];
//            hits.isohit   = false;
//            wtrack.hit.push_back(hits);
//
//            wtrack.diff_pos  = -1000.;
//            wtrack.diff_time = -1000.;
//            wtrack.diff_ang  = -1000.;
//
//            wmtrk.push_back(wtrack);
//          }
//          hasingtrk = true;
//
//#ifdef DEBUG_THREEDIMRECON
//          cout << " Jointed."  << endl;
//#endif
//        }
//
//      }
//    }
//  }
//
//  return hasingtrk;
//};

//bool fPMHitWMJoint( vector<TrackIng> &itrk, vector<TrackWM> &wmtrk, vector<TrackPM> &pmtrk, bool vertical, int mod){
//
//  if(MODE_DET==0){return false;}
//
//#ifdef DEBUG_THREEDIMRECON
//  cout << " -------------- fPMHitWMJoint " << vertical << " --------------- " << endl;
//#endif
//
//  float diff_pos;
//  //float joilik = -1e-5;
//  //int   joitra = -1;
//  bool  jointed;
//  bool  hasingtrk = false;
//  int view;
//  if(vertical) view = 1;
//  else         view = 0;
//
//
//  int used_nonreco[2][2][C_PMNumCh];
//  for(int ii=0;ii<2;ii++){
//    for(int jj=0;jj<2;jj++){
//      for(int kk=0;kk<C_PMNumCh;kk++){
//        used_nonreco[ii][jj][kk] = 0;
//      }
//    }
//  }
//
//  for(int i=0; i<(int)wmtrk.size(); i++){
//#ifdef DEBUG_THREEDIMRECON
//      cout 
//        << " wmtrk" << i
//        << " ipln=" << wmtrk[i].ipln
//        << endl;
//#endif
//
//    //if( wmtrk[i].ipln>7  ){continue;}
//    if( wmtrk[i].pm_match){continue;}
//    //if(!wmtrk[i].ing_trk ){continue;}
//
//    jointed = false;
//    for(int pln=C_PMNumPln-1; pln>=16; pln--){
//      for(int ch=0; ch<C_PMNumCh; ch++){
//
//        if(nonrechits[mod][view][pln][ch]<PEth_extrahits){ continue;}
//        if(used_nonreco[view][pln-16][ch]>2){continue;}
//        double zpos  = zposi (mod,view,pln,ch);
//        double xypos = xyposi(mod,view,pln,ch);
//        zpos += C_B2CHPosZ -C_B2WMPosZ;
//        if(view==1){ xypos += C_B2CHPosX-C_B2WMPosX; }
//        else       { xypos += C_B2CHPosY-C_B2WMPosY; }
//
//        diff_pos = (wmtrk[i].intcpt+wmtrk[i].slope*zpos) - xypos;
//
//#ifdef DEBUG_THREEDIMRECON
//      cout 
//        << " pmhit:"
//        << " pln=" << pln
//        << " ch="  << ch
//        << " diff_pos=" << diff_pos
//        << endl;
//#endif
//
//
//        if(fabs(diff_pos)<pos_th_hit3){
//          used_nonreco[view][pln-16][ch]++;
//          if(jointed){
//
//            Hits hits;
//            hits.clear();
//            hits.cyc    = cyc;
//            hits.mod    = mod;
//            hits.view   = view;
//            hits.pln    = pln;
//            hits.ch     = ch;
//            hits.pe     = nonrechits     [mod][view][pln][ch];
//            hits.lope   = nonrechits_lope[mod][view][pln][ch];
//            hits.pdg    = nonrechits_pdg [mod][view][pln][ch];
//            hits.hit_id = nonrechits_id  [mod][view][pln][ch];
//            hits.isohit = false;
//            int lptrk = (int)pmtrk.size()-1;
//            pmtrk[lptrk].hit.push_back(hits);
//
//            int nphits = pmtrk[lptrk].hit.size();
//            if(pmtrk[lptrk].fpln<pln){
//              pmtrk[lptrk].fpln = pln;
//              pmtrk[lptrk].fxy  = xyposi(mod,view,pln,ch);
//              pmtrk[lptrk].fz   = zposi (mod,view,pln,ch);
//              if(ch==C_PMNumCh-1||ch==0) pmtrk[lptrk].edge = true;
//            }else if(pmtrk[lptrk].fpln==pln){
//              pmtrk[lptrk].fxy  = (pmtrk[lptrk].fxy*(nphits-1)+xyposi(mod,view,pln,ch))/nphits;
//            }
//            if(pmtrk[lptrk].ipln>pln){
//              pmtrk[lptrk].ipln = pln;
//              pmtrk[lptrk].ixy  = xyposi(mod,view,pln,ch);
//              pmtrk[lptrk].iz   = zposi (mod,view,pln,ch);
//            }else if(pmtrk[lptrk].ipln==pln){
//              pmtrk[lptrk].ixy  = (pmtrk[lptrk].ixy*(nphits-1)+xyposi(mod,view,pln,ch))/nphits;
//            }
//
//          }
//          else{
//            jointed   = true;
//
//            TrackPM ptrack;
//            ptrack.clear();
//            ptrack.view      = view;
//            ptrack.ipln      = pln;
//            ptrack.fpln      = pln;
//            ptrack.ixy       = xyposi(mod,view,pln,ch);
//            ptrack.fxy       = xyposi(mod,view,pln,ch);
//            ptrack.iz        = zposi (mod,view,pln,ch);
//            ptrack.fz        = zposi (mod,view,pln,ch);
//            ptrack.slope     = wmtrk[i].slope;
//            ptrack.intcpt    = ptrack.ixy - ptrack.slope*ptrack.iz;
//            ptrack.ang       = wmtrk[i].ang;
//            if(wmtrk[i].ing_trk) ptrack.clstime = itrk[wmtrk[i].ing_num].clstime;
//            else                 ptrack.clstime = wmtrk[i].clstime;
//            ptrack.veto      = false;
//            if(ch==C_PMNumCh-1||ch==0) ptrack.edge = true;
//            else                       ptrack.edge = false;
//            ptrack.stop      = false;
//            ptrack.ing_imod  = wmtrk[i].ing_imod;
//            ptrack.ing_fmod  = wmtrk[i].ing_fmod;
//            ptrack.ing_ipln  = wmtrk[i].ing_ipln;
//            ptrack.ing_fpln  = wmtrk[i].ing_fpln;
//            ptrack.wm_ipln   = wmtrk[i].ipln;
//            ptrack.wm_fpln   = wmtrk[i].fpln;
//            ptrack.ing_trk   = wmtrk[i].ing_trk;
//            ptrack.ing_num   = wmtrk[i].ing_num;
//            ptrack.wm_trk    = true;
//            ptrack.wm_num    = i;
//            ptrack.pm_stop   = false;
//            ptrack.wm_stop   = wmtrk[i].stop;
//            ptrack.ing_stop  = wmtrk[i].ing_stop;
//            ptrack.iron_pene = wmtrk[i].iron_pene;
//
//            Hits hits;
//            hits.clear();
//            hits.cyc      = cyc;
//            hits.mod      = mod;
//            hits.view     = view;
//            hits.pln      = pln;
//            hits.ch       = ch;
//            hits.pe       = nonrechits     [mod][view][pln][ch];
//            hits.lope     = nonrechits_lope[mod][view][pln][ch];
//            hits.pdg      = nonrechits_pdg [mod][view][pln][ch];
//            hits.hit_id   = nonrechits_id  [mod][view][pln][ch];
//            hits.isohit   = false;
//            ptrack.hit.push_back(hits);
//
//            ptrack.diff_pos  = -1000.;
//            ptrack.diff_time = -1000.;
//            ptrack.diff_ang  = -1000.;
//
//            pmtrk.push_back(ptrack);
//          }
//          wmtrk[i].pm_match = true;
//          hasingtrk = true;
//
//#ifdef DEBUG_THREEDIMRECON
//          cout << " Jointed."  << endl;
//#endif
//        }
//
//      }
//    }
//  }
//
//  return hasingtrk;
//};

bool fPMWMJoint( vector<TrackIng> &itrk, vector<TrackWM> &wmtrk, vector<TrackPM> &pmtrk, bool vertical, int mod){

  if(MODE_DET==0){return false;}

#ifdef DEBUG_THREEDIMRECON
  cout << " -------------- fPMWMJoint " << vertical << "  --------------- " << endl;
#endif

  float diff_ang , diff_pos;
  float ang_th_pm, pos_th_pm;
  float joilik = -1e-5;
  int   joitra = -1;
  bool  jointed   = false;
  bool  hasingtrk = false;

  float ang_th = ang_th1;
  float pos_th = pos_th1;

  for(int j=0; j<(int)wmtrk.size(); j++){
    //if(wmtrk[j].ipln>7) continue;

    jointed=false;
    for(int i=0; i<(int)pmtrk.size(); i++){

      double zpos_at_joint = (  C_B2WMPosZ + zposi(WMMODNUM,0,0,0)
			      + C_B2CHPosZ + zposi(PMMODNUM,1,C_PMNumPln-1,0) ) /2.; 
      diff_pos = ( (pmtrk[i].intcpt + pmtrk[i].slope*(zpos_at_joint-C_B2CHPosZ))
		  -(wmtrk[j].intcpt + wmtrk[j].slope*(zpos_at_joint-C_B2WMPosZ)) );

      diff_ang = pmtrk[i].ang-wmtrk[j].ang;

      ang_th_pm = ang_th;
      pos_th_pm = pos_th;

#ifdef DEBUG_THREEDIMRECON
      cout 
        << " wmtrk" << j
        << " pmtrk" << i
        << " diff_ang=" << diff_ang
        << " diff_pos=" << diff_pos
        << " zpos_at_joint=" << zpos_at_joint
        << endl;
#endif

      if(fabs(diff_ang)<ang_th_pm && fabs(diff_pos)<pos_th_pm){
        if(jointed){
          if(joilik>sqrt( fabs(diff_ang)*fabs(diff_ang)/ang_th_pm/ang_th_pm
			 +fabs(diff_pos)*fabs(diff_pos)/pos_th_pm/pos_th_pm))
          {           
#ifdef DEBUG_THREEDIMRECON
            cout << " pmtrk:" << joitra << " is reapleced" << endl;
#endif
            pmtrk[joitra].wm_trk = false;
          }
          else continue;
        }

	//If pmtrk[i] already has jointed wmtrk
        if(pmtrk[i].wm_trk)
        {
          bool replace = false;
          if(pmtrk[i].ing_trk){
            if(wmtrk[j].ing_trk){
              if(pmtrk[i].ing_fpln<wmtrk[j].ing_fpln){ replace=true; }
              else{ replace=false; }
            }
            else{ replace=false; }
          }else{
            if(wmtrk[j].ing_trk){replace=true;}
            else{
              if(pmtrk[i].wm_fpln<wmtrk[j].fpln){ replace=true; }
              else{ replace=false; }
            }
          }
#ifdef DEBUG_THREEDIMRECON
          if(replace){ 
            cout << "wmtrk:" << pmtrk[i].wm_num << " is replaced." << endl; 
          }
#endif
          if(!replace){continue;}

          int l = pmtrk[i].wm_num;
          wmtrk[l].pm_match = false;
          wmtrk[j].pm_match = true;

        }

        wmtrk[j].pm_match = true;

        pmtrk[i].ing_imod  = wmtrk[j].ing_imod;
        pmtrk[i].ing_fmod  = wmtrk[j].ing_fmod;
        pmtrk[i].ing_ipln  = wmtrk[j].ing_ipln;
        pmtrk[i].ing_fpln  = wmtrk[j].ing_fpln;
        pmtrk[i].wm_ipln   = wmtrk[j].ipln;
        pmtrk[i].wm_fpln   = wmtrk[j].fpln;
        pmtrk[i].ing_trk   = wmtrk[j].ing_trk;
        pmtrk[i].ing_num   = wmtrk[j].ing_num;
        pmtrk[i].wm_trk    = true;
        pmtrk[i].wm_num    = j;
        pmtrk[i].ing_stop  = wmtrk[j].ing_stop;	  
        pmtrk[i].iron_pene = wmtrk[j].iron_pene;

        pmtrk[i].diff_pos  = diff_pos;
        pmtrk[i].diff_ang  = diff_ang;
        pmtrk[i].diff_time = pmtrk[i].clstime - wmtrk[j].clstime;


        jointed = true;
        joitra  = i;
        joilik  = sqrt( fabs(diff_ang)*fabs(diff_ang)/ang_th_pm/ang_th_pm
		       +fabs(diff_pos)*fabs(diff_pos)/pos_th_pm/pos_th_pm);
        hasingtrk = true;

#ifdef DEBUG_THREEDIMRECON
        cout << "  Jointed " << endl;
#endif

      }
    }
  }

//#ifdef USEINGTRK 
//
//  for(int i=0; i<(int)pmtrk.size(); i++){
//    if(pmtrk[i].wm_trk){
//      int wmnum  = pmtrk[i].wm_num;
//      int ingnum = wmtrk[wmnum].ing_num;
//      float xip[Cpln], yip[Cpln], xeip[Cpln], yeip[Cpln];
//      int ntp    = 0;
//
//      for(int hitnum=0; hitnum<(int)pmtrk[i].hit.size(); hitnum++){
//        int view  = pmtrk[i].hit[hitnum].view;
//        int pln   = pmtrk[i].hit[hitnum].pln;
//        int ch    = pmtrk[i].hit[hitnum].ch;
//        xip[ntp]  = zposi(mod,view,pln,ch,0);
//        yip[ntp]  = xyposi(mod,view,pln,ch,0);
//        xeip[ntp] = scithick(mod,view,pln,yip[ntp],0);
//        yeip[ntp] = sciwidth(mod,view,pln,ch,0);
//
//        xip[ntp] += C_B2CHPosZ;
//        if(vertical) yip[ntp] += C_B2CHPosX;
//        else         yip[ntp] += C_B2CHPosY;
//
//        ntp++;
//      }
//
//      for(int hitnum=0; hitnum<(int)wmtrk[wmnum].hit.size(); hitnum++){
//        int view  = wmtrk[wmnum].hit[hitnum].view;
//        int pln   = wmtrk[wmnum].hit[hitnum].pln;
//        int ch    = wmtrk[wmnum].hit[hitnum].ch;
//        xip[ntp]  = zposi(   WMMODNUM,view,pln,ch,0);
//        yip[ntp]  = xyposi(  WMMODNUM,view,pln,ch,0);
//        xeip[ntp] = scithick(WMMODNUM,view,pln,yip[ntp],0);
//        yeip[ntp] = sciwidth(WMMODNUM,view,pln,ch,0);
//
//        xip[ntp] += C_B2WMPosZ;
//        if(vertical) yip[ntp] += C_B2WMPosX;
//        else         yip[ntp] += C_B2WMPosY;
//
//        ntp++;
//      }
//
//      if(wmtrk[wmnum].ing_trk){
//        for(int hitnum=0; hitnum<(int)itrk[ingnum].hit.size(); hitnum++){
//          int ingmod = itrk[ingnum].mod ;
//          int view   = itrk[ingnum].hit[hitnum].view;
//          int pln    = itrk[ingnum].hit[hitnum].pln;
//          int ch     = itrk[ingnum].hit[hitnum].ch;
//          xip[ntp]  = zposi (ingmod,view,pln,ch,0);
//          yip[ntp]  = xyposi(ingmod,view,pln,ch,0);
//          xeip[ntp] = scithick(ingmod,view,pln,yip[ntp],0);
//          yeip[ntp] = sciwidth(ingmod,view,pln,ch,0);
//
//          xip[ntp] += C_B2INGPosZ;
//          if(vertical) yip[ntp] += C_B2INGPosX;
//          else         yip[ntp] += C_B2INGPosY;
//
//          ntp++;
//        }
//      }
//
//      SortArray(ntp,xip,yip,xeip,yeip);
//      int cnt=0;
//      float last_x=-1e+5;
//      float max_y =-1e+5;
//      float min_y = 1e+5;
//      for(int itp=0;itp<ntp;itp++){
//        if(last_x==xip[itp]){
//          if(max_y<yip[itp]){max_y=yip[itp];}
//          if(min_y>yip[itp]){min_y=yip[itp];}
//          cnt++;
//          yip [itp-1] *= cnt;
//          yip [itp-1] += yip[itp];
//          yip [itp-1] /= cnt+1;
//          yeip[itp-1] = max_y-min_y+yeip[itp];
//          for(int jtp=itp;jtp<ntp-1;jtp++){
//            xip [jtp] = xip [jtp+1];
//            yip [jtp] = yip [jtp+1];
//            xeip[jtp] = xeip[jtp+1];
//            yeip[jtp] = yeip[jtp+1];
//          }
//          ntp--;
//          itp--;
//        }
//        else{
//          cnt=0;
//          last_x=-1e+5;
//          max_y =-1e+5;
//          min_y = 1e+5;
//          last_x = xip[itp];
//        }
//      }
//
//      TGraphErrors *gip = new TGraphErrors(ntp,xip,yip,xeip,yeip);
//      TF1 *fip = new TF1("fip","[0]+[1]*x");
//      fip->SetParameters( yip[0]-xip[0]*(yip[ntp-1]-yip[0])/(xip[ntp-1]-xip[0]),
//			 (yip[ntp-1]-yip[0])/(xip[ntp-1]-xip[0]));
//      gip->Fit("fip","Q");
//
//      TF1 *funcip    = gip->GetFunction("fip");
//      float intcptip = funcip->GetParameter(0);
//      float slopeip  = funcip->GetParameter(1);
//
//#ifdef DEBUG_THREEDIMRECON
//      cout << " Fit:" 
//        << " slope=" << slopeip
//        << " ang="   << atan(slopeip)*180./PI
//        << endl;
//#endif
//      gip->Delete();
//      fip->Delete();
//
//      pmtrk[i].ang    = atan(slopeip)*180./PI;
//      pmtrk[i].slope  = slopeip;
//      pmtrk[i].intcpt = pmtrk[i].ixy-pmtrk[i].slope*pmtrk[i].iz;
//      //pmtrk[i].intcpt = intcptip + slopeip*C_B2CHPosZ;
//
//      wmtrk[wmnum].ang    = atan(slopeip)*180./PI;
//      wmtrk[wmnum].slope  = slopeip;
//      wmtrk[wmnum].intcpt = wmtrk[wmnum].ixy-wmtrk[wmnum].slope*wmtrk[wmnum].iz;
//      //wmtrk[wmnum].intcpt = intcptip + slopeip*C_B2WMPosZ;
//
//      if(pmtrk[i].ing_trk){
//        itrk[ingnum].ang    = atan(slopeip)*180./PI;
//        itrk[ingnum].slope  = slopeip;
//        itrk[ingnum].intcpt = itrk[ingnum].ixy - itrk[ingnum].slope*itrk[ingnum].iz;
//        //itrk[ingnum].intcpt = intcptip + slopeip*C_B2INGPosZ;
//      }
//
//    }
//  }
//#endif

  return hasingtrk;
};


float recalcMuCL(TrackWM &htrk, TrackWM &vtrk, int hstart, int vstart){
  float remucl;

  float trkang = 180./PI
    *atan(sqrt( pow(tan(htrk.ang*PI/180.),2)
	       +pow(tan(vtrk.ang*PI/180.),2)));
  initMuCL();
  addMuCLRe(htrk.hit,trkang,htrk.ipln,vtrk.intcpt,vtrk.slope,hstart);
  addMuCLRe(vtrk.hit,trkang,vtrk.ipln,htrk.intcpt,htrk.slope,vstart);

#ifdef USEINGPID
  if(htrk.ing_trk&&vtrk.ing_trk){
    addIngMuCL(hingtrack[htrk.ing_num].hit,
	       trkang,
	       vtrk.intcpt-(vingtrack[vtrk.ing_num].mod-INGMODNUM_mid)*C_INGSpace,
	       vtrk.slope);
    addIngMuCL(vingtrack[vtrk.ing_num].hit, 
	       trkang, 
	       htrk.intcpt, 
	       htrk.slope);
  }
  else{
    if(!(htrk.stop&&vtrk.stop)){
      int ingmod = Ingmod(htrk.intcpt,htrk.slope,vtrk.intcpt,vtrk.slope);
      if(ingmod>=INGMODNUM_start && ingmod<=INGMODNUM_end){
        addIngHitMuCL(ingmod,0,trkang,htrk.intcpt,htrk.slope,
		      vtrk.intcpt-(ingmod-INGMODNUM_mid)*C_INGSpace,vtrk.slope);
        addIngHitMuCL(ingmod,1,trkang,vtrk.intcpt-(ingmod-INGMODNUM_mid)*C_INGSpace,vtrk.slope,
		      htrk.intcpt,htrk.slope);
      }
    }
  }
#endif

  remucl=calcMuCL();

  return remucl;
};

float recalcMuCL(TrackPM &htrk, TrackPM &vtrk, int hstart, int vstart){
  float remucl;

  float trkang = 180./PI
    *atan(sqrt(pow(tan(htrk.ang*PI/180.),2)
	       +pow(tan(vtrk.ang*PI/180.),2)));
  initMuCL();
  addMuCLRe(htrk.hit,trkang,htrk.ipln,vtrk.intcpt,vtrk.slope, hstart);
  addMuCLRe(vtrk.hit,trkang,vtrk.ipln,htrk.intcpt,htrk.slope, vstart);

#ifdef USEINGPID
  if(htrk.ing_trk&&vtrk.ing_trk){
    addIngMuCL(hingtrack[htrk.ing_num].hit,
	       trkang,
	       vtrk.intcpt-(vingtrack[vtrk.ing_num].mod-INGMODNUM_mid)*C_INGSpace,
	       vtrk.slope);
    addIngMuCL(vingtrack[vtrk.ing_num].hit, 
	       trkang, 
	       htrk.intcpt, 
	       htrk.slope);
  }
  else{
    if(!(htrk.stop&&vtrk.stop)){
      int ingmod = Ingmod(htrk.intcpt,htrk.slope,vtrk.intcpt,vtrk.slope);
      if(ingmod>=INGMODNUM_start && ingmod<=INGMODNUM_end){
        addIngHitMuCL(ingmod,0,trkang,htrk.intcpt,htrk.slope,
		      vtrk.intcpt-(ingmod-INGMODNUM_mid)*C_INGSpace,vtrk.slope);
        addIngHitMuCL(ingmod,1,trkang,vtrk.intcpt-(ingmod-INGMODNUM_mid)*C_INGSpace,vtrk.slope,
		      htrk.intcpt,htrk.slope);
      }
    }
  }
#endif

  remucl=calcMuCL();

  return remucl;
};


void fTrackMatch(Trk &trk, TrackIng &htrk, TrackIng &vtrk){
#ifdef DEBUG_THREEDIMRECON
  cout << "  >>> fTrackMatch ING" << endl;
#endif

  trk.oneview  = 0;

  trk.ing_trk  = true;
  trk.pm_stop  = false;
  trk.ing_stop = (htrk.stop && vtrk.stop);

  trk.startmod  = INGMODNUM_mid;
  trk.stopmodx  = INGMODNUM_mid;
  trk.stopmody  = INGMODNUM_mid;

  trk.startxpln = htrk.ipln;
  trk.startypln = vtrk.ipln;
  trk.startxch  = htrk.ixy;
  trk.startych  = vtrk.ixy;
  trk.x         = htrk.ixy;
  trk.y         = vtrk.ixy;
  trk.endxpln   = htrk.fpln;
  trk.endypln   = vtrk.fpln;
  trk.endxch    = htrk.fxy;
  trk.endych    = vtrk.fxy;	
  trk.thetax    = htrk.ang;
  trk.thetay    = vtrk.ang;
  trk.intcptx   = htrk.intcpt;
  trk.intcpty   = vtrk.intcpt;
  trk.slopex    = htrk.slope;
  trk.slopey    = vtrk.slope;

  float trkang = 180./PI
    *atan(sqrt(pow(tan(htrk.ang*PI/180.),2)
          +pow(tan(vtrk.ang*PI/180.),2)));

  trk.angle         = trkang;
  trk.vetowtracking = htrk.veto||vtrk.veto;
  trk.edgewtracking = htrk.edge||vtrk.edge;

  trk.ing_startmod = INGMODNUM_mid;
  trk.ing_endmod   = INGMODNUM_mid;

  if(vtrk.ipln<htrk.ipln){
    trk.ing_startpln = vtrk.ipln;
  }
  else{
    trk.ing_startpln = htrk.ipln;
  }

  if(vtrk.fpln > htrk.fpln){
    trk.ing_endpln = vtrk.fpln;
  }
  else{
    trk.ing_endpln = htrk.fpln;
  }

  int h_iron_pene = htrk.fpln-htrk.ipln;
  int v_iron_pene = vtrk.fpln-vtrk.ipln;
  if(v_iron_pene>h_iron_pene){
    trk.iron_pene  = v_iron_pene;
    trk.iron_range = v_iron_pene/cos(trkang*PI/180.);
  }
  else{
    trk.iron_pene  = h_iron_pene;
    trk.iron_range = h_iron_pene/cos(trkang*PI/180.);
  }

  if((vtrk.fpln-vtrk.ipln)>(htrk.fpln-htrk.ipln)){
    trk.sci_range = (vtrk.fpln-vtrk.ipln)/cos(trkang*PI/180.);
  }
  else{
    trk.sci_range = (htrk.fpln-htrk.ipln)/cos(trkang*PI/180.);
  }

  float totalpe  = 0.;
  int   totalhit = 0;
  int   trkpdg[7];
  memset(trkpdg,0,sizeof(trkpdg));
  fAddPE(htrk.hit,trkang,totalpe,totalhit,trkpdg,htrk.ipln,vtrk.intcpt,vtrk.slope);
  fAddPE(vtrk.hit,trkang,totalpe,totalhit,trkpdg,vtrk.ipln,htrk.intcpt,htrk.slope);

  if(totalhit>0) trk.trkpe = totalpe/totalhit;
  else           trk.trkpe = 0;
  trk.pdg=num2pdg(trkpdg);

  initMuCL();
  addMuCL(htrk.hit,trkang,htrk.ipln,vtrk.intcpt,vtrk.slope);
  addMuCL(vtrk.hit,trkang,vtrk.ipln,htrk.intcpt,htrk.slope);

#ifdef USEINGPID
  if(htrk.ing_trk&&vtrk.ing_trk){
    addIngMuCL(hingtrack[htrk.ing_num].hit,
	       trkang, 
	       vtrk.intcpt-(vingtrack[vtrk.ing_num].mod-INGMODNUM_mid)*C_INGSpace,
	       vtrk.slope);
    addIngMuCL(vingtrack[vtrk.ing_num].hit,
	       trkang,
	       htrk.intcpt,
	       htrk.slope);
  }
  else{
    if(!trk.pm_stop){
      int ingmod = Ingmod(htrk.intcpt,htrk.slope,vtrk.intcpt,vtrk.slope);
      if(ingmod>=0 && ingmod<MOD_INGRID_H){
        addIngHitAng(ingmod,trk,htrk,vtrk, 0 ,1);
        addIngHitMuCL(ingmod,0,trkang,htrk.intcpt,htrk.slope,
		      vtrk.intcpt-(ingmod-INGMODNUM_mid)*C_INGSpace,vtrk.slope);
        addIngHitMuCL(ingmod,1,trkang,vtrk.intcpt-(ingmod-INGMODNUM_mid)*C_INGSpace,vtrk.slope,
		      htrk.intcpt,htrk.slope);
      }
    }
  }
#endif

  trk.mucl = calcMuCL();
};

void fTrackMatch(Trk &trk, TrackWM &htrk, TrackWM &vtrk){
#ifdef DEBUG_THREEDIMRECON
  cout << "  >>> fTrackMatch WM" << endl;
#endif

  trk.oneview  = 0;

  trk.ing_trk  = (htrk.ing_trk && vtrk.ing_trk);
  trk.pm_stop  = (htrk.stop    && vtrk.stop);
  trk.ing_stop = (htrk.ing_stop&& vtrk.ing_stop);

  if(htrk.ing_trk){
    trk.diff_posx  = htrk.diff_pos;
    trk.diff_timex = htrk.diff_time;
    trk.diff_angx  = htrk.diff_ang;
  }
  if(vtrk.ing_trk){
    trk.diff_posy  = vtrk.diff_pos;
    trk.diff_timey = vtrk.diff_time;
    trk.diff_angy  = vtrk.diff_ang;
  }
  


  htrk.ing_trk = trk.ing_trk;
  vtrk.ing_trk = trk.ing_trk;

  trk.startmod  = WMMODNUM;
  if(htrk.ing_trk){ trk.stopmodx = INGMODNUM_mid;}
  else            { trk.stopmodx = WMMODNUM     ;}
  if(vtrk.ing_trk){ trk.stopmody = INGMODNUM_mid;}
  else            { trk.stopmody = WMMODNUM     ;}

  trk.startxpln = htrk.ipln;
  trk.startypln = vtrk.ipln;
  trk.startxch  = htrk.ixy;
  trk.startych  = vtrk.ixy;
  trk.x         = htrk.ixy;
  trk.y         = vtrk.ixy;
  if(htrk.ing_trk){
    trk.endxpln   = hingtrack[htrk.ing_num].fpln;
    trk.endxch    = hingtrack[htrk.ing_num].fxy;
  }else{
    trk.endxpln   = htrk.fpln;
    trk.endxch    = htrk.fxy;
  }
  if(vtrk.ing_trk){
    trk.endypln   = vingtrack[vtrk.ing_num].fpln;
    trk.endych    = vingtrack[vtrk.ing_num].fxy;
  }else{
    trk.endypln   = vtrk.fpln;
    trk.endych    = vtrk.fxy;
  }

  trk.thetax    = htrk.ang;
  trk.thetay    = vtrk.ang;
  trk.intcptx   = htrk.intcpt;
  trk.intcpty   = vtrk.intcpt;
  trk.slopex    = htrk.slope;
  trk.slopey    = vtrk.slope;

  float trkang = 180./PI
    *atan(sqrt(pow(tan(htrk.ang*PI/180.),2)
          +pow(tan(vtrk.ang*PI/180.),2)));

  trk.angle         = trkang;
  trk.vetowtracking = htrk.veto||vtrk.veto;
  trk.edgewtracking = htrk.edge||vtrk.edge;

  if(abs(vtrk.ing_imod-INGMODNUM_mid)<abs(htrk.ing_imod-INGMODNUM_mid)){
    trk.ing_startmod = vtrk.ing_imod;
  }
  else{
    trk.ing_startmod = htrk.ing_imod;
  }

  if(abs(vtrk.ing_fmod-INGMODNUM_mid)>abs(htrk.ing_fmod-INGMODNUM_mid)){
    trk.ing_endmod = vtrk.ing_fmod;
  }
  else{
    trk.ing_endmod = htrk.ing_fmod;
  }

  if(vtrk.ing_ipln<htrk.ing_ipln){
    trk.ing_startpln = vtrk.ing_ipln;
  }
  else{
    trk.ing_startpln = htrk.ing_ipln;
  }

  if(vtrk.ing_fpln > htrk.ing_fpln){
    trk.ing_endpln = vtrk.ing_fpln;
  }
  else{
    trk.ing_endpln = htrk.ing_fpln;
  }


  if(!trk.ing_trk){
    trk.iron_pene  = 0;
    trk.iron_range = 0;
  }
  else if(vtrk.iron_pene>htrk.iron_pene){
    trk.iron_pene  = vtrk.iron_pene;
    trk.iron_range = vtrk.iron_pene/cos(trkang*PI/180.);
  }
  else{
    trk.iron_pene  = htrk.iron_pene;
    trk.iron_range = htrk.iron_pene/cos(trkang*PI/180.);
  }

  if((vtrk.fpln-vtrk.ipln)>(htrk.fpln-htrk.ipln)){
    trk.sci_range = (vtrk.fpln-vtrk.ipln)/cos(trkang*PI/180.);
  }
  else{
    trk.sci_range = (htrk.fpln-htrk.ipln)/cos(trkang*PI/180.);
  }


  float totalpe  = 0.;
  int   totalhit = 0;
  int   trkpdg[7];
  memset(trkpdg,0,sizeof(trkpdg));
  fAddPE(htrk.hit,trkang,totalpe,totalhit,trkpdg,htrk.ipln,vtrk.intcpt,vtrk.slope);
  fAddPE(vtrk.hit,trkang,totalpe,totalhit,trkpdg,vtrk.ipln,htrk.intcpt,htrk.slope);

  if(totalhit>0) trk.trkpe = totalpe/totalhit;
  else           trk.trkpe = 0;
  trk.pdg=num2pdg(trkpdg);

  initMuCL();
  addMuCL(htrk.hit,trkang,htrk.ipln,vtrk.intcpt,vtrk.slope);
  addMuCL(vtrk.hit,trkang,vtrk.ipln,htrk.intcpt,htrk.slope);

#ifdef USEINGPID
  if(htrk.ing_trk&&vtrk.ing_trk){
    addIngMuCL(hingtrack[htrk.ing_num].hit,
	       trkang, 
	       vtrk.intcpt-(vingtrack[vtrk.ing_num].mod-INGMODNUM_mid)*C_INGSpace,
	       vtrk.slope);
    addIngMuCL(vingtrack[vtrk.ing_num].hit,
	       trkang,
	       htrk.intcpt,
	       htrk.slope);
  }
  else{
    if(!trk.pm_stop){
      int ingmod = Ingmod(htrk.intcpt,htrk.slope,vtrk.intcpt,vtrk.slope);
      if(ingmod>=0 && ingmod<MOD_INGRID_H){
        addIngHitAng(ingmod,trk,htrk,vtrk, 0 ,1);
        addIngHitMuCL(ingmod,0,trkang,htrk.intcpt,htrk.slope,
		      vtrk.intcpt-(ingmod-INGMODNUM_mid)*C_INGSpace,vtrk.slope);
        addIngHitMuCL(ingmod,1,trkang,vtrk.intcpt-(ingmod-INGMODNUM_mid)*C_INGSpace,vtrk.slope,
		      htrk.intcpt,htrk.slope);
      }
    }
  }
#endif

  trk.mucl = calcMuCL();
};

void fTrackMatch(Trk &trk, TrackPM &htrk, TrackPM &vtrk){
#ifdef DEBUG_THREEDIMRECON
  cout << "  >>> fTrackMatch PM" << endl;
#endif

  trk.oneview  = 0;

  trk.ing_trk  = (htrk.ing_trk && vtrk.ing_trk );
  trk.pm_stop  = (htrk.stop    && vtrk.stop    );
  trk.ing_stop = (htrk.ing_stop&& vtrk.ing_stop);
  bool wg_trk  = (htrk.wm_trk  && vtrk.wm_trk  );


  if(htrk.ing_trk&&htrk.wm_trk){
    trk.diff_posx   = hwmtrack[htrk.wm_num].diff_pos;
    trk.diff_timex  = hwmtrack[htrk.wm_num].diff_time;
    trk.diff_angx   = hwmtrack[htrk.wm_num].diff_ang;
    trk.diff_posx2  = htrk.diff_pos;
    trk.diff_timex2 = htrk.diff_time;
    trk.diff_angx2  = htrk.diff_ang;
  }
  if(htrk.ing_trk&&!htrk.wm_trk){
    trk.diff_posx3  = htrk.diff_pos;
    trk.diff_timex3 = htrk.diff_time;
    trk.diff_angx3  = htrk.diff_ang;
  }
  if(vtrk.ing_trk&&vtrk.wm_trk){
    trk.diff_posy   = vwmtrack[vtrk.wm_num].diff_pos;
    trk.diff_timey  = vwmtrack[vtrk.wm_num].diff_time;
    trk.diff_angy   = vwmtrack[vtrk.wm_num].diff_ang;
    trk.diff_posy2  = vtrk.diff_pos;
    trk.diff_timey2 = vtrk.diff_time;
    trk.diff_angy2  = vtrk.diff_ang;
  }
  if(vtrk.ing_trk&&!vtrk.wm_trk){
    trk.diff_posy3  = vtrk.diff_pos;
    trk.diff_timey3 = vtrk.diff_time;
    trk.diff_angy3  = vtrk.diff_ang;
  }

  htrk.ing_trk = trk.ing_trk;
  vtrk.ing_trk = trk.ing_trk;
  if(!trk.ing_trk){
    htrk.wm_trk  = wg_trk;
    vtrk.wm_trk  = wg_trk;
  }

  trk.startmod  = PMMODNUM;
  if     ( htrk.ing_trk){ trk.stopmodx = INGMODNUM_mid;}
  else if( htrk.wm_trk ){ trk.stopmodx = WMMODNUM     ;}
  else                  { trk.stopmodx = PMMODNUM     ;}
  if     ( vtrk.ing_trk){ trk.stopmody = INGMODNUM_mid;}
  else if( vtrk.wm_trk ){ trk.stopmody = WMMODNUM     ;}
  else                  { trk.stopmody = PMMODNUM     ;}

  trk.startxpln = htrk.ipln;
  trk.startypln = vtrk.ipln;
  trk.startxch  = htrk.ixy;
  trk.startych  = vtrk.ixy;
  trk.x         = htrk.ixy;
  trk.y         = vtrk.ixy;
  if(htrk.ing_trk){
    trk.endxpln   = hingtrack[htrk.ing_num].fpln;
    trk.endxch    = hingtrack[htrk.ing_num].fxy;
  }else if(htrk.wm_trk){
    trk.endxpln   = hwmtrack[htrk.wm_num].fpln;
    trk.endxch    = hwmtrack[htrk.wm_num].fxy;
  }else{
    trk.endxpln   = htrk.fpln;
    trk.endxch    = htrk.fxy;
  }
  if(vtrk.ing_trk){
    trk.endypln   = vingtrack[vtrk.ing_num].fpln;
    trk.endych    = vingtrack[vtrk.ing_num].fxy;
  }else if(vtrk.wm_trk){
    trk.endypln   = vwmtrack[vtrk.wm_num].fpln;
    trk.endych    = vwmtrack[vtrk.wm_num].fxy;
  }else{
    trk.endypln   = vtrk.fpln;
    trk.endych    = vtrk.fxy;
  }


  trk.thetax    = htrk.ang;
  trk.thetay    = vtrk.ang;
  trk.intcptx   = htrk.intcpt;
  trk.intcpty   = vtrk.intcpt;
  trk.slopex    = htrk.slope;
  trk.slopey    = vtrk.slope;

  float trkang=180./PI
    *atan(sqrt( pow(tan(htrk.ang*PI/180.),2)
	       +pow(tan(vtrk.ang*PI/180.),2)));

  trk.angle         = trkang;
  trk.vetowtracking = htrk.veto||vtrk.veto;
  trk.edgewtracking = htrk.edge||vtrk.edge;

  if(abs(vtrk.ing_imod-INGMODNUM_mid)<abs(htrk.ing_imod-INGMODNUM_mid)){
    trk.ing_startmod = vtrk.ing_imod;
  }
  else{
    trk.ing_startmod = htrk.ing_imod;
  }

  if(vtrk.ing_trk&&htrk.ing_trk){
    trk.ing_endmod = INGMODNUM_mid;
  }
  else if(vtrk.wm_trk&&htrk.wm_trk){
    trk.ing_endmod = WMMODNUM;
  }
  else{
    trk.ing_endmod = -1;
  }


  if(vtrk.ing_trk&&htrk.ing_trk){
    if(vtrk.ing_ipln < htrk.ing_ipln){
      trk.ing_startpln = vtrk.ing_ipln;
    }
    else{
      trk.ing_startpln = htrk.ing_ipln;
    }

    if(vtrk.ing_fpln > htrk.ing_fpln){
      trk.ing_endpln = vtrk.ing_fpln;
    }
    else{
      trk.ing_endpln = htrk.ing_fpln;
    }
  }
  else if(vtrk.wm_trk&&htrk.wm_trk){
    if(vtrk.wm_ipln < htrk.wm_ipln){
      trk.ing_startpln = vtrk.wm_ipln;
    }
    else{
      trk.ing_startpln = htrk.wm_ipln;
    }

    if(vtrk.wm_fpln > htrk.wm_fpln){
      trk.ing_endpln = vtrk.wm_fpln;
    }
    else{
      trk.ing_endpln = htrk.wm_fpln;
    }
  }


  if(!trk.ing_trk){
    trk.iron_pene  = 0;
    trk.iron_range = 0;
  }
  else if(vtrk.iron_pene>htrk.iron_pene){
    trk.iron_pene  = vtrk.iron_pene;
    trk.iron_range = vtrk.iron_pene/cos(trkang*PI/180.);
  }
  else{
    trk.iron_pene  = htrk.iron_pene;
    trk.iron_range = htrk.iron_pene/cos(trkang*PI/180.);
  }

  if((vtrk.fpln-vtrk.ipln)>(htrk.fpln-htrk.ipln))
    trk.sci_range = (vtrk.fpln-vtrk.ipln)/cos(trkang*PI/180.);
  else
    trk.sci_range = (htrk.fpln-htrk.ipln)/cos(trkang*PI/180.);



  float totalpe  = 0.;
  int   totalhit = 0;
  int   trkpdg[7];
  memset(trkpdg,0,sizeof(trkpdg));
  fAddPE(htrk.hit,trkang,totalpe,totalhit,trkpdg,htrk.ipln,vtrk.intcpt,vtrk.slope);
  fAddPE(vtrk.hit,trkang,totalpe,totalhit,trkpdg,vtrk.ipln,htrk.intcpt,htrk.slope);
  if(totalhit>0) trk.trkpe=totalpe/totalhit;
  else           trk.trkpe=0;
  trk.pdg=num2pdg(trkpdg);
  
  initMuCL();
  addMuCL(htrk.hit,trkang,htrk.ipln,vtrk.intcpt,vtrk.slope);
  addMuCL(vtrk.hit,trkang,vtrk.ipln,htrk.intcpt,htrk.slope);
  
  trk.mucl = calcMuCL();
};


void fTrackMatchBack(Trk &trk, TrackWM &htrk, TrackWM &vtrk){
#ifdef DEBUG_THREEDIMRECON
  cout << "  >>> fTrackMatchBack wM" << endl;
#endif

  trk.oneview  = 0;

  trk.startxpln=htrk.fpln;
  trk.startypln=vtrk.fpln;
  trk.startxch=htrk.fxy;
  trk.startych=vtrk.fxy;
  trk.x=htrk.fxy;
  trk.y=vtrk.fxy;
  trk.endxpln=htrk.ipln;
  trk.endypln=vtrk.ipln;
  trk.endxch=htrk.ixy;
  trk.endych=vtrk.ixy;	

  if(htrk.ang>0) trk.thetax = -180.+htrk.ang;
  else           trk.thetax =  180.+htrk.ang;

  if(vtrk.ang>0) trk.thetay = -180.+vtrk.ang;
  else           trk.thetay =  180.+vtrk.ang;

  trk.intcptx=htrk.intcpt;
  trk.intcpty=vtrk.intcpt;
  trk.slopex=htrk.slope;
  trk.slopey=vtrk.slope;

  float trkang = 180.-180./PI
    *atan(sqrt( pow(tan(htrk.ang*PI/180.),2)
	       +pow(tan(vtrk.ang*PI/180.),2)));

  trk.angle         = trkang;
  trk.vetowtracking = false;
  trk.edgewtracking = false;

  if(abs(vtrk.ing_imod-INGMODNUM_mid)<abs(htrk.ing_imod-INGMODNUM_mid)){
    trk.ing_startmod = vtrk.ing_imod;
  }
  else{
    trk.ing_startmod = htrk.ing_imod;
  }

  if(abs(vtrk.ing_fmod-INGMODNUM_mid)>abs(htrk.ing_fmod-INGMODNUM_mid)){
    trk.ing_endmod = vtrk.ing_fmod;
  }
  else{
    trk.ing_endmod = htrk.ing_fmod;
  }

  if(vtrk.ing_ipln < htrk.ing_ipln){
    trk.ing_startpln = vtrk.ing_ipln;
  }
  else{
    trk.ing_startpln = htrk.ing_ipln;
  }

  if(vtrk.ing_fpln > htrk.ing_fpln){
    trk.ing_endpln = vtrk.ing_fpln;
  }
  else{
    trk.ing_endpln = htrk.ing_fpln;
  }

  trk.ing_trk  = (htrk.ing_trk  || vtrk.ing_trk);
  trk.pm_stop  = (htrk.stop     && vtrk.stop);
  trk.ing_stop = (htrk.ing_stop && vtrk.ing_stop);

  if(!trk.ing_trk){
    trk.iron_pene  = 0;
    trk.iron_range = 0;
  }
  else if(vtrk.iron_pene > htrk.iron_pene){
    trk.iron_pene  = vtrk.iron_pene;
    trk.iron_range = vtrk.iron_pene/cos(trkang*PI/180.);
  }
  else{
    trk.iron_pene  = htrk.iron_pene;
    trk.iron_range = htrk.iron_pene/cos(trkang*PI/180.);
  }

  if((vtrk.fpln-vtrk.ipln)>(htrk.fpln-htrk.ipln)){
    trk.sci_range = (vtrk.fpln-vtrk.ipln)/cos(trkang*PI/180.);
  }
  else{
    trk.sci_range = (htrk.fpln-htrk.ipln)/cos(trkang*PI/180.);
  }

  float totalpe  = 0.;
  int   totalhit = 0;
  int   trkpdg[7];
  memset(trkpdg,0,sizeof(trkpdg));
  fAddPE(htrk.hit,trkang,totalpe,totalhit,trkpdg,htrk.fpln-1,vtrk.intcpt,vtrk.slope);
  fAddPE(vtrk.hit,trkang,totalpe,totalhit,trkpdg,vtrk.fpln-1,htrk.intcpt,htrk.slope);

  if(totalhit>0) trk.trkpe = totalpe/totalhit;
  else           trk.trkpe = 0;
  trk.pdg = num2pdg(trkpdg);

  initMuCL();
  addMuCL(htrk.hit,trkang,htrk.fpln-1,vtrk.intcpt,vtrk.slope);
  addMuCL(vtrk.hit,trkang,vtrk.fpln-1,htrk.intcpt,htrk.slope);
  trk.mucl=calcMuCL();
};

void fTrackMatchX(Trk &trk,AnaTrack &pmtrk, TrackWM &htrk){
#ifdef DEBUG_THREEDIMRECON
  cout << "  >>> fTrackMatchX WM" << endl;
#endif

  trk.oneview  = 1;

  trk.startmod = WMMODNUM;
  trk.stopmodx = WMMODNUM;

  trk.startxpln = htrk.ipln;
  trk.startypln = pmtrk.trk[0].startypln;
  trk.startxch  = htrk.ixy;
  trk.startych  = pmtrk.trk[0].startych;
  trk.x         = htrk.ixy;
  trk.y         = pmtrk.trk[0].y;
  trk.endxpln   = htrk.fpln;
  trk.endypln   = htrk.fpln;
  trk.endxch    = htrk.fxy;
  trk.endych    = pmtrk.trk[0].endych;
  trk.thetax    = htrk.ang;
  trk.thetay    = pmtrk.trk[0].thetay;
  trk.intcptx   = htrk.intcpt;
  trk.intcpty   = pmtrk.trk[0].intcpty;
  trk.slopex    = htrk.slope;
  trk.slopey    = pmtrk.trk[0].slopey;

  float trkang = 180./PI
    *atan(sqrt( pow(tan(htrk.ang*PI/180.),2)
	       +pow(tan(pmtrk.trk[0].thetay*PI/180.),2)));

  trk.angle         = trkang;
  trk.vetowtracking = htrk.veto||pmtrk.trk[0].vetowtracking;
  trk.edgewtracking = htrk.edge||pmtrk.trk[0].edgewtracking;
  trk.ing_trk       = false;
  trk.pm_stop       = htrk.stop;

  trk.iron_pene  = 0;
  trk.iron_range = 0;
  trk.sci_range  = (htrk.fpln-htrk.ipln)/cos(trkang*PI/180.);

  float totalpe = 0;
  int totalhit  = 0;
  int trkpdg[7];
  memset(trkpdg,0,sizeof(trkpdg));
  fAddPE(htrk.hit,trkang,totalpe,totalhit,trkpdg,htrk.ipln,
      pmtrk.trk[0].intcpty,pmtrk.trk[0].slopey);

  if(totalhit>0) trk.trkpe=totalpe/totalhit;
  else           trk.trkpe=0;
  trk.pdg=num2pdg(trkpdg);

  initMuCL();
  addMuCL(htrk.hit,trkang,htrk.ipln,
      pmtrk.trk[0].intcpty,pmtrk.trk[0].slopey);

  trk.mucl=calcMuCL();
};

void fTrackMatchX(Trk &trk, AnaTrack &pmtrk, TrackPM &htrk){
#ifdef DEBUG_THREEDIMRECON
  cout << "  >>> fTrackMatchX PM" << endl;
#endif

  trk.oneview  = 1;

  trk.startmod = PMMODNUM;
  trk.stopmodx = PMMODNUM;


  trk.startxpln = htrk.ipln;
  trk.startypln = pmtrk.trk[0].startypln;
  trk.startxch  = htrk.ixy;
  trk.startych  = pmtrk.trk[0].startych;
  trk.x         = htrk.ixy;
  trk.y         = pmtrk.trk[0].y;
  trk.endxpln   = htrk.fpln;
  trk.endypln   = htrk.fpln;
  trk.endxch    = htrk.fxy;
  trk.endych    = pmtrk.trk[0].endych;
  trk.thetax    = htrk.ang;
  trk.thetay    = pmtrk.trk[0].thetay;
  trk.intcptx   = htrk.intcpt;
  trk.intcpty   = pmtrk.trk[0].intcpty;
  trk.slopex    = htrk.slope;
  trk.slopey    = pmtrk.trk[0].slopey;

  float trkang = 180./PI
    *atan(sqrt(pow(tan(htrk.ang*PI/180.),2)
          +pow(tan(pmtrk.trk[0].thetay*PI/180.),2)));

  trk.angle         = trkang;
  trk.vetowtracking = htrk.veto||pmtrk.trk[0].vetowtracking;
  trk.edgewtracking = htrk.edge||pmtrk.trk[0].edgewtracking;
  trk.ing_trk       = false;
  trk.pm_stop       = htrk.stop;

  trk.iron_pene  = 0;
  trk.iron_range = 0;
  trk.sci_range  = (htrk.fpln-htrk.ipln)/cos(trkang*PI/180.);

  float totalpe = 0;
  int totalhit  = 0;
  int trkpdg[7];
  memset(trkpdg,0,sizeof(trkpdg));
  fAddPE(htrk.hit,trkang,totalpe,totalhit,trkpdg,htrk.ipln,
      pmtrk.trk[0].intcpty,pmtrk.trk[0].slopey);

  if(totalhit>0) trk.trkpe=totalpe/totalhit;
  else           trk.trkpe=0;
  trk.pdg=num2pdg(trkpdg);

  initMuCL();
  addMuCL(htrk.hit,trkang,htrk.ipln,
	  pmtrk.trk[0].intcpty,pmtrk.trk[0].slopey);

  trk.mucl=calcMuCL();
};

//Currently not used.
void fTrackMatchBackX(Trk &trk,AnaTrack &pmtrk, TrackWM &htrk){
#ifdef DEBUG_THREEDIMRECON
  cout << "  >>> fTrackMatchBackX WM" << endl;
#endif

  trk.oneview  = 1;

  trk.startmod  = WMMODNUM;
  trk.stopmodx  = WMMODNUM;

  trk.startxpln = htrk.fpln;
  trk.startypln = pmtrk.trk[0].endypln;
  trk.startxch  = htrk.fxy;
  trk.startych  = pmtrk.trk[0].endych;
  trk.x         = htrk.fxy;
  trk.y         = pmtrk.trk[0].endych;
  trk.endxpln   = htrk.ipln;
  trk.endypln   = htrk.ipln;
  trk.endxch    = htrk.ixy;
  trk.endych    = pmtrk.trk[0].startych;

  if(htrk.ang>0) trk.thetax = -180.+htrk.ang;
  else           trk.thetax =  180.+htrk.ang;

  if(pmtrk.trk[0].thetay>0){
    trk.thetay = -180.+pmtrk.trk[0].thetay;
  }
  else{
    trk.thetay = 180.+pmtrk.trk[0].thetay;
  }

  trk.intcptx = htrk.intcpt;
  trk.intcpty = pmtrk.trk[0].intcpty;
  trk.slopex  = htrk.slope;
  trk.slopey  = pmtrk.trk[0].slopey;

  float trkang = 180.-180./PI
    *atan(sqrt( pow(tan(htrk.ang*PI/180.),2)
	       +pow(tan(pmtrk.trk[0].thetay*PI/180.),2)));

  trk.angle         = trkang;
  trk.vetowtracking = false;
  trk.edgewtracking = false;
  trk.ing_trk       = false;
  trk.pm_stop       = htrk.stop;
  trk.iron_pene     = 0;
  trk.iron_range    = 0;
  trk.sci_range     = (htrk.fpln-htrk.ipln)/cos(trkang*PI/180.);

  float totalpe = 0;
  int totalhit  = 0;
  int trkpdg[7];
  memset(trkpdg,0,sizeof(trkpdg));
  fAddPE(htrk.hit,trkang,totalpe,totalhit,trkpdg,htrk.fpln-1,
      pmtrk.trk[0].intcpty,pmtrk.trk[0].slopey);

  if(totalhit>0) trk.trkpe = totalpe/totalhit;
  else           trk.trkpe = 0;
  trk.pdg=num2pdg(trkpdg);

  initMuCL();
  addMuCL(htrk.hit,trkang,htrk.fpln-1,pmtrk.trk[0].intcpty,pmtrk.trk[0].slopey);
  trk.mucl = calcMuCL();
};

void fTrackMatchY(Trk &trk,AnaTrack &pmtrk, TrackWM &vtrk){
#ifdef DEBUG_THREEDIMRECON
  cout << "  >>> fTrackMatchY WM" << endl;
#endif

  trk.oneview  = 2;

  trk.startmod = WMMODNUM;
  trk.stopmody = WMMODNUM;

  trk.startxpln = pmtrk.trk[0].startxpln;
  trk.startypln = vtrk.ipln;
  trk.startxch  = pmtrk.trk[0].startxch;
  trk.startych  = vtrk.ixy;
  trk.x         = pmtrk.trk[0].x;
  trk.y         = vtrk.ixy;
  trk.endxpln   = vtrk.fpln;
  trk.endypln   = vtrk.fpln;
  trk.endxch    = pmtrk.trk[0].endxch;
  trk.endych    = vtrk.fxy;
  trk.thetax    = pmtrk.trk[0].thetax;
  trk.thetay    = vtrk.ang;
  trk.intcptx   = pmtrk.trk[0].intcptx;
  trk.intcpty   = vtrk.intcpt;
  trk.slopex    = pmtrk.trk[0].slopex;
  trk.slopey    = vtrk.slope;

  float trkang = 180./PI
    *atan(sqrt( pow(tan(vtrk.ang*PI/180.),2)
	       +pow(tan(pmtrk.trk[0].thetax*PI/180.),2)));
  trk.angle         = trkang;
  trk.vetowtracking = vtrk.veto||pmtrk.trk[0].vetowtracking;
  trk.edgewtracking = vtrk.edge||pmtrk.trk[0].edgewtracking;
  trk.ing_trk       = false;
  trk.pm_stop       = vtrk.stop;
  trk.iron_pene     = 0;
  trk.iron_range    = 0;
  trk.sci_range     = (vtrk.fpln-vtrk.ipln)/cos(trkang*PI/180.);

  float totalpe = 0;
  int totalhit  = 0;
  int trkpdg[7];
  memset(trkpdg,0,sizeof(trkpdg));
  fAddPE(vtrk.hit,trkang,totalpe,totalhit,trkpdg,vtrk.ipln,
      pmtrk.trk[0].intcptx,pmtrk.trk[0].slopex);

  if(totalhit>0) trk.trkpe = totalpe/totalhit;
  else           trk.trkpe = 0;
  trk.pdg=num2pdg(trkpdg);

  initMuCL();
  addMuCL(vtrk.hit,trkang,vtrk.ipln,pmtrk.trk[0].intcptx,pmtrk.trk[0].slopex);

  trk.mucl=calcMuCL();
};

void fTrackMatchY(Trk &trk, AnaTrack &pmtrk, TrackPM &vtrk){
#ifdef DEBUG_THREEDIMRECON
  cout << "  >>> fTrackMatchY PM" << endl;
#endif

  trk.oneview  = 2;

  trk.startmod = PMMODNUM;
  trk.stopmody = PMMODNUM;

  trk.startxpln = pmtrk.trk[0].startxpln;
  trk.startypln = vtrk.ipln;
  trk.startxch  = pmtrk.trk[0].startxch;
  trk.startych  = vtrk.ixy;
  trk.x         = pmtrk.trk[0].x;
  trk.y         = vtrk.ixy;
  trk.endxpln   = vtrk.fpln;
  trk.endypln   = vtrk.fpln;
  trk.endxch    = pmtrk.trk[0].endxch;
  trk.endych    = vtrk.fxy;
  trk.thetax    = pmtrk.trk[0].thetax;
  trk.thetay    = vtrk.ang;
  trk.intcptx   = pmtrk.trk[0].intcptx;
  trk.intcpty   = vtrk.intcpt;
  trk.slopex    = pmtrk.trk[0].slopex;
  trk.slopey    = vtrk.slope;

  float trkang = 180./PI
    *atan(sqrt(pow(tan(vtrk.ang*PI/180.),2)
          +pow(tan(pmtrk.trk[0].thetax*PI/180.),2)));
  trk.angle         = trkang;
  trk.vetowtracking = vtrk.veto||pmtrk.trk[0].vetowtracking;
  trk.edgewtracking = vtrk.edge||pmtrk.trk[0].edgewtracking;
  trk.ing_trk       = false;
  trk.pm_stop       = vtrk.stop;
  trk.iron_pene     = 0;
  trk.iron_range    = 0;
  trk.sci_range     = (vtrk.fpln-vtrk.ipln)/cos(trkang*PI/180.);

  float totalpe = 0;
  int totalhit  = 0;
  int trkpdg[7];
  memset(trkpdg,0,sizeof(trkpdg));
  fAddPE(vtrk.hit,trkang,totalpe,totalhit,trkpdg,vtrk.ipln,
      pmtrk.trk[0].intcptx,pmtrk.trk[0].slopex);

  if(totalhit>0) trk.trkpe = totalpe/totalhit;
  else           trk.trkpe = 0;
  trk.pdg=num2pdg(trkpdg);

  initMuCL();
  addMuCL(vtrk.hit,trkang,vtrk.ipln,pmtrk.trk[0].intcptx,pmtrk.trk[0].slopex);

  trk.mucl=calcMuCL();
};


//Currently not used.
void fTrackMatchBackY(Trk &trk,AnaTrack &pmtrk, TrackWM &vtrk){
#ifdef DEBUG_THREEDIMRECON
  cout << "  >>> fTrackMatchBackY WM" << endl;
#endif

  trk.oneview  = 2;

  trk.startmod  = WMMODNUM;
  trk.stopmody  = WMMODNUM;

  trk.startxpln = pmtrk.trk[0].endxpln;
  trk.startypln = vtrk.fpln;
  trk.startxch  = pmtrk.trk[0].endxch;
  trk.startych  = vtrk.fxy;
  trk.x         = pmtrk.trk[0].endxch;
  trk.y         = vtrk.fxy;
  trk.endxpln   = vtrk.ipln;
  trk.endypln   = vtrk.ipln;
  trk.endxch    = pmtrk.trk[0].startxch;
  trk.endych    = vtrk.ixy;

  if(pmtrk.trk[0].thetax>0){
    trk.thetax = -180.+pmtrk.trk[0].thetax;
  }
  else{
    trk.thetax =  180.+pmtrk.trk[0].thetax;
  }

  if(vtrk.ang>0) trk.thetay = -180.+vtrk.ang;
  else           trk.thetay =  180.+vtrk.ang;

  trk.intcptx = pmtrk.trk[0].intcptx;
  trk.intcpty = vtrk.intcpt;
  trk.slopex  = pmtrk.trk[0].slopex;
  trk.slopey  = vtrk.slope;

  float trkang = 180.-180./PI
    *atan(sqrt( pow(tan(vtrk.ang*PI/180.),2)
	       +pow(tan(pmtrk.trk[0].thetax*PI/180.),2)));

  trk.angle         = trkang;
  trk.vetowtracking = false;
  trk.edgewtracking = false;
  trk.ing_trk       = false;
  trk.pm_stop       = vtrk.stop;
  trk.iron_pene     = 0;
  trk.iron_range    = 0;
  trk.sci_range     = (vtrk.fpln-vtrk.ipln)/cos(trkang*PI/180.);

  float totalpe = 0;
  int totalhit  = 0;
  int trkpdg[7];
  memset(trkpdg,0,sizeof(trkpdg));
  fAddPE(vtrk.hit,trkang,totalpe,totalhit,trkpdg,vtrk.fpln-1,
      pmtrk.trk[0].intcptx,pmtrk.trk[0].slopex);
  if(totalhit>0) trk.trkpe = totalpe/totalhit;
  else           trk.trkpe = 0;
  trk.pdg=num2pdg(trkpdg);

  initMuCL();
  addMuCL(vtrk.hit,trkang,vtrk.fpln-1,pmtrk.trk[0].intcptx,pmtrk.trk[0].slopex);
  trk.mucl=calcMuCL();
};



//====================================================
//==========     Main Analysis Function     ==========
//====================================================

bool fAnalyze_new(){

  AnaTrack track;
  Trk      trk;
  analyzed_trk.clear();
  bool vertical;

  //if(hpmtrack.size()==0||vpmtrack.size()==0) return false;

  fIngSortTrack(hingtrack);
  fIngSortTrack(vingtrack);
  fWMSortTrack (hwmtrack );
  fWMSortTrack (vwmtrack );
  fPMSortTrack (hpmtrack );
  fPMSortTrack (vpmtrack );


  if(!NOTJOINT){
    vertical = false;
    fIngWMJoint     (hingtrack,hwmtrack,vertical,WMMODNUM);
    if(PMMODNUM>0){
      fPMWMJoint    (hingtrack,hwmtrack,hpmtrack,vertical,PMMODNUM);
      fIngPMJoint   (hingtrack,hwmtrack,hpmtrack,vertical,PMMODNUM);
    }
    vertical = true;
    fIngWMJoint     (vingtrack,vwmtrack,vertical,WMMODNUM);
    if(PMMODNUM>0){
      fPMWMJoint (vingtrack,vwmtrack,vpmtrack,vertical,PMMODNUM);
      fIngPMJoint   (vingtrack,vwmtrack,vpmtrack,vertical,PMMODNUM);
    }
  }

  vector<int>  id_h, id_v;
  vector<bool> used_h, used_v;
  vector<int>  track_h, track_v;
  vector<bool> tracked_h , tracked_v ; // vertex in PM
  vector<bool> tracked_h2, tracked_v2; // vetex in WM
  vector<bool> tracked_h3, tracked_v3; // vetex in INGRID

  tracked_h.clear();tracked_v.clear();
  //bool ingtrk_tmp = false;
  int hpmtrksize  = (int)hpmtrack .size();
  int vpmtrksize  = (int)vpmtrack .size();
  int hwmtrksize  = (int)hwmtrack .size();
  int vwmtrksize  = (int)vwmtrack .size();
  int hingtrksize = (int)hingtrack.size();
  int vingtrksize = (int)vingtrack.size();

#ifdef DEBUG_THREEDIMRECON
  cout << " hpmtrksize  = " << hpmtrksize   << endl;
  cout << " vpmtrksize  = " << vpmtrksize   << endl;
  cout << " hwmtrksize  = " << hwmtrksize   << endl;
  cout << " vwmtrksize  = " << vwmtrksize   << endl;
  cout << " hingtrksize = " << hingtrksize  << endl;
  cout << " vingtrksize = " << vingtrksize  << endl;
#endif


  for(int i=0; i<hpmtrksize ; i++){ tracked_h .push_back(false); }
  for(int i=0; i<vpmtrksize ; i++){ tracked_v .push_back(false); }
  for(int i=0; i<hwmtrksize ; i++){ tracked_h2.push_back(false); }
  for(int i=0; i<vwmtrksize ; i++){ tracked_v2.push_back(false); }
  for(int i=0; i<hingtrksize; i++){ tracked_h3.push_back(false); }
  for(int i=0; i<vingtrksize; i++){ tracked_v3.push_back(false); }


  // =============================================
  //             VERTEX in Proton Module
  // =============================================


  if(hpmtrksize>0&&vpmtrksize>0){

    // --------------------------------------------------------
    //Looking for tracks stopping in INGRID, with a vertex in PM
    for(int dif=0; dif<diff_th1; dif++){
      for(int pln=0; pln<plnmax(PMMODNUM,0,0,0); pln++){

#ifdef DEBUG_THREEDIMRECON
        cout
          << " ------ Looking for tracks by initial points : Start at PM, Stop at ING------" << endl
          << " dif=" << dif
          << " pln=" << pln
          << endl;
#endif
        //Search tracks with initial points close to each other
        id_h.clear();
        id_v.clear();
        used_h.clear();
        used_v.clear();
        for(int i=0; i<(int)hpmtrack.size(); i++){
          if(tracked_h[i]){ continue; }
          if((hpmtrack[i].ipln-pln)>dif||(hpmtrack[i].ipln-pln)<0){continue;}
          id_h.push_back(i);
          used_h.push_back(false);
#ifdef DEBUG_THREEDIMRECON
          cout 
            << " pmtrk(h):" << i
            << " ipln:"     << hpmtrack[i].ipln
            << " wm_trk:"   << hpmtrack[i].wm_trk
            << " ing_trk:"  << hpmtrack[i].ing_trk;
          if(hpmtrack[i].wm_trk ){
            cout << " wm_num:"   << hpmtrack[i].wm_num;
            cout << " wm_fpln:"  << hpmtrack[i].wm_fpln;
          }
          if(hpmtrack[i].ing_trk){
            cout << " ing_num:"  << hpmtrack[i].ing_num;
            cout << " ing_fpln:" << hpmtrack[i].ing_fpln;
          }
          cout << " id_h.size=" << (int)id_h.size();
          cout << " id_v.size=" << (int)id_v.size();
          cout << endl;
#endif

        }
        for(int j=0; j<(int)vpmtrack.size(); j++){
          if(tracked_v[j]){ continue; }
          if((vpmtrack[j].ipln-pln)>dif||(vpmtrack[j].ipln-pln)<0){continue;}
          id_v.push_back(j);
          used_v.push_back(false);
#ifdef DEBUG_THREEDIMRECON
          cout 
            << " pmtrk(v):" << j
            << " ipln:"     << vpmtrack[j].ipln
            << " wm_trk:"   << vpmtrack[j].wm_trk
            << " ing_trk:"  << vpmtrack[j].ing_trk;
          if(vpmtrack[j].wm_trk ){
            cout << " wm_num:"   << vpmtrack[j].wm_num;
            cout << " wm_fpln:"  << vpmtrack[j].wm_fpln;
          }
          if(vpmtrack[j].ing_trk){
            cout << " ing_num:"  << vpmtrack[j].ing_num;
            cout << " ing_fpln:" << vpmtrack[j].ing_fpln;
          }
          cout << " id_h.size=" << (int)id_h.size();
          cout << " id_v.size=" << (int)id_v.size();
          cout << endl;
#endif
        }

#ifdef DEBUG_THREEDIMRECON
        cout << " id_h.size=" << (int)id_h.size() << endl;
        cout << " id_v.size=" << (int)id_v.size() << endl;
#endif

        //Search tracks with final points in INGRID close to each other
        track_h.clear();track_v.clear();
        for(int ddif=0; ddif<C_INGNumPln-1; ddif++){
          for(int dpln=C_INGNumPln-1; dpln>=0; dpln--){
            for(int i=0; i<(int)id_h.size(); i++){
              if(used_h[i]) continue;
              if((hpmtrack[id_h[i]].ing_fpln-dpln)>ddif||(hpmtrack[id_h[i]].ing_fpln-dpln)<0){continue;}
              for(int j=0; j<(int)id_v.size(); j++){
                if(used_v[j]) continue;
                if((vpmtrack[id_v[j]].ing_fpln-dpln)>ddif||(vpmtrack[id_v[j]].ing_fpln-dpln)<0){continue;}
                if(hpmtrack[id_h[i]].clstime!=vpmtrack[id_v[j]].clstime){continue;}
                //if(used_h[i]) continue;
                //if(used_v[j]) continue;
                track_h.push_back(id_h[i]);
                track_v.push_back(id_v[j]);
                used_h[i]=true;
                used_v[j]=true;
#ifdef DEBUG_THREEDIMRECON
                cout <<"   OK    "<<endl;
                cout 
                  << " hpmtrack:" << id_h[i]
                  << " ing_fpln:" << hpmtrack[id_h[i]].ing_fpln
                  << " dpln:" << dpln
                  << " ddif:" << ddif
                  << endl;
                cout 
                  << " vpmtrack:" << id_v[j]
                  << " ing_fpln:" << vpmtrack[id_v[j]].ing_fpln
                  << " dpln:" << dpln
                  << " ddif:" << ddif
                  << endl;
                cout << " clstime (" << i<<","<<j<<")"
                 << " h=" << hpmtrack[id_h[i]].clstime
                 << " v=" << vpmtrack[id_v[j]].clstime
                 << endl;
#endif

              }
            }
          }
        }

        for(int k=0; k<(int)track_h.size(); k++){
          int h = track_h[k];
          int v = track_v[k];
          track.clear();
          tracked_h [h]=true;
          tracked_v [v]=true;
          if(hpmtrack[h].wm_trk )tracked_h2[hpmtrack[h].wm_num ]=true;
          if(vpmtrack[v].wm_trk )tracked_v2[vpmtrack[v].wm_num ]=true;
          if(hpmtrack[h].ing_trk)tracked_h3[hpmtrack[h].ing_num]=true;
          if(vpmtrack[v].ing_trk)tracked_v3[vpmtrack[v].ing_num]=true;
#ifdef DEBUG_THREEDIMRECON
          cout  << " Vertex:PM Stop:ING";
          cout  << " pmtrk=" << h                  ;
          cout  << " pmtrk=" << v                  ;
          if(hpmtrack[h].wm_trk)cout  << " wmtrk=" << hpmtrack[h].wm_num ;
          if(vpmtrack[v].wm_trk)cout  << " wmtrk=" << vpmtrack[v].wm_num ;
          if(hpmtrack[h].wm_trk)cout  << " ingtrk="<< hpmtrack[h].ing_num;
          if(vpmtrack[v].wm_trk)cout  << " ingtrk="<< vpmtrack[v].ing_num;
          cout  << endl;
#endif
          //Pushing the track info of both H/V tracks, including Muon CL.
          trk.clear();
          fTrackMatch(trk,hpmtrack[h],vpmtrack[v]);
          trk.hnum=h;
          trk.vnum=v;
          track.trk.push_back(trk);
          track.clstime       = (hpmtrack[h].clstime+vpmtrack[v].clstime)/2;
          track.vetowtracking = hpmtrack[h].veto||vpmtrack[v].veto;
          track.edgewtracking = hpmtrack[h].edge||vpmtrack[v].edge;
          track.Ntrack    = 1;
          track.Ningtrack = 1;
          //ingtrk_tmp = true;	
          analyzed_trk.push_back(track);	
        }
      }
    }

    // --------------------------------------------------------
    // Looking for tracks stopping in WAGASCI, with a vertex in PM
    for(int dif=0; dif<diff_th1; dif++){
      for(int pln=0; pln<plnmax(PMMODNUM,0,0,0); pln++){

        //Search tracks with initial points close to each other
        id_h.clear();
        id_v.clear();
        used_h.clear();
        used_v.clear();
        for(int i=0; i<(int)hpmtrack.size(); i++){
          if(tracked_h[i]){ continue;}
          if((hpmtrack[i].ipln-pln)>dif||(hpmtrack[i].ipln-pln)<0){ continue;}
          id_h.push_back(i);
          used_h.push_back(false);
        }
        for(int j=0; j<(int)vpmtrack.size(); j++){
          if(tracked_v[j]) continue;
          if((vpmtrack[j].ipln-pln)>dif||(vpmtrack[j].ipln-pln)<0){ continue;}
          id_v.push_back(j);
          used_v.push_back(false);
        }

        //Search tracks with final points in WaterModule close to each other
        track_h.clear();
        track_v.clear();
        for(int ddif=0; ddif<C_WMNumPln*3-1; ddif++){
          for(int dpln=C_WMNumPln*3-1;dpln>=0;dpln--){
            for(int i=0; i<(int)id_h.size(); i++){
              if((hpmtrack[id_h[i]].wm_fpln-dpln)>ddif||(hpmtrack[id_h[i]].wm_fpln-dpln)<0){continue;}
              for(int j=0; j<(int)id_v.size(); j++){
                if(hpmtrack[id_h[i]].ing_trk&&vpmtrack[id_v[j]].ing_trk){continue;}
                if((vpmtrack[id_v[j]].wm_fpln-dpln)>ddif || (vpmtrack[id_v[j]].wm_fpln-dpln)<0){continue;}
                if(hpmtrack[id_h[i]].clstime!=vpmtrack[id_v[j]].clstime){continue;}
                if(used_h[i]){continue;}
                if(used_v[j]){continue;}
                track_h.push_back(id_h[i]);
                track_v.push_back(id_v[j]);
                used_h[i] = true;
                used_v[j] = true;
              }
            }
          }
        }

        for(int k=0; k<(int)track_h.size(); k++){
          int h = track_h[k];
          int v = track_v[k];
          track.clear();
          tracked_h[h] = true;
          tracked_v[v] = true;
          tracked_h2[hpmtrack[h].wm_num]=true;
          tracked_v2[vpmtrack[v].wm_num]=true;

#ifdef DEBUG_THREEDIMRECON
          cout
            << " Vertex:PM Stop:WM"
            << " pmtrk="<< h                  
            << " pmtrk="<< v                  
            << " wmtrk="<< hpmtrack[h].wm_num 
            << " wmtrk="<< vpmtrack[v].wm_num 
            << endl;
#endif

          //Pushing the track info of both H/V tracks, including Muon CL.
          trk.clear();
          fTrackMatch(trk,hpmtrack[h],vpmtrack[v]);
          trk.hnum=h;
          trk.vnum=v;
          track.trk.push_back(trk);
          track.clstime=(hpmtrack[h].clstime+vpmtrack[v].clstime)/2;
          track.vetowtracking=hpmtrack[h].veto||vpmtrack[v].veto;
          track.edgewtracking=hpmtrack[h].edge||vpmtrack[v].edge;
          //if(ingtrk_tmp) track.Ntrack += 1;
          //else           track.Ntrack  = 1;
          track.Ntrack    = 1;
          track.Ningtrack = 0;
          analyzed_trk.push_back(track);	
        }
      }
    }

    // --------------------------------------------------------
    // Looking for tracks stopping in PM, with a vertex in PM
    for(int dif=0; dif<diff_th2; dif++){
      for(int pln=0; pln<plnmax(PMMODNUM,0,0,0); pln++){

        //Search tracks with initial points close to each other
        id_h.clear();
        id_v.clear();
        used_h.clear();
        used_v.clear();
        for(int i=0; i<(int)hpmtrack.size(); i++){
          if(tracked_h[i]){ continue;}
          if((hpmtrack[i].ipln-pln)>dif||(hpmtrack[i].ipln-pln)<0){ continue;}
          id_h.push_back(i);
          used_h.push_back(false);
        }
        for(int j=0; j<(int)vpmtrack.size(); j++){
          if(tracked_v[j]) continue;
          if((vpmtrack[j].ipln-pln)>dif||(vpmtrack[j].ipln-pln)<0){ continue;}
          id_v.push_back(j);
          used_v.push_back(false);
        }

        //Search tracks with final points in PM close to each other
        track_h.clear();
        track_v.clear();
        for(int ddif=0; ddif<C_PMNumPln-1; ddif++){
          for(int dpln=C_PMNumPln-1;dpln>=0;dpln--){
            for(int i=0; i<(int)id_h.size(); i++){
              if((hpmtrack[id_h[i]].fpln-dpln)>ddif||(hpmtrack[id_h[i]].fpln-dpln)<0){continue;}
              for(int j=0; j<(int)id_v.size(); j++){
                if(hpmtrack[id_h[i]].ing_trk&&vpmtrack[id_v[j]].ing_trk){continue;}
                if(hpmtrack[id_h[i]].wm_trk &&vpmtrack[id_v[j]].wm_trk ){continue;}
                if((vpmtrack[id_v[j]].fpln-dpln)>ddif||(vpmtrack[id_v[j]].fpln-dpln)<0){continue;}
                if(hpmtrack[id_h[i]].clstime!=vpmtrack[id_v[j]].clstime){continue;}
                if(used_h[i]){continue;}
                if(used_v[j]){continue;}
                track_h.push_back(id_h[i]);
                track_v.push_back(id_v[j]);
                used_h[i] = true;
                used_v[j] = true;
              }
            }
          }
        }

        for(int k=0; k<(int)track_h.size(); k++){
          int h = track_h[k];
          int v = track_v[k];
          track.clear();
          tracked_h[h] = true;
          tracked_v[v] = true;

#ifdef DEBUG_THREEDIMRECON
          cout 
            << " Vertex:PM Stop:PM"
            << " pmtrk="<< h                  
            << " pmtrk="<< v                  
            << endl;
#endif

          //Pushing the track info of both H/V tracks, including Muon CL.
          trk.clear();
          fTrackMatch(trk,hpmtrack[h],vpmtrack[v]);
          trk.hnum=h;
          trk.vnum=v;
          track.trk.push_back(trk);
          track.clstime=(hpmtrack[h].clstime+vpmtrack[v].clstime)/2;
          track.vetowtracking=hpmtrack[h].veto||vpmtrack[v].veto;
          track.edgewtracking=hpmtrack[h].edge||vpmtrack[v].edge;
          //if(ingtrk_tmp) track.Ntrack += 1;
          //else           track.Ntrack  = 1;
          track.Ntrack    = 1;
          track.Ningtrack = 0;
          analyzed_trk.push_back(track);	
        }
      }
    }
  }


  // =============================================
  //             VERTEX in WAGASCI
  // =============================================

  if(hwmtrksize>0&&vwmtrksize>0){

    // --------------------------------------------------------
    //Looking for tracks stopping in INGRID, with a vertex in WM
    for(int dif=0; dif<diff_th2; dif++){
      for(int pln=0; pln<plnmax(WMMODNUM,0,0,0)-1; pln++){

        //Search tracks with initial points close to each other
        id_h.clear();
        id_v.clear();
        used_h.clear();
        used_v.clear();
        for(int i=0; i<(int)hwmtrack.size(); i++){
          if(tracked_h2[i]){ continue; }
          if((hwmtrack[i].ipln-pln)>dif||(hwmtrack[i].ipln-pln)<0){continue;}
          id_h.push_back(i);
          used_h.push_back(false);
        }
        for(int j=0; j<(int)vwmtrack.size(); j++){
          if(tracked_v2[j]){ continue; }
          if((vwmtrack[j].ipln-pln)>dif||(vwmtrack[j].ipln-pln)<0){continue;}
          id_v.push_back(j);
          used_v.push_back(false);
        }

        //Search tracks with final points in INGRID close to each other
        track_h.clear();track_v.clear();
        for(int ddif=0; ddif<C_INGNumPln-1; ddif++){
          for(int dpln=C_INGNumPln-1; dpln>=0; dpln--){
            for(int i=0; i<(int)id_h.size(); i++){
              if((hwmtrack[id_h[i]].ing_fpln-dpln)>ddif||(hwmtrack[id_h[i]].ing_fpln-dpln)<0){continue;}
              for(int j=0; j<(int)id_v.size(); j++){
                if((vwmtrack[id_v[j]].ing_fpln-dpln)>ddif||(vwmtrack[id_v[j]].ing_fpln-dpln)<0){continue;}
                if(used_h[i]) continue;
                if(used_v[j]) continue;
                track_h.push_back(id_h[i]);
                track_v.push_back(id_v[j]);
                used_h[i]=true;
                used_v[j]=true;
              }
            }
          }
        }

        for(int k=0; k<(int)track_h.size(); k++){
          int h = track_h[k];
          int v = track_v[k];
          track.clear();
          tracked_h2[h]=true;
          tracked_v2[v]=true;
          tracked_h3[hwmtrack[h].ing_num]=true;
          tracked_v3[vwmtrack[v].ing_num]=true;

          //Pushing the track info of both H/V tracks, including Muon CL.
          trk.clear();
          fTrackMatch(trk,hwmtrack[h],vwmtrack[v]);
          trk.hnum=h;
          trk.vnum=v;
          track.trk.push_back(trk);
          track.clstime = 
            (hingtrack[hwmtrack[h].ing_num].clstime+vingtrack[vwmtrack[v].ing_num].clstime)/2;
          track.vetowtracking = hwmtrack[h].veto||vwmtrack[v].veto;
          track.edgewtracking = hwmtrack[h].edge||vwmtrack[v].edge;
          track.Ntrack    = 1;
          track.Ningtrack = 1;
          //ingtrk_tmp = true;	
          analyzed_trk.push_back(track);	
        }
      }
    }

    // --------------------------------------------------------
    // Looking for tracks stopping in WAGASCI, with a vertex in WM
    for(int dif=0; dif<diff_th2; dif++){
      for(int pln=0; pln<plnmax(WMMODNUM,0,0,0)-1; pln++){

        //Search tracks with initial points close to each other
        id_h.clear();
        id_v.clear();
        used_h.clear();
        used_v.clear();
        for(int i=0; i<(int)hwmtrack.size(); i++){
          if(tracked_h2[i]){ continue;}
          if((hwmtrack[i].ipln-pln)>dif||(hwmtrack[i].ipln-pln)<0){ continue;}
          id_h.push_back(i);
          used_h.push_back(false);
        }
        for(int j=0; j<(int)vwmtrack.size(); j++){
          if(tracked_v2[j]) continue;
          if((vwmtrack[j].ipln-pln)>dif||(vwmtrack[j].ipln-pln)<0){ continue;}
          id_v.push_back(j);
          used_v.push_back(false);
        }

        //Search tracks with final points in WaterModule close to each other
        track_h.clear();
        track_v.clear();
        for(int ddif=0; ddif<C_WMNumPln*3-1; ddif++){
          for(int dpln=C_WMNumPln*3-1;dpln>=0;dpln--){
            for(int i=0; i<(int)id_h.size(); i++){
              if((hwmtrack[id_h[i]].fpln-dpln)>ddif||(hwmtrack[id_h[i]].fpln-dpln)<0){continue;}
              for(int j=0; j<(int)id_v.size(); j++){
                if(hwmtrack[id_h[i]].ing_trk&&vwmtrack[id_v[j]].ing_trk){continue;}
                if((vwmtrack[id_v[j]].fpln-dpln)>ddif || (vwmtrack[id_v[j]].fpln-dpln)<0){continue;}
                //if(hwmtrack[id_h[i]].clstime!=vwmtrack[id_v[j]].clstime){continue;}
                if(used_h[i]){continue;}
                if(used_v[j]){continue;}
                track_h.push_back(id_h[i]);
                track_v.push_back(id_v[j]);
                used_h[i] = true;
                used_v[j] = true;
              }
            }
          }
        }

        for(int k=0; k<(int)track_h.size(); k++){
          int h = track_h[k];
          int v = track_v[k];
          track.clear();
          tracked_h2[h] = true;
          tracked_v2[v] = true;

          //Pushing the track info of both H/V tracks, including Muon CL.
          trk.clear();
          fTrackMatch(trk,hwmtrack[h],vwmtrack[v]);
          trk.hnum=h;
          trk.vnum=v;
          track.trk.push_back(trk);
          track.clstime=(hwmtrack[h].clstime+vwmtrack[v].clstime)/2;
          track.vetowtracking=hwmtrack[h].veto||vwmtrack[v].veto;
          track.edgewtracking=hwmtrack[h].edge||vwmtrack[v].edge;
          //if(ingtrk_tmp) track.Ntrack += 1;
          //else           track.Ntrack  = 1;
          track.Ntrack    = 1;
          track.Ningtrack = 0;
          analyzed_trk.push_back(track);	
        }
      }
    }
  }

  // =============================================
  //             VERTEX in INGRID
  // =============================================

  if(hingtrksize>0&&vingtrksize>0){

    // --------------------------------------------------------
    //Looking for tracks stopping in INGRID, with a vertex in ING
    for(int dif=0; dif<diff_th0; dif++){
      //for(int pln=0; pln<plnmax(INGMODNUM_mid,0,0,0)-1; pln++){
      for(int pln=0; pln<plnmax(INGMODNUM_mid,0,0,0); pln++){

        //Search tracks with initial points close to each other
        id_h.clear();
        id_v.clear();
        used_h.clear();
        used_v.clear();
        for(int i=0; i<(int)hingtrack.size(); i++){
          if(tracked_h3[i]){ continue; }
          if((hingtrack[i].ipln-pln)>dif||(hingtrack[i].ipln-pln)<0){continue;}
          id_h.push_back(i);
          used_h.push_back(false);
        }
        for(int j=0; j<(int)vingtrack.size(); j++){
          if(tracked_v3[j]){ continue; }
          if((vingtrack[j].ipln-pln)>dif||(vingtrack[j].ipln-pln)<0){continue;}
          id_v.push_back(j);
          used_v.push_back(false);
        }

        //Search tracks with final points in INGRID close to each other
        track_h.clear();track_v.clear();
        for(int ddif=0; ddif<C_INGNumPln-1; ddif++){
          for(int dpln=C_INGNumPln-1; dpln>=2; dpln--){
            for(int i=0; i<(int)id_h.size(); i++){
              if((hingtrack[id_h[i]].fpln-dpln)>ddif||(hingtrack[id_h[i]].fpln-dpln)<0){continue;}
              for(int j=0; j<(int)id_v.size(); j++){
                if((vingtrack[id_v[j]].fpln-dpln)>ddif||(vingtrack[id_v[j]].fpln-dpln)<0){continue;}
                if(hingtrack[id_h[i]].clstime!=vingtrack[id_v[j]].clstime){continue;}
                if(used_h[i]) continue;
                if(used_v[j]) continue;
                track_h.push_back(id_h[i]);
                track_v.push_back(id_v[j]);
                used_h[i]=true;
                used_v[j]=true;
              }
            }
          }
        }

        for(int k=0; k<(int)track_h.size(); k++){
          int h = track_h[k];
          int v = track_v[k];
          track.clear();
          tracked_h3[h]=true;
          tracked_v3[v]=true;

          //Pushing the track info of both H/V tracks, including Muon CL.
          trk.clear();
          fTrackMatch(trk,hingtrack[h],vingtrack[v]);
          trk.hnum=h;
          trk.vnum=v;
          track.trk.push_back(trk);
          track.clstime       = (hingtrack[h].clstime+vingtrack[v].clstime)/2;
          track.vetowtracking = hingtrack[h].veto||vingtrack[v].veto;
          track.edgewtracking = hingtrack[h].edge||vingtrack[v].edge;
          track.Ntrack    = 1;
          track.Ningtrack = 1;
          //ingtrk_tmp = true;	
          analyzed_trk.push_back(track);	
        }
      }
    }
  }


  // =============================================
  //             Search Vertex
  // =============================================
#ifdef DEBUG_THREEDIMRECON
  cout << "----- Searching Vertex -----" << endl;
  cout << "   num track = " << analyzed_trk.size() << endl;
#endif 
    

  for(int i=0; i<(int)analyzed_trk.size(); i++){
    for(int j=i+1; j<(int)analyzed_trk.size(); j++){

      if(analyzed_trk[i].Ntrack==0 || analyzed_trk[j].Ntrack==0){continue;}

      int vertexmod = analyzed_trk[i].trk[0].startmod;
      int diff_mod  = analyzed_trk[i].trk[0].startmod - analyzed_trk[j].trk[0].startmod;
      int diff_pln  = 
          abs((analyzed_trk[i].trk[0].startxpln)-(analyzed_trk[j].trk[0].startxpln))
         +abs((analyzed_trk[i].trk[0].startypln)-(analyzed_trk[j].trk[0].startypln));
      double diff_ch = 
         fabs((analyzed_trk[i].trk[0].y)-(analyzed_trk[j].trk[0].y))
        +fabs((analyzed_trk[i].trk[0].x)-(analyzed_trk[j].trk[0].x));

#ifdef DEBUG_THREEDIMRECON
      cout 
        << " vertexmod=" << vertexmod
        << " track: " << i << "," << j
        << " diff_mod:" << diff_mod
        << " diff_pln:" << diff_pln
        << " diff_ch:"  << diff_ch
        << endl;
#endif

      int pln_th;
      double ch_th;
      if     (vertexmod==WMMODNUM     ){pln_th=pln_th2;ch_th=ch_th2; }
      else if(vertexmod==PMMODNUM     ){pln_th=pln_th1;ch_th=ch_th1; }
      else if(vertexmod==INGMODNUM_mid){pln_th=pln_th0;ch_th=ch_th0; }
      else{continue;}

#ifdef DEBUG_THREEDIMRECON
      cout
        << " pln_th:" << pln_th
        << " ch_th:"  << ch_th
        << endl;
#endif

      if(diff_mod!=0      ){continue;}
      if(diff_pln>= pln_th){continue;}
      if(diff_ch >= ch_th ){continue;}

      bool former = false;

      int i_stopmod = analyzed_trk[i].trk[0].stopmodx; //this should be same as stopmody here.
      int j_stopmod = analyzed_trk[j].trk[0].stopmodx;
      if     (i_stopmod==INGMODNUM_mid&&j_stopmod!=INGMODNUM_mid){ former=true; }
      else if(i_stopmod!=INGMODNUM_mid&&j_stopmod==INGMODNUM_mid){ former=false;}
      else if(i_stopmod==WMMODNUM     &&j_stopmod==PMMODNUM     ){ former=true; }
      else if(i_stopmod==PMMODNUM     &&j_stopmod==WMMODNUM     ){ former=true; }
      else{ //stopmod is the same as each other
        int i_endpln = analyzed_trk[i].trk[0].endxpln+analyzed_trk[i].trk[0].endypln;
        int j_endpln = analyzed_trk[j].trk[0].endxpln+analyzed_trk[j].trk[0].endypln;
        if     (i_endpln>j_endpln){ former=true; } 
        else if(i_endpln<j_endpln){ former=false;}
        else{
          bool i_vetoedge = (analyzed_trk[i].vetowtracking||analyzed_trk[i].edgewtracking);
          bool j_vetoedge = (analyzed_trk[j].vetowtracking||analyzed_trk[j].edgewtracking);
          if     (!i_vetoedge && j_vetoedge){ former=true; }
          else if( i_vetoedge &&!j_vetoedge){ former=false;} 
          else{
            former = true;
          }
        }
      }

#ifdef DEBUG_THREEDIMRECON
      cout << "  >> Found: " << i << " " << j  << " former:" << former << endl;
#endif
      int ii,jj;
      if(former){ii=i;jj=j;}else{ii=j;jj=i;}

      analyzed_trk[ii].vetowtracking = analyzed_trk[ii].vetowtracking;
      analyzed_trk[ii].edgewtracking = analyzed_trk[ii].edgewtracking;
      analyzed_trk[ii].Ntrack += analyzed_trk[jj].Ntrack;
      analyzed_trk[jj].Ntrack = 0;

      if(analyzed_trk[jj].Ningtrack>0){
        analyzed_trk[ii].Ningtrack += analyzed_trk[jj].Ningtrack;
      }
      analyzed_trk[jj].Ningtrack =0;

      for(int t=0; t<(int)analyzed_trk[jj].trk.size(); t++){
        analyzed_trk[ii].trk.push_back(analyzed_trk[jj].trk[t]);
      }
      analyzed_trk[jj].trk.clear();
    }
  }



  int maxpdif; 
  if(false){
    // ------------------------------------------------------------
    // Track matching, for non-jointed track around the vertex in PM
    maxpdif = 10;
    for(int k=0; k<(int)analyzed_trk.size(); k++){
      if(analyzed_trk[k].Ntrack==0) continue;
      if(analyzed_trk[k].trk[0].startmod!=PMMODNUM){continue;}

      for(int h=0; h<(int)hpmtrack.size(); h++){
        if(tracked_h[h]){ continue; }
        for(int v=0; v<(int)vpmtrack.size(); v++){
          if(tracked_v[v]){ continue; }
          for(int pdif=0; pdif<maxpdif; pdif++){
            if((hpmtrack[h].fpln-vpmtrack[v].fpln>pdif+1)||(vpmtrack[v].fpln-hpmtrack[h].fpln>pdif)){ continue;}
            if( abs((analyzed_trk[k].trk[0].startxpln)-(hpmtrack[h].ipln))
                +abs((analyzed_trk[k].trk[0].startypln)-(vpmtrack[v].ipln)) 
                >= pln_th1 ){ continue;}
            if( fabs((analyzed_trk[k].trk[0].y)-(vpmtrack[v].ixy))
                +fabs((analyzed_trk[k].trk[0].x)-(hpmtrack[h].ixy))
                >= ch_th1 ){ continue;}
            tracked_h[h] = true;
            tracked_v[v] = true;
            trk.clear();
            fTrackMatch(trk,hpmtrack[h],vpmtrack[v]);
            trk.hnum = h;
            trk.vnum = v;
            analyzed_trk[k].trk.push_back(trk);
            analyzed_trk[k].vetowtracking = analyzed_trk[k].vetowtracking;
            analyzed_trk[k].edgewtracking = analyzed_trk[k].edgewtracking;
            if(hpmtrack[h].ing_trk && vpmtrack[v].ing_trk){
              analyzed_trk[k].Ningtrack++;
            }
            analyzed_trk[k].Ntrack++;
          }
        }
      }

      //Track matching (in X horizontal view)
      for(int h=0; h<(int)hpmtrack.size(); h++){
        if(tracked_h[h]){continue;}
        if(abs((analyzed_trk[k].trk[0].startxpln)-(hpmtrack[h].ipln))>=pln_th1){continue;}
        if(fabs((analyzed_trk[k].trk[0].x)-(hpmtrack[h].ixy))        >=ch_th1 ){continue;}
        tracked_h[h] = true;
        trk.clear();
        fTrackMatchX(trk,analyzed_trk[k],hpmtrack[h]);
        trk.hnum = h;
        trk.vnum = analyzed_trk[k].trk[0].vnum;
        analyzed_trk[k].trk.push_back(trk);
        analyzed_trk[k].vetowtracking = analyzed_trk[k].vetowtracking;
        analyzed_trk[k].edgewtracking = analyzed_trk[k].edgewtracking;
        analyzed_trk[k].Ntrack++;
        analyzed_trk[k].trk[0].mucl = recalcMuCL(hpmtrack[analyzed_trk[k].trk[0].hnum],
            vpmtrack[analyzed_trk[k].trk[0].vnum],
            0,
            trk.endxpln);
      }

      //Track matching (in Y vertical view)
      for(int v=0; v<(int)vpmtrack.size(); v++){
        if(tracked_v[v]) continue;
        if(abs((analyzed_trk[k].trk[0].startypln)-(vpmtrack[v].ipln))>=pln_th1){continue;}
        if(fabs((analyzed_trk[k].trk[0].y)-(vpmtrack[v].ixy))        >=ch_th1 ){continue;}
        tracked_v[v] = true;
        trk.clear();
        fTrackMatchY(trk,analyzed_trk[k],vpmtrack[v]);
        trk.hnum = analyzed_trk[k].trk[0].hnum;
        trk.vnum = v;
        analyzed_trk[k].trk.push_back(trk);
        analyzed_trk[k].vetowtracking = analyzed_trk[k].vetowtracking;
        analyzed_trk[k].edgewtracking = analyzed_trk[k].edgewtracking;
        analyzed_trk[k].Ntrack++;
        analyzed_trk[k].trk[0].mucl = recalcMuCL(hpmtrack[analyzed_trk[k].trk[0].hnum],
            vpmtrack[analyzed_trk[k].trk[0].vnum],
            trk.endypln,
            0);
      }
    }
  }

  if(USE_NON_3DTRK){
    // ------------------------------------------------------------
    // Track matching, for non-jointed track around the vertex in WM
    maxpdif = 6;
    for(int k=0; k<(int)analyzed_trk.size(); k++){
      if(analyzed_trk[k].Ntrack==0) continue;
      if(analyzed_trk[k].trk[0].startmod!=WMMODNUM){continue;}

      for(int h=0; h<(int)hwmtrack.size(); h++){
        if(tracked_h2[h]){ continue; }
        for(int v=0; v<(int)vwmtrack.size(); v++){
          if(tracked_v2[v]){ continue; }
          for(int pdif=0; pdif<maxpdif; pdif++){
            if((hwmtrack[h].fpln-vwmtrack[v].fpln>pdif+1)||(vwmtrack[v].fpln-hwmtrack[h].fpln>pdif)){ continue;}
            if( abs((analyzed_trk[k].trk[0].startxpln)-(hwmtrack[h].ipln))
                +abs((analyzed_trk[k].trk[0].startypln)-(vwmtrack[v].ipln)) 
                >= pln_th2 ){ continue;}
            if( fabs((analyzed_trk[k].trk[0].y)-(vwmtrack[v].ixy))
                +fabs((analyzed_trk[k].trk[0].x)-(hwmtrack[h].ixy))
                >= ch_th2 ){ continue;}
            tracked_h2[h] = true;
            tracked_v2[v] = true;
            trk.clear();
            fTrackMatch(trk,hwmtrack[h],vwmtrack[v]);
            trk.hnum = h;
            trk.vnum = v;
            analyzed_trk[k].trk.push_back(trk);
            analyzed_trk[k].vetowtracking = analyzed_trk[k].vetowtracking;
            analyzed_trk[k].edgewtracking = analyzed_trk[k].edgewtracking;
            if(hwmtrack[h].ing_trk && vwmtrack[v].ing_trk){
              analyzed_trk[k].Ningtrack++;
            }
            analyzed_trk[k].Ntrack++;
          }
        }
      }

      //Track matching (in X horizontal view)
      for(int h=0; h<(int)hwmtrack.size(); h++){
        if(tracked_h2[h]){continue;}
        if(abs((analyzed_trk[k].trk[0].startxpln)-(hwmtrack[h].ipln))>=pln_th2){continue;}
        if(fabs((analyzed_trk[k].trk[0].x)-(hwmtrack[h].ixy))        >=ch_th2 ){continue;}

        tracked_h2[h] = true;
        trk.clear();
        fTrackMatchX(trk,analyzed_trk[k],hwmtrack[h]);
        trk.hnum = h;
        trk.vnum = analyzed_trk[k].trk[0].vnum;
        analyzed_trk[k].trk.push_back(trk);
        analyzed_trk[k].vetowtracking = analyzed_trk[k].vetowtracking;
        analyzed_trk[k].edgewtracking = analyzed_trk[k].edgewtracking;
        analyzed_trk[k].Ntrack++;
        analyzed_trk[k].trk[0].mucl = recalcMuCL(hwmtrack[analyzed_trk[k].trk[0].hnum],
            vwmtrack[analyzed_trk[k].trk[0].vnum],
            0,
            trk.endxpln);
      }

      //Track matching (in Y vertical view)
      for(int v=0; v<(int)vwmtrack.size(); v++){
        if(tracked_v2[v]) continue;
        if(abs((analyzed_trk[k].trk[0].startypln)-(vwmtrack[v].ipln))>=pln_th2){continue;}
        if(fabs((analyzed_trk[k].trk[0].y)-(vwmtrack[v].ixy))        >=ch_th2 ){continue;}

        tracked_v2[v] = true;
        trk.clear();
        fTrackMatchY(trk,analyzed_trk[k],vwmtrack[v]);
        trk.hnum = analyzed_trk[k].trk[0].hnum;
        trk.vnum = v;
        analyzed_trk[k].trk.push_back(trk);
        analyzed_trk[k].vetowtracking = analyzed_trk[k].vetowtracking;
        analyzed_trk[k].edgewtracking = analyzed_trk[k].edgewtracking;
        analyzed_trk[k].Ntrack++;
        analyzed_trk[k].trk[0].mucl = recalcMuCL(hwmtrack[analyzed_trk[k].trk[0].hnum],
            vwmtrack[analyzed_trk[k].trk[0].vnum],
            trk.endypln,
            0);
      }
    }
  }
  


  return true;
};


/// ===================================================================================================
/// ===================================================================================================
/// ===================================================================================================
/// ===================================================================================================


#endif
