#include "HitEfficiency_2.hxx"

FileStat_t fs;
void PrintUsage(){
  std::cout << "./HitEfficiency [options]" << std::endl;
  std::cout << "Either a pair of run/srun number or "
    << "a pair of input/output filename must be set." << endl;
  std::cout << "A mode must be selected between sand/cosmic muons." << std::endl;
  std::cout << "-r <runid>         : Set Run ID." << std::endl;
  std::cout << "-s <srunid/acqid>  : Set Sub-Run/Acq ID."  << std::endl;
  std::cout << "-i <filename.root> : Input ROOT file, after 2d reconstruction." << std::endl;
  std::cout << "-m <mode>          : Select mode; " << std::endl;
  std::cout << "                       0->Sand muon, 1->Cosmic muon." << std::endl;
  std::cout << "                       2->MC Sand muon, 3->MC Cosmic muon." << std::endl;
  std::cout << "-o <filename.root> : Outpu ROOT file." << std::endl;
  std::cout << "-w                 : To analyze B2 WAGASCI cosmic trigger data." << endl;
  std::cout << "-b                 : To open all bad channels. Default: masked." << std::endl;
  std::cout << "-t <(int) pe>      : To set threshold for P.E."  << std::endl;
  std::cout << std::endl;
}

int main(int argc, char** argv){

  //------- Arguments -------
  int c = -1;
  std::string ifilename;
  std::string ofilename = "tmp.root";
  int  mode      = -1;
  bool BadCh_Ana = true;
  bool irename = false;
  bool orename = false;
  bool wg_cosmic = false;
  int    THRES_PE_int = 0;
  double THRES_PE = 0.;

  int runid  = -1;
  int srunid = -1;

  while ((c = getopt(argc, argv, "r:s:i:o:m:t:e:wbh")) != -1){
    switch(c){
      case 'r':
        runid = atoi(optarg);
        break;
      case 's':
        srunid = atoi(optarg);
        break;
      case 'i':
        ifilename = optarg;
        irename = true;
        break;
      case 'o':
        ofilename = optarg;
        orename = true;
        break;
      case 'm':
        mode = atoi(optarg);
        break;
      case 't':
        THRES_PE_int = atoi(optarg);
        break;
      case 'b':
        BadCh_Ana = false;
        break;
      case 'w':
        wg_cosmic = true;
        break;
      case 'h':
        PrintUsage();
        exit(0);
      default:
        PrintUsage();
        exit(0);
    }
  }
  THRES_PE = ((double)THRES_PE_int)-0.5;

  if(
      (mode!=0&&mode!=1&&mode!=2&&mode!=3)||
      ((runid<0||srunid<0)&&(!irename||!orename))||
      (mode==0&&wg_cosmic)
      )
  {
    PrintUsage();
    exit(1);
  }

  if(!irename){
    if(mode==0){
      ifilename = Form("%s/ingrid_%08d_%04d_recon.root",
          dst_ingrid,runid,srunid);
    }
    else if(mode==1){
      if(wg_cosmic){
        ifilename = Form("%s/run_%05d_%03d_recon.root",
            cosmic_wagasci,runid,srunid);
      }
      else{
        ifilename = Form("%s/ingrid_%08d_%04d_recon.root",
            cosmic_data,runid,srunid);
      }
    }
    else if(mode==2){
      ifilename = Form("%s/ingrid_%08d_%04d_recon.root",
          mc_sandmu,runid,srunid);
    }
    else if(mode==3){
      ifilename = Form("%s/ingrid_%08d_%04d_recon.root",
          mc_cosmic,runid,srunid);
    }
  }

  if(!orename){
    if(mode==0){
      ofilename = Form("%s/ingrid_%08d_%04d_hiteff_2_%d.root",
          dst_ingrid,runid,srunid,THRES_PE_int);
    }
    else if(mode==1){
      if(wg_cosmic){
        ofilename = Form("%s/run_%05d_%03d_hiteff_2_%d.root",
            cosmic_wagasci,runid,srunid,THRES_PE_int);
      }
      else{
        ofilename = Form("%s/ingrid_%08d_%04d_hiteff_2_%d.root",
            cosmic_data,runid,srunid,THRES_PE_int);
      }
    }
    else if(mode==2){
      ofilename = Form("%s/ingrid_%08d_%04d_hiteff_2_%d.root",
          mc_sandmu,runid,srunid,THRES_PE_int);
    }
    else if(mode==3){
      ofilename = Form("%s/ingrid_%08d_%04d_hiteff_2_%d.root",
          mc_cosmic,runid,srunid,THRES_PE_int);
    }

  }


  if(gSystem->GetPathInfo(ifilename.c_str(),fs)){
    std::cout << "No such a file: "<< ifilename << std::endl;
    exit(1);
  }

  // ----------------------------------------------------
  // ------- Open root file and get tree information -------
  TFile* ifile = new TFile(ifilename.c_str(),"read");
  if(ifile->IsZombie()){
    std::cout << "Cannot open file : " << ifilename << std::endl;
    exit(1);
  }
  std::cout << "Open file : " << ifilename << std::endl;
  TTree   *tree  = (TTree*)ifile ->Get("tree");
  TBranch *evtbr = tree->GetBranch("fDefaultReco.");
  EventSummary* evtsum = new EventSummary();
  evtbr -> SetAddress(&evtsum);
  tree  -> SetBranchAddress("fDefaultReco.",&evtsum);

  int nevt = (int)tree -> GetEntries();
  std::cout << "Total # of events = " << nevt << std::endl; 


  //------- Define histograms ------
  TFile *ofile = new TFile(ofilename.c_str(),"recreate");
  std::cout << "Output file name : " << ofilename << std::endl;
  std::string name, xtitle, modname, axisname, axtitle;

  const int NMOD   =  5;
  const int NVIEW  =  2;
  const int NAXIS  =  2; //0: z-axis, 1: xy-axis
  const int NBLOCK = 20; //equivalent to plane in this analysis
  TH1F *h0 [NMOD][NVIEW]; 
  TH1F *h1 [NMOD][NVIEW][NAXIS][NBLOCK];
  TH1F *h2 [NMOD][NVIEW][NAXIS][NBLOCK];
  TH1F *h3 [NMOD][NVIEW][NAXIS];
  TH1F *h4 [NMOD][NVIEW][NAXIS][NBLOCK];
  TH1F *h5 [NMOD][NVIEW][NAXIS][NBLOCK];
  TH1F *h6 [NMOD][NVIEW][NAXIS];
  TH1F *h7 [NMOD][NVIEW][NAXIS];
  TH1F *h8 [NMOD][NVIEW][NAXIS];
  TH1F *h9 [NMOD][NVIEW][NAXIS];
  TH1F *h10[NMOD][NVIEW][NAXIS];
  TH1F *h11[NMOD][NVIEW][NAXIS];
  TH1F *h12[NMOD][NVIEW][NAXIS];
  for(int i=0;i<NMOD;i++){
    if     (i==0){modname="INGmod3";} //INGRID mod3
    else if(i==1){modname="INGWM"  ;} //INGRID Water Module
    else if(i==2){modname="WM"     ;} //WAGASCI
    else if(i==3){modname="PM"     ;} //Proton Module
    else if(i==4){modname="B2ING"  ;} //B2 INGRID

    for(int j=0;j<NVIEW;j++){

      //Track angle
      name  =  Form("h0_%s_view%d",modname.c_str(),j);
      xtitle="Track angle from the  Z axis [deg]";
      h0[i][j] =  new TH1F(name.c_str(),name.c_str(),180,-90.,90.);
      h0[i][j] -> GetXaxis()->SetTitle(xtitle.c_str());

      for(int k=0;k<NAXIS;k++){
        if     (k==0){ axisname="Z" ;}
        else if(k==1){ axisname="XY";}


        //P.E.
        name  =  Form("h3_%s_view%d_axis%s",modname.c_str(),j,axisname.c_str());
        xtitle="Hit light yield [p.e.]";
        h3[i][j][k] =  new TH1F(name.c_str(),name.c_str(),2000,0.,200.);
        h3[i][j][k] -> GetXaxis()->SetTitle(xtitle.c_str());
        
        //Penetrated Layer, for track angle
        name  =  Form("h6_%s_view%d_axis%s",modname.c_str(),j,axisname.c_str());
        xtitle="Track angle from the  Z axis [deg]";
        h6[i][j][k] =  new TH1F(name.c_str(),name.c_str(),180,-90.,90.);
        h6[i][j][k] -> GetXaxis()->SetTitle(xtitle.c_str());
        //Penetrated Layer, for track angle
        name  =  Form("h7_%s_view%d_axis%s",modname.c_str(),j,axisname.c_str());
        xtitle="Track angle from the  Z axis [deg]";
        h7[i][j][k] =  new TH1F(name.c_str(),name.c_str(),180,-90.,90.);
        h7[i][j][k] -> GetXaxis()->SetTitle(xtitle.c_str());
        //Penetrated Layer, for posistion in scinti
        name  =  Form("h8_%s_view%d_axis%s",modname.c_str(),j,axisname.c_str());
        xtitle="Hit position in scintillator bar [mm]";
        h8[i][j][k] =  new TH1F(name.c_str(),name.c_str(),2000,-1000.,1000.);
        h8[i][j][k] -> GetXaxis()->SetTitle(xtitle.c_str());
        //Penetrated Layer, for posistion in scinti
        name  =  Form("h9_%s_view%d_axis%s",modname.c_str(),j,axisname.c_str());
        xtitle="Hit position in scintillator bar [mm]";
        h9[i][j][k] =  new TH1F(name.c_str(),name.c_str(),2000,-1000.,1000.);
        h9[i][j][k] -> GetXaxis()->SetTitle(xtitle.c_str());


        //Hit PE, for posistion in scinti
        name  =  Form("h10_%s_view%d_axis%s",modname.c_str(),j,axisname.c_str());
        xtitle="Hit position in scintillator bar [mm]";
        h10[i][j][k] =  new TH1F(name.c_str(),name.c_str(),2000,-1000.,1000.);
        h10[i][j][k] -> GetXaxis()->SetTitle(xtitle.c_str());
        //Hit PE/mm, for posistion in scinti
        name  =  Form("h11_%s_view%d_axis%s",modname.c_str(),j,axisname.c_str());
        xtitle="Hit position in scintillator bar [mm]";
        h11[i][j][k] =  new TH1F(name.c_str(),name.c_str(),2000,-1000.,1000.);
        h11[i][j][k] -> GetXaxis()->SetTitle(xtitle.c_str());
        //Num Hit, for posistion in scinti
        name  =  Form("h12_%s_view%d_axis%s",modname.c_str(),j,axisname.c_str());
        xtitle="Hit position in scintillator bar [mm]";
        h12[i][j][k] =  new TH1F(name.c_str(),name.c_str(),2000,-1000.,1000.);
        h12[i][j][k] -> GetXaxis()->SetTitle(xtitle.c_str());

        for(int l=0;l<NBLOCK;l++){

          //Penetrated Layer, for track angle
          name  =  Form("h1_%s_view%d_axis%s_%02d",modname.c_str(),j,axisname.c_str(),l);
          xtitle="Track angle from the  Z axis [deg]";
          h1[i][j][k][l] =  new TH1F(name.c_str(),name.c_str(),180,-90.,90.);
          h1[i][j][k][l] -> GetXaxis()->SetTitle(xtitle.c_str());
          //Penetrated && Hit Layer, for track angle
          name  =  Form("h2_%s_view%d_axis%s_%02d",modname.c_str(),j,axisname.c_str(),l);
          xtitle="Track angle from the  Z axis [deg]";
          h2[i][j][k][l] =  new TH1F(name.c_str(),name.c_str(),180,-90.,90.);
          h2[i][j][k][l] -> GetXaxis()->SetTitle(xtitle.c_str());


          //Penetrated Layer, for position in scinti
          name  =  Form("h4_%s_view%d_axis%s_%02d",modname.c_str(),j,axisname.c_str(),l);
          xtitle="Hit position in scintillator bar [mm]";
          h4[i][j][k][l] =  new TH1F(name.c_str(),name.c_str(),2000,-1000.,1000.);
          h4[i][j][k][l] -> GetXaxis()->SetTitle(xtitle.c_str());
          //Penetrated && Hit Layer, for position in scinti
          name  =  Form("h5_%s_view%d_axis%s_%02d",modname.c_str(),j,axisname.c_str(),l);
          xtitle="Hit position in scintillator bar [mm]";
          h5[i][j][k][l] =  new TH1F(name.c_str(),name.c_str(),2000,-1000.,1000.);
          h5[i][j][k][l] -> GetXaxis()->SetTitle(xtitle.c_str());

        }
      }
    }
  }

  //------- Event loop -------
  detdim   = new DetectorDimension();
  badch  = new INGRID_BadCh_mapping();
  badch -> set_BadCh(BadCh_Ana);
  //bool newcanv = true;
  int start_cyc,stop_cyc;
  if     (mode==0){start_cyc= 4;stop_cyc=11;}
  else if(mode==1){start_cyc=14;stop_cyc=15;}
  else if(mode==2){start_cyc= 4;stop_cyc= 4;}
  else if(mode==3){start_cyc= 4;stop_cyc= 4;}

  //int n_count = 0;
  for(int ievt=0;ievt<nevt;ievt++){

    if(ievt%100==0){
      std::cout << "Event # is " << ievt << std::endl;
    }

    tree  -> GetEntry(ievt);	
    for(int icyc=start_cyc;icyc<=stop_cyc;icyc++){
      for(int imod=0;imod<NMOD;imod++){
        int c_mod = -1;
        double THRES_DIST = -1.; 
        if     (imod==0){c_mod=3 ;THRES_DIST=75.;}//100.;} //THRES_DIST= //INGRID mod3
        else if(imod==1){c_mod=15;THRES_DIST=38.;}// 50.;} //THRES_DIST= //INGRID Water Module
        else if(imod==2){c_mod=21;THRES_DIST=38.;}// 50.;} //THRES_DIST= //WAGASCI
        else if(imod==3){c_mod=16;THRES_DIST=38.;}// 50.;} //THRES_DIST= //Proton Module
        else if(imod==4){c_mod=14;THRES_DIST=75.;}//100.;} //THRES_DIST= //B2 INGRID

  
        int n2drec[2] = {0,0};
        for(int iview=0;iview<2;iview++){
          n2drec[iview] = evtsum->NModTwoDimRecons(c_mod,icyc,iview);
        }
        if(n2drec[0]!=1||n2drec[1]!=1){continue;}

        double slopex,slopey;
        double intcptx,intcpty;
        double anglex,angley;
        double iX,iY,iZ,iZx,iZy;
        double fX,fY,fZ,fZx,fZy;
        int    nhitx,nhity;
        for(int iview=0;iview<2;iview++){
          TwoDimReconSummary* reconsum = evtsum->GetModTwoDimRecon(0,c_mod,icyc,iview);
          double slope  = reconsum->slope  ;
          double angle  = reconsum->angle  ;
          double intcpt = reconsum->intcpt ;
          double startxy= reconsum->startxy;
          double startz = reconsum->startz ;
          double endxy  = reconsum->endxy  ;
          double endz   = reconsum->endz   ;
          int    nhit   = reconsum->Nhits();
          if(iview==0){
            slopey = slope;
            intcpty= intcpt;
            angley = angle;
            nhity  = nhit;
            iY     = startxy;
            iZy    = startz;
            fY     = endxy;
            fZy    = endz;
          }
          else{
            slopex = slope;
            intcptx= intcpt;
            anglex = angle;
            nhitx  = nhit;
            iX     = startxy;
            iZx    = startz;
            fX     = endxy;
            fZx    = endz;
          }
        }
        if(nhitx < 4 ||nhity < 4 ){continue;}
        if(slopex==0.||slopey==0.){continue;}
        //if(fabs(iZx-iZy) >  200. ){continue;}
        //if(fabs(fZx-fZy) >  200. ){continue;}
        if(fabs(anglex)  >  80.  ){continue;}
        if(fabs(angley)  >  80.  ){continue;}

        // Define start/stop position of the track
        if(iZx>iZy){iZ=iZx;} else{iZ=iZy;}
        if(fZx<fZy){fZ=fZx;} else{fZ=fZy;}
        iX = intcptx + slopex*iZ;
        iY = intcpty + slopey*iZ;
        fX = intcptx + slopex*fZ;
        fY = intcpty + slopey*fZ;


        h0[imod][0]->Fill(angley);
        h0[imod][1]->Fill(anglex);


        bool hitpln1[2][NAXIS][NBLOCK];
        bool hitpln2[2][NAXIS][NBLOCK];
        for(int view=0;view<NVIEW;view++){
          for(int axis=0;axis<NAXIS;axis++){
            for(int blc=0;blc<NBLOCK;blc++){
              hitpln1[view][axis][blc]=false;
              hitpln2[view][axis][blc]=false;
            }
          }
        }

        //Get penetrated layers (==blocks)
        for(int view=0;view<NVIEW;view++){
          for(int axis=0;axis<NAXIS;axis++){
            for(int blc=0;blc<NBLOCK;blc++){
              int pln,ch;
              double posx,posy,posz;
              if(axis==0){  //along z
                if(c_mod==21||c_mod==15){
                  if(blc<8){ pln=blc; ch=0; }
                  else     {continue;}
                }
                else if(c_mod==3||c_mod==14){
                  if(blc<11){ pln=blc; ch=0; }
                  else      {continue;}
                }
                else if(c_mod==16){
                  if(blc<18){ pln=blc; ch=0; }
                  else      {continue;}
                }
                else{continue;}
                detdim->GetPosInMod(c_mod,pln,view,ch,&posx,&posy,&posz);
                if(posz>=iZ&&posz<=fZ){
                  hitpln1[view][axis][blc]=true;
                }
              }
              else if(axis==1){  //along x/y
                if(c_mod==MOD_B2_WM||c_mod==MOD_ONAXIS_WM){
                  if(blc<20){ pln = 0; ch  = blc+40; }
                  else      {continue;}
                  detdim->GetPosInMod(c_mod,pln,view,ch,&posx,&posy,&posz);
                  if(
                      (view==0&&((posy>=iY&&posy<=fY)||(posy>=fY&&posy<=iY)))||
                      (view==1&&((posx>=iX&&posx<=fX)||(posx>=fX&&posx<=iX))))
                  {
                    hitpln1[view][axis][blc]=true;
                  }
                }
                else{continue;}
              }
              else{continue;}
            }
          }
        }

        //Check if there is a hit or more around the reconstrcuted track
        int nmodhits = evtsum->NModHits(c_mod,icyc);
        double sciposx[NBLOCK];
        double sciposy[NBLOCK];
        for(int ihits =0; ihits< nmodhits; ihits++){
          HitSummary* hitsum = evtsum->GetModHit(ihits,c_mod,icyc);
          int    hitmod  = hitsum->mod;
          int    hitview = hitsum->view;
          int    hitpln  = hitsum->pln;
          int    hitch   = hitsum->ch;
          double hitpe   = hitsum->pe;
          if(hitpe<THRES_PE){continue;}

          double path = calc_pathlength_wg(slopex,slopey,hitmod,hitview,hitpln,hitch); 
          double hitpemm = 0.;
          if(path>0.) hitpemm = hitpe/path;

          double hitx,hity,hitz,hitxy;
          double intcpt,slope;
          double posx,posy;

          detdim->GetPosInMod(hitmod,hitpln,hitview,hitch,&hitx,&hity,&hitz);
          if(hitview==0){
            hitxy  = hity;
            slope  = slopey;
            intcpt = intcpty;
          }else{
            hitxy  = hitx;
            slope  = slopex;
            intcpt = intcptx;
          }

          double dist = fabs(hitxy-intcpt-slope*hitz)/sqrt(1.+slope*slope);
          if(dist>THRES_DIST){continue;}

          get_hitpos_in_scinti(intcptx,slopex,intcpty,slopey,hitmod,hitview,hitpln,hitch,posx,posy);

          for(int axis=0;axis<NAXIS;axis++){
            if(axis!=0&&hitmod!=21&&hitmod!=15){continue;}
            if(axis!=0&&hitch< 40){continue;}
            if(axis==0&&hitch>=40){continue;}

            int blc;
            if(axis==0){blc=hitpln;}
            else       {blc=hitch%20;}
            if(!hitpln1[hitview][axis][blc]){continue;}
            if( hitpln2[hitview][axis][blc]){continue;}
            hitpln2[hitview][axis][blc]=true;

            h3 [imod][hitview][axis]->Fill(hitpe);
            h10[imod][hitview][axis]->Fill(posx,hitpe  );
            h11[imod][hitview][axis]->Fill(posx,hitpemm);
            h12[imod][hitview][axis]->Fill(posx,1.     );

            sciposx[blc] = posx;
            sciposy[blc] = posy;
          }

        }//ihits

        for(int view=0;view<NVIEW;view++){
          double angle;
          if(view==0){angle=angley;}else{angle=anglex;}
          for(int axis=0;axis<NAXIS;axis++){
            for(int blc=0;blc<NBLOCK;blc++){
              if(hitpln1[view][axis][blc]){h1[imod][view][axis][blc]->Fill(angle);}
              if(hitpln2[view][axis][blc]){h2[imod][view][axis][blc]->Fill(angle);}
              if(hitpln1[view][axis][blc]){h6[imod][view][axis]->Fill(angle);}
              if(hitpln2[view][axis][blc]){h7[imod][view][axis]->Fill(angle);}
              if(hitpln1[view][axis][blc]){h4[imod][view][axis][blc]->Fill(sciposx[blc]);}
              if(hitpln2[view][axis][blc]){h5[imod][view][axis][blc]->Fill(sciposx[blc]);}
              if(hitpln1[view][axis][blc]){h8[imod][view][axis]->Fill(sciposx[blc]);}
              if(hitpln2[view][axis][blc]){h9[imod][view][axis]->Fill(sciposx[blc]);}
            }
          }
        }
      } //for(imod)
    } //icyc
  } //ievt

  ofile->cd();
  for(int i=0;i<NMOD;i++){
    for(int j=0;j<NVIEW;j++){
      h0[i][j]->Write();
      for(int k=0;k<NAXIS;k++){
        h3 [i][j][k]->Write();
        h6 [i][j][k]->Write();
        h7 [i][j][k]->Write();
        h8 [i][j][k]->Write();
        h9 [i][j][k]->Write();
        h10[i][j][k]->Write();
        h11[i][j][k]->Write();
        h12[i][j][k]->Write();
        for(int l=0;l<NBLOCK;l++){
          h1[i][j][k][l]->Write();
          h2[i][j][k][l]->Write();
          h4[i][j][k][l]->Write();
          h5[i][j][k][l]->Write();
        }
      }
    }
  }
  ofile->Write();
  ofile->Close();

  delete badch;
  delete detdim;
  return 0;
}
