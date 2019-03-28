#include "HitEfficiency.hxx"

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
      ofilename = Form("%s/ingrid_%08d_%04d_hiteff%d.root",
          dst_ingrid,runid,srunid,THRES_PE_int);
    }
    else if(mode==1){
      if(wg_cosmic){
        ofilename = Form("%s/run_%05d_%03d_hiteff%d.root",
            cosmic_wagasci,runid,srunid,THRES_PE_int);
      }
      else{
        ofilename = Form("%s/ingrid_%08d_%04d_hiteff%d.root",
            cosmic_data,runid,srunid,THRES_PE_int);
      }
    }
    else if(mode==2){
      ofilename = Form("%s/ingrid_%08d_%04d_hiteff%d.root",
          mc_sandmu,runid,srunid,THRES_PE_int);
    }
    else if(mode==3){
      ofilename = Form("%s/ingrid_%08d_%04d_hiteff%d.root",
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
  const int NAXIS  =  3; //0: z-axis, 1: y-axis, 2: x-axis
  const int NBLOCK = 32;
  TH1F *h0[NMOD][NAXIS]; 
  TH1F *h1[NMOD][NAXIS][NBLOCK];
  TH1F *h2[NMOD][NAXIS][NBLOCK];
  TH1F *h3[NMOD][NAXIS][NBLOCK];
  TH2F *h4[NMOD][NAXIS][NBLOCK];
  TH2F *h5[NMOD][NAXIS][NBLOCK];
  for(int i=0;i<NMOD;i++){
    if     (i==0){modname="INGmod3";} //INGRID mod3
    else if(i==1){modname="INGWM"  ;} //INGRID Water Module
    else if(i==2){modname="WM"     ;} //WAGASCI
    else if(i==3){modname="PM"     ;} //Proton Module
    else if(i==4){modname="B2ING"  ;} //B2 INGRID

    for(int j=0;j<NAXIS;j++){
      if     (j==0){ axisname="Z";xtitle="Track angle from the Z axis [cos]"; }
      else if(j==1){ axisname="Y";xtitle="Track angle from the Y axis [cos]"; }
      else if(j==2){ axisname="X";xtitle="Track angle from the X axis [cos]"; }
      name  =  Form("angle_%s_axis%s",modname.c_str(),axisname.c_str());
      h0[i][j] =  new TH1F(name.c_str(),name.c_str(),100,-1.,1.);
      h0[i][j] -> GetXaxis()->SetTitle(xtitle.c_str());

      for(int k=0;k<NBLOCK;k++){
        name  =  Form("h1_%s_axis%s_%02d",modname.c_str(),axisname.c_str(),k);
        h1[i][j][k] =  new TH1F(name.c_str(),name.c_str(),100,-1.,1.);
        h1[i][j][k] -> GetXaxis()->SetTitle(xtitle.c_str());
        name  =  Form("h2_%s_axis%s_%02d",modname.c_str(),axisname.c_str(),k);
        h2[i][j][k] =  new TH1F(name.c_str(),name.c_str(),100,-1.,1.);
        h2[i][j][k] -> GetXaxis()->SetTitle(xtitle.c_str());
        name  =  Form("h3_%s_axis%s_%02d",modname.c_str(),axisname.c_str(),k);
        h3[i][j][k] =  new TH1F(name.c_str(),name.c_str(),2000,0.,200.);
        h3[i][j][k] -> GetXaxis()->SetTitle(xtitle.c_str());

        name  =  Form("h4_%s_axis%s_%02d",modname.c_str(),axisname.c_str(),k);
        h4[i][j][k] =  new TH2F(name.c_str(),name.c_str(),100,-1.,1.,100,-1.,1.);
        axtitle="Track angle from the X axis [cos]";
        h4[i][j][k] -> GetXaxis()->SetTitle(axtitle.c_str());
        axtitle="Track angle from the Y axis [cos]";
        h4[i][j][k] -> GetYaxis()->SetTitle(axtitle.c_str());
        name  =  Form("h5_%s_axis%s_%02d",modname.c_str(),axisname.c_str(),k);
        h5[i][j][k] =  new TH2F(name.c_str(),name.c_str(),100,-1.,1.,100,-1.,1.);
        axtitle="Track angle from the X axis [cos]";
        h5[i][j][k] -> GetXaxis()->SetTitle(axtitle.c_str());
        axtitle="Track angle from the Y axis [cos]";
        h5[i][j][k] -> GetYaxis()->SetTitle(axtitle.c_str());
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
    //int nhits = evtsum->NHits();
    for(int icyc=start_cyc;icyc<=stop_cyc;icyc++){
      for(int imod=0;imod<NMOD;imod++){
        int c_mod = -1;
        double THRES_DIST = -1.; 
        if     (imod==0){c_mod=3 ;THRES_DIST=75.;} //INGRID mod3
        else if(imod==1){c_mod=15;THRES_DIST=38.;} //INGRID Water Module
        else if(imod==2){c_mod=21;THRES_DIST=38.;} //WAGASCI
        else if(imod==3){c_mod=16;THRES_DIST=38.;} //Proton Module
        else if(imod==4){c_mod=14;THRES_DIST=75.;} //B2 INGRID

  
        int n2drec[2] = {0,0};
        for(int iview=0;iview<2;iview++){
          n2drec[iview] = evtsum->NModTwoDimRecons(c_mod,icyc,iview);
        }
        if(n2drec[0]!=1||n2drec[1]!=1){continue;}

        bool penetrate_all[2] = {false,false};
        for(int iview=0;iview<2;iview++){
          TwoDimReconSummary* reconsum = evtsum->GetModTwoDimRecon(0,c_mod,icyc,iview);
          int   hitmod   = reconsum->hitmod;
          int   startpln = reconsum->startpln;
          int   endpln   = reconsum->endpln;
          float startxy  = reconsum->startxy;
          float endxy    = reconsum->endxy;
          bool isWAGASCI = (hitmod==MOD_ONAXIS_WM||hitmod==MOD_B2_WM);
          bool isPM      = (hitmod==MOD_PM);
          bool isINGRID  = ((hitmod>=0&&hitmod<NUMINGMOD)||hitmod==MOD_B2_INGRID);
          int max_pln = -1;
          int scinti_wdith = -1;
          if     (isINGRID ){max_pln = C_INGNumPln; scinti_wdith = C_INGScintiWidth;}
          else if(isPM     ){max_pln = C_PMNumPln ; scinti_wdith = C_PMScintiWidth ;}
          else if(isWAGASCI){max_pln = C_WMNumPln ; scinti_wdith = C_WMScintiWidth ;}
          else{continue;}
          double surface_xy = fabs(xyposi(hitmod,iview,max_pln-1,0));

          if(
              (startpln>1        &&fabs(fabs(startxy)-surface_xy)>scinti_wdith*3)||
              (endpln  <max_pln-2&&fabs(fabs(endxy)  -surface_xy)>scinti_wdith*3)
            )
          {continue;}
          else{
            penetrate_all[iview] = true;
          }
        } //iview
        //if(!penetrate_all[0]||!penetrate_all[1]){continue;}
        if(!penetrate_all[0]&&!penetrate_all[1]){continue;}

        bool hitpln1[NAXIS][NBLOCK];
        bool hitpln2[NAXIS][NBLOCK];
        for(int i=0;i<NAXIS;i++){
          for(int j=0;j<NBLOCK;j++){
            hitpln1[i][j]=false;
            hitpln2[i][j]=false;
          }
        }
        double slopex=-1.,slopey=-1.;
        double intcptx=-1.,intcpty=-1.;
        for(int iview=0;iview<2;iview++){
          TwoDimReconSummary* reconsum = evtsum->GetModTwoDimRecon(0,c_mod,icyc,iview);
          //int   hitmod   = reconsum->hitmod;
          double slope   = reconsum->slope;
          double intcpt  = reconsum->intcpt;
          if(iview==0){slopey=slope;intcpty=intcpt;}
          else        {slopex=slope;intcptx=intcpt;}
          int   nrechits = reconsum->Nhits();
          for(int i=0;i<nrechits;i++){
            HitSummary* rechit = reconsum->GetHit(i);
            int    hitmod  = rechit->mod;
            int    hitview = rechit->view;
            int    hitpln  = rechit->pln;
            int    hitch   = rechit->ch;
            double hitpe   = rechit->pe;
            //double hitx,hity,hitz,hitxy;
            //detdim->GetPosInMod(hitmod,hitpln,hitview,hitch,&hitx,&hity,&hitz);
            //if(iview==0){hitxy=hity;}else{hitxy=hitx;}
            //cout << fabs(slope*hitz-hitxy+intcpt)/sqrt(1.+slope*slope) << endl;

            if(hitpe<THRES_PE){continue;}
            hitpln1[0][hitpln] = true;
            if(hitmod==MOD_ONAXIS_WM||hitmod==MOD_B2_WM){
              if(hitch<40){ hitch /=  2; }
              else        { hitch %= 20; }
            }
            hitpln1[hitview+1][hitch] = true;
          }
        }//iview

        double cos[NAXIS];
        double norm = sqrt(1.+slopex*slopex+slopey*slopey);
        for(int i=0;i<NAXIS;i++){
          if     (i==0){ cos[i] = 1./norm; }
          else if(i==1){ cos[i] = slopey/norm; }
          else if(i==2){ cos[i] = slopex/norm; }
        }

        int nmodhits = evtsum->NModHits(c_mod,icyc);
        for(int ihits =0; ihits< nmodhits; ihits++){
          HitSummary* hitsum = evtsum->GetModHit(ihits,c_mod,icyc);
          int    hitmod  = hitsum->mod;
          int    hitview = hitsum->view;
          int    hitpln  = hitsum->pln;
          int    hitch   = hitsum->ch;
          double hitpe   = hitsum->pe;
          if(hitpe<THRES_PE){continue;}
          double hitx,hity,hitz;
          detdim->GetPosInMod(hitmod,hitpln,hitview,hitch,&hitx,&hity,&hitz);
          if(hitview==0){
            if(fabs(slopey*hitz-hity+intcpty)/sqrt(1.+slopey*slopey)>THRES_DIST){
              continue;
            }
          }
          else{
            if(fabs(slopex*hitz-hitx+intcptx)/sqrt(1.+slopex*slopex)>THRES_DIST){
              continue;
            }
          }
          hitpln2[0][hitpln] = true;
          if(hitmod==MOD_ONAXIS_WM||hitmod==MOD_B2_WM){
            if(hitch<40){ hitch /=  2; }
            else        { hitch %= 20; }
          }
          hitpln2[hitview+1][hitch] = true;

          h3[imod][0][hitpln]->Fill(hitpe);
          if(hitview==0){
            h3[imod][1][hitch]->Fill(hitpe);
          }
          else{
            h3[imod][2][hitch]->Fill(hitpe);
          }

        }//ihits

        //cout << "DEBUG: "
        //  << cos[0]<<" "<<cos[1]<<" "<<cos[2]<<" "
        //  <<cos[0]*cos[0]+cos[1]*cos[1]+cos[2]*cos[2]
        //  <<endl;

        for(int i=0;i<NAXIS;i++){
          h0[imod][i] -> Fill(cos[i]);
          for(int j=1;j<NBLOCK-1;j++){
            if(hitpln1[i][j-1]&&hitpln1[i][j+1]){
              h1[imod][i][j]->Fill(cos[i]);
              h4[imod][i][j]->Fill(cos[2],cos[1]);
              if(hitpln2[i][j]){
                h2[imod][i][j]->Fill(cos[i]);
                h5[imod][i][j]->Fill(cos[2],cos[1]);
              }
            }
          }
        }
      } //for(imod)
    } //icyc
  } //ievt

  ofile->cd();
  for(int i=0;i<NMOD;i++){
    for(int j=0;j<NAXIS;j++){
      h0[i][j]->Write();
      for(int k=0;k<NBLOCK;k++){
        h1[i][j][k]->Write();
        h2[i][j][k]->Write();
        h3[i][j][k]->Write();
        h4[i][j][k]->Write();
        h5[i][j][k]->Write();
      }
    }
  }
  ofile->Write();
  ofile->Close();

  delete badch;
  delete detdim;
  return 0;
}
