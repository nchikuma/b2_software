#include "B2pe_Scinti.hxx"

FileStat_t fs;
void PrintUsage(){
  std::cout << "./B2pe [options]" << std::endl;
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
  bool BadCh_Ana = false;
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
        BadCh_Ana = true;
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
      ofilename = Form("%s/ingrid_%08d_%04d_pe_s_%d.root",
          dst_ingrid,runid,srunid,THRES_PE_int);
    }
    else if(mode==1){
      if(wg_cosmic){
        ofilename = Form("%s/run_%05d_%03d_pe_s_%d.root",
            cosmic_wagasci,runid,srunid,THRES_PE_int);
      }
      else{
        ofilename = Form("%s/ingrid_%08d_%04d_pe_s_%d.root",
            cosmic_data,runid,srunid,THRES_PE_int);
      }
    }
    else if(mode==2){
      ofilename = Form("%s/ingrid_%08d_%04d_pe_s_%d.root",
          mc_sandmu,runid,srunid,THRES_PE_int);
    }
    else if(mode==3){
      ofilename = Form("%s/ingrid_%08d_%04d_pe_s_%d.root",
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

  cout << "START_UNIXT = " << START_UNIXT << endl;
  cout << "STOP_UNIXT  = " << STOP_UNIXT  << endl;
  cout << "RANGE_UNIXT = " << RANGE_UNIXT << endl;
  cout << "NBIN_UNIXT  = " << NBIN_UNIXT  << endl;


  //------- Define histograms ------
  TFile *ofile = new TFile(ofilename.c_str(),"recreate");
  std::cout << "Output file name : " << ofilename << std::endl;
  std::string name, xtitle, ytitle, modname, type,axisname;

  const int NMOD   =  5;
  const int NAXIS  =  3; //0: z-axis, 1: y-axis, 2: x-axis
  const int NVIEW  =  2;
  const int NPLN   = 18;
  const int NCH    = 80;
  int maxpln[NMOD] = {NPLN,NPLN,NPLN,NPLN,NPLN};
  int maxch [NMOD] = {NCH,NCH,NCH,NCH,NCH};
  //maxpln[0]=11;maxch[0]=24; //INGRID mod3
  //maxpln[1]= 8;maxch[1]=80; //INGRID Water Module
  //maxpln[2]= 8;maxch[2]=80; //WAGASCI
  //maxpln[3]=18;maxch[3]=32; //Proton Module
  //maxpln[4]=11;maxch[4]=24; //B2 INGRID

  TH1F *h1[NMOD][NVIEW][2];
  TH1F *h2[NMOD][NVIEW][2];
  TH1F *h3[NMOD][NVIEW][2];
  TH1F *h4[NMOD][NVIEW][2];
  TH1F *h5[NMOD][NVIEW][2];


  int    binX =  2000;
  int    binY =   100;
  double maxX =  1000.;
  double minX = -1000.;
  double maxY =    50.;
  double minY =   -50.;

  for(int i=0;i<NMOD;i++){
    if     (i==0){modname="INGmod3";} //INGRID mod3
    else if(i==1){modname="INGWM"  ;} //INGRID Water Module
    else if(i==2){modname="WM"     ;} //WAGASCI
    else if(i==3){modname="PM"     ;} //Proton Module
    else if(i==4){modname="B2ING"  ;} //B2 INGRID
    xtitle = "Position in Scinti [mm]";
    ytitle = "Position in Scinti [mm]";

    for(int j=0;j<NVIEW;j++){
      for(int k=0;k<2;k++){
        if(i==0||i==4){
          if(k==0){type="INGRID";}
          else{continue;}
        }else if(i==1||i==2){
          if(k==0){type="XY";}
          else    {type="Grid";}
        }else if(i==3){
          if(k==0){type="INGRID";}
          else    {type="SciBar";}
        }
        name  =  Form("h1_%s_%s_%01d",modname.c_str(),type.c_str(),j);     
        h1[i][j][k] =  new TH1F(name.c_str(),name.c_str(),binX,minX,maxX);         
        h1[i][j][k] -> GetXaxis()->SetTitle(xtitle.c_str());
        h1[i][j][k] -> GetYaxis()->SetTitle("Light Yield [p.e./mm]");
        name  =  Form("h2_%s_%s_%01d",modname.c_str(),type.c_str(),j);     
        h2[i][j][k] =  new TH1F(name.c_str(),name.c_str(),binX,minX,maxX); 
        h2[i][j][k] -> GetXaxis()->SetTitle(ytitle.c_str());
        h2[i][j][k] -> GetYaxis()->SetTitle("Light Yield [p.e.]");

        name  =  Form("h3_%s_%s_%01d",modname.c_str(),type.c_str(),j);     
        h3[i][j][k] =  new TH1F(name.c_str(),name.c_str(),binX,minX,maxX);         
        h3[i][j][k] -> GetXaxis()->SetTitle(xtitle.c_str());
        h3[i][j][k] -> GetYaxis()->SetTitle("Num of hits");

        name  =  Form("h4_%s_%s_%01d",modname.c_str(),type.c_str(),j);     
        h4[i][j][k] =  new TH1F(name.c_str(),name.c_str(),3000,0.,300.); 
        h4[i][j][k] -> GetXaxis()->SetTitle(ytitle.c_str());
        h4[i][j][k] -> GetYaxis()->SetTitle("Light Yield [p.e./mm]");
        name  =  Form("h5_%s_%s_%01d",modname.c_str(),type.c_str(),j);     
        h5[i][j][k] =  new TH1F(name.c_str(),name.c_str(),3000,0.,300.); 
        h5[i][j][k] -> GetXaxis()->SetTitle(ytitle.c_str());
        h5[i][j][k] -> GetYaxis()->SetTitle("Light Yield [p.e.]");

      }
    }
  }

  //------- Event loop -------
  detdim  = new DetectorDimension();
  badch   = new INGRID_BadCh_mapping();
  badch   -> set_BadCh(BadCh_Ana);
  int start_cyc,stop_cyc;
  if     (mode==0){start_cyc= 4;stop_cyc=11;}
  else if(mode==1){start_cyc=14;stop_cyc=15;}
  else if(mode==2){start_cyc= 4;stop_cyc= 4;}
  else if(mode==3){start_cyc= 4;stop_cyc= 4;}

  for(int ievt=0;ievt<nevt;ievt++){

    if(ievt%100==0){
      std::cout << "Event # is " << ievt << std::endl;
    }

    tree  -> GetEntry(ievt);	
    for(int icyc=start_cyc;icyc<=stop_cyc;icyc++){
      for(int imod=0;imod<NMOD;imod++){
        int c_mod = -1;
        if     (imod==0){c_mod=3 ;} //INGRID mod3
        else if(imod==1){c_mod=15;} //INGRID Water Module
        else if(imod==2){c_mod=21;} //WAGASCI
        else if(imod==3){c_mod=16;} //Proton Module
        else if(imod==4){c_mod=14;} //B2 INGRID


        bool nrechit = true;

        int n2drec[2] = {0,0};
        for(int iview=0;iview<2;iview++){
          n2drec[iview] = evtsum->NModTwoDimRecons(c_mod,icyc,iview);
        }
        if(n2drec[0]!=1||n2drec[1]!=1){continue;}

        bool isWAGASCI = (c_mod==MOD_ONAXIS_WM||c_mod==MOD_B2_WM);
        bool isPM      = (c_mod==MOD_PM);
        bool isINGRID  = ((c_mod>=0&&c_mod<NUMINGMOD)||c_mod==MOD_B2_INGRID);
        int max_pln = -1;
        int scinti_width = -1;
        if     (isINGRID ){max_pln = C_INGNumPln; scinti_width = C_INGScintiWidth;}
        else if(isPM     ){max_pln = C_PMNumPln ; scinti_width = C_PMScintiWidth ;}
        else if(isWAGASCI){max_pln = C_WMNumPln ; scinti_width = C_WMScintiWidth ;}
        else{continue;}

        double slopex=-1.,slopey=-1.;
        double intcptx=-1.,intcpty=-1.;
        double anglex=-1.,angley=-1.;
        for(int iview=0;iview<2;iview++){
          TwoDimReconSummary* reconsum = evtsum->GetModTwoDimRecon(0,c_mod,icyc,iview);
          double slope   = reconsum->slope;
          double angle   = reconsum->angle;
          double intcpt  = reconsum->intcpt;
          if(iview==0){slopey=slope;intcpty=intcpt;anglex=angle;}
          else        {slopex=slope;intcptx=intcpt;angley=angle;}

          if(reconsum->Nhits()<6){nrechit = false;}

        }//iview

        if(!nrechit){continue;}
        if(slopex==0.||slopey==0.){continue;}
        //if(fabs(slopex)>100.||fabs(slopey)>100.){continue;}
        if(fabs(anglex)>50.||fabs(angley)>50.){continue;}

        for(int iview=0;iview<2;iview++){
          TwoDimReconSummary* reconsum = evtsum->GetModTwoDimRecon(0,c_mod,icyc,iview);
          int nhit = reconsum->Nhits();
          for(int ihit=0;ihit<nhit;ihit++){
            HitSummary *hitsum = reconsum->GetHit(ihit);
            int    hitmod  = hitsum->mod;
            int    hitview = hitsum->view;
            int    hitpln  = hitsum->pln;
            int    hitch   = hitsum->ch;
            double hitpe   = hitsum->pe;

            if(hitmod!=c_mod){continue;}
            if(badch->is_BadCh(hitmod,hitview,hitpln,hitch)){continue;}
            if(hitpe<THRES_PE){continue;}

            double path = calc_pathlength_wg(slopex,slopey,hitmod,hitview,hitpln,hitch); 
            double hitpemm = 0.;
            if(path>0.) hitpemm = hitpe/path;

            double posx=0.,posy=0.;
            if(!get_hitpos_in_scinti(intcptx,slopex,intcpty,slopey,hitmod,hitview,hitpln,hitch,posx,posy))
            {continue;}

            int itype = -1;
            if(hitmod==15||hitmod==21){
              if(hitch<40){itype=0;}else{itype=1;}
            }else if(hitmod==16){
              if(hitpln==0||hitch<8||hitch>=24){itype=0;}
              else                             {itype=1;}
            }
            else{ itype = 0; }

            h1[imod][hitview][itype] ->Fill(posx,hitpemm);
            h2[imod][hitview][itype] ->Fill(posx,hitpe);
            h3[imod][hitview][itype] ->Fill(posx,1.);

            h4[imod][hitview][itype] ->Fill(hitpemm);
            h5[imod][hitview][itype] ->Fill(hitpe);

          }
        }
      } //for(imod)
    } //icyc
  } //ievt


  ofile->cd();
  for(int i=0;i<NMOD;i++){
    for(int j=0;j<NVIEW;j++){
      for(int k=0;k<2;k++){
        if(i==0||i==4){ if(k==1){continue;} }
        h1[i][j][k] -> Write();
        h2[i][j][k] -> Write();
        h3[i][j][k] -> Write();
        h4[i][j][k] -> Write();
        h5[i][j][k] -> Write();
      }
    }
  }

  ofile->Write();
  ofile->Close();

  delete badch;
  delete detdim;
  return 0;
}
