#include "B2pe.hxx"
#include "B2Plot.hxx"

FileStat_t fs;
void PrintUsage(){
  std::cout << "./B2pe [options]" << std::endl;
  std::cout << "Either a pair of run/srun number or "
    << "a pair of input/output filename must be set." << endl;
  std::cout << "A mode must be selected between sand/cosmic muons." << std::endl;
  std::cout << "-i <filename.root> : Input ROOT file, after 2d reconstruction." << std::endl;
  std::cout << "-o <filename.root> : Outpu ROOT file." << std::endl;
  //std::cout << "-m <mode>          : Select mode; " << std::endl;
  //std::cout << "                       0->Sand muon, 1->Cosmic muon." << std::endl;
  //std::cout << "                       2->MC Sand muon, 3->MC Cosmic muon." << std::endl;
  std::cout << "-m <fdid>            : Select fdid;  0-> data, 7->WM, 8->PM, 9->ING, ... " << std::endl;
  std::cout << "-b                 : To open all bad channels. Default: masked." << std::endl;
  std::cout << "-t <(int) pe>      : To set threshold for P.E. Default: 3pe (>2.5pe)"  << std::endl;
  std::cout << std::endl;
}


int main(int argc, char** argv){

  //------- Arguments -------
  int c = -1;
  std::string ifilename;
  std::string ofilename = "tmp.root";
  int fdid = -1;
  int sim_mod = -1;
  bool BadCh_Ana = false;
  bool irename = false;
  bool orename = false;
  int    THRES_PE_int = 3;
  double THRES_PE = 0.;

  while ((c = getopt(argc, argv, "i:o:m:t:bh")) != -1){
    switch(c){
      case 'i':
        ifilename = optarg;
        irename = true;
        break;
      case 'o':
        ofilename = optarg;
        orename = true;
        break;
      case 'm':
        fdid = atoi(optarg);
        if(fdid==7){sim_mod=2;}
        if(fdid==8){sim_mod=1;}
        if(fdid==9){sim_mod=0;}
        break;
      case 't':
        THRES_PE_int = atoi(optarg);
        break;
      case 'b':
        BadCh_Ana = true;
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

  if( (fdid<0) || (!irename||!orename))
  {
    PrintUsage();
    exit(1);
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

  // ------------------ Flux Tuning File ----------------
  char tunefile[300];
  string tunefilename = "~/b2_data/jnubeam/tunefile/tune_nd7_8_9.root";
  sprintf(tunefile,"%s",tunefilename.c_str());
  FluxTuning* fluxtune;
  if(fdid==7||fdid==8||fdid==9) fluxtune = new FluxTuning(fdid,tunefile);

  //------- Define histograms ------
  TFile *ofile = new TFile(ofilename.c_str(),"recreate");
  std::cout << "Output file name : " << ofilename << std::endl;
  std::string name, xtitle, modname, axisname;

  const int NPEMOD   =  3;
  const int NAXIS  =  3; //0: z-axis, 1: y-axis, 2: x-axis
  const int NVIEW  =  2;
  const int NPLN   = 18;
  const int NCH    = 80;
  int maxpln = NPLN;
  int maxch  = NCH ;
  //maxpln[0]=11;maxch[0]=24; //INGRID mod3
  //maxpln[1]= 8;maxch[1]=80; //INGRID Water Module
  //maxpln[2]= 8;maxch[2]=80; //WAGASCI
  //maxpln[3]=18;maxch[3]=32; //Proton Module
  //maxpln[4]=11;maxch[4]=24; //B2 INGRID

  TH1F *h0[NPEMOD][NAXIS]; 
  TH1F *h1[NPEMOD][NVIEW][NPLN][NCH];
  TH1F *h2[NPEMOD][NVIEW][NPLN][NCH];
  TH1F *h3[NPEMOD]; 
  TH1F *h4[NPEMOD]; 
  TTree *wt;
  int date;
  double mean_pe   [NPEMOD][NVIEW][NPLN][NCH];
  double mean_pemm [NPEMOD][NVIEW][NPLN][NCH];
  double pathlength[NPEMOD][NVIEW][NPLN][NCH];
  int    num_hit   [NPEMOD][NVIEW][NPLN][NCH];
  if(fdid==0)
  {
    wt = new TTree("tree","tree");
    wt->Branch("date"      ,&date,"date/I");
    wt->Branch("mean_pe"   ,mean_pe   ,Form("mean_pe[%d][%d][%d][%d]/D"   ,NPEMOD,NVIEW,NPLN,NCH));
    wt->Branch("mean_pemm" ,mean_pemm ,Form("mean_pemm[%d][%d][%d][%d]/D" ,NPEMOD,NVIEW,NPLN,NCH));
    wt->Branch("pathlength",pathlength,Form("pathlength[%d][%d][%d][%d]/D",NPEMOD,NVIEW,NPLN,NCH));
    wt->Branch("num_hit"   ,num_hit   ,Form("num_hit[%d][%d][%d][%d]/I"   ,NPEMOD,NVIEW,NPLN,NCH));
  }
  for(int i=0;i<NPEMOD;i++){
    //if     (i==0){modname="INGmod3";} //INGRID mod3
    //else if(i==1){modname="INGWM"  ;} //INGRID Water Module
    if     (i==0){modname="WM"     ;} //WAGASCI
    else if(i==1){modname="PM"     ;} //Proton Module
    else if(i==2){modname="B2ING"  ;} //B2 INGRID

    for(int j=0;j<NAXIS;j++){
      if     (j==0){ axisname="Z";xtitle="Track angle from the Z axis [cos]"; }
      else if(j==1){ axisname="Y";xtitle="Track angle from the Y axis [cos]"; }
      else if(j==2){ axisname="X";xtitle="Track angle from the X axis [cos]"; }
      name  =  Form("angle_%s_axis%s",modname.c_str(),axisname.c_str());
      h0[i][j] =  new TH1F(name.c_str(),name.c_str(),100,-1.,1.);
      h0[i][j] -> GetXaxis()->SetTitle(xtitle.c_str());
    }
    for(int j=0;j<NVIEW;j++){
      for(int k=0;k<maxpln;k++){
        for(int l=0;l<maxch;l++){
          xtitle = "Light Yield [p.e.]";
          name  =  Form("h1_%s_%01d%02d%02d",modname.c_str(),j,k,l);
          h1[i][j][k][l] =  new TH1F(name.c_str(),name.c_str(),4000,-10.,300.);
          h1[i][j][k][l] -> GetXaxis()->SetTitle(xtitle.c_str());
          xtitle = "Light Yield [p.e./mm]";
          name  =  Form("h2_%s_%01d%02d%02d",modname.c_str(),j,k,l);
          h2[i][j][k][l] =  new TH1F(name.c_str(),name.c_str(),4000,-10.,300.);
          h2[i][j][k][l] -> GetXaxis()->SetTitle(xtitle.c_str());

          mean_pe   [i][j][k][l] = 0.;
          mean_pemm [i][j][k][l] = 0.;
          num_hit   [i][j][k][l] = 0 ;
          pathlength[i][j][k][l] = 0 ;
        }
      }
    }
    xtitle = "Light Yield [p.e.]";
    name  =  Form("h3_%s",modname.c_str());
    h3[i] =  new TH1F(name.c_str(),name.c_str(),4000,-10.,300.);
    h3[i] -> GetXaxis()->SetTitle(xtitle.c_str());

    xtitle = "Light Yield [p.e./mm]";
    name  =  Form("h4_%s",modname.c_str());
    h4[i] =  new TH1F(name.c_str(),name.c_str(),4000,-10.,300.);
    h4[i] -> GetXaxis()->SetTitle(xtitle.c_str());

  }

  //------- Event loop -------
  detdim   = new DetectorDimension();
  badch  = new INGRID_BadCh_mapping();
  badch -> set_BadCh(BadCh_Ana);
  int start_cyc,stop_cyc;
  if(fdid==0){start_cyc= 4;stop_cyc=11;}
  else       {start_cyc= 4;stop_cyc= 4;}

  int last_date = -1;
  for(int ievt=0;ievt<nevt;ievt++){
    if(ievt%100==0){
      std::cout << "Event # is " << ievt << std::endl;
    }
    tree  -> GetEntry(ievt);	

    double norm_tot=1.,reweight=1.;
    double norm=1.,totcrsne=1.;
    double nuE=-1.;
    int    nutype=-1;

    int nsimver = evtsum->NSimVertexes();
    if(nsimver>0){
      SimVertexSummary* simver = evtsum->GetSimVertex(0);
      norm      = simver->norm;
      totcrsne  = simver->totcrsne;
      nuE       = simver->nuE;
      nutype    = simver->nutype;

      if(fdid==7||fdid==8||fdid==9){
        reweight = fluxtune->reweight(nuE,nutype);
      }
      if(sim_mod>=0){
        norm_tot = norm*totcrsne*1e-38*Avogadro*density[sim_mod]*thickness[sim_mod]*reweight;
      }
      else{
        float thick = 0.;
        if     (fdid==10){thick=thickness_wallbg;}
        else if(fdid==11){thick=thickness_ceiling;}
        else if(fdid==12){thick=thickness_floor;}
        else if(fdid==13){thick=thickness_pillar_r;}
        else if(fdid==14){thick=thickness_pillar_l;}
        norm_tot = norm*totcrsne*1e-38*Avogadro*density_wallbg*thick;
      }
    }


    double evttime = evtsum->time;
    evttime -= START_UNIXT;
    date     = evttime/86400;
    for(int icyc=start_cyc;icyc<=stop_cyc;icyc++){
      for(int imod=0;imod<NPEMOD;imod++){
        int c_mod = -1;
        //if     (imod==0){c_mod=3 ;} //INGRID mod3
        //else if(imod==1){c_mod=15;} //INGRID Water Module
        if     (imod==0){c_mod=21;} //WAGASCI
        else if(imod==1){c_mod=16;} //Proton Module
        else if(imod==2){c_mod=14;} //B2 INGRID

        int cyc_offset[2] = {0,0};

        bool muon_trk=true;
        for(int iview=0;iview<2;iview++){
          if(fdid==0&&c_mod==21){
            int numreco_tmp = 0;
            for(int iicyc=0;iicyc<2;iicyc++){
              numreco_tmp += evtsum->NModTwoDimRecons(c_mod,icyc+iicyc,iview);
            }
            if(numreco_tmp!=1){
              muon_trk=false;
              break;
            }
            int numrecohit = 0;
            for(int iicyc=0;iicyc<2;iicyc++){
              if(evtsum->NModTwoDimRecons(c_mod,icyc+iicyc,iview)==1){
                numrecohit = evtsum->GetModTwoDimRecon(0,c_mod,icyc+iicyc,iview)->Nhits();
                cyc_offset[iview] = iicyc;
              }
            }
            if(numrecohit<5){
              muon_trk=false;
              break;
            }
          }
          else{
            if(evtsum->NModTwoDimRecons(c_mod,icyc,iview)!=1){
              muon_trk=false;
              break;
            }
            if(evtsum->GetModTwoDimRecon(0,c_mod,icyc,iview)->Nhits()<5){
              muon_trk=false;
              break;
            }
          }
        }
        if(!muon_trk){continue;}

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
        for(int iview=0;iview<2;iview++){
          TwoDimReconSummary* reconsum = evtsum->GetModTwoDimRecon(0,c_mod,icyc+cyc_offset[iview],iview);
          double slope   = reconsum->slope;
          double intcpt  = reconsum->intcpt;
          if(iview==0){slopey=slope;intcpty=intcpt;}
          else        {slopex=slope;intcptx=intcpt;}
        }//iview

        for(int iview=0;iview<2;iview++){
          TwoDimReconSummary* reconsum = evtsum->GetModTwoDimRecon(0,c_mod,icyc+cyc_offset[iview],iview);
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

            h1[imod][hitview][hitpln][hitch]->Fill(hitpe  ,norm_tot);
            h2[imod][hitview][hitpln][hitch]->Fill(hitpemm,norm_tot);

            h3[imod]->Fill(hitpe  ,norm_tot);
            h4[imod]->Fill(hitpemm,norm_tot);

            num_hit   [imod][hitview][hitpln][hitch] += 1;
            pathlength[imod][hitview][hitpln][hitch] += path;
            mean_pe   [imod][hitview][hitpln][hitch] += hitpe;
            mean_pemm [imod][hitview][hitpln][hitch] += hitpemm;
          }
        }

        double cos[NAXIS];
        double norm = sqrt(1.+slopex*slopex+slopey*slopey);
        for(int i=0;i<NAXIS;i++){
          if     (i==0){ cos[i] = 1./norm; }
          else if(i==1){ cos[i] = slopey/norm; }
          else if(i==2){ cos[i] = slopex/norm; }
          h0[imod][i] -> Fill(cos[i],norm_tot);
        }

      } //for(imod)
    } //icyc

    if(last_date!=-1&&(date!=last_date || ievt==nevt-1)){
      for(int i=0;i<NPEMOD;i++){
        for(int j=0;j<NVIEW;j++){
          for(int k=0;k<maxpln;k++){
            for(int l=0;l<maxch;l++){
              if(num_hit  [i][j][k][l] != 0){
                mean_pe   [i][j][k][l] /= num_hit[i][j][k][l];
                mean_pemm [i][j][k][l] /= num_hit[i][j][k][l];
                pathlength[i][j][k][l] /= num_hit[i][j][k][l];
              }
            }
          }
        }
      }
      if(fdid==0) wt->Fill();
      for(int i=0;i<NPEMOD;i++){
        for(int j=0;j<NVIEW;j++){
          for(int k=0;k<maxpln;k++){
            for(int l=0;l<maxch;l++){
              mean_pe   [i][j][k][l] = 0.; 
              mean_pemm [i][j][k][l] = 0.; 
              pathlength[i][j][k][l] = 0.;
              num_hit   [i][j][k][l] = 0 ;
            }
          }
        }
      }
    }
    last_date = date;
  } //ievt

  ofile->cd();
  for(int i=0;i<NPEMOD;i++){
    for(int j=0;j<NAXIS;j++){
      h0[i][j] -> Write();
    }
    for(int j=0;j<NVIEW;j++){
      for(int k=0;k<maxpln;k++){
        for(int l=0;l<maxch;l++){
          h1[i][j][k][l] -> Write();
          h2[i][j][k][l] -> Write();
        }
      }
    }
    h3[i]-> Write();
    h4[i]-> Write();
  }

  if(fdid==0) wt   ->Write();
  ofile->Write();
  ofile->Close();

  delete badch;
  delete detdim;
  return 0;
}
