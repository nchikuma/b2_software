#include "ThreeDimRecon.hxx"

void fMode(int mode);

long fcTime;

void PrintUsage(){
  cout << "Usage:\n";
  cout << "./ThreeDimRecon [options]\n";
  cout << "Either a pair of run/srun number or a pair of input/output filename must be set.\n";
  cout << "Analysis mode must be set from 0 to 2.\n";
  cout << " -r <int>  : Run Number\n";
  cout << " -s <int>  : Sub Run Number\n";
  cout << " -f <char> : Input filename, if you directly indicate it.\n" ;
  cout << " -o <char> : Output filename, if you directly indicate it.\n";
  cout << " -m <int>  : Select analysis mode. {0: SS floor, 1: B2 floor}\n";

  cout << " -i <int>  : To set an event number for starting analysis in the middle.\n";
  cout << " -a        : To analyze all cycles. Only 9 cycles of 4 - 12 are anlaysis by defualt.\n";
  cout << " -c        : To analyze 14-15 cycles for cosmic trigger.\n";
  cout << " -w        : To analyze B2 WAGASCI cosmic trigger data." << endl;
  cout << " -b        : To open all bad channel masked. Bad channels are masked by defualt.\n";
  cout << " -x        : To required one INGRID track or more in a vertex\n";
  cout << " -d        : To display.\n";
  cout << " -v        : For simulation data.\n";
  cout << " -q        : For not jointing b/w modules.\n";
  cout << " -n        : To estimate detector responce.\n";
  cout << " -t        : To use non-3d track.\n";
  cout << " -h        : To show this help.\n";
  cout << endl;
}

int main(int argc,char *argv[]){

  Int_t run_number = -1;
  Int_t sub_run_number = -1;
  Int_t c = -1;
  char  Output[300];
  char  buff1[200];
  bool  rename  = false;
  bool  renamef = false;
  Int_t Scyc    =  4;
  Int_t Ncyc    = 12;
  Int_t Nini    =  0;
  bool  OneEvt    = false;
  bool  disp      = false; 
  bool  cosmic    = false;
  bool  wg_cosmic = false;
  bool  BadChAna  = true;
  bool  simulation = false;
  bool  BeamAna    = false;
  int   detres_para = -1;

  while((c = getopt(argc, argv, "r:s:f:cwbavdxtqo:i:m:n:")) != -1) {
    switch(c){
    case 'r':
      run_number = atoi(optarg);
      break;
    case 's':
      sub_run_number = atoi(optarg);
      break;
    case 'f':
      sprintf(buff1,"%s",optarg);
      renamef    = true;
      break;
    case 'o':
      sprintf(Output,"%s",optarg);
      rename  = true;
      break;
    case 'b':
      BadChAna = false;
      break;
    case 'x':
      BeamAna = true;
      break;
    case 'a':
      Scyc = 0;
      Ncyc = 23;
      break;
    case 'c':
      Scyc   = 14;
      Ncyc   = 16;
      cosmic = true;
      break;
    case 'w':
      wg_cosmic = true;
      break;
    case 'v':
      Scyc   = 4;
      Ncyc   = 5;
      simulation = true;
      break;
    case 'd':
      disp = true;
      break;
    case 'i':
      Nini   = atoi(optarg);
      OneEvt = true;
      break;
    case 'm':
      MODE_DET = atoi(optarg);
      break;
    case 'q':
      NOTJOINT = true;
      break;
    case 't':
      USE_NON_3DTRK = true;
      break;
    case 'n':
      detres_para = atoi(optarg);
      if     (detres_para== 0){ cout << "Nominal" << endl;}
      else if(detres_para== 1){ ang_th0 = 30; } 
      else if(detres_para== 2){ ang_th0 = 40; } 
      else if(detres_para== 3){ pos_th0 = 140; } 
      else if(detres_para== 4){ pos_th0 = 160; } 
      else if(detres_para== 5){ ang_th1 = 30; } 
      else if(detres_para== 6){ ang_th1 = 40; } 
      else if(detres_para== 7){ pos_th1 = 140; } 
      else if(detres_para== 8){ pos_th1 = 160; } 
      else if(detres_para== 9){ ang_th2 = 30; } 
      else if(detres_para==10){ ang_th2 = 40; } 
      else if(detres_para==11){ pos_th2 = 140; } 
      else if(detres_para==12){ pos_th2 = 160; }
      else if(detres_para==13){ hitpos_th2 = 65; }
      else if(detres_para==14){ hitpos_th2 = 85; }
      else if(detres_para==15){ diff_th2 =  6; }
      else if(detres_para==16){ diff_th2 = 12; }
      else if(detres_para==17){ diff_th1 =  2; }
      else if(detres_para==18){ diff_th1 =  4; }
      else if(detres_para==19){ pln_th2  =  6; }
      else if(detres_para==20){ pln_th2  = 12; }
      else if(detres_para==21){ pln_th1  =  2; }
      else if(detres_para==22){ pln_th1  =  4; }
      else if(detres_para==23){ ch_th2   = 100; }
      else if(detres_para==24){ ch_th2   = 200; }
      else if(detres_para==25){ ch_th1   = 100; }
      else if(detres_para==26){ ch_th1   = 200; }
      else{
        cout << "Unexpected parameter." << endl;
        exit(0);
      }
      cout << "ang_th0   =" <<  ang_th0   <<endl; 
      cout << "pos_th0   =" <<  pos_th0   <<endl; 
      cout << "ang_th1   =" <<  ang_th1   <<endl; 
      cout << "pos_th1   =" <<  pos_th1   <<endl; 
      cout << "ang_th2   =" <<  ang_th2   <<endl; 
      cout << "pos_th2   =" <<  pos_th2   <<endl; 
      cout << "hitpos_th2=" <<  hitpos_th2<<endl; 
      cout << "diff_th2  =" <<  diff_th2  <<endl; 
      cout << "diff_th1  =" <<  diff_th1  <<endl; 
      cout << "pln_th2   =" <<  pln_th2   <<endl; 
      cout << "pln_th1   =" <<  pln_th1   <<endl; 
      cout << "ch_th2    =" <<  ch_th2    <<endl; 
      cout << "ch_th1    =" <<  ch_th1    <<endl; 
      break;
    default:
      PrintUsage();
      exit(0);
    }
  }

  if(
      ((run_number<0||sub_run_number<0)&&(!rename||!renamef))||
      (wg_cosmic&&!cosmic)||
      (wg_cosmic&&MODE_DET!=1))
  {
    PrintUsage(); exit(0);
  }

  //Modules @SS
  if(MODE_DET==0){
    INGMODNUM_start = MOD_INGRID_C;
    INGMODNUM_end   = MOD_INGRID_C;
    INGMODNUM_mid   = MOD_INGRID_C;
    WMMODNUM        = MOD_ONAXIS_WM;
    PMMODNUM        = -1;//MOD_PM;
  }
  //Modules @B2
  else if(MODE_DET==1){
    INGMODNUM_start = MOD_B2_INGRID;
    INGMODNUM_end   = MOD_B2_INGRID;
    INGMODNUM_mid   = MOD_B2_INGRID;
    WMMODNUM        = MOD_B2_WM;
    PMMODNUM        = MOD_PM;   //MOD_B2_CH;
  }
  else{
    PrintUsage();
    exit(1);
  }

  if(wg_cosmic&&!renamef){
    sprintf(buff1,"%s/run_%05d_%03d_recon.root",
	    cosmic_wagasci,run_number,sub_run_number);
  }
  else if(cosmic&&!renamef){
    sprintf(buff1,"%s/ingrid_%08d_%04d_recon.root",
	    cosmic_data, run_number, sub_run_number);
  }
  else if(!renamef){
    sprintf(buff1,"%s/ingrid_%08d_%04d_recon.root",
	    //dst_ingrid, run_number, sub_run_number);
	    dst_ingrid2, run_number, sub_run_number);
  }

  //============================================
  //==========     Read root file     ==========
  //============================================
  FileStat_t fs;
  if(gSystem->GetPathInfo(buff1,fs)){
    cout << "Cannot open file " << buff1 << endl;
    exit(1);
  }
  cout << "Input filename : " << buff1 << endl;
  EventSummary* evt   =  new EventSummary();
  TFile*        rfile =  new TFile(buff1, "read");
  TTree*        tree  =  (TTree*)rfile -> Get("tree");
  TBranch*      EvtBr =  tree->GetBranch("fDefaultReco.");
  EvtBr               -> SetAddress(&evt);
  tree                -> SetBranchAddress("fDefaultReco.", &evt);
  int            nevt =  (int)tree -> GetEntries();

  //==========================================================
  //==========     Make rootfile after analysis     ==========
  //==========================================================
  if(wg_cosmic&&!rename){
    sprintf(Output, "%s/run_%05d_%03d_anas%d.root", 
	    cosmic_wagasci, run_number, sub_run_number,MODE_DET); 
  }
  else if(cosmic&&!rename){
    sprintf(Output, "%s/ingrid_%08d_%04d_anas%d.root", 
	    cosmic_data, run_number, sub_run_number,MODE_DET); 
  }
  else if(!rename){
    sprintf(Output, "%s/ingrid_%08d_%04d_anas%d.root", 
	    //dst_ingrid, run_number, sub_run_number,MODE_DET); 
	    dst_ingrid2, run_number, sub_run_number,MODE_DET); 
  }
  cout << "Output filename : " << Output << endl;
  TFile*            wfile =  new TFile(Output, "recreate");
  TTree*            wtree =  new TTree("tree","tree");
  wtree                   -> SetMaxTreeSize(5000000000);
  EventSummary* wsummary  =  new EventSummary(); 
  wtree                   -> Branch   ("fDefaultReco.","EventSummary", 
				       &wsummary, 64000, 99);
  HitSummary*           hitsum;
  TwoDimReconSummary*   recon;
  ThreeDimReconSummary* anasum = new ThreeDimReconSummary();
  Hits                  hits;
  TrackPM               pmtrack;
  TrackWM               wmtrack;
  TrackIng              ingtrack;

  detdim  = new DetectorDimension();
  badch = new INGRID_BadCh_mapping();
  badch->set_BadCh(BadChAna);

  int EvtStop = -1;
  if(OneEvt) EvtStop = Nini+1;
  else       EvtStop = nevt;

  cout << "Total # of events = " << nevt <<endl;
  cout << "Analyze : event# " << Nini << " to " << EvtStop << endl;
  for(int ievt=Nini; ievt<EvtStop; ievt++){

#ifdef DEBUG_THREEDIMRECON
    cout << "\n\n ====== analyze event# " << ievt << " ====== " <<endl;
#else
    if(ievt%100==0){
      cout << "====== analyze event# " << ievt << " ====== " <<endl;
    }
#endif 

    wsummary -> Clear();
    evt      -> Clear();
    tree     -> GetEntry(ievt);
    
    wsummary = evt;

    for(cyc=Scyc; cyc<Ncyc; cyc++){
#ifdef DEBUG_THREEDIMRECON
      cout << " ----- cyc = " << cyc << " ------------------ " << endl;
#endif

      //------- Fill tracks on INGRID Horizontal Modules & B2 INGRID -------
      hingtrack.clear();
      vingtrack.clear();
      for(int mod=INGMODNUM_start; mod<=INGMODNUM_end; mod++){
        for(int view=0; view<2; view++){
          int ningtrack = evt -> NModTwoDimRecons(mod, cyc, view);
          for(int i=0; i<ningtrack; i++){

            recon   = (TwoDimReconSummary*)(evt -> GetModTwoDimRecon(i, mod, cyc, view));
            ingtrack.clear();
            ingtrack.mod     = mod;
            ingtrack.view    = recon -> view;
            ingtrack.ipln    = recon -> startpln;
            ingtrack.fpln    = recon -> endpln;
            ingtrack.iz      = recon -> startz;
            ingtrack.fz      = recon -> endz;
            ingtrack.ang     = recon -> angle;
            ingtrack.slope   = recon -> slope;
            ingtrack.veto    = recon -> vetowtracking;
            ingtrack.edge    = recon -> edgewtracking;
            ingtrack.stop    = recon -> modfc;
            ingtrack.clstime = recon -> clstime;
            ingtrack.ixy     = recon -> startxy;
            ingtrack.fxy     = recon -> endxy;
            ingtrack.intcpt  = recon -> intcpt;

            if(view==0){
              ingtrack.ixy    += (mod-INGMODNUM_mid)*C_INGSpace;
              ingtrack.fxy    += (mod-INGMODNUM_mid)*C_INGSpace;
              ingtrack.intcpt += (mod-INGMODNUM_mid)*C_INGSpace;
            }

            for(int NHIT=0; NHIT<(recon->Nhits()); NHIT++){
              hitsum = recon -> GetHit(NHIT);
              hits.clear();
              hits.mod      = hitsum->mod;
              if(simulation) hits.cyc = 4;
              else           hits.cyc = hitsum->cyc;
              hits.view     = hitsum->view;
              hits.pln      = hitsum->pln;
              hits.ch       = hitsum->ch;
              hits.pe       = hitsum->pe;
              hits.lope     = hitsum->lope;
              hits.isohit   = hitsum->isohit;
              hits.recon_id = i;
              hits.hit_id   = NHIT;

              if((hitsum -> NSimHits()) > 0){
                hits.pdg = hitsum -> GetSimHit(0)->pdg;
	      }
              else{
                hits.pdg = 0;
	      }
              ingtrack.hit.push_back(hits);
            }

            if(view==0) hingtrack.push_back(ingtrack);
            else        vingtrack.push_back(ingtrack);
          }
        }
      }

      //------- Fill tracks on Water Module -------
      hwmtrack.clear();
      vwmtrack.clear();
      if((WMMODNUM==MOD_B2_WM)&&(cyc<22)&&(!simulation)){max_dif_cyc=1;} else{max_dif_cyc=0;}
      for(int icyc=0;icyc<=max_dif_cyc;icyc++){
        for(int view=0; view<2; view++){
          int nwmtrack = evt->NModTwoDimRecons(WMMODNUM, cyc+icyc, view);
          for(int i=0; i<nwmtrack; i++){
            recon   = (TwoDimReconSummary*)(evt -> GetModTwoDimRecon(i, WMMODNUM, cyc+icyc, view) );
            wmtrack.clear();
            wmtrack.view    = recon -> view;
            wmtrack.ipln    = recon -> startpln;
            wmtrack.fpln    = recon -> endpln;
            wmtrack.ixy     = recon -> startxy;
            wmtrack.fxy     = recon -> endxy;
            wmtrack.iz      = recon -> startz;
            wmtrack.fz      = recon -> endz;
            wmtrack.ang     = recon -> angle;
            wmtrack.slope   = recon -> slope;
            wmtrack.intcpt  = recon -> intcpt;
            wmtrack.veto    = recon -> vetowtracking;
            wmtrack.edge    = recon -> edgewtracking;
            wmtrack.stop    = recon -> modfc;
            if(simulation) wmtrack.clstime = recon -> clstime;
            else           wmtrack.clstime = recon -> clstime - wg_hittime_corr;
            wmtrack.ing_trk = false;

            for(int NHIT=0; NHIT<(recon->Nhits()); NHIT++){
              hitsum = recon -> GetHit(NHIT);
              hits.clear();
              hits.mod      = hitsum->mod;
              if(simulation) hits.cyc = 4;
              else           hits.cyc = hitsum->cyc;
              hits.view     = hitsum->view;
              hits.pln      = hitsum->pln;
              hits.ch       = hitsum->ch;
              hits.pe       = hitsum->pe;
              hits.lope     = hitsum->lope;
              hits.isohit   = hitsum->isohit;
              hits.recon_id = i;
              hits.hit_id   = NHIT;

              int reconpln, reconch;
              int axis = 0;
              detdim -> GetReconPlnCh(hits.mod, hits.view, hits.pln, hits.ch, axis, &reconpln, &reconch);
              hits.pln = reconpln;
              hits.ch  = reconch;

              if((hitsum->NSimHits())>0){
                hits.pdg = hitsum->GetSimHit(0)->pdg;
              }
              else{
                hits.pdg = 0;
              }
              wmtrack.hit.push_back(hits);
            }

            if(view==0) hwmtrack.push_back(wmtrack);
            else        vwmtrack.push_back(wmtrack);
          }
        }
      }

      //------- Fill tracks on Proton Module -------
      hpmtrack.clear();
      vpmtrack.clear();
      for(int view=0; view<2; view++){
        if(PMMODNUM<0){continue;}
        int npmtrack = evt -> NModTwoDimRecons(PMMODNUM, cyc, view);
        for(int i=0; i<npmtrack; i++){
          recon   = (TwoDimReconSummary*)(evt -> GetModTwoDimRecon(i, PMMODNUM, cyc, view) );
          pmtrack.clear();
          pmtrack.view    = recon -> view;
          pmtrack.ipln    = recon -> startpln;
          pmtrack.fpln    = recon -> endpln;
          pmtrack.ixy     = recon -> startxy;
          pmtrack.fxy     = recon -> endxy;
          pmtrack.iz      = recon -> startz;
          pmtrack.fz      = recon -> endz;
          pmtrack.ang     = recon -> angle;
          pmtrack.slope   = recon -> slope;
          pmtrack.intcpt  = recon -> intcpt;
          pmtrack.veto    = recon -> vetowtracking;
          pmtrack.edge    = recon -> edgewtracking;
          pmtrack.stop    = recon -> modfc;
          pmtrack.clstime = recon -> clstime;
          pmtrack.ing_trk = false;

          for(int NHIT=0; NHIT<(recon->Nhits()); NHIT++){
            hitsum = recon->GetHit(NHIT);
            hits.clear();
            hits.mod      = hitsum->mod;
            if(simulation) hits.cyc = 4;
            else           hits.cyc = hitsum->cyc;
            hits.view     = hitsum->view;
            hits.pln      = hitsum->pln;
            hits.ch       = hitsum->ch;
            hits.pe       = hitsum->pe;
            hits.lope     = hitsum->lope;
            hits.isohit   = hitsum->isohit;
            hits.recon_id = i;
            hits.hit_id   = NHIT;

            int reconpln, reconch;
            int axis = 0;
            detdim -> GetReconPlnCh( hits.mod, hits.view, hits.pln, hits.ch, axis, &reconpln, &reconch );
            hits.pln = reconpln;
            hits.ch  = reconch;

            if((hitsum->NSimHits())>0){
              hits.pdg = hitsum -> GetSimHit(0)->pdg;
	    }
            else{
              hits.pdg = 0;
	    }

            pmtrack.hit.push_back(hits);
          }

          if(view==0) hpmtrack.push_back(pmtrack);
          else        vpmtrack.push_back(pmtrack);
        }
      }

      GetNonRecHits(WMMODNUM,INGMODNUM_start,INGMODNUM_end,PMMODNUM,evt,cyc);

      if(!fAnalyze_new()){ 
#ifdef DEBUG_THREEDIMRECON
        cout << "      No vertex." << endl;
#endif
        continue;
      }
#ifdef DEBUG_THREEDIMRECON
      cout << "------- Start Filling the vertex info. --------------" <<endl;
#endif

      for(int i=0; i<(int)analyzed_trk.size(); i++){
        if(analyzed_trk[i].Ntrack==0) continue;
        if(BeamAna&&!NOTJOINT&&analyzed_trk[i].Ningtrack==0) continue; //To remove duplicate track since random jitter.

        anasum -> Ntrack        = analyzed_trk[i].Ntrack;
        anasum -> Ningtrack     = analyzed_trk[i].Ningtrack;
        if(anasum->Ntrack>RECON_MAXTRACKS){
	  anasum->Ntrack    = RECON_MAXTRACKS;
	}
        if(anasum->Ningtrack>RECON_MAXTRACKS){
	  anasum->Ningtrack = RECON_MAXTRACKS;
	}
        anasum -> vetowtracking = analyzed_trk[i].vetowtracking;
        anasum -> edgewtracking = analyzed_trk[i].edgewtracking;
        anasum -> hitcyc        = cyc;
        anasum -> clstime       = analyzed_trk[i].clstime;
        anasum -> exptime       = ((int)anasum -> clstime - TDCOffset )% GapbwBunch - fExpBeamTime;
        if(fabs(anasum->exptime)<beamontime){
          anasum->ontime = true;
	}
        else{
          anasum->ontime = false;
	}

        anasum -> startmod    .clear();
        anasum -> stopmodx    .clear();
        anasum -> stopmody    .clear();

        anasum -> x           .clear();
        anasum -> y           .clear();
        anasum -> z           .clear();
        anasum -> zx          .clear();
        anasum -> zy          .clear();

        anasum -> startxpln   .clear();
        anasum -> startypln   .clear();
        anasum -> startxch    .clear();
        anasum -> startych    .clear();
        anasum -> endxpln     .clear();
        anasum -> endypln     .clear();
        anasum -> endxch      .clear();
        anasum -> endych      .clear();
        anasum -> thetax      .clear();
        anasum -> thetay      .clear();
        anasum -> angle       .clear();

        anasum -> ing_startmod.clear();
        anasum -> ing_endmod  .clear();
        anasum -> ing_startpln.clear();
        anasum -> ing_endpln  .clear();
        anasum -> ing_trk     .clear();
        anasum -> pm_stop     .clear();
        anasum -> ing_stop    .clear();
        anasum -> sci_range   .clear();
        anasum -> iron_range  .clear();
        anasum -> iron_pene   .clear();

        anasum -> veto        .clear();
        anasum -> edge        .clear();
        anasum -> pdg         .clear();
        anasum -> trkpe       .clear();
        anasum -> mucl        .clear();
        anasum -> oneview     .clear();
        anasum -> diff_posx   .clear();
        anasum -> diff_posy   .clear();
        anasum -> diff_angx   .clear();
        anasum -> diff_angy   .clear();
        anasum -> diff_timex  .clear();
        anasum -> diff_timey  .clear();
        anasum -> diff_posx2  .clear();
        anasum -> diff_posy2  .clear();
        anasum -> diff_angx2  .clear();
        anasum -> diff_angy2  .clear();
        anasum -> diff_timex2 .clear();
        anasum -> diff_timey2 .clear();
        anasum -> diff_posx3  .clear();
        anasum -> diff_posy3  .clear();
        anasum -> diff_angx3  .clear();
        anasum -> diff_angy3  .clear();
        anasum -> diff_timex3 .clear();
        anasum -> diff_timey3 .clear();

        anasum->nhits = 0;
        for(int itrk=0; itrk<RECON_MAXTRACKS;itrk++){
	  anasum->nhitTs[itrk] = 0;
	}

        for(int t=0; t<(int)analyzed_trk[i].trk.size(); t++){
          if(t>=RECON_MAXTRACKS) continue;

          anasum -> startmod    .push_back(analyzed_trk[i].trk[t].startmod);
          anasum -> stopmodx    .push_back(analyzed_trk[i].trk[t].stopmodx);
          anasum -> stopmody    .push_back(analyzed_trk[i].trk[t].stopmody);

          anasum -> x           .push_back(analyzed_trk[i].trk[t].x);
          anasum -> y           .push_back(analyzed_trk[i].trk[t].y);
          anasum -> z           .push_back(analyzed_trk[i].trk[t].z);
          anasum -> zx          .push_back(analyzed_trk[i].trk[t].zx);
          anasum -> zy          .push_back(analyzed_trk[i].trk[t].zy);

          anasum -> startxpln   .push_back(analyzed_trk[i].trk[t].startxpln);
          anasum -> startypln   .push_back(analyzed_trk[i].trk[t].startypln);
          anasum -> startxch    .push_back(analyzed_trk[i].trk[t].startxch);
          anasum -> startych    .push_back(analyzed_trk[i].trk[t].startych);
          anasum -> endxpln     .push_back(analyzed_trk[i].trk[t].endxpln);
          anasum -> endypln     .push_back(analyzed_trk[i].trk[t].endypln);
          anasum -> endxch      .push_back(analyzed_trk[i].trk[t].endxch);
          anasum -> endych      .push_back(analyzed_trk[i].trk[t].endych);
          anasum -> thetax      .push_back(analyzed_trk[i].trk[t].thetax);
          anasum -> thetay      .push_back(analyzed_trk[i].trk[t].thetay);
          anasum -> angle       .push_back(analyzed_trk[i].trk[t].angle);

          anasum -> ing_startmod.push_back(analyzed_trk[i].trk[t].ing_startmod);
          anasum -> ing_endmod  .push_back(analyzed_trk[i].trk[t].ing_endmod);
          anasum -> ing_startpln.push_back(analyzed_trk[i].trk[t].ing_startpln);
          anasum -> ing_endpln  .push_back(analyzed_trk[i].trk[t].ing_endpln);
          anasum -> ing_trk     .push_back(analyzed_trk[i].trk[t].ing_trk);
          anasum -> pm_stop     .push_back(analyzed_trk[i].trk[t].pm_stop);
          anasum -> ing_stop    .push_back(analyzed_trk[i].trk[t].ing_stop);
          anasum -> sci_range   .push_back(analyzed_trk[i].trk[t].sci_range);
          anasum -> iron_range  .push_back(analyzed_trk[i].trk[t].iron_range);
          anasum -> iron_pene   .push_back(analyzed_trk[i].trk[t].iron_pene);

          anasum -> veto        .push_back(analyzed_trk[i].trk[t].vetowtracking);
          anasum -> edge        .push_back(analyzed_trk[i].trk[t].edgewtracking);
          anasum -> pdg         .push_back(analyzed_trk[i].trk[t].pdg);
          anasum -> trkpe       .push_back(analyzed_trk[i].trk[t].trkpe);
          anasum -> mucl        .push_back(analyzed_trk[i].trk[t].mucl);
          anasum -> oneview     .push_back(analyzed_trk[i].trk[t].oneview);
          anasum -> diff_posx   .push_back(analyzed_trk[i].trk[t].diff_posx);
          anasum -> diff_posy   .push_back(analyzed_trk[i].trk[t].diff_posy);
          anasum -> diff_angx   .push_back(analyzed_trk[i].trk[t].diff_angx);
          anasum -> diff_angy   .push_back(analyzed_trk[i].trk[t].diff_angy);
          anasum -> diff_timex  .push_back(analyzed_trk[i].trk[t].diff_timex);
          anasum -> diff_timey  .push_back(analyzed_trk[i].trk[t].diff_timey);
          anasum -> diff_posx2  .push_back(analyzed_trk[i].trk[t].diff_posx2);
          anasum -> diff_posy2  .push_back(analyzed_trk[i].trk[t].diff_posy2);
          anasum -> diff_angx2  .push_back(analyzed_trk[i].trk[t].diff_angx2);
          anasum -> diff_angy2  .push_back(analyzed_trk[i].trk[t].diff_angy2);
          anasum -> diff_timex2 .push_back(analyzed_trk[i].trk[t].diff_timex2);
          anasum -> diff_timey2 .push_back(analyzed_trk[i].trk[t].diff_timey2);
          anasum -> diff_posx3  .push_back(analyzed_trk[i].trk[t].diff_posx3);
          anasum -> diff_posy3  .push_back(analyzed_trk[i].trk[t].diff_posy3);
          anasum -> diff_angx3  .push_back(analyzed_trk[i].trk[t].diff_angx3);
          anasum -> diff_angy3  .push_back(analyzed_trk[i].trk[t].diff_angy3);
          anasum -> diff_timex3 .push_back(analyzed_trk[i].trk[t].diff_timex3);
          anasum -> diff_timey3 .push_back(analyzed_trk[i].trk[t].diff_timey3);

	  //double anglex = analyzed_trk[i].trk[t].thetax;
	  //double angley = analyzed_trk[i].trk[t].thetay;
	  double slopex = analyzed_trk[i].trk[t].slopex;
	  double slopey = analyzed_trk[i].trk[t].slopey;

          int hnum    = analyzed_trk[i].trk[t].hnum;
          int vnum    = analyzed_trk[i].trk[t].vnum;
          int vmod    = analyzed_trk[i].trk[t].startmod;
          int oneview = analyzed_trk[i].trk[t].oneview;
          bool hview = (oneview==0||oneview==1);
          bool vview = (oneview==0||oneview==2);

          vector<Hits>::iterator iter;

          if(vmod==MOD_ONAXIS_WM||vmod==MOD_B2_WM){
            if(hview){
              for(iter=hwmtrack[hnum].hit.begin(); iter<hwmtrack[hnum].hit.end(); iter++){
                Hits hit = (*iter);
#ifdef DEBUG_THREEDIMRECON
                cout
                  << " {" << hit.mod
                  << ","  << hit.view
                  << ","  << hit.pln
                  << ","  << hit.ch
                  << "}";
#endif

                TwoDimReconSummary* recsum = 
                  evt->GetModTwoDimRecon(hit.recon_id,hit.mod,hit.cyc,hit.view);
                HitSummary* hitsummary = recsum->GetHit(hit.hit_id);
                double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                hitsummary->pathlength = path;
                if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                else       { hitsummary->pe_permm = 0.;         }
                anasum->AddHitTrk(hitsummary,t);
              }
            }
            if(vview){
              for(iter=vwmtrack[vnum].hit.begin(); iter<vwmtrack[vnum].hit.end(); iter++){
                Hits hit = (*iter);
#ifdef DEBUG_THREEDIMRECON
                cout
                  << " {" << hit.mod
                  << ","  << hit.view
                  << ","  << hit.pln
                  << ","  << hit.ch
                  << "}";
#endif

                TwoDimReconSummary* recsum = 
                  evt->GetModTwoDimRecon(hit.recon_id,hit.mod,hit.cyc,hit.view);
                HitSummary* hitsummary = recsum->GetHit(hit.hit_id);
                double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                hitsummary->pathlength = path;
                if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                else       { hitsummary->pe_permm = 0.;         }
                anasum->AddHitTrk(hitsummary,t);
              }
            }

            int hingnum = hwmtrack[hnum].ing_num;
            int vingnum = vwmtrack[vnum].ing_num;
            if(hwmtrack[hnum].ing_trk&&hview){
              for(iter=hingtrack[hingnum].hit.begin(); iter<hingtrack[hingnum].hit.end(); iter++){
                Hits hit = (*iter);
                if(hit.recon_id!=-1){
                  TwoDimReconSummary* recsum = 
                    evt->GetModTwoDimRecon(hit.recon_id,hit.mod,hit.cyc,hit.view);
                  HitSummary* hitsummary = recsum->GetHit(hit.hit_id);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
                else{
                  HitSummary* hitsummary = evt->GetModHit(hit.hit_id,hit.mod,hit.cyc);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
              }
            }
            if(vwmtrack[vnum].ing_trk&&vview) {
              for(iter=vingtrack[vingnum].hit.begin(); iter<vingtrack[vingnum].hit.end(); iter++){
                Hits hit = (*iter);
                if(hit.recon_id!=-1){
                  TwoDimReconSummary* recsum = 
                    evt->GetModTwoDimRecon(hit.recon_id,hit.mod,hit.cyc,hit.view);
                  HitSummary* hitsummary = recsum->GetHit(hit.hit_id);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
                else{
                  HitSummary* hitsummary = evt->GetModHit(hit.hit_id,hit.mod,hit.cyc);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
              }
            }
          }
          else if(vmod==MOD_PM||vmod==MOD_B2_CH)
          {
            int  hnum    = analyzed_trk[i].trk[t].hnum;
            int  vnum    = analyzed_trk[i].trk[t].vnum;
            bool h_wmtrk = hpmtrack[hnum].wm_trk;
            bool v_wmtrk = vpmtrack[vnum].wm_trk;
            bool h_ingtrk= hpmtrack[hnum].ing_trk;
            bool v_ingtrk= vpmtrack[vnum].ing_trk;
#ifdef DEBUG_THREEDIMRECON
            cout
              << " hnum="     << hnum    
              << " vnum="     << vnum    
              << " h_wmtrk="  << h_wmtrk 
              << " v_wmtrk="  << v_wmtrk 
              << " h_ingtrk=" << h_ingtrk
              << " v_ingtrk=" << v_ingtrk
              << endl;
#endif
            if(hview){
              for(iter=hpmtrack[hnum].hit.begin(); iter<hpmtrack[hnum].hit.end(); iter++){
                Hits hit = (*iter);
#ifdef DEBUG_THREEDIMRECON
                cout
                  << " {" << hit.mod
                  << ","  << hit.view
                  << ","  << hit.pln
                  << ","  << hit.ch
                  << "}";
#endif
                if(hit.recon_id!=-1){
                  TwoDimReconSummary* recsum = 
                    evt->GetModTwoDimRecon(hit.recon_id,hit.mod,hit.cyc,hit.view);
                  HitSummary* hitsummary = recsum->GetHit(hit.hit_id);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
                else{
                  HitSummary* hitsummary = evt->GetModHit(hit.hit_id,hit.mod,hit.cyc);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
              }
            }
            if(vview){
              for(iter=vpmtrack[vnum].hit.begin(); iter<vpmtrack[vnum].hit.end(); iter++){
                Hits hit = (*iter);
#ifdef DEBUG_THREEDIMRECON
                cout
                  << " {" << hit.mod
                  << ","  << hit.view
                  << ","  << hit.pln
                  << ","  << hit.ch
                  << "}";
#endif
                if(hit.recon_id!=-1){
                  TwoDimReconSummary* recsum = 
                    evt->GetModTwoDimRecon(hit.recon_id,hit.mod,hit.cyc,hit.view);
                  HitSummary* hitsummary = recsum->GetHit(hit.hit_id);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
                else{
                  HitSummary* hitsummary = evt->GetModHit(hit.hit_id,hit.mod,hit.cyc);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
              }
            }
            int hwmnum = hpmtrack[hnum].wm_num;
            int vwmnum = vpmtrack[vnum].wm_num;
            if(h_wmtrk&&hview){
              for(iter=hwmtrack[hwmnum].hit.begin(); iter<hwmtrack[hwmnum].hit.end(); iter++){
                Hits hit = (*iter);
#ifdef DEBUG_THREEDIMRECON
              cout
                << " {" << hit.mod
                << ","  << hit.view
                << ","  << hit.pln
                << ","  << hit.ch
                << "}";
#endif
                if(hit.recon_id!=-1){
                  TwoDimReconSummary* recsum = 
                    evt->GetModTwoDimRecon(hit.recon_id,hit.mod,hit.cyc,hit.view);
                  HitSummary* hitsummary = recsum->GetHit(hit.hit_id);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
                else{
                  HitSummary* hitsummary = evt->GetModHit(hit.hit_id,hit.mod,hit.cyc);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
              }
            }
            if(v_wmtrk&&vview){
              for(iter=vwmtrack[vwmnum].hit.begin(); iter<vwmtrack[vwmnum].hit.end(); iter++){
                Hits hit = (*iter); 
#ifdef DEBUG_THREEDIMRECON
              cout
                << " {" << hit.mod
                << ","  << hit.view
                << ","  << hit.pln
                << ","  << hit.ch
                << "}";
#endif
                if(hit.recon_id!=-1){
                  TwoDimReconSummary* recsum = 
                    evt->GetModTwoDimRecon(hit.recon_id,hit.mod,hit.cyc,hit.view);
                  HitSummary* hitsummary = recsum->GetHit(hit.hit_id);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
                else{
                  HitSummary* hitsummary = evt->GetModHit(hit.hit_id,hit.mod,hit.cyc);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
              }
            }
            int hingnum = hpmtrack[hnum].ing_num;
            int vingnum = vpmtrack[vnum].ing_num;
            if(h_ingtrk&&hview){
              for(iter=hingtrack[hingnum].hit.begin(); iter<hingtrack[hingnum].hit.end(); iter++){
                Hits hit = (*iter); 
#ifdef DEBUG_THREEDIMRECON
              cout
                << " {" << hit.mod
                << ","  << hit.view
                << ","  << hit.pln
                << ","  << hit.ch
                << "}";
#endif

                if(hit.recon_id!=-1){
                  TwoDimReconSummary* recsum = 
                    evt->GetModTwoDimRecon(hit.recon_id,hit.mod,hit.cyc,hit.view);
                  HitSummary* hitsummary = recsum->GetHit(hit.hit_id);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
                else{
                  HitSummary* hitsummary = evt->GetModHit(hit.hit_id,hit.mod,hit.cyc);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
              }
            }
            if(v_ingtrk&&vview) {
              for(iter=vingtrack[vingnum].hit.begin(); iter<vingtrack[vingnum].hit.end(); iter++){
                Hits hit = (*iter);
#ifdef DEBUG_THREEDIMRECON
              cout
                << " {" << hit.mod
                << ","  << hit.view
                << ","  << hit.pln
                << ","  << hit.ch
                << "}";
#endif

                if(hit.recon_id!=-1){
                  TwoDimReconSummary* recsum = 
                    evt->GetModTwoDimRecon(hit.recon_id,hit.mod,hit.cyc,hit.view);
                  HitSummary* hitsummary = recsum->GetHit(hit.hit_id);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
                else{
                  HitSummary* hitsummary = evt->GetModHit(hit.hit_id,hit.mod,hit.cyc);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
              }
            }
#ifdef DEBUG_THREEDIMRECON
            cout << endl;
#endif

          }
          else if(vmod==MOD_INGRID_C||vmod==MOD_B2_INGRID)
          {
            int hnum = analyzed_trk[i].trk[t].hnum;
            int vnum = analyzed_trk[i].trk[t].vnum;
            if(hview){
              for(iter=hingtrack[hnum].hit.begin(); iter<hingtrack[hnum].hit.end(); iter++){
                Hits hit = (*iter); 
                if(hit.recon_id!=-1){
                  TwoDimReconSummary* recsum = 
                    evt->GetModTwoDimRecon(hit.recon_id,hit.mod,hit.cyc,hit.view);
                  HitSummary* hitsummary = recsum->GetHit(hit.hit_id);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
                else{
                  HitSummary* hitsummary = evt->GetModHit(hit.hit_id,hit.mod,hit.cyc);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
              }
            }
            if(vview){
              for(iter=vingtrack[vnum].hit.begin(); iter<vingtrack[vnum].hit.end(); iter++){
                Hits hit = (*iter);
                if(hit.recon_id!=-1){
                  TwoDimReconSummary* recsum = 
                    evt->GetModTwoDimRecon(hit.recon_id,hit.mod,hit.cyc,hit.view);
                  HitSummary* hitsummary = recsum->GetHit(hit.hit_id);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
                else{
                  HitSummary* hitsummary = evt->GetModHit(hit.hit_id,hit.mod,hit.cyc);
                  double path = calc_pathlength_wg(slopex,slopey,hit.mod,hit.view,hit.pln,hit.ch);
                  hitsummary->pathlength = path;
                  if(path>0.){ hitsummary->pe_permm = hit.pe/path;}
                  else       { hitsummary->pe_permm = 0.;         }
                  anasum->AddHitTrk(hitsummary,t);
                }
              }
            }
          }
        }

#ifdef DEBUG_THREEDIMRECON
        std::cout << " \nNtrack       " << anasum->Ntrack;
        std::cout << " \nNingtrack    " << anasum->Ningtrack;
        std::cout << " \nclstime      " << anasum->clstime;
        std::cout << " \nclstimecorr  " << anasum->clstimecorr;
        std::cout << " \nexptime      " << anasum->exptime;
        std::cout << " \nhitcyc       " << anasum->hitcyc;
        std::cout << " \nontime       " << anasum->ontime;
        std::cout << " \nvetowtracking" << anasum->vetowtracking;
        std::cout << " \nedgewtracking" << anasum->edgewtracking;
        std::cout << " \nnhits        " << anasum->nhits;
        std::cout << " \nnhitTs[RECON_MAXTRACKS]";
        for(int inhitTs=0; inhitTs<anasum->Ntrack; inhitTs++){
          std::cout << " " << anasum->nhitTs[inhitTs];
	}
        std::cout << " \nvact[100]    ";
        for(int ivact=0; ivact<100; ivact++){
          std::cout << " " << anasum->vact[ivact];
	}
        for(int t=0; t<(int)analyzed_trk[i].trk.size(); t++){
          std::cout << " \nstartmod     " << anasum->startmod[t];
          std::cout << " \nstopmodx     " << anasum->stopmodx[t];
          std::cout << " \nstopmody     " << anasum->stopmody[t];
          std::cout << " \nx            " << anasum->x[t];
          std::cout << " \ny            " << anasum->y[t];
          std::cout << " \nz            " << anasum->z[t];
          std::cout << " \nzx           " << anasum->zx[t];
          std::cout << " \nzy           " << anasum->zy[t];
          std::cout << " \nstartxpln    " << anasum->startxpln[t];
          std::cout << " \nstartypln    " << anasum->startypln[t];
          std::cout << " \nstartxch     " << anasum->startxch[t];
          std::cout << " \nstartych     " << anasum->startych[t];
          std::cout << " \nendxpln      " << anasum->endxpln[t];
          std::cout << " \nendypln      " << anasum->endypln[t];
          std::cout << " \nendxch       " << anasum->endxch[t];
          std::cout << " \nendych       " << anasum->endych[t];
          std::cout << " \nthetax       " << anasum->thetax[t];
          std::cout << " \nthetay       " << anasum->thetay[t];
          std::cout << " \nangle        " << anasum->angle[t];
          std::cout << " \ning_startmod " << anasum->ing_startmod[t];
          std::cout << " \ning_endmod   " << anasum->ing_endmod[t];
          std::cout << " \ning_startpln " << anasum->ing_startpln[t];
          std::cout << " \ning_endpln   " << anasum->ing_endpln[t];
          std::cout << " \ning_trk      " << anasum->ing_trk[t];
          std::cout << " \npm_stop      " << anasum->pm_stop[t];
          std::cout << " \ning_stop     " << anasum->ing_stop[t];
          std::cout << " \nsci_range    " << anasum->sci_range[t];
          std::cout << " \niron_range   " << anasum->iron_range[t];
          std::cout << " \niron_pene    " << anasum->iron_pene[t];
          std::cout << " \nveto         " << anasum->veto[t];
          std::cout << " \nedge         " << anasum->edge[t];
          std::cout << " \npdg          " << anasum->pdg[t];
          std::cout << " \nmucl         " << anasum->mucl[t];
          std::cout << " \ntrkpe        " << anasum->trkpe[t];
          std::cout << " \noneview      " << anasum->oneview[t];
          std::cout << "\n";
        }
#endif
        //evt->AddThreeDimRecon(anasum);
        wsummary->AddThreeDimRecon(anasum);
      }
    }

    //wsummary = evt;
    wtree -> Fill();
  }

  wfile->cd();
  //------- Write and Close -------
  wtree->Write();
  wfile->Write();
  wfile->Close();

  delete detdim;
  delete badch;
}


void fMode(int mode){
  //CC interaction
  if(mode==1)       cout << "CCQE"      << endl;
  else if(mode==11) cout << "CC1pi+"    << endl;
  else if(mode==12) cout << "CC1pi0"    << endl;
  else if(mode==13) cout << "CC1pi-"    << endl;
  else if(mode==16) cout << "CCcohpi"   << endl;
  else if(mode==17) cout << "CC1gamma"  << endl;
  else if(mode==21) cout << "CCmultipi" << endl;
  else if(mode==22) cout << "CC1eta"    << endl;
  else if(mode==23) cout << "CC1K"      << endl;
  else if(mode==26) cout << "CCDIS"     << endl;
  else if(mode==31) cout << "NC1pi0"    << endl;
  //NC interaction
  else if(mode==32) cout << "NC1pi0"    << endl;
  else if(mode==33) cout << "NC1pi-"    << endl;
  else if(mode==34) cout << "NC1pi+"    << endl;
  else if(mode==36) cout << "NCcohpi"   << endl;
  else if(mode==37) cout << "NC1gamma"  << endl;
  else if(mode==38) cout << "NC1gamma"  << endl;
  else if(mode==41) cout << "NCmultipi" << endl;
  else if(mode==42) cout << "NC1eta"    << endl;
  else if(mode==43) cout << "NC1eta"    << endl;
  else if(mode==44) cout << "NC1K0"     << endl;
  else if(mode==45) cout << "NC1K+"     << endl;
  else if(mode==46) cout << "NCDIS"     << endl;
  else if(mode==51) cout << "NCE"       << endl;
  else if(mode==52) cout << "NCE"       << endl;
  else              cout << mode        << endl;
}
