#include "TwoDimRecon.hxx"


vector<int>  bad_mod;
vector<int>  bad_pln;
vector<int>  bad_view;
vector<int>  bad_ch;
void fMode(int mode);

void PrintUsage(){
  cout << "Usage:" << endl;
  cout << "./TwoDimRecon [options]" << endl;
  cout << "Either a pair of run/srun number or a pair of input/output filename must be set." << endl;
  cout << " -r <int>  : Run Number    " << endl;
  cout << " -s <int>  : Sub Run Number" << endl;
  cout << " -f <char> : Input filename, if you directly indicate it."  << endl;
  cout << " -o <char> : Output filename, if you directly indicate it." << endl;
  cout << " -i <int>  : To set an event number for starting analysis in the middle." << endl;
  cout << " -a        : To analyze all cycles. Only 9 cycles of 4 - 12 are anlaysis by defualt." << endl;
  cout << " -c        : To analyze 14-15 cycles for cosmic trigger." << endl;
  cout << " -w        : To analyze B2 WAGASCI cosmic trigger data." << endl;
  cout << " -b        : To open all bad channel masked. Bad channels are masked by defualt." << endl;
  cout << " -n <int>  : To implement random noise" << endl;
  cout << "      0 -> WAGASCI : 1 hits/cycle"      << endl;
  cout << "      1 -> PM :  5 hits/cycle"          << endl;
  cout << "      2 -> PM : 10 hits/cycle"          << endl;
  cout << "      3 -> PM : 15 hits/cycle"          << endl;
  cout << "      4 -> PM : 20 hits/cycle"          << endl;
  cout << "      5 -> PM : 25 hits/cycle"          << endl;
  cout << " -x        : To implement cross talk b/w scintillators" << endl;
  cout << " -t <int>  : To adjust the hit pe threshold." << endl;
  cout << "      0 -> WAGASCI : 1.5pe" << endl;
  cout << "      1 -> WAGASCI : 2.5pe" << endl;
  cout << "      2 -> WAGASCI : 3.5pe" << endl;
  cout << "      3 -> PM      : 2.5pe" << endl;
  cout << "      4 -> PM      : 3.5pe" << endl;
  cout << "      5 -> PM      : 4.5pe" << endl;
  cout << " -d        : To display."        << endl;
  cout << " -h        : To show this help." << endl;
}


long fcTime;
int main(int argc,char *argv[]){

  char buff1[200];
  char Output[300];

  Int_t run_number     = -1;
  Int_t sub_run_number = -1;

  Int_t c    = -1;
  Int_t Scyc =  4;
  Int_t Ncyc = 13;
  Int_t Nmod = 24;
  Int_t Nini = 0;

  bool rename   = false;
  bool renamef  = false;
  bool disp     = false;
  bool cosmic   = false;
  bool wg_cosmic= false;
  bool BadChAna = true ;
  bool OneEvent = false;

  int  Num_RandNoise = -1;
  bool Crosstalk = false;
  int  PE_THRES   = -1;

  while ( (c = getopt(argc, argv, "r:s:f:o:i:n:t:xcwadbh")) != -1) {
    switch(c){
      case 'r':
        run_number=atoi(optarg);
        break;
      case 's':
        sub_run_number=atoi(optarg);
        break;
      case 'f':
        sprintf(buff1,"%s",optarg);
        run_number=0;
        renamef = true;
        break;
      case 'o':
       sprintf(Output,"%s",optarg);
        rename  = true;
        break;
      case 'i':
        Nini=atoi(optarg);
        OneEvent = true;
        break;
      case 'n':
        Num_RandNoise=atoi(optarg);
        break;
      case 'x':
        Crosstalk = true;
        break;
      case 't':
        PE_THRES = atoi(optarg);
        break;
      case 'a':
        Scyc = 0;
        Ncyc = 23;
        break;
      case 'c':
        Scyc = 14;
        Ncyc = 16;
        cosmic = true;
        break;
      case 'w':
        Scyc = 14;
        Ncyc = 16;
        wg_cosmic = true;
        break;
      case 'd':
        disp = true;
        break;
      case 'b':
        BadChAna = false;
        break;
      case 'h':
        PrintUsage();
        exit(0);
    }
  }

  if(
      ((run_number<0||sub_run_number<0)&&(!rename||!renamef))||
      (!cosmic&&wg_cosmic))
  {
    PrintUsage(); exit(0);
  }

  if(!renamef && wg_cosmic){
    sprintf(buff1,"%s/run_%05d_%03d_inglib.root",cosmic_wagasci,run_number,sub_run_number);
  }
  else if(!renamef && cosmic){
    sprintf(buff1,"%s/ingrid_%08d_%04d_calib.root",cosmic_data,run_number,sub_run_number);
  }

  else if(!renamef){
    sprintf(buff1,"%s/ingrid_%08d_%04d_bsd.root",dst_ingrid,run_number,sub_run_number);
    //sprintf(buff1,"%s/ingrid_%08d_%04d_spill.root",dst_ingrid,run_number,sub_run_number);
    //sprintf(buff1,"%s/ingrid_%08d_%04d_spill.root",dst_ingrid2,run_number,sub_run_number);
  }


  //######## read root file #############################
  //#####################################################
  cout<<"reading "<<buff1<<"....."<<endl;
  FileStat_t fs;
  if(gSystem->GetPathInfo(buff1,fs)){
    cout<<"Cannot open file "<<buff1<<endl;
    exit(1);
  }

  EventSummary* evt       =  new EventSummary();
  TFile*            rfile =  new TFile(buff1,"read");
  TTree*             tree =  (TTree*)rfile -> Get("tree");
  TBranch*          EvtBr =  tree->GetBranch("fDefaultReco.");
  EvtBr                   -> SetAddress(&evt);
  tree                    -> SetBranchAddress("fDefaultReco.", &evt);

  int                nevt = (int)tree -> GetEntries();
  cout << "Total # of events = " << nevt <<endl;


  //#### make rootfile after analysis #####
  //#######################################
  if(!rename && wg_cosmic){
    sprintf(Output, "%s/run_%05d_%03d_recon.root",cosmic_wagasci,run_number,sub_run_number); 
  }
  else if(!rename && cosmic){
    sprintf(Output, "%s/ingrid_%08d_%04d_recon.root",cosmic_data,run_number,sub_run_number); 
  }
  else if( !rename ){
    //sprintf(Output, "%s/ingrid_%08d_%04d_recon.root",dst_ingrid,run_number,sub_run_number); 
    sprintf(Output, "%s/ingrid_%08d_%04d_recon.root",dst_ingrid2,run_number,sub_run_number); 
  }

  cout << "Output filename : "<< Output << endl;

  TFile*            wfile =  new TFile(Output, "recreate");
  TTree*            wtree =  new TTree("tree","tree");
  wtree                   -> SetMaxTreeSize(5000000000);
  EventSummary* wsummary  =  new EventSummary(); 
  wtree                   -> Branch   ("fDefaultReco.","EventSummary", 
      &wsummary,  64000,  99);


  //INGRID_BadCh_mapping *badch = new INGRID_BadCh_mapping();
  badch = new INGRID_BadCh_mapping();
  badch->set_BadCh(BadChAna);

  detdim = new DetectorDimension();

  HitSummary*         hitsum;
  //SimVertexSummary*   simver;
  TwoDimReconSummary* recon = new TwoDimReconSummary();
  Hit                 hit;
  if(OneEvent){nevt=Nini+1;}
  for(int ievt=Nini; ievt<nevt; ievt++){

#ifdef DEBUG_TWODIMRECON
    cout << "==============================================================================================" << endl;
    cout << "==============================================================================================" << endl;
    cout << "==============================================================================================" << endl;
    cout << "==============================================================================================" << endl;
    cout << "analyze event# " << ievt<<endl;
#else
    if(ievt%100==0) 
      cout << "analyze event# " << ievt<<endl;
#endif

    wsummary -> Clear();
    evt      -> Clear();
    tree     -> GetEntry(ievt);

    // ===============================================
    // ============ Threshold Adjust =================

    if(PE_THRES!=-1){
      if     (PE_THRES==0){hitpe_threshold_WM=1.5;}
      else if(PE_THRES==1){hitpe_threshold_WM=2.5;}
      else if(PE_THRES==2){hitpe_threshold_WM=3.5;}
      else if(PE_THRES==3){hitpe_threshold_PM=2.5;}
      else if(PE_THRES==4){hitpe_threshold_PM=3.5;}
      else if(PE_THRES==5){hitpe_threshold_PM=4.5;}
      else{
        cout << "Wrong input.. for PE Threshold" << endl;
        return 1;
      }
    }


    for( int cyc=Scyc; cyc<Ncyc; cyc++ ){  //### Cycle Loop
      for( int mod=0; mod<Nmod; mod++ ){   //### Module Loop

        bool isWAGASCI = (mod==MOD_ONAXIS_WM||mod==MOD_B2_WM);
        bool isPM      = (mod==MOD_PM);
        bool isINGRID  = ((mod>=0&&mod<NUMINGMOD)||mod==MOD_B2_INGRID);

        allhit.clear();
        int nhit = evt -> NModHits(mod, cyc);

#ifdef DEBUG_TWODIMRECON
        //std::cout << "-------------------------------------------------" << std::endl;
        //std::cout << "MOD=" << mod  << " nhits=" <<  nhit << std::endl;
#endif

        for(int i=0; i<nhit; i++){
          hitsum   = (HitSummary*) (evt -> GetModHit(i, mod, cyc) );
          if(hitsum->cyc==-2) continue;
          if(badch->is_BadCh(mod,hitsum->view,hitsum->pln,hitsum->ch)){continue;}

          if(isWAGASCI){
            if(hitsum->pe + hitsum->pe_cross < hitpe_threshold_WM){continue;}
#ifdef MASK_GRID_CHANNEL
            if(hitsum->ch>=40){continue;}
#endif
          }
          else if(isPM){
            if(hitsum->pe < hitpe_threshold_PM){continue;}
          }
          else if(isINGRID){
            if(hitsum -> pe < hitpe_threshold_ING){continue;}
          }
          else{
            continue;
          }

          hitsum -> addbasicrecon = false;

          hit.mod      = mod;
          hit.id       = i;
          hit.pe       = hitsum -> pe;
          hit.lope     = hitsum -> lope;
          hit.time     = (Long_t)hitsum -> time;
          hit.view     = hitsum -> view;
          hit.pln      = hitsum -> pln;
          hit.ch       = hitsum -> ch;
          hit.pe_cross = hitsum -> pe_cross;

#ifdef DEBUG_TWODIMRECON
          std::cout
            << " mod="  << hit.mod
            << " view=" << hit.view
            << " pln="  << hit.pln
            << " ch="   << hit.ch
            << " pe="   << hit.pe
            << " time=" << hit.time
            << std::endl;
#endif


          allhit.push_back(hit);
          allhit_for_disp.push_back(hit);
        }

        // ===============
        // Add Random Hits 
        //
        if(cyc==4&&Num_RandNoise>=0&&Num_RandNoise<=5){
          srand(time(NULL));
          Add_RandomNoise(Num_RandNoise,mod,cyc,evt,5.0);
        }
        
        // ==========================
        // Add Scintillator Crosstalk
        //
        if(cyc==4&&mod==MOD_B2_WM&&Crosstalk){ 
          srand(time(NULL));
          Add_Crosstalk(0.01,mod,cyc,evt); 
        }


        if(mod==MOD_B2_WM||mod==MOD_ONAXIS_WM){
          if((int)allhit.size()<nhit_threshold_WM){continue;}
        }else{
          if((int)allhit.size()<nhit_threshold){continue;}
        }

        fSortTime(allhit);
        while(fFindTimeClster(allhit, hitcls, fcTime, cosmic)){

#ifdef DEBUG_TWODIMRECON
          std::cout << "Found a time cluster" << std::endl;
          int hitclssize = hitcls.size();
          for(int itmp=0;itmp<hitclssize;itmp++){
            std::cout
              << " mod="  << hitcls[itmp].mod
              << " view=" << hitcls[itmp].view
              << " pln="  << hitcls[itmp].pln
              << " ch="   << hitcls[itmp].ch
              << " pe="   << hitcls[itmp].pe
              << " time=" << hitcls[itmp].time
              << std::endl;
          }
#endif

          int nactnum;
          recon->nactpln = fNactpln(mod);
          nactnum = fNactpln(mod);
#ifdef DEBUG_TWODIMRECON
          std::cout << "Nactpln=" << nactnum << std::endl;
#endif
          if(isINGRID){
            if(recon->nactpln < nactpln_threshold){continue;}
          }
          recon->layerpe = fLayerpe(mod);
#ifdef DEBUG_TWODIMRECON
          std::cout << "LayerPe=" << recon->layerpe << std::endl;
#endif
          if(!isWAGASCI){
            if(recon->layerpe < layerpe_threshold){continue;}
          }

          hitcls_for_joint=hitcls;

          int alltrack_size_bak=0;
          int Nconnect = 0;
          if(isWAGASCI){ fTracking_WM(mod); }
          else         { fTracking   (mod); }
          alltrack_size_bak=(int)alltrack.size();
          while(fConnectTracks(mod)){
            Nconnect++;
            if(Nconnect>alltrack_size_bak) break;
          }

          for(int i=0;i<(int)alltrack.size();i++){
            recon->Clear();
            recon -> view          = alltrack[i].view;
            recon -> startpln      = alltrack[i].ipln;
            recon -> endpln        = alltrack[i].fpln;
            recon -> startxy       = alltrack[i].ixy;
            recon -> endxy         = alltrack[i].fxy;
            recon -> startz        = alltrack[i].iz;
            recon -> endz          = alltrack[i].fz;
            recon -> angle         = alltrack[i].ang;
            recon -> slope         = alltrack[i].slope;
            recon -> intcpt        = alltrack[i].intcpt;
            recon -> vetowtracking = alltrack[i].veto;
            recon -> edgewtracking = alltrack[i].edge;
            recon -> modfc         = alltrack[i].stop;
            recon -> vetodist      = alltrack[i].vetodist;
            //recon -> totpe        = alltrack[i].totpe;
            //recon -> isohit       = alltrack[i].isohit;
            //for(TMP=0;TMP<4;TMP++)
            //recon -> pdg[TMP]       = alltrack[i].pdg[TMP];
            for(int NHIT=0;NHIT<(int)alltrack[i].hitid.size();NHIT++){
              hitsum   = (HitSummary*) (evt -> GetModHit(alltrack[i].hitid[NHIT], mod, cyc) );
              //hitsum -> isohit   = alltrack[i].isohit[NHIT];
              //hitsum -> gocosmic = true;
              recon -> AddHit(hitsum);
            }
            recon -> clstime      = fcTime;
            recon -> hitmod       = mod;
            recon -> hitcyc       = cyc;
            evt   -> AddModTwoDimRecon( recon , mod, cyc, alltrack[i].view);
          }
          hitcls.clear(); //Reset hit clster
        }//Find Time Clster
      }//mod
    }//cyc
    allhit_for_disp.clear();

    wsummary = evt;
    wtree -> Fill();
  }//Event Loop


  //######## Write and Close ####################
  wfile -> cd();
  wtree  -> Write();
  wfile  -> Write();
  wfile  -> Close();

  delete detdim;
  delete badch;
  delete recon;

}

void fMode(int mode){
  if(mode==1)
    cout<<"CCQE"<<endl;
  else if(mode==11)
    cout<<"CC1pi+"<<endl;
  else if(mode==12)
    cout<<"CC1pi0"<<endl;
  else if(mode==13)
    cout<<"CC1pi-"<<endl;
  else if(mode==16)
    cout<<"CCcohpi"<<endl;
  else if(mode==17)
    cout<<"CC1gamma"<<endl;
  else if(mode==21)
    cout<<"CCmultipi"<<endl;
  else if(mode==22)
    cout<<"CC1eta"<<endl;
  else if(mode==23)
    cout<<"CC1K"<<endl;
  else if(mode==26)
    cout<<"CCDIS"<<endl;

  else if(mode==31)
    cout<<"NC1pi0"<<endl;
  else if(mode==32)
    cout<<"NC1pi0"<<endl;
  else if(mode==33)
    cout<<"NC1pi-"<<endl;
  else if(mode==34)
    cout<<"NC1pi+"<<endl;
  else if(mode==36)
    cout<<"NCcohpi"<<endl;
  else if(mode==37)
    cout<<"NC1gamma"<<endl;
  else if(mode==38)
    cout<<"NC1gamma"<<endl;
  else if(mode==41)
    cout<<"NCmultipi"<<endl;
  else if(mode==42)
    cout<<"NC1eta"<<endl;
  else if(mode==43)
    cout<<"NC1eta"<<endl;
  else if(mode==44)
    cout<<"NC1K0"<<endl;
  else if(mode==45)
    cout<<"NC1K+"<<endl;
  else if(mode==46)
    cout<<"NCDIS"<<endl;
  else if(mode==51)
    cout<<"NCE"<<endl;
  else if(mode==52)
    cout<<"NCE"<<endl;


}
