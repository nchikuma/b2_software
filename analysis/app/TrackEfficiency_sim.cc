#include "TrackEfficiency_sim.hxx"

FileStat_t fs;
void PrintUsage(){
  std::cout << "./TrackEfficiency_sim [options]" << std::endl;
  std::cout << "Either a pair of run/srun number or "
    << "a pair of input/output filename must be set." << endl;
  std::cout << "A mode must be selected between sand/cosmic muons." << std::endl;
  std::cout << "-r <runid>           : Set Run ID." << std::endl;
  std::cout << "-s <srunid/acqid>    : Set Sub-Run/Acq ID."  << std::endl;
  std::cout << "-i <filename.root>   : Input ROOT file, after 2d reconstruction." << std::endl;
  std::cout << "-m <mode>            : Select mode; " << std::endl;
  std::cout << "                       0->Sand muon, 1->Cosmic muon." << std::endl;
  std::cout << "                       2->MC Sand muon, 3->MC Cosmic muon." << std::endl;
  std::cout << "                       4->MC NEUT." << std::endl;
  std::cout << "-t <Nhits Threshold> : To set the threshold for num of his." << std::endl;
  std::cout << "                       *Default = 3 " << std::endl;
  std::cout << "-o <filename.root>   : Outpu ROOT file." << std::endl;
  std::cout << "-w                   : To analyze B2 WAGASCI cosmic trigger data." << endl;
  std::cout << "-b                   : To open all bad channels. Default: masked." << std::endl;
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
  int  LIM_NUM_HITS = 3;

  int runid  = -1;
  int srunid = -1;

  while ((c = getopt(argc, argv, "r:s:i:o:m:t:wh")) != -1){
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
        LIM_NUM_HITS = atoi(optarg);
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


  if(
      (mode!=0&&mode!=1&&mode!=2&&mode!=3&&mode!=4)||
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
          mc_muon,runid,srunid);
    }
    else if(mode==4){
      ifilename = Form("%s/ingrid_%08d_%04d_recon.root",
          mc_neut,runid,srunid);
    }
  }

  if(!orename){
    if(mode==0){
      ofilename = Form("%s/ingrid_%08d_%04d_trkeff%d_sim.root",
          dst_ingrid,runid,srunid,LIM_NUM_HITS);
    }
    else if(mode==1){
      if(wg_cosmic){
        ofilename = Form("%s/run_%05d_%03d_trkeff%d_sim.root",
            cosmic_wagasci,runid,srunid,LIM_NUM_HITS);
      }
      else{
        ofilename = Form("%s/ingrid_%08d_%04d_trkeff%d_sim.root",
            cosmic_data,runid,srunid,LIM_NUM_HITS);
      }
    }
    else if(mode==2){
      ofilename = Form("%s/ingrid_%08d_%04d_trkeff%d_sim.root",
          mc_sandmu,runid,srunid,LIM_NUM_HITS);
    }
    else if(mode==3){
      ofilename = Form("%s/ingrid_%08d_%04d_trkeff%d_sim.root",
          mc_muon,runid,srunid,LIM_NUM_HITS);
    }
    else if(mode==4){
      ofilename = Form("%s/ingrid_%08d_%04d_trkeff%d_sim.root",
          mc_neut,runid,srunid,LIM_NUM_HITS);
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
  std::string name, xtitle, modname, axisname;

  const int NMOD   =  2;
  TH1F *h1[2][NMOD][8];
  TH1F *h2[2][NMOD][8];
  TH1F *h3[2][NMOD][8];
  TH1F *h4[2][NMOD][8];
  TH1F *h5[2][NMOD][8];
  TH1F *h6[2][NMOD][8];
  for(int i=0;i<NMOD;i++){
    for(int j=0;j<2;j++){
      for(int k=0;k<8;k++){
        //if     (i==0){modname="INGWM";} //INGRID Water Module
        //else if(i==1){modname="WM"   ;} //WAGASCI
        if(i==0){modname="WM"   ;} //WAGASCI
        else if(i==1){modname="PM"   ;} //Proton Module
        axisname="Z";xtitle="Track angle from the Z axis [cos]";
        name  =  Form("h1_%s_axis%s_view%d_hit%d",
            modname.c_str(),axisname.c_str(),j,k+3);
        h1[j][i][k] =  new TH1F(name.c_str(),name.c_str(),100,0.,1.);
        h1[j][i][k] -> GetXaxis()->SetTitle(xtitle.c_str());
        name  =  Form("h2_%s_axis%s_view%d_hit%d",
            modname.c_str(),axisname.c_str(),j,k+3);
        h2[j][i][k] =  new TH1F(name.c_str(),name.c_str(),100,0.,1.);
        h2[j][i][k] -> GetXaxis()->SetTitle(xtitle.c_str());
        name  =  Form("h3_%s_axis%s_view%d_hit%d",
            modname.c_str(),axisname.c_str(),j,k+3);
        h3[j][i][k] =  new TH1F(name.c_str(),name.c_str(),180,-90.,90.);
        h3[j][i][k] -> GetXaxis()->SetTitle("Angular Difference [deg]");
        name  =  Form("h4_%s_axis%s_view%d_hit%d",
            modname.c_str(),axisname.c_str(),j,k+3);
        h4[j][i][k] =  new TH1F(name.c_str(),name.c_str(),120,0.,1.2);
        h4[j][i][k] -> GetXaxis()->SetTitle("Ratio of hits used in recon.");
        name  =  Form("h5_%s_axis%s_view%d_hit%d",
            modname.c_str(),axisname.c_str(),j,k+3);
        h5[j][i][k] =  new TH1F(name.c_str(),name.c_str(),1000,0.,1000.);
        h5[j][i][k] -> GetXaxis()->SetTitle("Start position difference [mm]");
        name  =  Form("h6_%s_axis%s_view%d_hit%d",
            modname.c_str(),axisname.c_str(),j,k+3);
        h6[j][i][k] =  new TH1F(name.c_str(),name.c_str(),1000,0.,1000.);
        h6[j][i][k] -> GetXaxis()->SetTitle("Hit distance from track [mm]");
      }
    }
  }

  //------- Event loop -------
  detdim   = new DetectorDimension();
  badch  = new INGRID_BadCh_mapping();
  badch -> set_BadCh(BadCh_Ana);
  //bool newcanv = true;

  for(int ievt=1;ievt<nevt;ievt++){

#ifdef DEBUG_TRKEFF_SIM
    cout << "================================" << endl;
    cout << "Event : " << ievt << endl;
#else
    if(ievt%100==0){
      std::cout << "Event # is " << ievt << std::endl;
    }
#endif

    tree  -> GetEntry(ievt);	

    // ========================
    // Sim Particle Track

    double offset_simpar[3] = { C_B2MotherPosX,C_B2MotherPosY,C_B2MotherPosZ };
    int nsimpar = evtsum->NSimParticles();
    double sim_slope[2]  = {0.,0.};
    double sim_intcpt[2] = {0.,0.};
    double sim_ipos[3] = {0.,0.,0.};
    double sim_fpos[3] = {0.,0.,0.};
    bool muon = false;
    for(int isimpar=0;isimpar<nsimpar;isimpar++){
      SimParticleSummary* simpar = evtsum->GetSimParticle(isimpar);
      if(abs(simpar->pdg)==13){
        muon = true;
        double momx = simpar->momentum[0];
        double momy = simpar->momentum[1];
        double momz = simpar->momentum[2];
        if(momz==0.){
          sim_slope[0]=1e+8;
          sim_slope[1]=1e+8;
        }
        else{
          sim_slope[0]=momx/momz;
          sim_slope[1]=momy/momz;
        }
        for(int i=0;i<3;i++){
          sim_ipos[i] = simpar->ipos[i]*10.-offset_simpar[i];
          sim_fpos[i] = simpar->fpos[i]*10.-offset_simpar[i];
        }
        sim_intcpt[0] = sim_ipos[0]-sim_ipos[2]*sim_slope[0];
        sim_intcpt[1] = sim_ipos[1]-sim_ipos[2]*sim_slope[1];
        break;
      }
    }
    if(!muon){continue;}

#ifdef DEBUG_TRKEFF_SIM
    cout << " --- SimParticle --- " << endl;
    cout << " slope={" << sim_slope[0] <<","<<sim_slope[1] << "}" << endl;
#endif

    int icyc = 4;
    double offset[3];
    for(int imod=0;imod<1;imod++){
      int c_mod = -1;
      if     (imod==0){c_mod=21;} //WAGASCI
      else if(imod==1){c_mod=15;} //INGRID Water Module
      else if(imod==2){c_mod=16;} //Proton Module
      if     (imod==0){
        offset[0] = C_B2WMPosX;
        offset[1] = C_B2WMPosY;
        offset[2] = C_B2WMPosZ;
      }
  
      // ===================
      // Hits
      int v_hitmod [2][100];
      int v_hitview[2][100];
      int v_hitpln [2][100];
      int v_hitch  [2][100];
      int nhit_around_simtrk[2] = {0,0};
      int nhit_reconstructed[2] = {0,0};
      int nmodhits = evtsum->NModHits(c_mod,icyc);
      for(int ihits =0; ihits< nmodhits; ihits++){
        HitSummary* hitsum = evtsum->GetModHit(ihits,c_mod,icyc);
        int    hitmod  = hitsum->mod;
        int    hitview = hitsum->view;
        int    hitpln  = hitsum->pln;
        int    hitch   = hitsum->ch;
        double hitpe   = hitsum->pe;
        int    nsimhit = hitsum->NSimHits();

        double hitx,hity,hitz;
        detdim->GetPosInMod(hitmod,hitpln,hitview,hitch,&hitx,&hity,&hitz);
        hitx += offset[0];
        hity += offset[1];
        hitz += offset[2];
        double hitpos[3] = {hitx,hity,hitz};
        int r_view = 1- hitview;
        double dist = 
          fabs(sim_slope[r_view]*hitpos[2]-hitpos[r_view]+sim_intcpt[r_view])
            /sqrt(1.+sim_slope[r_view]*sim_slope[r_view]);
        bool muon_hit = false;
        int  hitpdg = -1;
        for(int isimhit=0;isimhit<nsimhit;isimhit++){
          SimHitSummary* simhitsum = hitsum->GetSimHit(isimhit);
          hitpdg = simhitsum->pdg; 
          if(abs(hitpdg)==13){muon_hit=true;}
        }

#ifdef DEBUG_TRKEFF_SIM
        cout
          << " mod="  << hitmod
          << " view=" << hitview
          << " pln="  << hitpln
          << " ch="   << hitch
          << " pe="   << hitpe
          << " nsim=" << nsimhit
          << " dist=" << dist
          << " pdg="  << hitpdg
          << endl;
#endif



        //if(!muon_hit){continue;}

        if(badch->is_BadCh(hitmod,hitview,hitpln,hitch)){
#ifdef DEBUG_TRKEFF_SIM
            cout << "              ///bad." << endl;
#endif
          continue;
        }
        if(hitpe<THRES_PE){
#ifdef DEBUG_TRKEFF_SIM
            cout << "              ///low pe." << endl;
#endif
          continue;
        }


        if(dist<MAX_HIT_DIST)
        {
          if(nhit_around_simtrk[hitview]<100){
            v_hitmod [hitview][nhit_around_simtrk[hitview]]=hitmod ;
            v_hitview[hitview][nhit_around_simtrk[hitview]]=hitview;
            v_hitpln [hitview][nhit_around_simtrk[hitview]]=hitpln ;
            v_hitch  [hitview][nhit_around_simtrk[hitview]]=hitch  ;
            nhit_around_simtrk[hitview]++;
#ifdef DEBUG_TRKEFF_SIM
            cout << "             --> added." << endl;
#endif
          }
        }

      }//ihits

      // ===================
      // Two Dim Recon
      double slope  [2] = { -1.,-1.};
      double intcpt [2] = { -1.,-1.};
      double startxy[2] = { -1.,-1.};
      double startz [2] = { -1.,-1.};
      double endxy  [2] = { -1.,-1.};
      double endz   [2] = { -1.,-1.};
      bool reconstructed [2] = {false,false};
      for(int iview=0;iview<2;iview++){
        if(nhit_around_simtrk[iview]<3){continue;}

        int n2drec = evtsum->NModTwoDimRecons(c_mod,icyc,iview);
#ifdef DEBUG_TRKEFF_SIM
        cout << " Num 2D recon = " << n2drec << endl;
        if(n2drec>1){
          cout <<"   -- track ang diff : "
            << (atan(evtsum->GetModTwoDimRecon(1,c_mod,icyc,iview)->slope)
                -atan(evtsum->GetModTwoDimRecon(0,c_mod,icyc,iview)->slope))*180./PI << endl;
        }
#endif 

        if(n2drec<1){continue;}
        //if(n2drec!=1){continue;}


        TwoDimReconSummary* reconsum = evtsum->GetModTwoDimRecon(0,c_mod,icyc,iview);

        slope  [1-iview] = reconsum->slope;
        intcpt [1-iview] = reconsum->intcpt;
        startxy[1-iview] = reconsum->startxy += offset[1-iview];
        startz [1-iview] = reconsum->startz  += offset[2];
        endxy  [1-iview] = reconsum->endxy   += offset[1-iview];
        endz   [1-iview] = reconsum->endz    += offset[2];


        if(fabs(atan(slope[1-iview])-atan(sim_slope[1-iview]))<20.){
          reconstructed[iview]=true;
        }

        int nrechit = reconsum->Nhits();
        for(int irechit=0;irechit<nrechit;irechit++){
          HitSummary* rechit = reconsum->GetHit(irechit);
          int rechitmod  = rechit->mod;
          int rechitview = rechit->view;
          int rechitpln  = rechit->pln;
          int rechitch   = rechit->ch;

          double rechitx,rechity,rechitz;
          detdim->GetPosInMod(rechitmod,rechitpln,rechitview,rechitch,&rechitx,&rechity,&rechitz);
          double rechitpos[3] = {rechitx,rechity,rechitz};
          double dist = 
            fabs(slope[1-iview]*rechitpos[2]-rechitpos[1-iview]+intcpt[1-iview])
            /sqrt(1.+slope[1-iview]*slope[1-iview]);
          for(int ithres=3;ithres<=nhit_around_simtrk[iview];ithres++){
            h6[iview][imod][ithres-3]->Fill(dist);
          }

#ifdef DEBUG_TRKEFF_SIM
              cout 
                << " mod="  << rechitmod
                << " view=" << rechitview 
                << " pln="  << rechitpln 
                << " ch="   << rechitch
                << endl;
#endif


          for(int isimhit=0;isimhit<nhit_around_simtrk[iview];isimhit++){
            if(
                v_hitmod [iview][isimhit]==rechitmod &&
                v_hitview[iview][isimhit]==rechitview&&
                v_hitpln [iview][isimhit]==rechitpln &&
                v_hitch  [iview][isimhit]==rechitch  )
            {
#ifdef DEBUG_TRKEFF_SIM
              cout << "   ---> recon.d" << endl;
#endif 
              nhit_reconstructed[iview]++;
              break;
            }
          }
        }
      } //iview

#ifdef DEBUG_TRKEFF_SIM
      cout << " --- Recon ---" << endl;
      cout << " slope = {" << slope[0]<<","<<slope[1]<<"}" << endl;
      cout << " --- Recon hits --- "<< endl;
      cout << " num hits : {" << nhit_around_simtrk[0] 
        <<","<< nhit_around_simtrk[1]<< "}" << endl;
      cout << " num recon hits : {" << nhit_reconstructed[0]
        <<","<< nhit_reconstructed[1] << "}" << endl;
#endif


      for(int iview=0;iview<2;iview++){
        double cos = 1./sqrt(1.+sim_slope[1-iview]*sim_slope[1-iview]);
        for(int thres_nhit=3;thres_nhit<11;thres_nhit++){
          if( (nhit_around_simtrk[iview]==thres_nhit ) ||
              (thres_nhit==10&&nhit_around_simtrk[iview]>=10) )
          {
            h1[iview][imod][thres_nhit-3]->Fill(cos);
            if(reconstructed[iview]){
              h2[iview][imod][thres_nhit-3]->Fill(cos);
              double diff_ang = (atan(slope[1-iview]) - atan(sim_slope[1-iview]))*180./PI;
              if     (diff_ang> 90.)diff_ang-=180.;
              else if(diff_ang<-90.)diff_ang+=180.;
              h3[iview][imod][thres_nhit-3]->Fill(diff_ang);

              double ratio_nhit = 
                ((double)nhit_reconstructed[iview])/((double)nhit_around_simtrk[iview]);

#ifdef DEBUG_TRKEFF_SIM
              cout << " --- Used Hit Ratio ---" << endl;
              cout 
                << " view=" << iview 
                << " thres=" << thres_nhit 
                << " ratio = " << ratio_nhit 
                << endl;
#endif
              h4[iview][imod][thres_nhit-3]->Fill(ratio_nhit);
              double dist1 = 
                sqrt(pow(sim_ipos[2]-startz[1-iview],2)+pow(sim_ipos[1-iview]-startxy[1-iview],2));
              double dist2 = 
                sqrt(pow(sim_ipos[2]-endz[1-iview],2)+pow(sim_ipos[1-iview]-endxy[1-iview],2));
              if(dist2<dist1){dist1=dist2;}
              h5[iview][imod][thres_nhit-3]->Fill(dist1);
            }
          }
        }
      }
    } //for(imod)
  } //ievt

  ofile->cd();
  for(int i=0;i<NMOD;i++){
    for(int j=0;j<2;j++){
      for(int k=0;k<8;k++){
        h1[j][i][k] -> Write();
        h2[j][i][k] -> Write();
        h3[j][i][k] -> Write();
        h4[j][i][k] -> Write();
        h5[j][i][k] -> Write();
        h6[j][i][k] -> Write();
      }
    }
  }
  ofile->Write();
  ofile->Close();

  delete badch;
  delete detdim;
  return 0;
}
