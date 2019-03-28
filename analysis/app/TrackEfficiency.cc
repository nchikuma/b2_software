#include "TrackEfficiency.hxx"

FileStat_t fs;
void PrintUsage(){
  std::cout << "./TrackEfficiency [options]" << std::endl;
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
          mc_cosmic,runid,srunid);
    }
    else if(mode==4){
      ifilename = Form("%s/ingrid_%08d_%04d_recon.root",
          mc_neut,runid,srunid);
    }
  }

  if(!orename){
    if(mode==0){
      ofilename = Form("%s/ingrid_%08d_%04d_trkeff%d.root",
          dst_ingrid,runid,srunid,LIM_NUM_HITS);
    }
    else if(mode==1){
      if(wg_cosmic){
        ofilename = Form("%s/run_%05d_%03d_trkeff%d.root",
            cosmic_wagasci,runid,srunid,LIM_NUM_HITS);
      }
      else{
        ofilename = Form("%s/ingrid_%08d_%04d_trkeff%d.root",
            cosmic_data,runid,srunid,LIM_NUM_HITS);
      }
    }
    else if(mode==2){
      ofilename = Form("%s/ingrid_%08d_%04d_trkeff%d.root",
          mc_sandmu,runid,srunid,LIM_NUM_HITS);
    }
    else if(mode==3){
      ofilename = Form("%s/ingrid_%08d_%04d_trkeff%d.root",
          mc_cosmic,runid,srunid,LIM_NUM_HITS);
    }
    else if(mode==4){
      ofilename = Form("%s/ingrid_%08d_%04d_trkeff%d.root",
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

  const int NMOD   =  3; //
  TH1F *h1[NMOD];
  TH1F *h2[NMOD];
  for(int i=0;i<NMOD;i++){
    if     (i==0){modname="INGWM";} //INGRID Water Module
    else if(i==1){modname="WM"   ;} //WAGASCI
    else if(i==2){modname="PM"   ;} //Proton Module
    axisname="Z";xtitle="Track angle from the Z axis [cos]";
    name  =  Form("h1_%s_axis%s",modname.c_str(),axisname.c_str());
    h1[i] =  new TH1F(name.c_str(),name.c_str(),100,-1.,1.);
    h1[i] -> GetXaxis()->SetTitle(xtitle.c_str());
    name  =  Form("h2_%s_axis%s",modname.c_str(),axisname.c_str());
    h2[i] =  new TH1F(name.c_str(),name.c_str(),100,-1.,1.);
    h2[i] -> GetXaxis()->SetTitle(xtitle.c_str());
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
  else if(mode==4){start_cyc= 4;stop_cyc= 4;}

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
        int t_mod = -1;
        int max_diff_cyc = 0;
        //double THRES_DIST = -1.; 
        if     (imod==0){c_mod=15;t_mod= 3;max_diff_cyc=0;} //INGRID Water Module
        else if(imod==1){c_mod=21;t_mod=14;max_diff_cyc=1;} //WAGASCI
        else if(imod==2){c_mod=16;t_mod=14;max_diff_cyc=0;} //Proton Module
        double offset[3] = {0.,0.,0.};
        if     (imod==0){
          offset[0] = -C_PMMotherPosX;
          offset[1] = -C_PMMotherPosY;
          offset[2] = -C_PMMotherPosZ;
        }
        else if(imod==1){
          offset[0] = C_B2INGPosX-C_B2WMPosX;
          offset[1] = C_B2INGPosY-C_B2WMPosY;
          offset[2] = C_B2INGPosZ-C_B2WMPosZ;
        }
        else if(imod==2){
          offset[0] = C_B2INGPosX-C_B2CHPosX;
          offset[1] = C_B2INGPosY-C_B2CHPosY;
          offset[2] = C_B2INGPosZ-C_B2CHPosZ;
        }
  
        int n2drec[2] = {0,0};
        for(int iview=0;iview<2;iview++){
          n2drec[iview] = evtsum->NModTwoDimRecons(t_mod,icyc,iview);
        }
        if(n2drec[0]!=1||n2drec[1]!=1){continue;}


        double slopex=-1.,slopey=-1.;
        double intcptx=-1.,intcpty=-1.;
        bool ingrid_recon = false;
        int endpln_tmp = -1;
        for(int iview=0;iview<2;iview++){
          TwoDimReconSummary* reconsum = evtsum->GetModTwoDimRecon(0,t_mod,icyc,iview);
          int   startpln = reconsum->startpln;
          int   endpln   = reconsum->endpln;
          double slope   = reconsum->slope;
          double intcpt  = reconsum->intcpt;
          if(iview==0){slopey=slope;intcpty=intcpt;}
          else        {slopex=slope;intcptx=intcpt;}

          if(startpln!=0){ingrid_recon=false;break;}
          else{
            if(iview==0){endpln_tmp=endpln;}
            else{
              if(abs(endpln-endpln_tmp)>1){ continue; }
              else{
                ingrid_recon = true;
              }
            }
          }
        } //iview
        if(!ingrid_recon){continue;}


        int tmp             [2] = {0,0};
        int nhit_around_reco[2] = {0,0};
        int tmp_hitmod [2][100]; int v_hitmod [2][100];
        int tmp_hitview[2][100]; int v_hitview[2][100];
        int tmp_hitpln [2][100]; int v_hitpln [2][100];
        int tmp_hitch  [2][100]; int v_hitch  [2][100];
        double tmp_hitpe  [2][100];
        double tmp_hittime[2][100];
        int used_cyc[2] = {-1,-1};
        for(int diff_cyc=0;diff_cyc<=max_diff_cyc;diff_cyc++){
          int cyc = icyc+diff_cyc;
          int nmodhits = evtsum->NModHits(c_mod,cyc);
          for(int ihits =0; ihits< nmodhits; ihits++){
            HitSummary* hitsum = evtsum->GetModHit(ihits,c_mod,cyc);
            int    hitmod  = hitsum->mod;
            int    hitview = hitsum->view;
            int    hitpln  = hitsum->pln;
            int    hitch   = hitsum->ch;
            double hitpe   = hitsum->pe;
            double hittime = hitsum->time;
            
            if(badch->is_BadCh(hitmod,hitview,hitpln,hitch)){continue;}
            if(hitpe<hitpe_threshold_WM){ continue; }

            double hitx,hity,hitz;
            detdim->GetPosInMod(hitmod,hitpln,hitview,hitch,&hitx,&hity,&hitz);
            if(hitview==0){
              if(
                  fabs(slopey*(hitz-offset[2])-(hity-offset[1])+intcpty)
                  /sqrt(1.+slopey*slopey)<MAX_HIT_DIST)
              {
                if(tmp[0]<100){
                  tmp_hitmod [0][tmp[0]]=hitmod ;
                  tmp_hitview[0][tmp[0]]=hitview;
                  tmp_hitpln [0][tmp[0]]=hitpln ;
                  tmp_hitch  [0][tmp[0]]=hitch  ;
                  tmp_hittime[0][tmp[0]]=hittime;
                  tmp_hitpe  [0][tmp[0]]=hitpe  ;
                  tmp[0]++;
                }
              }
            }
            else{
              if(
                  fabs(slopex*(hitz-offset[2])-(hitx-offset[0])+intcptx)
                  /sqrt(1.+slopex*slopex)<MAX_HIT_DIST)
              {
                if(tmp[1]<100){
                  tmp_hitmod [1][tmp[1]]=hitmod ;
                  tmp_hitview[1][tmp[1]]=hitview;
                  tmp_hitpln [1][tmp[1]]=hitpln ;
                  tmp_hitch  [1][tmp[1]]=hitch  ;
                  tmp_hittime[1][tmp[1]]=hittime;
                  tmp_hitpe  [1][tmp[1]]=hitpe  ;
                  tmp[1]++;
                }
              }
            }
          }//ihits

#ifdef USE_TIME_CLUSTER
          //Find time cluster
          double TIME_WIDTH = 50.;
          if(imod==1){TIME_WIDTH=580.;} //WAGASCI
          for(int iview=0;iview<2;iview++){
            int    cls_id  [100][100];
            double cls_time[100];
            double cls_pe  [100];
            int    cls_size[100];
            int    usedhit [100];
            for(int i=0;i<100;i++){usedhit[i]=0;cls_pe[i]=0.;cls_size[i]=0;}
            int    cls    = 0;
            int    clshit = 0;            
            while(true){
              for(int i=0;i<tmp[iview];i++){
                if(usedhit[i]==1){continue;}
                if(cls_pe[cls]<tmp_hitpe[iview][i]){
                  cls_pe  [cls] = tmp_hitpe  [iview][i];
                  cls_time[cls] = tmp_hittime[iview][i];
                }
              }
              for(int i=0;i<tmp[iview];i++){
                if(usedhit[i]==1){continue;}
                if(fabs(tmp_hittime[iview][i]-cls_time[cls])<TIME_WIDTH){
                  usedhit[i]=1;
                  cls_id[cls][clshit] = i;
                  cls_size[cls]++;
                  clshit++;
                }
              }
              bool complete=true;
              for(int j=0;j<tmp[iview];j++){
                if(usedhit[j]==0){complete=false;break;}
              }
              if(complete){break;}
              cls++;
              clshit=0;
            }
            int max_clssize = 0;
            int max_cls_id  = -1;
            for(int i=0;i<=cls;i++){
              if(max_clssize<cls_size[i]){
                max_clssize = cls_size[i];
                max_cls_id=i;
              }
            }
            for(int i=0;i<tmp[iview];i++){
              bool clustered = false;
              for(int j=0;j<max_clssize;j++){
                if(i==cls_id[max_cls_id][j]){clustered=true;break;}
              }
              if(!clustered){
                for(int k=i;k<tmp[iview]-1;k++){
                  tmp_hitmod [iview][k]=tmp_hitmod [iview][k+1];
                  tmp_hitview[iview][k]=tmp_hitview[iview][k+1];
                  tmp_hitpln [iview][k]=tmp_hitpln [iview][k+1];
                  tmp_hitch  [iview][k]=tmp_hitch  [iview][k+1];
                }
                tmp[iview]--;
              }
            }
          }
          // Completed finding time cluster
#endif

          if(tmp[0]>nhit_around_reco[0]){
            nhit_around_reco[0]=tmp[0];
            used_cyc[0] = cyc;
            for(int i=0;i<100;i++){
              v_hitmod [0][i] = tmp_hitmod [0][i];
              v_hitview[0][i] = tmp_hitview[0][i];
              v_hitpln [0][i] = tmp_hitpln [0][i];
              v_hitch  [0][i] = tmp_hitch  [0][i];
            }
          }
          if(tmp[1]>nhit_around_reco[1]){
            nhit_around_reco[1]=tmp[1];
            used_cyc[1] = cyc;
            for(int i=0;i<100;i++){
              v_hitmod [1][i] = tmp_hitmod [1][i];
              v_hitview[1][i] = tmp_hitview[1][i];
              v_hitpln [1][i] = tmp_hitpln [1][i];
              v_hitch  [1][i] = tmp_hitch  [1][i];
            }
          }
          tmp[0] = 0;
          tmp[1] = 0;
        }

#ifdef DEBUG_TRKEFF
        if(c_mod==21)
        std::cout 
          << "ievt=" << ievt 
          << " mod=" << c_mod
          << " ingrid_recon : "
          << " nhit[0]=" << nhit_around_reco[0]
          << " nhit[1]=" << nhit_around_reco[1]
          << std::endl;
#endif

        double norm = sqrt(1.+slopex*slopex+slopey*slopey);
        //double cos  = 1./norm;

        for(int iview=0;iview<2;iview++){
          if(nhit_around_reco[iview]<LIM_NUM_HITS){continue;}
          bool target_recon = false;
          for(int diff_cyc=0;diff_cyc<=max_diff_cyc;diff_cyc++){
            int cyc = icyc+diff_cyc;
            if(cyc!=used_cyc[iview]){continue;}
            int n2dreco = evtsum->NModTwoDimRecons(c_mod,cyc,iview);
#ifdef DEBUG_TRKEFF
            std::cout << "   -n2dreco=" << n2dreco << std::endl;
#endif
            for(int irec=0;irec<n2dreco;irec++){
              TwoDimReconSummary* reconsum = evtsum->GetModTwoDimRecon(irec,c_mod,cyc,iview);
              double endxy  = reconsum->endxy;
              double endpln = reconsum->endpln;
              double slope  = reconsum->slope;
              double endx,endy,endz;
              detdim->GetPosInMod(c_mod,endpln,iview,0,&endx,&endy,&endz);
              double ang_diff,pos_diff;
              if(iview==0){
                ang_diff = fabs(atan(slope)-atan(slopey))*180./PI;
                pos_diff = 
                  fabs(slopey*(endz-offset[2])-(endxy-offset[1])+intcpty)
                  /sqrt(1.+slopey*slopey);
              }
              else{
                ang_diff = fabs(atan(slope)-atan(slopex))*180./PI;
                pos_diff =
                  fabs(slopex*(endz-offset[2])-(endxy-offset[0])+intcptx)
                  /sqrt(1.+slopex*slopex);
              }
              int nhits_ontrk = 0;
              //if((ang_diff<LIM_ANGLE_DIF)&&(pos_diff<LIM_POS_DIF))
              if(ang_diff<LIM_ANGLE_DIF)
              {
                int nrechits = reconsum->Nhits();
                for(int irhit=0;irhit<nrechits;irhit++){
                  HitSummary* rechits = reconsum->GetHit(irhit);
                  int hitmod  = rechits->mod;
                  int hitview = rechits->view;
                  int hitpln  = rechits->pln;
                  int hitch   = rechits->ch;
                  for(int i=0;i<nhit_around_reco[iview];i++){
                    if(
                        (hitmod ==v_hitmod [iview][i])&&
                        (hitview==v_hitview[iview][i])&&
                        (hitpln ==v_hitpln [iview][i])&&
                        (hitch  ==v_hitch  [iview][i]))
                    {
                      nhits_ontrk++;
                      break;
                    }
                  }
                }
                //if(nhits_ontrk>=LIM_NUM_HITS){
                if(nhits_ontrk>=LIM_NUM_HITS*0.7){
                  target_recon = true;
                }
              }
#ifdef DEBUG_TRKEFF
              if(c_mod==21)
              std::cout
                << " ang:" << ang_diff
                << " pos:" << pos_diff
                << " nhits:" << nhits_ontrk
                << std::endl;
#endif
            } //for(irec)
          } //for(diff_cyc)
          //h1[imod]->Fill(cos);
          if(iview==0){ h1[imod]->Fill(1./sqrt(1.+slopey*slopey)); }
          else        { h1[imod]->Fill(1./sqrt(1.+slopex*slopex)); }

          if(target_recon){
            //h2[imod]->Fill(cos);
            if(iview==0){ h2[imod]->Fill(1./sqrt(1.+slopey*slopey)); }
            else        { h2[imod]->Fill(1./sqrt(1.+slopex*slopex)); }

#ifdef DEBUG_TRKEFF
            if(c_mod==21)
            std::cout << "  >> wagasci recon OK, view=" << iview << std::endl;
#endif
          }
#ifdef DEBUG_TRKEFF
          else{
            if(c_mod==21)
            std::cout << "  >> wagasci recon NG, view=" << iview << std::endl;
          }
#endif
        }//iview

      } //for(imod)
    } //icyc
  } //ievt

  ofile->cd();
  for(int i=0;i<NMOD;i++){
    h1[i]->Write();
    h2[i]->Write();
  }
  ofile->Write();
  ofile->Close();

  delete badch;
  delete detdim;
  return 0;
}
