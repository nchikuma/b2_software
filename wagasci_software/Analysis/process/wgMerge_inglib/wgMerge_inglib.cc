/* ***********************************************************************
 * Reconstruction program for wagasci 
 * Program : wagasci_recon_exe.cc
 * Name: Naruhiro Chikuma
 * Date: 2017-05-09 02:45:48
 * ********************************************************************** */
#include <iostream>
#include <iomanip>

#include "wgReconClass.h"
#include "wgGetTree.h"
#include "wgGetCalibData.h"
#include "wgErrorCode.h"
#include "wgTools.h"
#include "Const.h"

#include "wgDetectorDimension.h"
#include "EVENTSUMMARY.h"

using namespace std;

int main(int argc, char* argv[]){

  int opt;
  wgConst *con = new wgConst();
  con->GetENV();
  string inputDir     = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/decode";
  string outputDir    = "";
  string outputDir_b  = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_dst_wagasci";
  string outputDir_c  = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_cosmic_wagasci";
  string logoutputDir = con->LOG_DIRECTORY;
  CheckExist *check = new CheckExist();
  

  int runid = -1;
  int acqid = -1;

  bool cosmic = false;
  bool rename = false;
  while((opt = getopt(argc,argv, "r:s:f:o:ch")) !=-1 ){
    switch(opt){
      case 'r':
        runid = atoi(optarg);
        break;
      case 's':
        acqid = atoi(optarg);
        break;
      case 'f':
        inputDir = optarg;
        break;
      case 'o': 
        rename = true;
        outputDir=optarg;
        if(!check->Dir(outputDir)){
          cout << "!!Error!! output directory:" << outputDir<< " don't exist!" << endl;
          return 1;
        }
        break;
      case 'c':
        cosmic = true;
        break;
      case 'h':
        cout <<"This program is for reconstrancting track from tree file."  << endl;
        cout <<"You can take several option..."                             << endl;
        cout <<"Usage : ./wgRecon [-h] [-f input] ..."                      << endl;
        cout <<"  -h         : help"                                        << endl;
        cout <<"  -f (char*) : choose input directory. *Default: "<<inputDir<< endl;
        cout <<"  -r (int)   : choose WAGASCI Run ID"                       << endl;
        cout <<"  -s (int)   : choose WAGASCI Acq ID"                       << endl;
        cout <<"  -o (char*) : choose output directory "                    << endl;
        cout <<"               *Default (beam trigger)   "   << outputDir_b << endl;
        cout <<"               *Default (cosmic trigger) "   << outputDir_c << endl;
        cout <<"  -c         : run for cosmic trigger mode."                << endl;
        cout <<" (ex.) ./wgRecon -f test"                                   << endl;
        exit(0);
      default:
        cout <<"See the help. : " << argv[0] << " -h"<<endl;
        exit(0);
    }
  }

 
  if(runid<0||acqid<0){
    cout <<"See the help. : " << argv[0] << " -h"<<endl;
    exit(1);
  }

  if(!rename){
    if(!cosmic){ outputDir = outputDir_b; }
    else       { outputDir = outputDir_c; }
  }

  // =====================
  // initialize reconstruction Class
  wgRecon wg_rec;
  
  // =====================
  // open ROOT Tree files

  string inputFileName = Form("%s/run_%05d_%03d",inputDir.c_str(),runid,acqid);
  string filename1 = Form("%s_dif_1_1_1_tree.root", inputFileName.c_str());
  string filename2 = Form("%s_dif_1_1_2_tree.root", inputFileName.c_str());
  
  TFile* datafile1 = new TFile(filename1.c_str(),"read");
  if(!datafile1){
    cerr << "Cannot open file : " << filename1.c_str() << endl;
    exit(1);
  }
  TFile* datafile2 = new TFile(filename2.c_str(),"read");
  if(!datafile2){
    cerr << "Cannot open file : " << filename2.c_str() << endl;
    exit(1);
  }

  TTree* tree[NumDif];
  TTree* info_out;
  Long64_t starttime = 0;
  Long64_t stoptime  = 0;
  for(int i=0;i<NumDif;i++){
    if(i==0){
      datafile1->cd();
      tree[i] = (TTree *)datafile1->Get("tree");
      TTree *tmp_info = (TTree*) datafile1->Get("info");
      tmp_info->SetBranchAddress("start_time",&starttime);
      tmp_info->SetBranchAddress("stop_time" ,&stoptime );
      tmp_info->GetEntry(0);
      info_out = tmp_info->CloneTree(tmp_info->GetEntries());
      delete tmp_info; 
    } 
    if(i==1){
      datafile2->cd();
      tree[i] = (TTree *)datafile2->Get("tree");
    } 
    tree[i]->SetBranchAddress("charge"     , wg_rec.type_raw[i].charge     );
    tree[i]->SetBranchAddress("time"       , wg_rec.type_raw[i].time       );
    tree[i]->SetBranchAddress("gs"         , wg_rec.type_raw[i].gs         );
    tree[i]->SetBranchAddress("hit"        , wg_rec.type_raw[i].hit        );
    tree[i]->SetBranchAddress("bcid"       , wg_rec.type_raw[i].bcid       );
    tree[i]->SetBranchAddress("col"        , wg_rec.type_raw[i].col        );
    tree[i]->SetBranchAddress("ch"         , wg_rec.type_raw[i].ch         );
    tree[i]->SetBranchAddress("chip"       , wg_rec.type_raw[i].chip       );
    tree[i]->SetBranchAddress("chipid"     , wg_rec.type_raw[i].chipid     );
    tree[i]->SetBranchAddress("spill"      ,&wg_rec.type_raw[i].spill      );
    tree[i]->SetBranchAddress("spill_flag" ,&wg_rec.type_raw[i].spill_flag );
    tree[i]->SetBranchAddress("spill_mode" ,&wg_rec.type_raw[i].spill_mode );
    tree[i]->SetBranchAddress("spill_count",&wg_rec.type_raw[i].spill_count);
    tree[i]->SetBranchAddress("debug"      , wg_rec.type_raw[i].debug      );
    tree[i]->SetBranchAddress("view"       ,&wg_rec.type_raw[i].view       );
    tree[i]->SetBranchAddress("pln"        , wg_rec.type_raw[i].pln        );
    tree[i]->SetBranchAddress("ch"         , wg_rec.type_raw[i].ch         );
    tree[i]->SetBranchAddress("grid"       , wg_rec.type_raw[i].grid       );
    tree[i]->SetBranchAddress("x"          , wg_rec.type_raw[i].x          );
    tree[i]->SetBranchAddress("y"          , wg_rec.type_raw[i].y          );
    tree[i]->SetBranchAddress("z"          , wg_rec.type_raw[i].z          );
    tree[i]->SetBranchAddress("pe"         , wg_rec.type_raw[i].pe         );
    tree[i]->SetBranchAddress("gain"       , wg_rec.type_raw[i].gain       );
    tree[i]->SetBranchAddress("pedestal"   , wg_rec.type_raw[i].pedestal   );
    tree[i]->SetBranchAddress("ped_nohit"  , wg_rec.type_raw[i].ped_nohit  );
    tree[i]->SetBranchAddress("time"       , wg_rec.type_raw[i].time       );
    tree[i]->SetBranchAddress("time_ns"    , wg_rec.type_raw[i].time_ns    );
  }

  int nevt1 = tree[0]->GetEntries();
  int nevt2 = tree[1]->GetEntries();
  
  // ====================
  // set output ROOT file and TTree
  OperateString *OpStr = new OperateString();
  string outputname = Form("%s/%s_inglib.root",outputDir.c_str(),OpStr->GetName(inputFileName).c_str());
  int spill,spill_mode,spill_count;
  int debug[NumDif][NumChip];

  TFile* foutput = new TFile(outputname.c_str(),"recreate");
  TTree*           wtree = new TTree("tree","tree");
  EventSummary *wsummary = new EventSummary();
  HitSummary   *hitsum   = new HitSummary  ();
  wtree->Branch ("fDefaultReco.","EventSummary", &wsummary,64000,99);
  
  wgDetectorDimension *detdim   = new wgDetectorDimension();
  wgGetCalibData      *getcalib = new wgGetCalibData     ();


  // =====================
  // start event loop

  cout << "** inputFile1 : " << filename1 << " **" << endl;
  cout << "** inputFile2 : " << filename2 << " **" << endl;
  cout << "Number of event: " << nevt1 << " (dif1), " << nevt2 << " (dif2)" << endl;
  int startevt = 0;
  int maxevent= 999999;
  int ievt = 0;
  int ievt1=startevt;
  int ievt2=startevt;

  int unixtime = starttime;
  while(ievt1<nevt1&&ievt2<nevt2&&ievt1<maxevent&&ievt2<maxevent){
    
    wsummary->Clear("C");
    wg_rec.clear();
#ifdef DEBUG_RECON
    cout << "===================================================="
      << endl
      << "event=" << ievt << endl;
#else
    if(ievt%1000==0) cout << "event: " << ievt1 << "(dif1), "  << ievt2 << "(dif2)" << endl;
#endif
    tree[0]->GetEntry(ievt1);
    tree[1]->GetEntry(ievt2);
    int diff_spillcount = wg_rec.type_raw[0].spill_count - wg_rec.type_raw[1].spill_count;
    int diff_spill      = wg_rec.type_raw[0].spill       - wg_rec.type_raw[1].spill;
    int diff_spillmode  = wg_rec.type_raw[0].spill_mode  - wg_rec.type_raw[1].spill_mode;
    int num_rep = 0;
    int repstart1 = ievt1;
    int repstart2 = ievt2;
    while(diff_spillcount!=0||diff_spill!=0||diff_spillmode!=0)
    {
#ifdef DEBUG_RECON
      cout 
        << "{spill_count,spill,spill_mode} are different between dif1 and dif2 :" << endl
        << "event: " << ievt1 << "(dif1) " << ievt2 << "(dif2)"
        << " dif1={"<<wg_rec.type_raw[0].spill_count<<","<<wg_rec.type_raw[0].spill<<","<<wg_rec.type_raw[0].spill_mode<<"}"
        << " dif2={"<<wg_rec.type_raw[1].spill_count<<","<<wg_rec.type_raw[1].spill<<","<<wg_rec.type_raw[1].spill_mode<<"}" 
        << endl;
#endif
      if(num_rep<3){
        if(diff_spillcount>0){ievt2++;} else{ievt1++;}
      }
      else{
        ievt1=repstart1+1;
        ievt2=repstart2+1;
        num_rep=0;
        repstart1 = ievt1;
        repstart2 = ievt2;
      }
      if(ievt1<nevt1&&ievt2<nevt2&&ievt1<maxevent&&ievt2<maxevent){
        tree[0]->GetEntry(ievt1);
        tree[1]->GetEntry(ievt2);
        diff_spillcount = wg_rec.type_raw[0].spill_count - wg_rec.type_raw[1].spill_count;
        diff_spill      = wg_rec.type_raw[0].spill       - wg_rec.type_raw[1].spill;
        diff_spillmode  = wg_rec.type_raw[0].spill_mode  - wg_rec.type_raw[1].spill_mode;
      }
      else{
        goto ENDOFEVT;
      }
      num_rep++;
    }

    // ===============
    // Push hits into vector
    // Note: This is required to be set with "chipid", not with "chip".

    spill       = wg_rec.type_raw[0].spill;
    spill_mode  = wg_rec.type_raw[0].spill_mode;
    spill_count = wg_rec.type_raw[0].spill_count;


    int trgid = -1;
    if     ( cosmic&&spill_mode==0){trgid=128;}
    else if(!cosmic&&spill_mode==1){trgid=  1;}
    else{
      ievt1++;
      ievt2++;
      ievt++;
      continue;
    }
    if(spill_mode==1){ unixtime+=2.48; }

    wsummary->time        = unixtime; 
    wsummary->event       = spill_count&0xfff;
    wsummary->trgtime     = 0;
    wsummary->trgid       = trgid;
    wsummary->nd280nspill = spill%0x8000; //wagasci spill has the MSB fixed. 
                                          //1-bit inversion might happen.

    for(int dif=0;dif<NumDif;dif++){
      for(int chip=0;chip<NumChip;chip++){
        debug[dif][chip] = wg_rec.type_raw[dif].debug[chip];
      }
    }

    for(int dif=0;dif<NumDif;dif++){
      for(int chip=0;chip<NumChip;chip++){
        int chipid = wg_rec.type_raw[dif].chipid[chip];
        for(int chipch=0;chipch<NumChipCh;chipch++){
          for(int sca=0;sca<NumSca;sca++){
            if( wg_rec.type_raw[dif].hit[chip][chipch][sca]==1 ){              
              int    view   = wg_rec.type_map.view[dif][chipid][chipch];
              int    pln    = wg_rec.type_map.pln [dif][chipid][chipch];
              int    ch     = wg_rec.type_map.ch  [dif][chipid][chipch];
              int    charge = wg_rec.type_raw[dif].charge[chip][chipch][sca]; 
              int    tdc    = wg_rec.type_raw[dif].time  [chip][chipch][sca]; 
              int    gs     = wg_rec.type_raw[dif].gs  [chip][chipch][sca]; 
              int    bcid   = wg_rec.type_raw[dif].bcid[chip][sca];
              double pe     = wg_rec.type_raw[dif].pe[chip][chipch][sca];
              double time   = wg_rec.type_raw[dif].time_ns[chip][chipch][sca];
              wg_rec.pushHitInfo(bcid,view,pln,ch,dif,chipid,chipch,sca,charge,gs,tdc,pe,time);
            }
          }
        }
      }
    }

    wg_rec.sort_byBCIDnMAP();

    getcalib->Set_Timewalk();
    for(int i=0;i<wg_rec.get_num_hit();i++){
      int    bcid   = wg_rec.get_hitbcid   (i);
      int    view   = wg_rec.get_hitview   (i);
      int    pln    = wg_rec.get_hitpln    (i);
      int    ch     = wg_rec.get_hitch     (i);
      int    dif    = wg_rec.get_hitdif    (i);
      int    chip   = wg_rec.get_hitchip   (i);
      int    chipch = wg_rec.get_hitchipch (i);
      int    sca    = wg_rec.get_hitsca    (i);
      int    adc    = wg_rec.get_hitadc    (i);
      int    gs     = wg_rec.get_hitgs     (i);
      int    tdc    = wg_rec.get_hittdc    (i);
      double pe     = wg_rec.get_hitpe     (i);
      double hittime= wg_rec.get_hittime   (i);
      if(tdc==0) tdc =4096;
      int tdc_mod=0;
      if(bcid%2==0) tdc_mod = tdc + 10*getcalib->Get_Timewalk(adc-wg_rec.type_raw[dif].pedestal[chip][chipch][sca]);
      if(bcid%2==1) tdc_mod = tdc - 10*getcalib->Get_Timewalk(adc-wg_rec.type_raw[dif].pedestal[chip][chipch][sca]);
      if(tdc_mod<0) tdc_mod=0;
      else if(tdc_mod>4096) tdc_mod=4096;
      
      int mod = MOD_B2_WM;
      int cyc = bcid -22;           //This offset could randomly be 22 or 23.

      if(!cosmic){
        if(cyc<0||cyc>=23){ continue; } //INGRID cycle = 0-22.
      }
      else{
        if(!wg_rec.pushHitStruct(bcid,view,pln,ch,dif,chip,chipch,sca,adc,gs,tdc,pe,hittime,tdc_mod)){
          continue;
        }
      }

      //For Beam trigger, all hits during beam coiming are filled.
      if(!cosmic){
        double x,y,z,xy;
        int h_adc = -1;
        int l_adc = -1;
        detdim->GetPosInMod(mod,pln,view,ch,&x,&y,&z);
        if(view==SideView){xy = y;}else{xy = x;}
        if(gs==1){h_adc=adc;}else{l_adc=adc;}

        hitsum->Clear("C");
        hitsum->cyc   = cyc;
        hitsum->mod   = mod;
        hitsum->view  = 1-view;
        hitsum->pln   = pln;
        hitsum->ch    = ch;
        hitsum->xy    = xy;
        hitsum->z     = z;
        hitsum->adc   = h_adc;
        hitsum->loadc = l_adc;
        hitsum->tdc   = tdc;
        hitsum->pe    = pe;
        hitsum->lope  = pe;
        hitsum->time  = hittime;
        wsummary -> AddModHit(hitsum,mod,cyc);
      }
    }

    //For cosmic trigger, hits are filled only when they make a BCID cluster
    if(cosmic){
      if(wg_rec.findBCIDCluster()){
        for(int i=0;i<wg_rec.type_hit.num_bcid_cluster;i++){
          for(int j=0;j<wg_rec.type_hit.num_bcid_hits[i];j++){
            int hitid = wg_rec.type_recon.bcid_cluster_hitid[i][j];
            int    view   = wg_rec.get_hitview   (hitid);
            int    pln    = wg_rec.get_hitpln    (hitid);
            int    ch     = wg_rec.get_hitch     (hitid);
            int    adc    = wg_rec.get_hitadc    (hitid);
            int    gs     = wg_rec.get_hitgs     (hitid);
            int    tdc    = wg_rec.get_hittdc    (hitid);
            double pe     = wg_rec.get_hitpe     (hitid);
            double hittime= wg_rec.get_hittime   (hitid);
            if(tdc==0) tdc =4096;

            int mod = MOD_B2_WM;
            double x,y,z,xy;
            int h_adc = -1;
            int l_adc = -1;
            int cyc = 14+i;
            if(cyc>=23){cyc=0;}
            detdim->GetPosInMod(mod,pln,view,ch,&x,&y,&z);
            if(view==SideView){xy = y;}else{xy = x;}
            if(gs==1){h_adc=adc;}else{l_adc=adc;}

            hitsum->Clear("C");
            hitsum->cyc   = cyc;
            hitsum->mod   = mod;
            hitsum->view  = 1-view;
            hitsum->pln   = pln;
            hitsum->ch    = ch;
            hitsum->xy    = xy;
            hitsum->z     = z;
            hitsum->adc   = h_adc;
            hitsum->loadc = l_adc;
            hitsum->tdc   = tdc;
            hitsum->pe    = pe;
            hitsum->lope  = pe;
            hitsum->time  = hittime;
            wsummary -> AddModHit(hitsum,mod,cyc);

          }
        }
      }
    }

    wtree->Fill();
    ievt1++;
    ievt2++;
    ievt++;
  }// ievt loop end
ENDOFEVT:
  cout << "Number of events filled: " << ievt << endl;
  cout << "Total number of decoded events: "  << nevt1 << "(dif1) " << nevt2 << "(dif2)" << endl;
  
  wtree    ->Write();
  info_out ->Write();
  foutput  ->Write();
  foutput  ->Close();
    
  delete getcalib;
  delete detdim;

  cout << "File : " << outputname << " is written!" << endl;

  return 0;
}
