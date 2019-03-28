#include <string>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <dirent.h>
#include <math.h>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <TCanvas.h>
#include <TLegend.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TBox.h>
#include <TPaveText.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TArrow.h>
#include <TFrame.h>
#include <TGaxis.h>
#include <TStyle.h>
#include <TSystem.h>

#include "wgTools.h"
#include "wgErrorCode.h"
#include "wgEditXML.h"
#include "wgColor.h"
#include "wgGetTree.h"
#include "wgGetCalibData.h"
#include "wgChannelMap.h"
#include "wgDetectorDimension.h"

#include "EVENTSUMMARY.h"

using namespace std;


FileStat_t fs;

void SpillCheck(string inputDirName, string outputDirName, 
    int runid, int srunid,
    string wagasciDirName,int start_wgrun, int stop_wgrun);

int main(int argc, char** argv){
  int opt;
  int runid  = -1;
  int srunid = -1;
  int start_wgrun = -1;
  int stop_wgrun  = -1;
  string inputDirName   = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_dst";
  string outputDirName  = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_dst";
  string wagasciDirName = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_dst_wagasci";

  while((opt = getopt(argc,argv, "f:o:w:r:s:i:j:h")) !=-1 ){
    switch(opt){
      case 'f':
        inputDirName = optarg;
        break;
      case 'o':
        outputDirName = optarg; 
        break;
      case 'w':
        wagasciDirName = optarg;
        break;
      case 'r':
        runid=atoi(optarg); 
        break;
      case 's':
        srunid=atoi(optarg); 
        break;
      case 'i':
        start_wgrun=atoi(optarg);
        break;
      case 'j':
        stop_wgrun=atoi(optarg);
        break;
      case 'h':
        cout <<"This program is for data quality check. "<<endl;
        cout <<"You can take several option..."<<endl;
        cout <<"  -h               : help"<<endl;
        cout <<"  -f <inputDir>    : Input Directory. *Default: "      <<inputDirName   << endl;
        cout <<"  -o <outputDir>   : Output Directory *Defualt: "      <<outputDirName  << endl;
        cout <<"  -w <wagasciDir>  : WAGASCI Data Directory *Defualt: "<<wagasciDirName << endl;
        cout <<"  -r <runid>       : Run ID"     << endl;
        cout <<"  -s <srunid>      : Sub Run ID" << endl;
        cout <<"  -i <start_wgrun> : Start WAGASCI Run ID" << endl;
        cout <<"  -j <stop_wgrun>  : Stop WAGASCI Run ID"  << endl;
        exit(0);
    }   
  }

  if(runid<0||srunid<0||start_wgrun<0||stop_wgrun<0){
    cout << "See the usage: " << argv[0] << " -h" << endl;
    exit(0);
  }

  cout << " *****  READ   DIRECTORY     :" << inputDirName  << "  *****" << endl;
  cout << " *****  OUTPUT DIRECTORY     :" << outputDirName << "  *****" << endl;
  cout << " *****  RUN ID = " << runid << " Sub-Run ID = " << srunid << endl;

  SpillCheck(inputDirName,outputDirName,runid,srunid,
      wagasciDirName,start_wgrun,stop_wgrun);
  cout << "End of " << argv[0] << endl;
  return 0;
}

//******************************************************************
void SpillCheck(string inputDirName, string outputDirName, 
    int runid, int srunid,
    string wagasciDirName,int start_wgrun, int stop_wgrun)
{
  // ======================
  // Open input INGRID file 
  //
  string filename = Form("%s/ingrid_%08d_%04d_calib.root",inputDirName.c_str(),runid,srunid);
  if(gSystem->GetPathInfo(filename.c_str(),fs)){
    cout << "Cannot open file: "<< filename << endl;
    return;
  }
  TFile* ingfile = new TFile(filename.c_str(),"read");
  cout << "Reading ... runid=" << runid << " srunid=" << srunid << endl;
  cout << "Opend a INGRID file : " << filename << endl;

  EventSummary *summary = new EventSummary();
  TTree        *tree    = (TTree*)ingfile->Get("tree"); 
  TBranch      *evtbr   = tree->GetBranch("fDefaultReco.");
  evtbr                 ->SetAddress(&summary);
  tree                  ->SetBranchAddress("fDefaultReco.",&summary);
  int stoptime  = tree->GetMaximum("fDefaultReco.time");
  int starttime = tree->GetMinimum("fDefaultReco.time");
  int nevt = (int)tree->GetEntries();

  int ievt     =  0;
  int ingspill = -1;

  // ==================
  // Open output file
  //
  string outfilename = Form("%s/ingrid_%08d_%04d_spill.root",outputDirName.c_str(),runid,srunid);
  cout << "Output file : "<< outfilename << endl;
  TFile        *ofile    = new TFile(outfilename.c_str(),"recreate");
  TTree        *otree    = new TTree("tree","tree");
  EventSummary *osummary = new EventSummary();
  otree                  -> Branch("fDefaultReco.","EventSummary",&osummary,64000,99);

  // =========================================
  // Open WAGASCI files, and start event loop
  //
  TFile   *wgfile;
  TTree   *wgtree;
  TBranch *wgevtbr;
  bool wg_newfile  = true;
  int  wgevt       = 0;
  int  wgrun       = start_wgrun;
  int  wgacq       = 0;
  int  wgspill     = -1;
  Long64_t  wgstarttime = -1;
  Long64_t  wgstoptime  = -1;
  int  wg_nevt     = -1;
  EventSummary *wgsummary;
  while(ievt<nevt && wgrun<=stop_wgrun && wgacq<=600){

    if(wg_newfile){
      string wgfilename = Form("%s/run_%05d_%03d_spillcorr.root",wagasciDirName.c_str(),wgrun,wgacq);
      cout << "Checking WAGASCI file : " << wgfilename << " ....";
      if(gSystem->GetPathInfo(wgfilename.c_str(),fs)){ 
        cout << " Not exist" << endl; 
        if(wgacq<600){wgacq++;} 
        else{wgrun++;wgacq=0;}       
        continue; 
      }
      wgfile    = new TFile(wgfilename.c_str(),"read");
      wgtree    = (TTree*)wgfile->Get("tree"); 
      wgsummary = new EventSummary();
      wgevtbr = wgtree->GetBranch("fDefaultReco.");
      wgevtbr          ->SetAddress(&wgsummary);
      wgtree           ->SetBranchAddress("fDefaultReco.",&wgsummary);
      TTree *info_tree = (TTree*)wgfile->Get("info");
      info_tree->SetBranchAddress("start_time",&wgstarttime);
      info_tree->SetBranchAddress("stop_time" ,&wgstoptime );
      info_tree->GetEntry(0);
      if(wgstoptime<starttime||wgstarttime>stoptime){ 
        cout << " NG." << endl;
        if(wgacq<600){wgacq++;}
        else{wgrun++;wgacq=0;}
        wgfile->Close();
        continue;
      }
      else{
        cout << " OK." << endl;
        wg_newfile = false;
        wg_nevt = wgtree->GetEntries();
      }
    }

    tree  ->GetEntry(ievt );
    wgtree->GetEntry(wgevt);
    ingspill = summary  ->nd280nspill%0x8000;
    wgspill  = wgsummary->nd280nspill%0x8000;

    int spilldiff = ingspill - wgspill;
    if(spilldiff==0){
      osummary = summary;
      for(int icyc=0;icyc<23;icyc++){
        int nwghit = wgsummary->NModHits(MOD_B2_WM,icyc);
        for(int ihit=0;ihit<nwghit;ihit++){
          HitSummary* hitsum = wgsummary->GetModHit(ihit,MOD_B2_WM,icyc);
          osummary->AddModHit(hitsum,MOD_B2_WM,icyc);
        }
      }
      otree->Fill();
      wgevt++;
      ievt++;
    }
    else{
      if(spilldiff<0){
        ievt++;
      }
      else{
        wgevt++;
      }
    }

    if(wgevt>=wg_nevt){
      wgfile->Close();
      if(wgacq<600){wgacq++;}
      else{wgrun++;wgacq=0;}
      wgevt=0;
      wg_newfile=true;
    }

  }

  ofile->cd();
  otree   ->Write();
  ofile   ->Write();
  ofile   ->Close();
  ingfile ->Close();

  return;
}
