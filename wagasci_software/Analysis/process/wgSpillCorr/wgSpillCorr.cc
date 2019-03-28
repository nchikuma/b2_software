#include <string>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <dirent.h>
#include <ctime>
#include <fstream>
#include <stdlib.h>
#include <math.h>

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

void SpillCheck(string inputDir, string outputDir, int runid, int acqid);
bool CorrSpillNb(int dspill,int ddspill);


int main(int argc, char** argv){
  int opt;
  int runid = -1;
  int acqid = -1;
  string inputDirName  = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_dst_wagasci";
  string outputDirName = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_dst_wagasci";

  while((opt = getopt(argc,argv, "f:o:r:s:h")) !=-1 ){
    switch(opt){
      case 'f':
        inputDirName = optarg;
        break;
      case 'o':
        outputDirName = optarg; 
        break;
      case 'r':
        runid=atoi(optarg); 
        break;
      case 's':
        acqid=atoi(optarg); 
        break;
      case 'h':
        cout <<"Usage: " << argv[0] << " [options] " <<endl;
        cout <<"  -h               : help"<<endl;
        cout <<"  -f <inputDir>    : Input Directory  (default: "<<inputDirName <<" )" << endl;
        cout <<"  -o <outputDir>   : Output Directory (defualt: "<<outputDirName<<" )" << endl;
        cout <<"  -r <runid>       : Run ID" << endl;
        cout <<"  -s <acqid>       : Acq ID" << endl;
        exit(0);
    }   
  }

  if(runid<0||acqid<0){
    cout << "See the usage: " << argv[0] << " -h" << endl;
    exit(0);
  }

  cout << " *****  READ   DIRECTORY     :" << inputDirName  << "  *****" << endl;
  cout << " *****  OUTPUT DIRECTORY     :" << outputDirName << "  *****" << endl;
  cout << " *****  RUN ID = " << runid << " ACQ ID = " << acqid << endl;

  SpillCheck(inputDirName,outputDirName,runid,acqid);
  cout << "End of wgSpillCorr...." << endl;
  return 0;
}

//******************************************************************
void SpillCheck(string inputDir, string outputDir, int runid, int acqid)
{
  // ================
  // Open input file 
  string filename = Form("%s/run_%05d_%03d_inglib.root",inputDir.c_str(),runid,acqid);
  if(gSystem->GetPathInfo(filename.c_str(),fs)){
    cout << "No such a file." << filename << endl;
    return;
  }
  cout << "Opend WAGASCI file : " << filename << endl;

  TFile        *ifile     = new TFile(filename.c_str(),"read");
  EventSummary *isummary  = new EventSummary();
  TTree        *itree     = (TTree*)ifile->Get("tree");
  TBranch      *ibranch   = itree->GetBranch("fDefaultReco.");
  ibranch                 ->SetAddress(&isummary);
  itree                   ->SetBranchAddress("fDefaultReco.",&isummary);
  TTree        *info_tree = (TTree*)ifile->Get("info");

  // =======================
  // Open output file
  string outfilename = Form("%s/run_%05d_%03d_spillcorr.root" ,outputDir.c_str(),runid,acqid);
  cout << "Output file : " << outfilename << endl;
  TFile        *ofile    = new TFile(outfilename.c_str(),"recreate");
  TTree        *otree    = new TTree("tree","tree");
  EventSummary *osummary = new EventSummary();
  otree                  ->Branch("fDefaultReco.","EventSummary",&osummary,64000,99);
  TTree* info_out = info_tree->CloneTree(info_tree->GetEntries());
  info_out->Write();

  // =======================
  // Event loop starts
  int nevt = itree->GetEntries();
  int ievt = 0;
  cout << "Total number of events = " << nevt << endl;

  int spill     = -1;
  int spill_out = -1;
  int last_spill[2] = {-1,-1};
  int dspill    [2] = {-1,-1};
  int ddspill       = -1;

  while(ievt<nevt+1){
    
    if(ievt<nevt){ 
      itree->GetEntry(ievt); 
      spill = isummary->nd280nspill;
    }
    else{ spill = spill+1; }


    bool isGoodSpill = false;
    if(last_spill[0]!=-1){
      dspill[0] = spill-last_spill[0];
      dspill[1] = last_spill[0]-last_spill[1];
      ddspill   = dspill[0]-dspill[1];

      if(last_spill[1]==-1){ isGoodSpill = true; }
      else{
        if(dspill[0]==1||dspill[1]==1){ isGoodSpill = true; }
        else if(spill_out==0&&dspill[1]==-0x7fff){ isGoodSpill = true; }
        else if(CorrSpillNb(dspill[0],ddspill)){
          spill_out = spill-1;
          last_spill[0] = spill_out;
          isGoodSpill = true;
        }
        else if(dspill[1]>1&&dspill[0]>0){
          isGoodSpill = true;
        }
      }
    }

    if(isGoodSpill){
      itree->GetEntry(ievt-1);
      osummary = isummary;
      osummary->nd280nspill = spill_out;
      otree->Fill();
    }

    // ======================
    // Set the current info for the next loop
    // and fille the info into the tree.
    last_spill[1] = last_spill[0];
    last_spill[0] = spill;
    spill_out     = spill;

    ievt++;
  }
  cout << "Event loop was done." << endl;

  otree->Write();
  ofile->Write();
  ofile->Close();
  ifile->Close();
}

bool CorrSpillNb(int dspill,int ddspill){
  if(
      (dspill==0x0001+1&&ddspill==0x0001*2)||
      (dspill==0x0002+1&&ddspill==0x0002*2)||
      (dspill==0x0004+1&&ddspill==0x0004*2)||
      (dspill==0x0008+1&&ddspill==0x0008*2)||
      (dspill==0x0010+1&&ddspill==0x0010*2)||
      (dspill==0x0020+1&&ddspill==0x0020*2)||
      (dspill==0x0040+1&&ddspill==0x0040*2)||
      (dspill==0x0080+1&&ddspill==0x0080*2)||
      (dspill==0x0100+1&&ddspill==0x0100*2)||
      (dspill==0x0200+1&&ddspill==0x0200*2)||
      (dspill==0x0400+1&&ddspill==0x0400*2)||
      (dspill==0x0800+1&&ddspill==0x0800*2)||
      (dspill==0x1000+1&&ddspill==0x1000*2)||
      (dspill==0x2000+1&&ddspill==0x2000*2)||
      (dspill==0x4000+1&&ddspill==0x4000*2)||
      (dspill==-0x0001+1&&ddspill==-0x0001*2)||
      (dspill==-0x0002+1&&ddspill==-0x0002*2)||
      (dspill==-0x0004+1&&ddspill==-0x0004*2)||
      (dspill==-0x0008+1&&ddspill==-0x0008*2)||
      (dspill==-0x0010+1&&ddspill==-0x0010*2)||
      (dspill==-0x0020+1&&ddspill==-0x0020*2)||
      (dspill==-0x0040+1&&ddspill==-0x0040*2)||
      (dspill==-0x0080+1&&ddspill==-0x0080*2)||
      (dspill==-0x0100+1&&ddspill==-0x0100*2)||
      (dspill==-0x0200+1&&ddspill==-0x0200*2)||
      (dspill==-0x0400+1&&ddspill==-0x0400*2)||
      (dspill==-0x0800+1&&ddspill==-0x0800*2)||
      (dspill==-0x1000+1&&ddspill==-0x1000*2)||
      (dspill==-0x2000+1&&ddspill==-0x2000*2)||
      (dspill==-0x4000+1&&ddspill==-0x4000*2))
  {
    return true;
  }
  else return false;

}
