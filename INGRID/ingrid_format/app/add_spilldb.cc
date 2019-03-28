////##### Standard C++ lib. ######
//#include <iostream>
//#include <sstream>
//#include <fstream>
//using namespace std;
//#include <iomanip>
//#include <stdlib.h>
//#include <sys/stat.h>
//#include <math.h> 
////##### Root Library ###########
//#include <TROOT.h>
//#include <TStyle.h>
//#include <TApplication.h>
//#include <TFile.h>
//#include <TCanvas.h>
//#include <TTree.h>
//#include <TClonesArray.h>
//#include <TObject.h>
//#include <TEventList.h>
//#include <TBranch.h>
//#include <TSystem.h>
//
//const char* data_dir = "/home/t2k/nchikuma/b2_data/bsd";

#include "header.hh"

void PrintUsage()
{
  cout << " Usage: " << endl;
  cout << "  ./add_spilldb -m <MR Run#> -t 't2krun<#>' -v <v01/p06>" << endl;
  cout << "  -i <sub-run id> : option for starting in the middle." << endl;
  cout << endl;
}

FileStat_t fs;
int main(int argc, char *argv[]){

  int   c  = -1;
  int   run_number;      //### nu DAQ run number
  int   sub_run_number;  //### nu DAQ sub run number
  int   mrun = -1;       //### MR run number
  char  trun   [300];    //### T2K run number
  char  version[300];    //### bsd version 
  sprintf(version,"p06");//### default
  char FileName[300];
  int offset = 0;


  while ((c = getopt(argc, argv, "m:t:v:i:h")) != -1) {
    switch(c){
      case 'm':
        mrun = atoi(optarg);
        break;
      case 't':
        sprintf(trun, "%s", optarg);
        break;
      case 'v':
        sprintf(version, "%s", optarg);
        break;
      case 'i':
        offset = atoi(optarg);
        break;
      case 'h':
        PrintUsage();
        return 0;
      others:
        PrintUsage();
        return 0;
    }
  }
  if(mrun<0){
    PrintUsage();
    return 0;
  }
  int start_num   =       mrun  * 10000;
  int   end_num   = ( mrun + 1 )* 10000;


  cout << start_num << "\t" << end_num<< endl;
  sprintf(FileName, 
      "%s/spilldb/run%07d.db",
      data_bsd_dir.c_str(),
      mrun);
  cout<<FileName<<endl;
  ofstream wf(FileName);
  start_num += offset;
  for(run_number = start_num; run_number < end_num; run_number++){
    for(sub_run_number = 0; sub_run_number < 100; sub_run_number++){
      //#### read BSD ##########
      //########################
      sprintf(FileName, 
          "%s/%s/%s/bsd_run0%06d_%02d%s.root",
          data_bsd_dir.c_str(),
          version, trun,
          run_number, sub_run_number,
          version);

      if(gSystem -> GetPathInfo(FileName, fs)) //### not exist, skip ####
        continue;   

      cout << "read " << FileName << "..." << endl;

      wf << Form("%s/%s/bsd_run0%06d_%02d%s.root ",
          version, trun, run_number, sub_run_number, version)
        << run_number << " ";

      Int_t  trg_sec[3];
      Int_t  spillnum;
      TFile* rfile  = new TFile( FileName, "read" );
      TTree* tree   = (TTree*)rfile->Get("bsd");
      tree -> SetBranchAddress("trg_sec" , trg_sec );
      tree -> SetBranchAddress("spillnum",&spillnum);
      int        nevt   = tree -> GetEntries();
      //#### write down first and last time and spill# to the file #####
      //################################################################      
      tree -> GetEntry(        0 );
      wf << trg_sec[0] << " ";
      tree -> GetEntry( nevt - 1 );      
      wf << trg_sec[0] << " ";
      tree -> GetEntry(        0 );
      wf << spillnum   << " ";
      tree -> GetEntry( nevt - 1 );      
      wf << spillnum   << " ";
      wf << endl;

    }///sub_run
  }//run
}
