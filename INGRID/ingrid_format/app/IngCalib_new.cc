////##### Standard C++ lib. ######
//#include <iostream>
//#include <sstream>
//#include <fstream>
//#include <iomanip>
//#include <sys/stat.h>
//using namespace std;
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
////##### INGRID Library #########
//#include "EVENTSUMMARY.h"
//
//string data_calib_dir  = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_calib" ;
//string data_dst_dir    = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_dst"   ;
//string data_cosmic_dir = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_cosmic";


#include "header.hh"

const int nMod  = 17;
const int nView =  2;
const int nPln  = 24;
const int nCh   = 80;

float       gain[ nMod ][ nView ][ nPln ][ nCh ];
float        ped[ nMod ][ nView ][ nPln ][ nCh ];
float     logain[ nMod ][ nView ][ nPln ][ nCh ];
float      loped[ nMod ][ nView ][ nPln ][ nCh ];


void InitilizeCalibConst(){
  for(int mod=0; mod<nMod; mod++){
    for(int view=0; view<nView; view++){
      for(int pln=0; pln<nPln; pln++){
	for(int ch=0; ch<nCh; ch++){
	  gain    [mod][view][pln][ch] =  1e5;
	  logain  [mod][view][pln][ch] =  1e5;
	  ped     [mod][view][pln][ch] = -1e-5;
	  loped   [mod][view][pln][ch] = -1e-5;
	}
      }
    }
  }
}


float HighGain[ nMod ][ nView ][ nPln ][ nCh ];
float HighPed [ nMod ][ nView ][ nPln ][ nCh ];
int   LowPed  [ nMod ][ nView ][ nPln ][ nCh ];
bool  GoodCh  [ nMod ][ nView ][ nPln ][ nCh ];


bool readonline_calib_run(int *online_calib_run){
  string cardfile = "/home/t2kingrid/shell/card.txt";
  std::ifstream ifs(cardfile.c_str());
  if(ifs.fail()){
    std::cout << "There is no such a file : " << cardfile << std::endl;
    return false;
  }
  std::cout << "Calibration file : " << cardfile << std::endl;
  char temp[1000];
  char calib[1000]="#calib";
  
  while(ifs >> temp){
    if(strcmp(temp,calib)==0){
      ifs >> *online_calib_run;
      std::cout << "calib_run_num " << *online_calib_run << "\n";
    }
  }
  ifs.close();
  return true;
}


bool GetCalibConst(int irun,int isrun){
  string calibfile = Form("%s/ingrid_%08d_%04d_Calib00.root",
      data_calib_dir.c_str(),irun,isrun);
  FileStat_t info;
  if(gSystem->GetPathInfo(calibfile.c_str(),info)){
    std::cout << "There is no calib file: " << calibfile << std::endl;
    return false;
  }
  std::cout << "Calibration file : " << calibfile << std::endl;

  TFile* fTCalib = new TFile(calibfile.c_str());
  TTree* calibtree = (TTree*)fTCalib->Get("calibtree");
  if(!calibtree){
    std::cout << "There is no calib tree. File:" << calibfile << std::endl;
    return false;
  }
  calibtree -> SetBranchAddress("HighGain", HighGain);
  calibtree -> SetBranchAddress("HighPed" , HighPed );
  calibtree -> SetBranchAddress("LowPed"  , LowPed  );
  calibtree -> SetBranchAddress("GoodCh"  , GoodCh  );
  calibtree -> GetEntry(0);
  for(int imod=0; imod<nMod; imod++){
    for(int iview=0; iview<nView; iview++){
      for(int ipln=0; ipln<nPln; ipln++){
	for(int ich=0; ich<nCh ; ich++){
	  if( GoodCh[imod][iview][ipln][ich] ){
	    gain  [imod][iview][ipln][ich] = HighGain[imod][iview][ipln][ich];
	    logain[imod][iview][ipln][ich] = 0.1 * HighGain[imod][iview][ipln][ich];
	    ped   [imod][iview][ipln][ich] = HighPed[imod][iview][ipln][ich];
	    loped [imod][iview][ipln][ich] = LowPed[imod][iview][ipln][ich];
	  }
	  else{
	    gain  [imod][iview][ipln][ich] = 1000;
	    logain[imod][iview][ipln][ich] = 1000;
	    ped   [imod][iview][ipln][ich] = 1000;
	    loped [imod][iview][ipln][ich] = 1000;
	  }
	}
      }
    }
  }
    
  //### double use VETO 
  for(int i=1; i<=6; i++){
    for(int j=0; j<22; j++){
      gain  [i][1][11][j] = gain  [i-1][1][12][j];
      ped   [i][1][11][j] = ped   [i-1][1][12][j];
      logain[i][1][11][j] = logain[i-1][1][12][j];
      loped [i][1][11][j] = loped [i-1][1][12][j];
    }
  }
  for(int i=8; i<=13; i++){
    for(int j=0; j<22; j++){
      gain[i][0][13][j]   = gain[i-1][0][14][j];
      ped [i][0][13][j]   = ped [i-1][0][14][j];
      logain[i][0][13][j] = gain[i-1][0][14][j];
      loped [i][0][13][j] = ped [i-1][0][14][j];
    }
  }

  fTCalib->Close();

  return true;

}

const static Int_t   cTimeCorrectionBase = 24;
//Time correction for difference of cable length

int main(int argc,char *argv[]){

  cout << "===============" << endl;
  cout << " IngCalib_new  " << endl;
  cout << "===============" << endl;

  char FileName[300];
  char Output  [300];
  int  run_number     = 0;
  int  sub_run_number = 0;
  int  c=-1;

  FileStat_t fs;
  bool cosmic        = false;
  bool rename_input  = false;
  bool rename_output = false;
  bool semioffline   = false; //if true, use the reference table
  while ((c = getopt(argc, argv, "r:s:i:o:cq")) != -1) {
    switch(c){
    case 'r':
      run_number=atoi(optarg);
      break;
    case 's':
      sub_run_number=atoi(optarg);
      break;
    case 'i':
      sprintf(FileName,"%s",optarg);
      rename_input = true;
      break;
    case 'o':
      sprintf(Output,"%s",optarg);
      rename_output = true;
      break;
    case 'c':
      cosmic = true;
      break;
    case 'q':
      semioffline = true;
      break;
    }
  }

  InitilizeCalibConst(); //### Initialize gain, pedestal, and so on

  //#### Read file before calibration ######
  //########################################
  EventSummary* summary = new EventSummary();

  if(!rename_input){
    if(cosmic){
      sprintf(FileName,"%s/ingrid_%08d_%04d.root",
          data_cosmic_dir.c_str(),run_number,sub_run_number);
    }
    else{
      sprintf(FileName,"%s/ingrid_%08d_%04d.root",
          data_dst_dir.c_str(),run_number,sub_run_number);
    }
  }
  
  if(gSystem->GetPathInfo(FileName,fs)){
    cout<<"Cannot open file "<<FileName<<endl;
    exit(1);
  }
  cout << "Input file : " << FileName << endl;

  TFile*   rfile = new TFile(FileName,"read");
  TTree*   tree  = (TTree*)rfile->Get("tree");
  TBranch* EvtBr = tree->GetBranch("fDefaultReco.");
  EvtBr          ->SetAddress(&summary);
  tree           ->SetBranchAddress("fDefaultReco.", &summary);
  int      nevt  = (int)tree -> GetEntries();


  //###### Get calibration info ########
  //####################################
  int calib_reference   = run_number;
  int calib_reference_s = sub_run_number;
  if(semioffline){
    int online_calib_run;
    if(!readonline_calib_run(&online_calib_run)){
      std::cout << "ERROR: Failed to open reference file." << std::endl;
      rfile->Close();
      return 1;
    }
    calib_reference   = online_calib_run;
    calib_reference_s = 0;
  }
  //if(!GetCalibConst(calib_reference,calib_reference_s)){
  //  std::cout << "ERROR: Failed to open calibration file." << std::endl;
  //  if(calib_reference_s!=0){
  //    std::cout << "       Trying to open another file with one previous sub-run." << std::endl;
  //    if(!GetCalibConst(calib_reference,calib_reference_s-1)){
  //      std::cout << "ERROR: Failed to open calibration file." << std::endl;
  //      rfile->Close();
  //      return 1;
  //    }
  //  }
  //  else{
  //    rfile->Close();
  //    return 1;
  //  }
  //}
  int ini_srun = calib_reference_s;
  bool up = false;
  while(!GetCalibConst(calib_reference,calib_reference_s)){
    if(calib_reference_s>ini_srun+3){
      std::cout << "ERROR: Failed to open calibration file. Run="
        << calib_reference << endl;
      rfile->Close();
      return 1;
    }
    else if(calib_reference_s==0){ up=true; }

    if(up){ calib_reference_s++; }
    else  { calib_reference_s--; }
  }

  //#### Create file after calibration ######
  //#########################################
  if(!rename_output){
    if(cosmic){
      sprintf(Output,"%s/ingrid_%08d_%04d_calib.root",
          data_cosmic_dir.c_str(),run_number,sub_run_number);
    }
    else{
      sprintf(Output,"%s/ingrid_%08d_%04d_calib.root",
          data_dst_dir.c_str(),run_number,sub_run_number);
    }
  }

  std::cout << "Output file : " << Output << std::endl;

  TFile*           wfile = new TFile(Output,"recreate");
  TTree*           wtree = new TTree("tree","tree");
  EventSummary* wsummary = new EventSummary(); 
  wtree              -> Branch   ("fDefaultReco.","EventSummary", 
				 &wsummary,  64000,   99);

  // For calibration information 
  TTree* calibtree = new TTree("calibtree","calibtree");
  int c_mod,c_view,c_pln,c_ch;
  double c_gain,c_ped,c_logain,c_loped;
  calibtree->Branch("mod"   ,   &c_mod ,"c_mod/I"   );
  calibtree->Branch("view"  ,  &c_view ,"c_view/I"  );
  calibtree->Branch("pln"   ,   &c_pln ,"c_pln/I"   );
  calibtree->Branch("ch"    ,    &c_ch ,"c_ch/I"    );
  calibtree->Branch("gain"  ,  &c_gain ,"c_gain/D"  );
  calibtree->Branch("ped"   ,   &c_ped ,"c_ped/D"   );
  calibtree->Branch("logain",&c_logain ,"c_logain/D");
  calibtree->Branch("loped" , &c_loped ,"c_loped/D" );
  for(c_mod=0; c_mod<nMod; c_mod++){
    if(c_mod!=3&&c_mod!=14&&c_mod!=15&&c_mod!=16){continue;}
    for(c_view=0; c_view<nView; c_view++){
      for(c_pln=0; c_pln<nPln; c_pln++){
        for(c_ch=0; c_ch<nCh ; c_ch++){
          if(GoodCh[c_mod][c_view][c_pln][c_ch]){
            c_gain   = gain  [c_mod][c_view][c_pln][c_ch];
            c_logain = logain[c_mod][c_view][c_pln][c_ch];
            c_ped    = ped   [c_mod][c_view][c_pln][c_ch];
            c_loped  = loped [c_mod][c_view][c_pln][c_ch];
            calibtree->Fill();
          }
        }
      }
    }
  }
  calibtree->Write();

  cout << "Start event loop : Total# of event = " << nevt << endl;
  for(int ievt = 0; ievt < nevt; ++ievt){
    if(ievt%100==0)cout<<"event:"<<ievt<<endl;
    summary -> Clear();
    wsummary-> Clear();
    tree    -> GetEntry(ievt);
    for(int mod=0; mod<nMod; mod++){
      for(int cyc=0; cyc<23; cyc++){
        int ninghit = summary -> NModHits(mod, cyc);
        for(int i=0; i<ninghit; i++){

          HitSummary *inghitsum;
          inghitsum   = (HitSummary*) (summary -> GetModHit(i, mod, cyc) );

          int view = inghitsum -> view;
          int pln  = inghitsum -> pln;
          int ch   = inghitsum -> ch;

          //##### Conversion from ADC to #p.e. ##############
          inghitsum -> pe   = 1.0 * ( inghitsum ->  adc - ped[mod][view][pln][ch] ) / gain[mod][view][pln][ch];

          //##### Conversion from ADC to #p.e.(Logain)#######
          inghitsum -> lope = 1.0 * ( inghitsum ->  loadc - loped[mod][view][pln][ch] ) / logain[mod][view][pln][ch] ;

          //##### Conversiont from TDC to time[nsec] ########

          long time = 2.5 * ( inghitsum ->  tdc ) - 10.0 * ( summary -> trgtime ); 

          //##### If we don't have Pulse Per Second, ########
          //##### time is larger than one second     ########
          //#################################################
          while(time>1000000000){
            time -= 1000000000;
          }
          while(cosmic&&time>1000000000-100000000){
            time -= 1000000000;
          }

          //##### We have to do  correction because of ##########
          //##### difference of cable length b/w  ###############
          //##### Back end board and front end board  ###########
          //##### but some VETO channels should be careful to do 
          //#####################################################
          float cTimeCorrection;
          if(!cosmic){ cTimeCorrection = cTimeCorrectionBase;     }
          else       { cTimeCorrection = 0.5*cTimeCorrectionBase; }

          //H and V TOF correcion is added from 2010 summer
          if(0<=mod&&mod<=6 ){ time -= 7; }
          if(7<=mod&&mod<=13){ time += 7; }

          switch ( mod ) {
            case  3:
              if( pln != 11 ) //#### Because pln 11 at mod 3 is pln 12 at mod 2        
                time -= cTimeCorrection;
              break;
            case  4:
              time -= cTimeCorrection;
              break;
            case  5:
              time -= cTimeCorrection;
              break;
            case  6:
              if( pln == 11 ) //#### Because pln 11 at mod 3 is pln 12 at mod 2        
                time -= cTimeCorrection;
              break;
            case  7:
              time += cTimeCorrection;
              break;
            case  8:
              time += cTimeCorrection;
              break;
            case  9:
              time += cTimeCorrection;
              break;
            case  10:
              if( pln == 13)
                time += cTimeCorrection;
              break;
            case 13:
              if( pln != 13)
                time += cTimeCorrection;
              break;
            default:
              break;
          }
          inghitsum -> time = time;

        }//Hit Loop
      }//Cyc
    }//Mod
    wsummary = summary;
    wtree -> Fill();
    if(ievt%1000==0)wtree->AutoSave();
  }
  wtree -> Write();
  wfile -> Write();

  wfile -> Close();

  std::cout << "IngCalib_new has been done." << std::endl;

} //main
