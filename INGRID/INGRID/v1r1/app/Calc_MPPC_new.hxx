// ND280 software includes
#include "TMidasBank.hxx"
#include "TMidasFile.hxx"
#include "TMidasBankProxy.hxx"
#include "TMidasBankProxyRegistry.hxx"
#include "TND280RawEvent.hxx"
#include "TRawDataHeader.hxx"
#include "TRunInfoBank.hxx"

// oaRawEvent includes
#include "TTripTDigitBank.hxx"
#include "TMidasTripTDigitItr.hxx"
#include "TMidasTripTDigit.hxx"
#include "TMCMBank.hxx"
#include "TTriggerBank.hxx"

// ROOT includes
#include "TApplication.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TStyle.h>
#include "TString.h"
#include "TSystem.h"
#include "TCanvas.h"
#include <stdlib.h>

#include <math.h>
#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <unistd.h> 

#include "setup.hxx"
#include "ana_MPPC.cxx"

using namespace std;

FileStat_t  fs;
ana_MPPC*   fana_MPPC;

Int_t       fMinAdcwTdcCut[NumMod][2][NumTFB][LayerNumCh];
Int_t       cAnaTrg;
Int_t       nCalib;
char        cCalibFile[300];


int   ANABREAK     = 900;
const static float ErrorMinGain =   5;
const static float ErrorMaxGain =  25;

TH1F* fH_sHighAdc[NumMod][2][NumTFB][LayerNumCh];
TH1F* fH_sLowAdc[NumMod][2][NumTFB][LayerNumCh];

TFile*   fTFile;
TTree*   tree;
//###############################################
//### following variables for calib. process ####

TH1F*    fH_HighAdc[NumMod][2][NumTFB][LayerNumCh]; //### High ADC dist. 
TH1F*    fH_LowAdc [NumMod][2][NumTFB][LayerNumCh]; //### Low  ADC dist.
Int_t    IngRun;       //### INGRID run number
Int_t    IngSubRun;    //### INGRID sub run number
Int_t    StartEvt;     //### Start event number of this calib. file
Int_t    StartTime;    //### Start time         of this calib. file
Int_t    EndEvt;       //### End   event number of this calib. file
Int_t    EndTime;      //### End   time         of this calib. file
Int_t    AnaEvt;       //###
Int_t    AnaTime;      //###
Int_t    Version;      //### version of anaMPPC
Int_t    MakeTime;     //### time this calib. file is made
Int_t    TrgMode;      //### trigger mode
Float_t  HighGain [NumMod][2][NumTFB][LayerNumCh]; //### MPPC gain[ADC counts] at HighGain ch
Float_t  HighPed  [NumMod][2][NumTFB][LayerNumCh]; //### Pedestal[ADC counts] at HighGain ch
Float_t  HighPedSigma  [NumMod][2][NumTFB][LayerNumCh]; //### Pedestal[ADC counts] at HighGain ch
Int_t    LowPed   [NumMod][2][NumTFB][LayerNumCh]; //### Pedestal[ADC counts] at LowGain ch
Float_t  MPPCNoise[NumMod][2][NumTFB][LayerNumCh]; //### MPPC noise[kHz] by HighGain ch
Float_t  MPPCCandA[NumMod][2][NumTFB][LayerNumCh]; //### MPPC CrossTalk&AfterPulse rate
Float_t  TdcThre  [NumMod][2][NumTFB][LayerNumCh];
Bool_t   GoodCh   [NumMod][2][NumTFB][LayerNumCh];


Bool_t   EvtCheck;
//###############################################
//### following variables for easy check ########
vector<int>    mod;
vector<int>    view;
vector<int>    pln;
vector<int>    ch;
vector<int>    rmm;
vector<int>    tfb;
vector<int>    trip;
vector<int>    trip_ch;
vector<double> gain;
vector<double> pedestal;
vector<double> refgain;
vector<double> refpedestal;
vector<double> pedestal_sigma;//added by koga
vector<double> noise;//added by koga
vector<double> canda;//added by koga
void CLEAR(){
  mod           .clear();
  view          .clear();
  pln           .clear();
  ch            .clear();
  rmm           .clear();
  tfb           .clear();
  trip          .clear();
  trip_ch.       clear();
  gain          .clear();
  pedestal      .clear();
  refgain       .clear();
  refpedestal   .clear();
  pedestal_sigma.clear();
  noise         .clear();
  canda         .clear();
}
Float_t  HighGainRef [NumMod][2][NumTFB][LayerNumCh]; //### MPPC gain[ADC counts] at HighGain ch
Float_t  HighPedRef [NumMod][2][NumTFB][LayerNumCh]; //### MPPC gain[ADC counts] at HighGain ch
void VClear(){
  mod.     clear();
  view.    clear();
  pln.     clear();
  ch.      clear();
  rmm.     clear();
  tfb.     clear();
  trip.    clear();
  trip_ch. clear();
  gain.    clear();
  pedestal.clear();
  refgain.    clear();
  refpedestal.clear();
  pedestal_sigma.clear();
  noise.clear();
  canda.clear();
}

void GetRef(){
  cout << "reading reference value" << endl;

  int tmod, tview, tpln, tch;
  for(tmod=0;tmod<NumMod;tmod++){
    for(tview=0;tview<2;tview++){
      for(tpln=0;tpln<NumTFB;tpln++){
        for(tch=0;tch<LayerNumCh;tch++){
          HighGainRef [tmod][tview][tpln][tch] = 0.;
          HighPedRef  [tmod][tview][tpln][tch] = 0.;
        }
      }
    }
  }

  string reffilename = Form("%s/ingrid_info/calib_ref_table.txt",ingrid_dir.c_str());
  ifstream f(reffilename.c_str());
  if(!f){
    cout << "There is no file for referential gains : " 
      << reffilename << endl;
  }
  else{
    double tgain, tped;
    while( f >> tmod >> tview >> tpln >> tch >> tgain >> tped ){
      HighGainRef [tmod][tview][tpln][tch] = tgain;
      HighPedRef  [tmod][tview][tpln][tch] = tped;
    }
  }
  f.close();

}
//__________________________________________________________ 
//
//
//__________________________________________________________ 
void Book(string fname){
  cout << "make calib file " << fname << endl;

  fTFile    = new TFile  (fname.c_str(), "recreate");
  tree      = new TTree  ("calibtree", "tree for calibration");


  tree -> Branch ("IngRun"   , &IngRun   , "IngRun/I"   );
  tree -> Branch ("IngSubRun", &IngSubRun, "IngSubRun/I");
  tree -> Branch ("StartTime", &StartTime, "StartTime/I");
  tree -> Branch ("StartEvt" , &StartEvt , "StartEvt/I" );
  tree -> Branch ("EndTime"  , &EndTime  , "EndTime/I"  );
  tree -> Branch ("EndEvt"   , &EndEvt   , "EndEvt/I"   );
  tree -> Branch ("AnaTime"  , &AnaTime  , "AnaTime/I"  );
  tree -> Branch ("AnaEvt"   , &AnaEvt   , "AnaEvt/I"   );
  tree -> Branch ("HighGain" , HighGain  , Form("HighGain[%d][%d][%d][%d]/F",NumMod,2,NumTFB,LayerNumCh));
  tree -> Branch ("HighPed"  , HighPed   , Form("HighPed[%d][%d][%d][%d]/F" ,NumMod,2,NumTFB,LayerNumCh));
  tree -> Branch ("LowPed"   , LowPed    , Form("LowPed[%d][%d][%d][%d]/I"  ,NumMod,2,NumTFB,LayerNumCh));
  tree -> Branch ("TdcThre"  , TdcThre   , Form("TdcThre[%d][%d][%d][%d]/F" ,NumMod,2,NumTFB,LayerNumCh));
  tree -> Branch ("GoodCh"   , GoodCh    , Form("GoodCh[%d][%d][%d][%d]/O"  ,NumMod,2,NumTFB,LayerNumCh));
  tree -> Branch ("EvtCheck" , &EvtCheck , "EvtCheck/O" );
  tree -> Branch ("mod"           ,&mod            );
  tree -> Branch ("view"          ,&view           );
  tree -> Branch ("pln"           ,&pln            );
  tree -> Branch ("ch"            ,&ch             );
  tree -> Branch ("gain"          ,&gain           );
  tree -> Branch ("pedestal"      ,&pedestal       );
  tree -> Branch ("refgain"       ,&refgain        );
  tree -> Branch ("refpedestal"   ,&refpedestal    );
  tree -> Branch ("pedestal_sigma",&pedestal_sigma );
  tree -> Branch ("noise"         ,&noise          );
  tree -> Branch ("canda"         ,&canda          );
  tree -> Branch ("rmm"           ,&rmm            );
  tree -> Branch ("tfb"           ,&tfb            );
  tree -> Branch ("trip"          ,&trip           );
  tree -> Branch ("trip_ch"       ,&trip_ch        );

  EvtCheck  = false;

  for(int imod = 0; imod < NumMod; imod++){
    for(int iview = 0; iview < 2; iview++){
      for(int itfb = 0; itfb < NumTFB; itfb++){
        for(int ich = 0; ich < LayerNumCh; ich++){
          HighGain       [imod][iview][itfb][ich] = -1.e-5;
          HighPed        [imod][iview][itfb][ich] = -1.e-5;
          LowPed         [imod][iview][itfb][ich] = -1;
          TdcThre        [imod][iview][itfb][ich] = -1.e-5;
          fMinAdcwTdcCut [imod][iview][itfb][ich] = 1000; //because this is used for minimum value
          GoodCh         [imod][iview][itfb][ich] = false;

          if(imod==15){
            fH_sHighAdc     [imod][iview][itfb][ich]
              = new TH1F( Form("fH_sHighAdc%02d%d%02d%02d",imod,iview,itfb,ich),
                  Form("High Gain ADC Mod.%02d view%d Plane%02d Ch%02d",imod,iview,itfb,ich),
                  150, 100, 250);
          }
          else{
            fH_sHighAdc     [imod][iview][itfb][ich]
              = new TH1F( Form("fH_sHighAdc%02d%d%02d%02d",imod,iview,itfb,ich),
                  Form("High Gain ADC Mod.%02d view%d Plane%02d Ch%02d",imod,iview,itfb,ich),
                  150, 100, 250);
          }
          fH_sLowAdc      [imod][iview][itfb][ich]
            = new TH1F( Form("fH_sLowAdc%02d%d%02d%02d",imod,iview,itfb,ich),
                Form("Low Gain ADC Mod.%02d view%d Plane%02d Ch%02d",imod,iview,itfb,ich),
                1024, 0, 1024);
        }//ich
      }//itfb
    }//iview
  }//imod
}
//___________________________________________________________________
//
//
//___________________________________________________________________
void MakeHist(){
  cout << "Making histograms" << endl;
  for(int imod = 0; imod < NumMod; imod++){
    for(int iview = 0; iview < 2; iview++){
      for(int itfb = 0; itfb < NumTFB; itfb++){
        for(int ich = 0; ich < LayerNumCh; ich++){
          if(imod==15){
            fH_HighAdc  [imod][iview][itfb][ich]
              = new TH1F( Form("fH_HighAdc%02d%d%02d%02d",imod,iview,itfb,ich),
                  Form("High Gain ADC Mod.%02d view%d Plane%02d Ch%02d",imod,iview,itfb,ich),
                  150, 100, 250);
          }
          else{
            fH_HighAdc  [imod][iview][itfb][ich]
              = new TH1F( Form("fH_HighAdc%02d%d%02d%02d",imod,iview,itfb,ich),
                  Form("High Gain ADC Mod.%02d view%d Plane%02d Ch%02d",imod,iview,itfb,ich),
                  150, 100, 250);
          }
          fH_LowAdc      [imod][iview][itfb][ich]
            = new TH1F( Form("fH_LowAdc%02d%d%02d%02d",imod,iview,itfb,ich),
                Form("Low Gain ADC Mod.%02d view%d Plane%02d Ch%02d",imod,iview,itfb,ich),
                1024, 0, 1024);
        }//ich
      }//itfb
    }//iview
  }//imod
}



//____________________ Start of Close ________________________________________ 
void Close(){
  cout << "close current calib. root file..." << endl;

  for(int imod = 0; imod < NumMod; imod++){
    for(int iview = 0; iview < 2; iview++){
      for(int itfb = 0; itfb < NumTFB; itfb++){
        for(int ich = 0; ich < LayerNumCh; ich++){
          fH_sHighAdc[imod][iview][itfb][ich]->Add(fH_HighAdc[imod][iview][itfb][ich]);
          fH_sLowAdc[imod][iview][itfb][ich] ->Add(fH_LowAdc [imod][iview][itfb][ich]);
        }//ich
      }//itfb
    }//iview
  }//imod

  tree      -> Fill();
  tree      -> Write();
  fTFile    -> Write();
  fTFile    -> Close();
}
//____________________ End of Close __________________________________________ 
//____________________________________________________________________________

//____________________ Start of Clear ________________________________________ 
void HistClear(){
  cout << "clear the variables and histograms...." << endl;
  for(int imod = 0; imod < NumMod; imod++){
    for(int iview = 0; iview < 2; iview++){
      for(int itfb = 0; itfb < NumTFB; itfb++){
        for(int ich = 0; ich < LayerNumCh; ich++){

          fH_HighAdc[imod][iview][itfb][ich] -> Reset();
          fH_HighAdc[imod][iview][itfb][ich] -> SetName ( Form("fH_HighAdc%02d%d%02d%02d",imod,iview,itfb,ich) );
          fH_HighAdc[imod][iview][itfb][ich] -> SetTitle( Form("High Gain ADC Mod.%02d view%d Plane%02d Ch%02d",imod,iview,itfb,ich) );

          fH_LowAdc [imod][iview][itfb][ich] -> Reset();
          fH_LowAdc [imod][iview][itfb][ich] -> SetName ( Form("fH_LowAdc%02d%d%02d%02d",imod,iview,itfb,ich) );
          fH_LowAdc [imod][iview][itfb][ich] -> SetTitle( Form("Low Gain ADC Mod.%02d view%d Plane%02d Ch%02d",imod,iview,itfb,ich) );
        }//ich
      }//itfb
    }//iview
  }//imod
}
//_______________________________________________________________ 
//
//
//_________________________________________________________________
void StartCalib(ND::TND280RawEvent* re){

  ND::TRawDataHeader header = re->GetHeader();
  StartTime  = header.GetTimeStamp();
  ND::THandle<ND::TRunInfoBank> RunInfoBank;
  while ( RunInfoBank = re->GetMidasBank<ND::TRunInfoBank>("XRUN",RunInfoBank) ) {
    ND::TRunInfoBank& runinfo = re->UseMidasBank<ND::TRunInfoBank>("XRUN");
    StartEvt                  = runinfo.GetSeqNumber();
  }
  AnaEvt     =    0;
  MakeTime   = time(NULL);
  Version    = fana_MPPC -> get_version();

}
//_______________________________________________________________ 
//
//
//_______________________________________________________________ 
void EndCalib(ND::TND280RawEvent* re){
  ND::TRawDataHeader header = re->GetHeader();
  EndTime  = header.GetTimeStamp();
  ND::THandle<ND::TRunInfoBank> RunInfoBank;
  while ( RunInfoBank = re->GetMidasBank<ND::TRunInfoBank>("XRUN",RunInfoBank) ) {
    ND::TRunInfoBank& runinfo = re->UseMidasBank<ND::TRunInfoBank>("XRUN");
    EndEvt                    = runinfo.GetSeqNumber();
  }
}
//_______________________________________________________________ 
//
//
//_______________________________________________________________ 

int evtnum(ND::TND280RawEvent* re){
  ND::THandle<ND::TRunInfoBank> RunInfoBank;
  while ( RunInfoBank = re->GetMidasBank<ND::TRunInfoBank>("XRUN",RunInfoBank) ) {
    ND::TRunInfoBank& runinfo = re->UseMidasBank<ND::TRunInfoBank>("XRUN");
    int num = runinfo.GetSeqNumber();
    return num;
  }
  return -1;
}

int trgid(ND::TND280RawEvent* re){

  ND::THandle<ND::TTriggerBank> triggerBank;
  while ( triggerBank = re->GetMidasBank<ND::TTriggerBank>("ITRI",triggerBank) ) {
    ND::TTriggerBank& trigger = re->UseMidasBank<ND::TTriggerBank>("ITRI");
    return (trigger.GetTriggerWord()>>48)&0xffff;
  }
  return -1;
}
