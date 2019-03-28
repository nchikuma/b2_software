#ifndef __DSTMAKER_HXX__
#define __DSTMAKER_HXX__

// ND280 software includes
#include "TMidasBank.hxx"
#include "TMidasFile.hxx"
#include "TMidasBankProxy.hxx"
#include "TMidasBankProxyRegistry.hxx"
#include "TND280RawEvent.hxx"
#include "TRawDataHeader.hxx"
// oaRawEvent includes
#include "TTripTDigitBank.hxx"
#include "TRunInfoBank.hxx"
#include "TMidasTripTDigitItr.hxx"
#include "TMidasTripTDigit.hxx"
#include "TMCMBank.hxx"
#include "TTriggerBank.hxx"
// ROOT includes
#include "TApplication.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TCanvas.h"
#include <TStyle.h>
#include "TString.h"
#include "TSystem.h"
#include "TSpectrum.h"
#include "TTree.h"
//C++ libraly includes
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <deque>
#include <vector>
#include <sys/stat.h>
#include <unistd.h> // using getopt      


#include "setup.hxx"

#include "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/lib/EVENTSUMMARY.h"
//#include "../../../../lib/HitSummary.h"

using namespace std;



FileStat_t    fs;

Int_t         NumEvt;
Int_t         TrgId;
Int_t         cAnaEvt;
Int_t         cAnaTrg;

TFile*        fTFile;
TTree*        tree;
EventSummary* summary;
HitSummary*   hitsum;
bool          cosmic;
bool          rename_input;
bool          rename_output;

class Hit{
  public:
    Int_t       mod;
    Int_t       pln;
    Int_t      view;
    Int_t        ch;
    Int_t       adc;
    Long_t      tdc;  
    Long_t       t0;
};


////__________________________________________________________

bool TDCHit(long tdc){
  if(tdc<16777201)
    return true;
  else
    return false;
}

void Book(Int_t run, Int_t subrun,char* output, Int_t qrun=0){

  char buff[1000];
  if(rename_output== false){ 
    if(cosmic){
      sprintf(buff,"%s/ingrid_%08d_%04d.root",data_cosmic_dir.c_str(),run,subrun);
    }
    else{
      sprintf(buff,"%s/ingrid_%08d_%04d.root",data_dst_dir.c_str(),run,subrun);
    }
  }
  else{
    sprintf(buff,"%s", output);
  }

  cout << "Make output ROOT files: " << buff << endl;
  fTFile   = new TFile(buff,"recreate");
  tree     = new TTree("tree","tree");
  summary  = new EventSummary();
  hitsum= new HitSummary();
  tree->Branch("fDefaultReco.", "EventSummary",
      &summary, 64000, 99);
}


void Save(){
  tree->AutoSave();
}

void WriteOnly(){
  cout << "treewrite" << endl;
  tree->Write();
  cout << "filewrite" << endl;
  fTFile->Write();

}

void Write(){
  cout << "treewrite" << endl;
  tree->Write();
  cout << "filewrite" << endl;
  fTFile->Write();
  fTFile->Close();
}


void Read(ND::TND280RawEvent* re){

  //Header
  ND::TRawDataHeader header = re->GetHeader();
  summary -> time           = header.GetTimeStamp();

  //TRunInfo
  ND::THandle<ND::TRunInfoBank> RunInfoBank;
  while ( RunInfoBank = re->GetMidasBank<ND::TRunInfoBank>("XRUN",RunInfoBank) ) {
    ND::TRunInfoBank& runinfo = re->UseMidasBank<ND::TRunInfoBank>("XRUN");
    NumEvt = runinfo.GetSeqNumber();
    summary -> event = NumEvt; //EVENT Number
  }

  //MCMBank
  ND::THandle<ND::TMCMBank> mcmBank;
  while ( mcmBank = re->GetMidasBank<ND::TMCMBank>("IMCM",mcmBank) ) {
    ND::TMCMBank& mcm = re->UseMidasBank<ND::TMCMBank>("IMCM");
    summary -> trgtime = mcm.GetUnixTimeSSecTrig();
    TrgId                  = (mcm.GetTriggerWord()>>48) & 0xffff;
    summary -> trgid       = TrgId; //Trigger ID, 1:beam, 2:Periodic, 128:Cosmic 
    summary -> nd280nspill = (mcm.GetTriggerWord()>>32) & 0xffff;
  }

  //###### Trip-t Bank ##### ########
  //################################# 
  if(TrgId==cAnaTrg){

    ND::THandle<ND::TTripTDigitBank> triptBank;
    while ( triptBank = re->GetMidasBank<ND::TTripTDigitBank>("",triptBank) ) {
      ND::TMidasTripTDigitItr itr(triptBank->GetMidasTripTDigitItr());
      while ( ! itr.EOD() ) {
        ND::TMidasTripTDigit digit(itr.Get());
        Int_t rmm     =  digit.GetRMMNum();
        Int_t tfb     =  digit.GetTFBNum();
        Int_t trip    =  digit.GetTripTNum();
        Int_t trip_ch =  digit.GetChannelNum();
        //Int_t cycle   =  digit.GetIntegrationNum();
        Int_t mod,view,plane,ch;

        bool  read = isWAGASCIconfig(&rmm,&tfb,&trip,&trip_ch,&mod,&view,&plane,&ch);
        if(!read)continue;

        long tdc = digit.GetTimeOffset();
        long t0  = digit.GetTimeOffsetT0();

        if( !TDCHit( tdc ) ){ continue; }

        read = read && !is_BadCh(mod,view,plane,ch);

        if(read)
        { 
          //const ND::TRawDataHeader& h = re->GetHeader();
          //int evno  = h.GetSeqNo();                 // Event Sequence Number
          int iint  = digit.GetIntegrationNum();    // = Capacitor number
          int icoff = triptBank->GetTFBStartCycle();
          int co    = iint - icoff;
          if (co<0) co += 23;
          if ( cosmic && ( co!=14&&co!= 15 ) )continue;

          double xy = -1.0e-5;
          double z  = -1.0e-5;
          if(mod==15){ fINGRID_Dimension->get_pos_loli_xy(mod,view,plane,ch,&xy,&z); }
          else       { fINGRID_Dimension->get_posXY      (mod,view,plane,ch,&xy,&z); }

          hitsum -> Clear("C");
          hitsum -> mod    = mod;
          hitsum -> view   = view;
          hitsum -> pln    = plane;
          hitsum -> ch     = ch;
          hitsum -> xy     = xy;
          hitsum -> z      = z;
          hitsum -> cyc    = co;
          hitsum -> adc    = digit.GetHighGainADC();
          hitsum -> loadc  = digit.GetLowGainADC();
          hitsum -> tdc    = tdc + t0;
          summary   -> AddModHit( hitsum, mod, co );

          //#### Speciall for VETO ######
          //#### Horizontal Right VETO #####
          if( 0<=mod && mod <=5 ){
            if( plane == 12 ){
              plane = 11;
              hitsum -> mod = mod+1;
              hitsum -> pln = plane;
              fINGRID_Dimension -> get_posXY( mod+1, view, plane, ch,
                  &xy, &z	 );
              hitsum -> xy     = xy;
              hitsum -> z      = z;
              summary   -> AddModHit( hitsum, mod+1, co );
            }
          }
          //#### Vertical Bottom VETO #####
          if( 7<=mod && mod <=12 ){
            if( plane == 14 ){
              hitsum -> mod = mod+1;
              plane = 13;
              hitsum -> pln = plane;
              fINGRID_Dimension -> get_posXY( mod+1, view, plane, ch,
                  &xy, &z	 );
              hitsum -> xy     = xy;
              hitsum -> z      = z;
              summary   -> AddModHit( hitsum, mod+1, co );
            }
          }
        }//read
      }//itr
    }//tribBank
  }//TrgId
}

//_________________________________________________________________
void ProcessFile(const char *FileName) {
  ND::TMidasFile mf;
  mf.Open(FileName);
  cout<<"loop all event..."<<endl;

  set_BadCh(BadCh_Ana);

  int count=0;
  while (ND::TND280RawEvent* re = mf.ReadRawEvent()) {
    re->PromoteMidasBanks(false);
    TrgId=-1;

    if(re->size()!=0){//eliminate header and laster
      summary -> Clear("C");
      count++;
      Read(re);
      if(NumEvt%10000==0){ cout<<"Event:"<<NumEvt<<endl; }
      if(cAnaTrg==TrgId){
        tree -> Fill();
      }
    }
    delete re;
    if(cAnaEvt!=-1&&NumEvt>cAnaEvt)break;
    if(NumEvt%10000==0)Save();
  }// End loop over events

}
//______________________________________________



#endif
