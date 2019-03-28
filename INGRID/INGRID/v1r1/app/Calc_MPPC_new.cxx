//###########################################################
//###     Make Calib. table and MPPC easy check file      ###
//###########################################################
//###                Made by N.Chikujma                   ###
//###########################################################

#include "Calc_MPPC_new.hxx"
Int_t NumEvt;

void Analysis(){
  cout<<"Analysis..."<<endl;
  set_BadCh(BadCh_Ana);
  for(int irmm=0;irmm<5;irmm++){
    for(int itfb=0; itfb<48; itfb++){  
      for(int itrip=0;itrip < 4; itrip++){
        for(int itrip_ch=0;itrip_ch<16;itrip_ch++){
          int imod, iview, ipln, numch;
          bool read = isWAGASCIconfig(&irmm,&itfb,&itrip,&itrip_ch,&imod,&iview,&ipln,&numch);
          read = read && !is_BadCh(imod,iview,ipln,numch);
          if(read){

            GoodCh   [imod][iview][ipln][numch] = true;
            if(fana_MPPC->analysis_old_version( fH_HighAdc[imod][iview][ipln][numch],imod))
            {
              fana_MPPC  -> analysis_logain     ( fH_LowAdc [imod][iview][ipln][numch] );
              HighPed      [imod][iview][ipln][numch] = fana_MPPC -> get_pedestal                ();
              LowPed       [imod][iview][ipln][numch] = fana_MPPC -> get_lowpedestal             ();
              HighPedSigma [imod][iview][ipln][numch] = fana_MPPC -> get_pedestal_sigma          ();
              HighGain     [imod][iview][ipln][numch] = fana_MPPC -> get_gain                    ();
              MPPCNoise    [imod][iview][ipln][numch] = fana_MPPC -> get_noise                   ();
              MPPCCandA    [imod][iview][ipln][numch] = fana_MPPC -> get_crosstalk_and_afterpulse();
              TdcThre      [imod][iview][ipln][numch] 
                = 1.0 * ( fMinAdcwTdcCut[imod][iview][ipln][numch] - HighPed[imod][iview][ipln][numch] )
                /HighGain[imod][iview][ipln][numch];
            }
            else{//no peak is found
              HighPed      [imod][iview][ipln][numch] = 0;
              LowPed       [imod][iview][ipln][numch] = 0;
              HighPedSigma [imod][iview][ipln][numch] = 0;
              HighGain     [imod][iview][ipln][numch] = -1;
              MPPCNoise    [imod][iview][ipln][numch] = 0;
              MPPCCandA    [imod][iview][ipln][numch] = 0;
              TdcThre      [imod][iview][ipln][numch] = 0;
            }

            mod           . push_back( imod     );
            view          . push_back( iview    );
            pln           . push_back( ipln     );
            ch            . push_back( numch    );
            rmm           . push_back( irmm     );
            tfb           . push_back( itfb     );
            trip          . push_back( itrip    );
            trip_ch       . push_back( itrip_ch );
            gain          . push_back( HighGain     [imod][iview][ipln][numch] );
            pedestal      . push_back( HighPed      [imod][iview][ipln][numch] );
            refgain       . push_back( HighGainRef  [imod][iview][ipln][numch] );
            refpedestal   . push_back( HighPedRef   [imod][iview][ipln][numch] );
            pedestal_sigma. push_back( HighPedSigma [imod][iview][ipln][numch] );
            noise         . push_back( MPPCNoise    [imod][iview][ipln][numch] );
            canda         . push_back( MPPCCandA    [imod][iview][ipln][numch] );
          }//active channel
        }//numch
      }//ipln
    }//view
  }//imod
}
//____________________________________________________________________
//
//
//____________________________________________________________________
void Event(ND::TND280RawEvent* re) 
{

  ND::THandle<ND::TRunInfoBank> RunInfoBank;
  while ( RunInfoBank = re->GetMidasBank<ND::TRunInfoBank>("XRUN",RunInfoBank) ) {
    ND::TRunInfoBank& runinfo = re->UseMidasBank<ND::TRunInfoBank>("XRUN");
    NumEvt = runinfo.GetSeqNumber();
  }

  // Loop over all banks of type TTripTHitBank
  ND::THandle<ND::TTripTDigitBank> triptBank;
  while ( triptBank = re->GetMidasBank<ND::TTripTDigitBank>("",triptBank) ) {

    // Create an iterator over digits
    ND::TMidasTripTDigitItr itr(triptBank->GetMidasTripTDigitItr());
    while ( ! itr.EOD() ) {

      ND::TMidasTripTDigit digit(itr.Get());
      Int_t imod, iview, ipln, ich;
      Int_t irmm      = digit.GetRMMNum();
      Int_t itrip     = digit.GetTripTNum();
      Int_t itrip_ch  = digit.GetChannelNum();
      Int_t itfb      = digit.GetTFBNum();

      bool read = isWAGASCIconfig(&irmm,&itfb,&itrip,&itrip_ch,&imod,&iview,&ipln,&ich);

      if(read){  // using channel
        Int_t highadc  =  digit.GetHighGainADC();
        Int_t loadc    =  digit.GetLowGainADC ();
        long  tdc      =  digit.GetTimeOffset ();
        fH_HighAdc[imod][iview][ipln][ich]  -> Fill(highadc);
        fH_LowAdc [imod][iview][ipln][ich]  -> Fill(loadc);

        if( NumEvt  > 0             && 
            tdc     < 16777201      && 
            highadc < fMinAdcwTdcCut[imod][iview][ipln][ich])
        {
          fMinAdcwTdcCut[imod][iview][ipln][ich] = highadc;
        }//Tdc threshold calibration	 
      }
    }   // End of loop over digits in this bank
  }     // End of loop over banks of digits in this event
}
//____________________________________________________________________
//
//
//_____________________________________________________________________
void ProcessFile(const char *FileName, char *Output, bool rename_output) 
{

  cout << "Start processing the file." << endl;

  ND::TMidasFile mf;
  mf.Open(FileName);

  MakeHist();


  cout << "Start event loop." << endl;
  bool  exist_file = false;
  while ( ND::TND280RawEvent* re = mf.ReadRawEvent() ) {
    re -> PromoteMidasBanks(false);
    NumEvt = evtnum(re); 
    if(NumEvt%10000==0){
      cout<< "Current event# = " << NumEvt;
      if(NumEvt==0) cout << " : The 1st events is removed.";
      cout << endl;
    }
    if (NumEvt   == 0      ||  //### first event or
        trgid(re)!= cAnaTrg )  //### trigger for not analysis 
    {
      delete re;
      continue;
    }

    if(!rename_output){
      sprintf(Output,
          "%s/ingrid_%08d_%04d_Calib%02d.root",
          data_calib_dir.c_str(),
          IngRun, IngSubRun, nCalib);
    }
    if( !exist_file ){//#### Start Calib. File #####
      Book(Output);
      StartCalib(re);
      exist_file = true;
      nCalib++;
    }

    if(trgid(re) == cAnaTrg){
      Event(re);
      AnaEvt++;
    }

    if( AnaEvt >= ANABREAK ){ //### End of Calib. File
      cout<<"start End Calib." 
        << " " << ANABREAK << " of events are analyzed." << endl;
      EvtCheck   = true; 
      if(!rename_output)exist_file = false; //### flag for start Calib. File
      EndCalib (re);      // Get End Time
      Analysis ();        // Analyze MPPC gain and so on
      Close    ();        // Write and Close ROOT file
      HistClear();        // Clear the contents of Histogram for Analysis and ROOT file
      CLEAR    ();
      delete re;
      return;             // Modified to avoid creating another Calib file.
    }
    delete re;
  }  // End loop over events
  if(!EvtCheck){
    cout << "This file do NOT have enough statistics to be analyzed." << endl;
    cout << " AnaEvt=" << AnaEvt << ", but " << ANABREAK << " of events are required."  << endl;
    if(fTFile){fTFile->Close();}
    HistClear();  
    CLEAR    ();
  }
  return;
}

//______________________________________________________________________


int main(int argc,char *argv[]){

  cout << "===============" << endl;
  cout << " Calc_MPPC_new " << endl;
  cout << "===============" << endl;

  int  c = -1;
  char FileName[300];
  char Output  [300];
  bool rename_input  = false;
  bool rename_output = false;
  cAnaTrg=1;

  NumEvt    = 0;
  IngRun    = -1;
  IngSubRun = -1;

  while ((c = getopt(argc, argv, "r:s:t:i:o:bwh")) != -1) {
    switch(c){

      case 'r':
        IngRun=atoi(optarg);
        break;
      case 's':
        IngSubRun=atoi(optarg);
        break;
      case 't':
        cAnaTrg=atoi(optarg);
        break;
      case 'i':
        sprintf(FileName,"%s",optarg);
        rename_input = true;
        break;
      case 'o':
        sprintf(Output,"%s",optarg);
        rename_output = true;
        break;
      case 'b':
        BadCh_Ana = true;
        break;
      case 'w':
        ANABREAK = 300;
        break;
      case 'h':
        cout<<"-r [run number]     "<<endl;
        cout<<"-s [sub run number] "<<endl;
        cout<<"-t [trigger ID]     "<<endl;
        cout<<"   1: bream trigger "<<endl;
        cout<<" 128: cosmic trigger"<<endl;
        cout<<"-b                  "<<endl;
        cout<<"  This option would be set for masking bad channels."<<endl;
        cout<<"-i [input file name]"<<endl;
        cout<<"  Input file can be directly set by this option"
          <<", instead of putting run/srun number."<<endl;
        cout<<"-o [output file name]"<<endl;
        cout<<"  If this option is not set"
          <<", output file name is automatically set."<<endl;
        cout<<"-w       : To set the number of analyzed events as 300"<<endl;
        cout<<"           *Default=900" << endl;
        exit(0);
        break;
      case '?':
        cout<<"Unknown option"<<endl;
        exit(1);
        break;
    }
  }//option end

  if( ((!rename_input)&&(IngRun==-1||IngSubRun==-1)) )
  {
    cout << "Put arguments; " << endl;
    cout << " - INGRID Run number/Sub-run number, otherwise input file name." << endl;
    cout << "Check 'help'" << endl;
    exit(0);
  }

  if(!rename_input){
    //int run_i = IngRun/1000;
    //string sub_dir = Form("%05d000_%05d999",run_i,run_i);
    sprintf(FileName,
        //"%s/%s/ingrid_%08d_%04d.daq.mid.gz",
        "%s/ingrid_%08d_%04d.daq.mid.gz",
        data_midas_dir.c_str(),
        //sub_dir.c_str(),
        IngRun,IngSubRun);
  }

  cout << "file name is :" << FileName << endl;
  if ( gSystem->GetPathInfo(FileName,fs) ) {
    cout << "Error: There is no such a file." << endl;
    exit(0);
  }

  fana_MPPC             = new ana_MPPC();
  fINGRID_Ch_config     = new INGRID_Ch_config();
  fPM_Ch_config         = new PM_Ch_config();
  fWM_Ch_config         = new WM_Ch_config();
  fB2ING_Ch_config      = new B2ING_Ch_config();

  GetRef();
  nCalib=0;
  ProcessFile(FileName,Output,rename_output);

  delete fana_MPPC        ;
  delete fINGRID_Ch_config;
  delete fPM_Ch_config    ;
  delete fWM_Ch_config    ;
  delete fB2ING_Ch_config ;

  return 0;
}


