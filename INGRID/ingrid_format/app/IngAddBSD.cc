////##### Standard C++ lib. ######
//#include<iostream>
//#include<sstream>
//#include<fstream>
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
////##### INGRID Library #########
//#include "INGRIDEVENTSUMMARY.h"
//#include "BeamInfoSummary.h"
//#include "IngridBasicReconSummary.h"
//#include "IngridConstants.h"
////##### BSD acc. Library #######/


#include "header.hh"

BEAM_DATA bsd;
BeamData& beam = BeamData::instance();
bool GoodSpillv01();
bool GoodSpill();
bool GoodSpilltemp();
bool RemoveEvtSkew();

void PrintUsage()
{
  cout << " Usage :" << endl;
  cout << " ./IngAddBSD -r <run#> -s <srun#> -m <mode>" << endl;
  cout << "    mode : 0 -> SS floor, 1 -> B2 floor."    << endl;
  cout << "  Option: -v -> For BSD analysis."           << endl;
}

int main(int argc, char *argv[]){

  beam.init_db();
  char FileName[300], Output[300];
  int    c  =     -1;
  int    run_number  = -1;
  int    srun_number = -1;
  int    mode        = -1;
  bool   v01  = false;

  while ((c = getopt(argc,argv,"r:s:m:v")) != -1) {
    switch(c){
    case 'r':
      run_number=atoi(optarg);
      break;
    case 's':
      srun_number=atoi(optarg);
      break;
    case 'm':
      mode=atoi(optarg);
      break;
    case 'v':
      v01 = true;
      break;
    }
  }
  if((mode!=0&&mode!=1)
      ||run_number<0
      ||srun_number<0
      )
  {
    PrintUsage();
    return 0;
  }

  //###### Read INGRID Data #######
  //###############################
  EventSummary* summary  = new EventSummary();

  FileStat_t fs;
  sprintf(FileName, "%s/ingrid_%08d_%04d_spill.root",
      data_dst_dir.c_str(), run_number, srun_number);
  if(gSystem->GetPathInfo(FileName, fs ) ){
    cout << "Cannot open input file : " << FileName << endl;
    exit(1);
  }
  TFile*            rfile      = new TFile( FileName, "read" );
  TTree*             tree      = (TTree*)rfile -> Get("tree");
  TBranch*          EvtBr      = tree->GetBranch("fDefaultReco.");
  EvtBr                        ->SetAddress(&summary);
  tree                         ->SetBranchAddress("fDefaultReco.", &summary);
  int                nevt      = (int)tree -> GetEntries();
  cout << "Total # of events = " << nevt <<endl;

  //###### New INGRID Data ########
  //###############################
  BeamInfoSummary*        bsum  = new BeamInfoSummary();
  sprintf(Output, "%s/ingrid_%08d_%04d_bsd.root",
        data_dst_dir.c_str(), run_number, srun_number);
  TFile*                 wfile  = new TFile( Output,"recreate");
  TTree*                 wtree  = new TTree("tree", "tree");
  EventSummary* wsummary        = new EventSummary(); 
  wtree                         -> Branch   ("fDefaultReco.",
					     "EventSummary", 
					     &wsummary,  64000,   99);
  wtree -> SetMaxTreeSize(10000000000);

  cout << "Processing... Run#" << run_number << "_" << srun_number << endl;
  cout << "Input file : " << FileName << endl;
  cout << "Output file : " << Output << endl;
  for(int ievt = 0; ievt < nevt; ++ievt){  
    if(ievt%1000==0) cout<<" "<<ievt<<flush;
    wsummary -> Clear("C");
    summary  -> Clear("C");
    bsum     -> Clear("C");
    tree     -> GetEntry(ievt);

    int nd280nspill = summary -> nd280nspill;
    int       utime = summary -> time;

    if( ! beam.get_spill( nd280nspill, utime, 10 )   ){ continue; }
    if( abs( beam.trigger_time_sec(0) - utime ) > 10 ){ continue; }

    bool ok = false;
    if     (  v01 && GoodSpillv01() ){ ok = true;}
    else if( !v01 && GoodSpill()    ){ ok = true;}
    
    if(!ok){ continue; }

    bsd = beam.data();
  
    //### Set bunch_flag ####
    for(int cyc=0; cyc<23; cyc++){
      if(4<=cyc&&cyc<=11){
	if( bsd.ct_np [4][ cyc - 4 + 1] > 1e11 )  //Use CT05, all bunches
        {
          summary -> bunch_flag[cyc] = true;
        }
        else{
          summary -> bunch_flag[cyc] = false;
        }
      }
      else{
        summary -> bunch_flag[cyc]=false;
      }
    }
    //#### Fill ########
    //##################
    bsum -> nurun           = bsd.nurun;
    bsum -> spill_flag      = bsd.spill_flag;
    bsum -> good_spill_flag = bsd.good_spill_flag;
    bsum -> trg_sec         = bsd.trg_sec[0];
    bsum -> spillnum        = bsd.spillnum;
    bsum -> nd280spillnum   = ( (bsd.spillnum) & 0xffff ) + 1;
    //bsum -> good_spill_flag = bsd.good_spill_flag;
    //bsum -> good_spill_flag = 1;
    bsum -> run_type        = bsd.run_type;
    for(int l=0; l<3; l++){
      bsum -> target_eff[l] = bsd.target_eff[l];
    }
    for(int l=0; l<5; l++){
      for(int k=0; k<9; k++){
	bsum -> beam_time[l][k] = bsd.beam_time[l][k];
      }
    }
    for(int l=0; l<5; l++){
      for(int k=0; k<9; k++){
	bsum -> ct_np    [l][k] = bsd.ct_np    [l][k];
      }
    }
    for(int l=0; l<3; l++){
      for(int k=0; k<5; k++){
	bsum -> hct[l][k]       = bsd.hct[l][k];
      }
    }
    for(int l=0; l<12; l++){
      bsum -> mumon[l]          = bsd.mumon[l];
    }
    for(int l=0; l<13; l++){
      bsum -> otr[l]            = bsd.otr[l];
    }

    for(Int_t ibas=0; ibas < summary->NBasicRecons(); ibas++ ){
      BasicReconSummary* ingbasrec = (BasicReconSummary*)( summary->GetBasicRecon(ibas) );
      float CTtime = bsum->beam_time[0][ingbasrec->hitcyc-bunch1st_cyc + 1] - GapbwBunch * 1e-9 * ( ingbasrec->hitcyc - bunch1st_cyc);
      ingbasrec -> clstimecorr = ingbasrec -> clstime - ( CTtime - CTtimeBase )*1e9; 
      ingbasrec -> exptime = 
	((int)ingbasrec -> clstimecorr - TDCOffset )% GapbwBunch - fExpBeamTime;
      if( fabs( ingbasrec->exptime ) < beamontime && summary->bunch_flag[ingbasrec->hitcyc] )
	ingbasrec -> ontime = true;
      else
	ingbasrec -> ontime = false;
    }//ibas


    summary -> AddBeamSummary(bsum);
    wsummary = summary;
    wtree   -> Fill();
  }
  wtree   -> Write();
  wfile   -> Write();
  wfile   -> Close();

  cout << "IngAddBSD finished" << endl;
}


bool GoodSpillv01(){
  if( beam.good_spill() )
    return true;
  else
    return false;
}

bool GoodSpill(){
  if( ! beam.spill_flag() )
    return false;

  RemoveEvtSkew();
   
  bool normal_cond = ( beam.run_type() == 1 ) &&
    ( fabs( fabs(beam.horn_current(0) ) - 250 ) < 5 ) && 
    ( fabs( fabs(beam.horn_current(1) ) - 250 ) < 5 ) && 
    ( fabs( fabs(beam.horn_current(2) ) - 250 ) < 5 );
  
  bool beamtrg     = ( beam.ct_np(4) > 1e11 );
  
  bool mucenter    = 
    ( sqrt( beam.mumon_si_x() * beam.mumon_si_x() + beam.mumon_si_y() * beam.mumon_si_y() ) < 10 );
  
  bool muct        =
    ( 31.7 < ( beam.mumon_si_totq() / beam.ct_np(4) * 1e9 ) ) &&
    ( ( beam.mumon_si_totq() / beam.ct_np(4) * 1e9 ) < 35.1 );
  if( normal_cond && beamtrg && mucenter && muct )
    return true;
  else
    return false;
}

bool RemoveEvtSkew()
{
  //if( 598833  <= beam.spill_number() && beam.spill_number() <= 598837  &&
  //    beam.trigger_time_sec(0) <= 1267455600 ) // till Mar. 2nd (this is bas spill at MR run 31) 
  //  return false;
  //if( 598937  <= beam.spill_number() && beam.spill_number() <= 598941  &&
  //    beam.trigger_time_sec(0) <= 1267455600 ) // till Mar. 2nd (this is bas spill at MR run 31) 
  //  return false;
  //if( 1201242 <= beam.spill_number() && beam.spill_number() <= 1201270 &&
  //    beam.trigger_time_sec(0) <= 1267455600 ) // till Mar. 2nd (this is bas spill at MR run 31) 
  //  return false;
  //if( 1202461 <= beam.spill_number() && beam.spill_number() <= 1202533 &&
  //    beam.trigger_time_sec(0) <= 1267455600 ) // till Mar. 2nd (this is bas spill at MR run 31) 
  //  return false;
  //if( 1234283 <= beam.spill_number() && beam.spill_number() <= 1234304 &&
  //    beam.trigger_time_sec(0) <= 1267455600 ) // till Mar. 2nd (this is bas spill at MR run 31) 
  //  return false;
  //if( 1235436 <= beam.spill_number() && beam.spill_number() <= 1235450 &&
  //    beam.trigger_time_sec(0) <= 1267455600 ) // till Mar. 2nd (this is bas spill at MR run 31) 
  //  return false;
  //if( 8022255 <= beam.spill_number() && beam.spill_number() <= 8022604 &&
  if( beam.trigger_time_sec(0) <= 1520843507 && beam.trigger_time_sec(0) >= 1520842641 ) // MR78 2018/03/13
    return false;
  if( beam.trigger_time_sec(0) <= 1525031280 && beam.trigger_time_sec(0) >= 1525030962 ) // MR79
    return false;

  return true;
}
