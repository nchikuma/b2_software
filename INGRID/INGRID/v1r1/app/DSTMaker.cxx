//###########################################################
//###            MIDAS File -> INGRID ROOT File           ###
//###########################################################
//###                Made by M.Otani                      ###
//###########################################################
#include"DSTMaker.hxx"

void Read       (ND::TND280RawEvent* re);
void ProcessFile(const char *FileName);

int run_number, subrun;
//_______________________________________________________________
int main(int argc,char *argv[]){

  cout << "===============" << endl;
  cout << " DSTMaker      " << endl;
  cout << "===============" << endl;

  int  c=-1;
  char FileName[1000] = "hoge.mid";  //Midas file name
  char Output  [1000] = "hoge.root"; //output ROOT file name
  fINGRID_Ch_config = new INGRID_Ch_config();
  fPM_Ch_config     = new PM_Ch_config();
  fWM_Ch_config     = new WM_Ch_config();
  fB2ING_Ch_config  = new B2ING_Ch_config();
  fINGRID_Dimension = new INGRID_Dimension();
  run_number = -1;
  subrun     = -1;
  cAnaEvt    = -1;
  cosmic        = false;  
  rename_input  = false;
  rename_output = false;
  while ((c = getopt(argc, argv, "hr:s:t:be:i:o:")) != -1) {
    switch(c){
    case 'r':
      run_number=atoi(optarg);
      break;
    case 's':
      subrun=atoi(optarg);
      break;
    case 't':
      cAnaTrg=atoi(optarg);
      break;
    case 'b':
      BadCh_Ana=true;
      break;
    case 'd':
      cAnaEvt=atoi(optarg);
      break;
    case 'i':
      sprintf(FileName,"%s",optarg);
      rename_input = true;
      break;
    case 'o':
      sprintf(Output,"%s",optarg);
      rename_output = true;
      break;
    case 'h':
      cout<<"-r [run number]     "<<endl;
      cout<<"-s [sub run number] "<<endl;
      cout<<"-t [trigger ID]     "<<endl;
      cout<<"   1: bream trigger "<<endl;
      cout<<" 128: cosmic trigger"<<endl;
      cout<<"-i [input file name]"<<endl;
      cout<<"  Input file can be directly set by this option"
        <<", instead of putting run/srun number."<<endl;
      cout<<"-o [output file name]"<<endl;
      cout<<"  If this option is not set"
        <<", output file name is automatically set."<<endl;
      cout<<"-e [end event number]"<<endl;
      cout<<"   This option would be set for breaking in the middle of file"<<endl;
      exit(0);
    case '?':
      cout<<"Unknown option"<<endl;
      exit(1);
      break;
    }
  }//option end

  if( (!rename_input)&&
      (run_number==-1||subrun==-1) )
  {
    cout << "Check the usage:" << endl;
    cout << "./DSTMaker -h"    << endl;
    exit(0);
  }


  if(cAnaTrg==128){ cosmic = true; } 


  if(!rename_input){
    //int run_i = run_number/1000;
    //string sub_dir = Form("%05d000_%05d999",run_i,run_i);
    sprintf(FileName,
        //"%s/%s/ingrid_%08d_%04d.daq.mid.gz",
        "%s/ingrid_%08d_%04d.daq.mid.gz",
        data_midas_dir.c_str(),
        //sub_dir.c_str(),
        run_number,subrun);
  }

  if ( gSystem->GetPathInfo(FileName,fs) ) {
    std::cout << "Cannot find file: " << FileName << std::endl;
    exit(0);
  }
  cout<<"Input MIDAS file: " << FileName << endl;
  cout<<"Start Book"<<endl;
  Book(run_number, subrun, Output);

  cout<<"Stat Process"<<endl;
  ProcessFile(FileName);
  cout<<"Start Write"<<endl;
  Write();

  delete fINGRID_Dimension;       
  delete fWM_Ch_config;           
  delete fPM_Ch_config;           
  delete fB2ING_Ch_config;           
  delete fINGRID_Ch_config;
  return 0;
}


//_________________________________________________________________


