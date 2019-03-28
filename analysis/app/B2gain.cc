#include "B2gain.hxx"
#define MAXCH 4000

int main(int argc, char** argv){

  int run_start = -1;
  int run_stop  = -1;
  bool BadCh_Ana = true;
  bool Only_plot = true;

  string calibfiledir = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_calib";
  string outputfile   = "B2gain.root";

  int c=-1;
  while((c=getopt(argc,argv,"i:l:d:o:bph")) != -1){
    switch(c){
      case 'i':
        run_start = std::atoi(optarg);
        break;
      case 'l':
        run_stop = std::atoi(optarg);
        break;
      case 'd':
        calibfiledir = optarg;
        break;
      case 'o':
        outputfile = optarg;
        break;
      case 'b':
        BadCh_Ana = false;
        break;
      case 'p':
        Only_plot = false;
        break;
      case 'h':
        std::cout << "Usage: " << argv[0] << " <option>" << std::endl;
        std::cout << " -h : To show this help.            "  << std::endl;
        std::cout << " -i <(int) run#>      : To set initial run number.  "  << std::endl;
        std::cout << " -l <(int) run#>      : To set last run number.     "  << std::endl;
        std::cout << " -d <(char) dir name> : To set input file directory."  << std::endl;
        std::cout << "      *Currently set as " << calibfiledir              << std::endl;
        std::cout << " -o <(char) file name>: To set output file name.    "  << std::endl;
        std::cout << "      *Currently set as " << outputfile                << std::endl;
        std::cout << " -b                   : To open bad channel masks.  "  << std::endl;
        std::cout << " -p                   : To remake the calib file.   "  << std::endl;
        exit(0);
      default:
        std::cout << "See the help: " << argv[0] << " -h" << std::endl;
        exit(0);
      }
  }

  if(run_start<0||run_stop<0){
    std::cout << "Put arguments for run numbers" << std::endl;
    std::cout << "See the help: " << argv[0] << " -h" << std::endl;
    exit(0);
  }
  

  std::cout << "CALIB FILE DIR : " << calibfiledir << std::endl;
  std::cout << "OUTPUT FILE    : " << outputfile   << std::endl;
  
  if(!Only_plot){

    TFile*  wfile  = new TFile (outputfile.c_str(), "recreate" );
    TChain* wchain = new TChain("calibtree");


    TFile* fTFile;
    TTree* fTTree;
    char   FileName[300];
    bool   EvtCheck;

    FileStat_t fs;
    for(int irun=run_start; irun <= run_stop; irun++ ){
      int tsrun=0;
      std::cout << "=== read RUN " << irun << " ===" << std::endl;
      for(int isrun=tsrun; isrun <= 200; isrun++ ){

        sprintf( FileName, "%s/ingrid_%08d_%04d_Calib00.root",calibfiledir.c_str(),irun,isrun);
        if( gSystem->GetPathInfo(FileName, fs) ){ continue; }

        fTFile = new TFile( FileName, "read");
        fTTree = (TTree*)fTFile->Get("calibtree");
        if(fTTree==0){ continue; }

        std::cout << "reading " << FileName << std::endl;
        fTTree ->SetBranchAddress("EvtCheck", &EvtCheck);
        fTTree ->GetEntry(0);
        if(EvtCheck){
          wfile  -> cd ();
          wchain -> Add( FileName );
        }

        fTFile ->Close();
      }//isrun
    }//irun
    wfile ->cd();

    wchain->Write();
    wfile ->Write();
    wfile ->Close();
    std::cout << "Output file was closed." << std::endl;

  }

  plot_gain_history(outputfile,BadCh_Ana);

  return 0;

}

void plot_gain_history(string filename,bool BadCh_Ana)
{
  std::cout << "Start plotting gain history." << std::endl;

  vector<int> v_bad_mod (0);
  vector<int> v_bad_view(0);
  vector<int> v_bad_pln (0);
  vector<int> v_bad_ch  (0);
  int num_v_bad = 0;
   
  FileStat_t fs;
  if( gSystem->GetPathInfo(filename.c_str(), fs) ){ return; }

  INGRID_BadCh_mapping* badch = new INGRID_BadCh_mapping();
  badch->set_BadCh(BadCh_Ana);


  cout << "Reading a ROOT file : " << filename << endl;
  TFile*  file  = new TFile(filename.c_str(),"read");
  TChain* chain = (TChain*)file->Get("calibtree");
  if(chain==0){ return; }
  int ntrees  = chain->GetNtrees();
  chain->GetEntry(0);
  TTree* starttree = chain->GetTree();
  int starttime = starttree->GetMinimum("StartTime");
  chain->GetEntry(ntrees-1);
  TTree* endtree = chain->GetTree();
  int endtime = endtree->GetMinimum("StartTime");
  int nbin = (int)((endtime-starttime)/3600*2.8*2);

  // ==== Canvas setup ==== //
  gStyle->SetOptStat(0);
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);

  //TFile *wfile = new TFile("B2gain_plot.root","recreate");
  TFile *wfile = new TFile("B2gain_tree.root","recreate");
  //TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,400);
  const double min_gain = -2.;
  const double max_gain = 13.;
  const int    nbin_gain = (max_gain-min_gain)/0.1;
  const int nummod = 4;
  const int maxid  = 1280;
  string modname = "";
  string name = "";
  //TH2F *h1[nummod][maxid];
  double gain[nummod][maxid];
  int    start_time;
  TTree *wtree = new TTree("calibtree","calibtree");
  wtree->Branch("gain",gain,Form("gain[%d][%d]/D",nummod,maxid));
  wtree->Branch("start_time",&start_time,"start_time/I");

  cout << "Setting histograms ..." << endl;
  for(int mod=0;mod<nummod;mod++){
    int im;
    if     (mod==0){im= 3;modname="INGmod3";}
    else if(mod==1){im=15;modname="INGWM"  ;}
    else if(mod==2){im=16;modname="PM"     ;}
    else if(mod==3){im=14;modname="B2ING"  ;}
    else           {im=99;modname="Unknown";}
    for(int id=0;id<maxid;id++){
      int iv,ip,ic;
      if(getchannel_id(im,id,iv,ip,ic)){
        //cout
        //  << modname
        //  << " id=" << id
        //  << " iv=" << iv
        //  << " ip=" << ip
        //  << " ic=" << ic
        //  << endl;
        //name = Form("%s_%02d%01d%02d%02d",modname.c_str(),im,iv,ip,ic);
        //h1[mod][id] = new TH2F(name.c_str(),name.c_str(),
        //    nbin,starttime,endtime,nbin_gain,min_gain,max_gain);
        //h1[mod][id]  -> GetXaxis() -> SetTimeDisplay(1);
        //h1[mod][id]  -> GetXaxis() -> SetTimeOffset(0,"jst");
        //h1[mod][id]  -> SetXTitle("time");
        //h1[mod][id]  -> SetYTitle("MPPC gain (ADC counts)");
        //h1[mod][id]  -> SetTitle (name.c_str());
        gain[mod][id] = 0;
      }
    }
  }

  for(int itree=0;itree<ntrees;itree++){
    chain->GetEntry(itree);
    TTree* tree = chain->GetTree();

    TBranch *branch = 0;
    int t_StartTime;
    std::vector<double> *t_gain = 0;
    std::vector<int>    *t_mod  = 0;
    std::vector<int>    *t_view = 0;
    std::vector<int>    *t_pln  = 0;
    std::vector<int>    *t_ch   = 0;
    tree->SetBranchAddress("StartTime",&t_StartTime);
    tree->SetBranchAddress("gain"     ,&t_gain     ,&branch);
    tree->SetBranchAddress("mod"      ,&t_mod      ,&branch);
    tree->SetBranchAddress("view"     ,&t_view     ,&branch);
    tree->SetBranchAddress("pln"      ,&t_pln      ,&branch);
    tree->SetBranchAddress("ch"       ,&t_ch       ,&branch);

    tree->GetEntry(0);
    Long64_t entry = tree->LoadTree(0);
    branch->GetEntry(entry);
    start_time = t_StartTime;
    std::cout << "StartTime=" << std::setw(10) << start_time << std::endl;

    int nchannels = (*t_gain).size();
    for(int ich=0;ich<nchannels;ich++){
      double c_gain = (*t_gain)[ich];
      int    c_mod  = (*t_mod )[ich];
      int    c_view = (*t_view)[ich];
      int    c_pln  = (*t_pln )[ich];
      int    c_ch   = (*t_ch  )[ich];

      if(badch->is_BadCh(c_mod,c_view,c_pln,c_ch)){continue;}
      //double lim_low = -1;
      //double lim_up  = -1;
      int c_id;
      if(getid_channel(c_mod,c_view,c_pln,c_ch,c_id)){
        int im;
        if     (c_mod== 3){im=0;}
        else if(c_mod==15){im=1;}
        else if(c_mod==16){im=2;}
        else if(c_mod==14){im=3;}
        gain[im][c_id] = c_gain;

        //cout
        //  << " mod="  << c_mod
        //  << " c_mod=" << c_view
        //  << " c_view="  << c_pln
        //  << " c_pln=" << c_ch
        //  << " c_ch="   << c_id
        //  << endl;
        //int imod;
        //if     (c_mod== 3){ imod=0; lim_low=9.2;lim_up=11.5;}
        //else if(c_mod==15){ imod=1; lim_low=8.2;lim_up=10.2;}
        //else if(c_mod==16){ imod=2; lim_low=9.0;lim_up=11.0;}
        //else if(c_mod==14){ imod=3; lim_low=9.0;lim_up=11.5;}
        //else{continue;}

        //h1[imod][c_id]->Fill(StartTime,c_gain);

        //if(c_gain<lim_low||c_gain>lim_up){
        //  bool registored=false;
        //  for(int ibad=0;ibad<num_v_bad;ibad++){
        //    if(
        //        v_bad_mod [ibad] == c_mod  &&
        //        v_bad_view[ibad] == c_view &&
        //        v_bad_pln [ibad] == c_pln  &&
        //        v_bad_ch  [ibad] == c_ch   )
        //    {
        //      registored = true;
        //      break;
        //    }
        //  }
        //  if(!registored){
        //    v_bad_mod .push_back(c_mod );
        //    v_bad_view.push_back(c_view);
        //    v_bad_pln .push_back(c_pln );
        //    v_bad_ch  .push_back(c_ch  );
        //    num_v_bad++;
        //  }
        //}
      }
    }
    wtree->Fill();
  }
  file->Close();
  delete badch;

  cout << "Writing histograms...." << endl;
  //for(int mod=0;mod<nummod;mod++){
  //  cout << "--- mod=" << mod << endl;
  //  for(int id=0;id<maxid;id++){
  //    if(id%100==0) cout << " " << id;
  //    int im = 0;
  //    if     (mod==0){im= 3;}
  //    else if(mod==1){im=15;}
  //    else if(mod==2){im=16;}
  //    else if(mod==3){im=14;}
  //    else           {continue;}
  //    for(int id=0;id<maxid;id++){
  //      int iv,ip,ic;
  //      if(getchannel_id(im,id,iv,ip,ic)){
  //        h1[mod][id]->Write();
  //      }
  //    }
  //  }
  //  cout << endl;
  //}
  wfile->cd();
  wtree->Write();
  wfile->Write();
  wfile->Close();

  //ofstream ofs("badchannel.txt");
  //ofs << " mod view pln ch" << endl;
  //for(int i=0;i<num_v_bad;i++){
  //  ofs
  //    << " " << v_bad_mod [i]
  //    << " " << v_bad_view[i]
  //    << " " << v_bad_pln [i]
  //    << " " << v_bad_ch  [i]
  //    << endl;
  //}
  //ofs.close();

  return;
}
