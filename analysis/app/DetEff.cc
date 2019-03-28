#include "DetEff.hxx"

void PrintUsage()
{
  cout << "-i Input ROOT file."    << endl;
  cout << "-o Output ROOT file."   << endl;
}


FileStat_t fs;
int main(int argc, char** argv){


  //*****arguments*******************//
  int c = -1;
  string readfilename, tunefilename;
  string outputfilename="";
  while ((c = getopt(argc, argv, "i:o:h")) != -1){
    switch(c){
    case 'i':
      readfilename = optarg;
      break;
    case 'o':
      outputfilename = optarg;
      break;
    case 'h':
      PrintUsage();
      exit(0);
    others:
      PrintUsage();
      exit(0);
    }
  }

  if(gSystem->GetPathInfo(readfilename.c_str(),fs)){
    cout << "No such a file (Input ROOT File) : " << readfilename << endl;
    PrintUsage();
    return 0;
  }
  if(outputfilename==""){
    cout << "Output filename is not set." << endl;
    PrintUsage();
    return 0;
  }

  //------- Open root file and get tree -------
  gROOT->SetStyle("Plain");
  TFile* readfile = new TFile(readfilename.c_str(),"read");
  if(readfile->IsZombie()){
    cout << "Cannot open file : " << readfilename << endl;
    exit(1);
  }
  cout << "Open file : " << readfilename << endl;

  EventSummary *evtsum = new EventSummary();
  TTree        *tree   = (TTree*)readfile ->Get("tree");
  TBranch      *evtbr  = tree->GetBranch("fDefaultReco.");
  evtbr -> SetAddress(&evtsum);
  tree  -> SetBranchAddress("fDefaultReco.",&evtsum);
  int nevt = tree->GetEntries();

  //------- Open output root file -------
  TFile* outfile = new TFile(outputfilename.c_str(),"recreate");

  // ---- Set up -----
  badch =  new INGRID_BadCh_mapping();
  badch -> set_BadCh(true);
  detdim  = new DetectorDimension();

  // ------ Set histograms ------
  TH1F *h[NumHist];
  for(int i=0; i<NumHist; i++){
    int nbin;
    double min,max;
    string title,xtitle,ytitle,name;
    SetHist(i,nbin,min,max,title,xtitle,ytitle,name);
    h[i] = new TH1F(name.c_str(),title.c_str(),nbin,min,max);
    h[i]->GetXaxis()->SetTitle(xtitle.c_str());
    h[i]->GetYaxis()->SetTitle(ytitle.c_str());
  }

  TH2F *h2[NumHist2];
  for(int i=0; i<NumHist2; i++){
    int    nbin1,nbin2;
    double min1, max1, min2, max2;
    string title,xtitle,ytitle,name;
    SetHist2(i,nbin1,min1,max1,nbin2,min2,max2,title,xtitle,ytitle,name);
    h2[i] = new TH2F(name.c_str(),title.c_str(),nbin1,min1,max1,nbin2,min2,max2);
    h2[i]->GetXaxis()->SetTitle(xtitle.c_str());
    h2[i]->GetYaxis()->SetTitle(ytitle.c_str());
  }


  const float FVCut[3][4] ={
    {FVcutX0,FVcutY0,zposi(14,1,VetoCutPln0,0)+5.0,zposi(14,1,FVcutPlnZ0,0)-5.0},
    {FVcutX1,FVcutY1,zposi(16,1,VetoCutPln1,0)+6.5,zposi(16,1,FVcutPlnZ1,0)-6.5},
    {FVcutX2,FVcutY2,zposi(21,1,VetoCutPln2,0)+1.5,zposi(21,0,FVcutPlnZ2,0)-1.5}
  };

  cout<< "------- FV Thickness ---------" << endl;
  cout<<"ING : "<<FVCut[0][3]-FVCut[0][2]<<" ("<<FVCut[0][2]<<","<<FVCut[0][3]<<")"<<endl;
  cout<<"PM  : "<<FVCut[1][3]-FVCut[1][2]<<" ("<<FVCut[1][2]<<","<<FVCut[1][3]<<")"<<endl;
  cout<<"WM  : "<<FVCut[2][3]-FVCut[2][2]<<" ("<<FVCut[2][2]<<","<<FVCut[2][3]<<")"<<endl;
  //cout<<"ING : "<<FVCut[0][3]-FVCut[0][2]<<" ("<<FVCut[0][2]+CoT[0][2]<<","<<FVCut[0][3]+CoT[0][2]<<")"<<endl;
  //cout<<"PM  : "<<FVCut[1][3]-FVCut[1][2]<<" ("<<FVCut[1][2]+CoT[1][2]<<","<<FVCut[1][3]+CoT[1][2]<<")"<<endl;
  //cout<<"WM  : "<<FVCut[2][3]-FVCut[2][2]<<" ("<<FVCut[2][2]+CoT[2][2]<<","<<FVCut[2][3]+CoT[2][2]<<")"<<endl;
  cout<< "------- Normalization Thickness ---------" << endl;
  cout<<"ING : "<<thickness[0]<<endl;
  cout<<"PM  : "<<thickness[1]<<endl;
  cout<<"WM  : "<<thickness[2]<<endl;

  //==============================================
  //==========     Event loop start     ==========
  //==============================================
  std::cout << "Total # of events = " << nevt << std::endl;
  for(int ievt = 0; ievt<nevt; ievt++){
    tree->GetEntry(ievt);
    if(ievt%1000==0){
      std::cout << ">> event " << ievt << std::endl;
    }
    
    // ----------------
    // SimVertex
    // ----------------
    double vertex[3]={0.,0.,0.};
    int nsimver = evtsum->NSimVertexes();
    int sim_mod = -1;
    if(nsimver>0){
      SimVertexSummary* simver = evtsum->GetSimVertex(0);
      vertex[0] = simver->xnu*10.;
      vertex[1] = simver->ynu*10.;
      vertex[2] = simver->znu*10.;
      for(int imod=0;imod<3;imod++){
        if(vertex[2]>=FVCut[imod][2]+CoT[imod][2]&&vertex[2]<=FVCut[imod][3]+CoT[imod][2]){
          sim_mod=imod;
          break;
        }
      }
    }
    if(sim_mod==-1){
      cout << "Out of fiducal volume" << endl;
      continue;
    }
    // ------------------
    // SimParticle
    // ------------------
    double simpar_ang=-1.;
    double simpar_mom=-1.;
    double simpar_pdg=-1.;
    int nsimpar = evtsum->NSimParticles();
    for(int isimpar=0;isimpar<nsimpar;isimpar++){
      if(isimpar!=0){continue;}
      SimParticleSummary* simpar = evtsum->GetSimParticle(isimpar);
      int   pdg        = simpar->pdg;
      float momentum   = simpar->momentum[3]; //GeV/c
      float momX       = simpar->momentum[0];
      float momY       = simpar->momentum[1];
      float momZ       = simpar->momentum[2];

      simpar_pdg = abs(pdg);
      //simpar_mom = momentum;
      simpar_mom = sqrt(momX*momX+momY*momY+momZ*momZ);
      if(momZ==0.        ){ simpar_ang = 90.; }
      else                { simpar_ang = atan( sqrt(momX*momX+momY*momY)/momZ )*180./PI;}
      if(simpar_ang < 0. ){ simpar_ang = 180. + simpar_ang; }
      if(simpar_ang==180.){ simpar_ang = 0.;}
    }

    h [0+sim_mod]->Fill(simpar_mom);
    h [3+sim_mod]->Fill(simpar_ang);
    h2[0+sim_mod]->Fill(simpar_mom,simpar_ang);

    // -------------------
    // SimHit
    // ------------------
    int nsimhit = evtsum->NSimHits();
    for(int isimhit=0;isimhit<nsimhit;isimhit++){
      SimHitSummary* simhit = evtsum->GetSimHit(isimhit);
      HitSummary   * hit    = evtsum->GetHit(isimhit);
      int simhitmod = hit->mod;
      int cmod = -1;
      if     (simhitmod==14){cmod=0;}
      else if(simhitmod==16){cmod=1;}
      else if(simhitmod==21){cmod=2;}
      else{continue;}
    }

    // -----------
    // Vertices
    // -----------
    bool first_vtx=true;
    bool simpar_found=false;
    int nvtx = evtsum->NThreeDimRecons();
    for(int ivtx=0;ivtx<nvtx;ivtx++){

      ThreeDimReconSummary *vtx = evtsum->GetThreeDimRecon(ivtx);
      int ntrk = vtx->Ntrack;
      for(int itrk=0;itrk<ntrk;itrk++){
        int   startmod    = vtx->startmod   [itrk];
        int   stopmodx    = vtx->stopmodx   [itrk]; //side view
        int   stopmody    = vtx->stopmody   [itrk]; //top view
        float startxch    = vtx->startxch   [itrk]; //y (side view)
        float startych    = vtx->startych   [itrk]; //x (top view)
        int   startxpln   = vtx->startxpln  [itrk];
        int   startypln   = vtx->startypln  [itrk];
        int   endxpln     = vtx->endxpln    [itrk];
        int   endypln     = vtx->endypln    [itrk];
        float endxch      = vtx->endxch     [itrk];
        float endych      = vtx->endych     [itrk];
        float diff_posx   = vtx->diff_posx  [itrk];
        float diff_posy   = vtx->diff_posy  [itrk];
        float diff_angx   = vtx->diff_angx  [itrk];
        float diff_angy   = vtx->diff_angy  [itrk];
        float diff_timex  = vtx->diff_timex [itrk];
        float diff_timey  = vtx->diff_timey [itrk];
        float diff_posx2  = vtx->diff_posx2 [itrk];
        float diff_posy2  = vtx->diff_posy2 [itrk];
        float diff_angx2  = vtx->diff_angx2 [itrk];
        float diff_angy2  = vtx->diff_angy2 [itrk];
        float diff_timex2 = vtx->diff_timex2[itrk];
        float diff_timey2 = vtx->diff_timey2[itrk];
        float diff_posx3  = vtx->diff_posx3 [itrk];
        float diff_posy3  = vtx->diff_posy3 [itrk];
        float diff_angx3  = vtx->diff_angx3 [itrk];
        float diff_angy3  = vtx->diff_angy3 [itrk];
        float diff_timex3 = vtx->diff_timex3[itrk];
        float diff_timey3 = vtx->diff_timey3[itrk];
        float oneview     = vtx->oneview    [itrk];
        float thetax      = vtx->thetax     [itrk]*PI/180.; //side
        float thetay      = vtx->thetay     [itrk]*PI/180.; //top

        double angle = fabs(atan(sqrt(pow(tan(thetax),2)+pow(tan(thetay),2))))*180./PI;

        int anamod = -1;
        if     (startmod==14){anamod=0;}
        else if(startmod==16){anamod=1;}
        else if(startmod==21){anamod=2;}
        else{
          cout << "Vertex out of the target module" << endl;
          continue;
        }

        if(anamod!=sim_mod){continue;}

        if(first_vtx){
          h[6+sim_mod]->Fill(simpar_mom);
          h[9+sim_mod]->Fill(simpar_ang);
          h2[3+sim_mod]->Fill(simpar_mom,simpar_ang);
          first_vtx=false;
        }


        int nanahit = vtx->NhitTs(itrk);
        // ===================
        // Particle Info
        int numhit_on_trk = 0;
        for(int ianahit=0;ianahit<nanahit;ianahit++){
          HitSummary* anahit = vtx->GetHitTrk(ianahit,itrk);
          if(anahit->NSimHits()<1){continue;}
          SimHitSummary* anasimhit = anahit->GetSimHit(0);
          int anapdg   = anasimhit->pdg;
          int anatrkid = anasimhit->trackid;
          anapdg = abs(anapdg);
          if(anatrkid==1&&anapdg==simpar_pdg){
            numhit_on_trk++;
          }
        }
        if(numhit_on_trk>2&&!simpar_found){
          h[12+sim_mod]->Fill(simpar_mom);
          h[15+sim_mod]->Fill(simpar_ang);
          h2[6+sim_mod]->Fill(simpar_mom,simpar_ang);
          simpar_found=true;
        }

      } //for(itrk)
    }//for(ivtx)

  }//for(ievt)

  if(evtsum) delete evtsum;

  outfile->cd();
  for(int i=0; i<NumHist; i++){
    if(h[i]!=NULL) h[i]->Write();
  }
  for(int i=0; i<NumHist2; i++){
    if(h2[i]!=NULL) h2[i]->Write();
  }
  outfile->Close();

  std::cout << "Done!!!" << std::endl;
  return 0;
}
