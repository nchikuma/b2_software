#include "B2Plot.hxx"

//#define DEBUG_B2PLOT

void PrintUsage()
{
  cout << "-i Input ROOT file."    << endl;
  cout << "-o Output ROOT file."   << endl;
  cout << "-t Flux tuning file."   << endl;
  cout << "-m <fdid>           "   << endl;
  cout << "  fdid : 0-> data, 7-> MC (WAGASCI), 8-> MC (PM), 9-> MC (ING)" << endl;
  cout << "-c <target>       "     << endl;
  cout << "  target : 0 -> CCQE, 1 -> MEC, 2-> CC1pi, 3-> CCcoh, 4-> CCDIF, 5->CC other, 6->NC" <<endl;
  cout << "           7 -> CC0pi0p, 8 -> CC0piNp, 9-> CC1pi, 10-> CCcother" <<endl;
  cout << "           11 -> CC0pi" <<endl;
  cout << "-n <dial> : For detector systematics" << endl;
  cout << "       27 -> Beam timing cut : OFF"                  << endl;
  cout << "       28 -> Upst Veto (WM)  : VetoCutPln2 = 7"      << endl;
  cout << "       29 -> Upst Veto (WM)  : VetoCutPln2 = 10"     << endl;
  cout << "       30 -> Upst Veto (PM)  : VetoCutPln1 =  2"     << endl;
  cout << "       31 -> Upst Veto (PM)  : VetoCutPln1 =  3"     << endl;
  cout << "       32 -> FV cut (WM)     : FVcutX2=FVcutY2=250." << endl;
  cout << "       33 -> FV cut (PM)     : FVcutX1=FVcutY1=250." << endl;
  cout << "       34 -> Acceptance cut  : ACCEPT_OFFSETZ  = +50." << endl;
  cout << "       35 -> Acceptance cut  : ACCEPT_OFFSETZ  = -50." << endl;
  cout << "       36 -> Acceptance cut  : ACCEPT_OFFSETY  = +50." << endl;
  cout << "       37 -> Acceptance cut  : ACCEPT_OFFSETY  = -50." << endl;
  cout << "       38 -> Acceptance cut  : ACCEPT_OFFSETX  = +50." << endl;
  cout << "       39 -> Acceptance cut  : ACCEPT_OFFSETX  = -50." << endl;
  cout << "-r <t2kreweight file>" <<endl;
  cout << "-f <flux weight ID>" << endl;
  cout << "-w <weight ID>" <<endl;
  cout << "-a : For reweight by an arbitrary function." << endl;
  cout << "-d : For unifying phase space between WM and PM" << endl;
  cout << "-p : For just extracting number of selected events." << endl;
}


FileStat_t fs;
int main(int argc, char** argv){


  // =================================================
  //    Set up the configuration
  //
  int c = -1;
  string readfilename, tunefilename;
  string outputfilename   = "";
  string t2krewfilename   = "";
  int fdid                = -1;
  int sim_mod             = -1;
  int target              = -1;
  int weight_id           = -1;
  int detsys_id           = -1;
  int fluxerr_id          = -1;
  bool t2kreweight        = false;
  bool xsec_ratio         = false;
  bool NselOnly           = false;
  bool ARBITRARY_REWEIGHT = false;
  bool FluxErr            = false;
  while ((c = getopt(argc, argv, "i:o:t:m:c:r:w:n:f:daph")) != -1){
    switch(c){
    case 'i':
      readfilename = optarg;
      break;
    case 'o':
      outputfilename = optarg;
      break;
    case 't':
      tunefilename = optarg;
      break;
    case 'm':
      fdid = atoi(optarg);
      if(fdid==7){sim_mod=2;}
      if(fdid==8){sim_mod=1;}
      if(fdid==9){sim_mod=0;}
      break;
    case 'c':
      target = atoi(optarg);
      break;
    case 'r':
      t2krewfilename = optarg;
      t2kreweight = true;
      break;
    case 'f':
      fluxerr_id = atoi(optarg);
      FluxErr = true;
      break;
    case 'w':
      weight_id = atoi(optarg);
      break;
    case 'n':
      detsys_id = atoi(optarg);
      break;
    case 'd':
      xsec_ratio = true;
      break;
    case 'p':
      NselOnly = true;
      break;
    case 'a':
      ARBITRARY_REWEIGHT = true;
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
  if((fdid==7||fdid==8||fdid==9)&&gSystem->GetPathInfo(tunefilename.c_str(),fs)){
    cout << "No such a file (JNUBEAM Tuning File) : " << tunefilename << endl;
    PrintUsage();
    return 0;
  }
  if(outputfilename==""){
    cout << "Output filename is not set." << endl;
    PrintUsage();
    return 0;
  }
  if(fdid<0){
    cout << "Flux id <fdid> must be set." << endl;
    PrintUsage();
    return 0;
  }
  if(t2kreweight&&gSystem->GetPathInfo(t2krewfilename.c_str(),fs)&&t2krewfilename!="local"){
    cout << "No such a file (T2K Reweight): " << t2krewfilename << endl;
    PrintUsage();
    return 0;
  }
  if(t2kreweight&&(weight_id<0&&weight_id>=num_reweight)){
    cout << "Set an ID for t2k reweight. 0 to " << num_reweight-1 << endl;
    PrintUsage();
    return 0;
  }
  if(xsec_ratio){
    if(PION_MOM_TH  [2]>PION_MOM_TH  [1]){ PION_MOM_TH  [1] = PION_MOM_TH  [2]; }
    else                                 { PION_MOM_TH  [2] = PION_MOM_TH  [1]; }
    if(PROTON_MOM_TH[2]>PROTON_MOM_TH[1]){ PROTON_MOM_TH[1] = PROTON_MOM_TH[2]; }
    else                                 { PROTON_MOM_TH[2] = PROTON_MOM_TH[1]; }
    if(PION_ANG_TH  [2]<PION_ANG_TH  [1]){ PION_ANG_TH  [1] = PION_ANG_TH  [2]; }
    else                                 { PION_ANG_TH  [2] = PION_ANG_TH  [1]; }
    if(PROTON_ANG_TH[2]<PROTON_ANG_TH[1]){ PROTON_ANG_TH[1] = PROTON_ANG_TH[2]; }
    else                                 { PROTON_ANG_TH[2] = PROTON_ANG_TH[1]; }
  }


  //------- Open root file and get tree -------
  gROOT->SetStyle("Plain");
  char tunefile[300];
  sprintf(tunefile,"%s",tunefilename.c_str());
  FluxTuning* fluxtune;
  if(fdid==7||fdid==8||fdid==9) fluxtune = new FluxTuning(fdid,tunefile);

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

  //------- T2K Reweight File -------
  TFile* t2krewfile;
  TTree* t2krewtree;
  TArrayF *t2krew = new TArrayF(num_reweight);
  if(t2kreweight){
    if(t2krewfilename=="local"){
      t2krewtree = (TTree*) readfile->Get("weightstree");
    }
    else{
      t2krewfile = new TFile(t2krewfilename.c_str(),"read");
      t2krewtree = (TTree*) t2krewfile->Get("weightstree");
    }
    if(t2krewtree->GetEntries()!=nevt){
      cout << "This t2kreweight file does not have the same number of entries" << endl;
      exit(1);
    }
    t2krewtree->SetBranchAddress("weights",&t2krew);
  }

  // ----------- Flux Error File --------------
  double *fluxwei;
  if(FluxErr){
    string fluxname = "/home/t2k/nchikuma/b2_data2/flux_err/tree/flux_err_weight.root";
    TFile *fluxfile = new TFile(fluxname.c_str(),"read");
    string fwei_name = Form("weight_%05d",fluxerr_id);
    TH1F *hfluxerr = (TH1F*) fluxfile->Get(fwei_name.c_str());
    int nbins = hfluxerr->GetNbinsX();
    fluxwei = (double*)malloc(sizeof(double)*nbins);
    for(int ibin=0;ibin<nbins;ibin++){
      fluxwei[ibin] = hfluxerr->GetBinContent(ibin+1);
    }
  }
  
  // ---- Set up -----
  badch =  new INGRID_BadCh_mapping();
  badch -> set_BadCh(true);
  detdim  = new DetectorDimension();


  bool beam_timing_mask = false;
  if(detsys_id>0){
    if     (detsys_id==27){beam_timing_mask = true  ;}
    else if(detsys_id==28){VetoCutPln2 = 7          ;}
    else if(detsys_id==29){VetoCutPln2 = 10         ;}
    else if(detsys_id==30){VetoCutPln1 =  2         ;}
    else if(detsys_id==31){VetoCutPln1 =  3         ;}
    else if(detsys_id==32){FVcutX2=250.;FVcutY2=250.;}
    else if(detsys_id==33){FVcutX1=250.;FVcutY1=250.;}
    else if(detsys_id==34){ACCEPT_OFFSETZ  =  50.   ;}
    else if(detsys_id==35){ACCEPT_OFFSETZ  = -50.   ;}
    else if(detsys_id==36){ACCEPT_OFFSETY  =  50.   ;}
    else if(detsys_id==37){ACCEPT_OFFSETY  = -50.   ;}
    else if(detsys_id==38){ACCEPT_OFFSETX  =  50.   ;}
    else if(detsys_id==39){ACCEPT_OFFSETX  = -50.   ;}
    else{
      PrintUsage();
      cout << "Unexpected dial input for detector systematics." << endl;
      return 0;
    }
  }



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

  TH1F *h3[NumHist3];
  for(int i=0;i<NumHist3;i++){
    int nbin;
    double min,max;
    string title,xtitle,ytitle,name;
    SetHist3(i,nbin,min,max,title,xtitle,ytitle,name);
    h3[i] = new TH1F(name.c_str(),title.c_str(),nbin,min,max);
    h3[i]->GetXaxis()->SetTitle(xtitle.c_str());
    h3[i]->GetYaxis()->SetTitle(ytitle.c_str());
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
  cout<< "------- Normalization Thickness ---------" << endl;
  cout<<"ING : "<<thickness[0]<<endl;
  cout<<"PM  : "<<thickness[1]<<endl;
  cout<<"WM  : "<<thickness[2]<<endl;
  cout<< "------- Phase Space for PION ---------" << endl;
  cout<<"ING : "<<PION_MOM_TH[0]<<", "<<PION_ANG_TH[0]<<endl;
  cout<<"PM  : "<<PION_MOM_TH[1]<<", "<<PION_ANG_TH[1]<<endl;
  cout<<"WM  : "<<PION_MOM_TH[2]<<", "<<PION_ANG_TH[2]<<endl;
  cout<< "------- Phase Space for PROTON ---------" << endl;
  cout<<"ING : "<<PROTON_MOM_TH[0]<<", "<<PROTON_ANG_TH[0]<<endl;
  cout<<"PM  : "<<PROTON_MOM_TH[1]<<", "<<PROTON_ANG_TH[1]<<endl;
  cout<<"WM  : "<<PROTON_MOM_TH[2]<<", "<<PROTON_ANG_TH[2]<<endl;


  //==============================================
  //==========     Event loop start     ==========
  //==============================================
  std::cout << "Total # of events = " << nevt << std::endl;
  int nutype = -1;
  for(int ievt = 0; ievt<nevt; ievt++){
    tree->GetEntry(ievt);
    if(ievt%1000==0){
      std::cout << ">> event " << ievt << std::endl;
    }

    double t2krewei = 1.;
    if(t2kreweight){ 
      t2krewtree->GetEntry(ievt);
      t2krewei = t2krew->At(weight_id);
    }

    // --------------------------------
    // Simulation : Vertex Information
    // --------------------------------
    double reweight=1.,flux_tot=1.,norm_tot=1.;
    double norm=1.,totcrsne=1.;
    double nuE=-1.,vertex[3]={0.,0.,0.};
    int    inttype=-1,targetz=-1;
    int    inttype_id=-1;
    bool   sim_inFV     = false;
    bool   sim_inFVarea = false;
    int nsimver = evtsum->NSimVertexes();
    double arb_rewei = 1.;
    if(nsimver>0){
      SimVertexSummary* simver = evtsum->GetSimVertex(0);
      norm      = simver->norm;
      totcrsne  = simver->totcrsne;
      nuE       = simver->nuE;
      inttype   = simver->inttype;
      nutype    = simver->nutype;
      targetz   = simver->targetz;
      vertex[0] = simver->xnu*10.;
      vertex[1] = simver->ynu*10.;
      vertex[2] = simver->znu*10.;

      inttype = abs(inttype);

      if(target==0){ //CCQE
        if(inttype!=1){continue;}
      }else if(target==1){ //MEC
        if(inttype!=2){continue;}
      }else if(target==2){ //CC1pi
        if(inttype!=11&&inttype!=12&&inttype!=13){continue;}
      }else if(target==3){ //CCcoh
        if(inttype!=16){continue;}
      }else if(target==4){ //CCdis
        if(inttype!=21&&inttype!=26){continue;}
      }else if(target==5){ //CC other
        if(inttype!=17&&inttype!=22&&inttype!=23){continue;}
      }else if(target==6){ //NC
        if(inttype!=31&&inttype!=32&&inttype!=33&&inttype!=34
            &&inttype!=36&&inttype!=38&&inttype!=39&&inttype!=41
            &&inttype!=42&&inttype!=43&&inttype!=44&&inttype!=45
            &&inttype!=46&&inttype!=51&&inttype!=52
          )
        {continue;}
      }

      if     (inttype==1){inttype_id=0;}
      else if(inttype==2){inttype_id=1;}
      else if(inttype==11||inttype==12||inttype==13){inttype_id=2;}
      else if(inttype==16){inttype_id=3;}
      else if(inttype==21||inttype==26){inttype_id=4;}
      else if(inttype==17||inttype==22||inttype==23){inttype_id=5;}
      else if(inttype==31||inttype==32||inttype==33||inttype==34
            ||inttype==36||inttype==38||inttype==39||inttype==41
            ||inttype==42||inttype==43||inttype==44||inttype==45
            ||inttype==46||inttype==51||inttype==52
          ){inttype_id=6;}

      bool tmp=true;
      if(sim_mod>=0){
        for(int i=0;i<3;i++){
          double pos_inmod = vertex[i]-CoT[sim_mod][i];
          double fvcut0    = FVCut[sim_mod][i];
          double fvcut1;
          if(i==2) fvcut1 = FVCut[sim_mod][i+1];
          if(
              (i< 2&&fabs(pos_inmod)>fvcut0)
              ||(i==2&&(pos_inmod<fvcut0||pos_inmod>fvcut1))
            )
          {
            tmp=false;break;
          }
        }
        sim_inFV = tmp;

        tmp=true;
        for(int i=0;i<2;i++){
          double pos_inmod = vertex[i]-CoT[sim_mod][i];
          double fvcut0    = FVCut[sim_mod][i];
          if(fabs(pos_inmod)>fvcut0) { tmp=false;break; }
        }
        sim_inFVarea = tmp;

      }

      if(ARBITRARY_REWEIGHT){
        if     (nuE<0.5){ arb_rewei = 1.+nuE; }
        else if(nuE<1.5){ arb_rewei = 2.-nuE; }
        else if(nuE<2.0){ arb_rewei = nuE-1.; }
        else            { arb_rewei = 1.;     } 
      }

      double fluxerr = 1.;
      if(FluxErr){
        int offset;
        if     (nutype==2){offset=20;}
        else if(nutype==1){offset= 0;}
        else{
          cout << "Unexpected neutrino type: " << nutype << endl;
          continue;
        }
        int id;
        if     (nuE< 3.){ id = (int)(nuE/0.2)      + offset; fluxerr = 1.+fluxwei[id];}
        else if(nuE< 4.){ id = (int)(nuE/1.0) + 15 + offset; fluxerr = 1.+fluxwei[id];}
        else if(nuE<10.){ id = (int)(nuE/2.0) + 16 + offset; fluxerr = 1.+fluxwei[id];}
        else if(nuE<30.){ id = (int)(nuE/20.) + 19 + offset; fluxerr = 1.+fluxwei[id];}
      }

      if(fdid==7||fdid==8||fdid==9){
        reweight = fluxtune->reweight(nuE,nutype);
        reweight *= fluxerr;
      }
      if(sim_mod>=0){
        norm_tot = norm*totcrsne*1e-38*Avogadro*density[sim_mod]*thickness[sim_mod]*reweight*t2krewei*arb_rewei;
        flux_tot = norm/area[sim_mod]*reweight;
      }else{
        float thick = 0.;
        if     (fdid==10){thick=thickness_wallbg;}
        else if(fdid==11){thick=thickness_ceiling;}
        else if(fdid==12){thick=thickness_floor;}
        else if(fdid==13){thick=thickness_pillar_r;}
        else if(fdid==14){thick=thickness_pillar_l;}
        norm_tot = norm*totcrsne*1e-38*Avogadro*density_wallbg*thick;
        flux_tot = norm;
      }
    }

    // ---------------------------------------------
    // Simulation : Injected Particles Information
    // ---------------------------------------------
    double mu_ang=-1., mu_mom=-1.;
    double p_ang=-1., p_mom=-1.; //with the largest momentum
    double pi_ang=-1., pi_mom=-1.; //with the largetst momentum
    bool   sim_muon      = false;
    bool   sim_proton    = false;
    bool   sim_pion      = false;
    bool   sim_mu_400MeV = false;
    bool   sim_mu_45deg  = false;
    bool   sim_mu_30deg  = false;
    bool   sim_hitING  = false;
    bool   sim_stopING = false;
    int    inttype2 = -1; //0: CC0pi0p, 1: CC0piNp, 2: CC1pi, 3: CCother
    int    num_sim_mu = 0;
    int    num_sim_pi = 0;
    int    num_sim_p  = 0;
    int nsimpar = evtsum->NSimParticles();
    for(int isimpar=0;isimpar<nsimpar;isimpar++){
      SimParticleSummary* simpar = evtsum->GetSimParticle(isimpar);
      int   pdg      = abs(simpar->pdg);
      float kineticE = simpar->momentum[3]; //GeV : Kinetic Energy
      float momX     = simpar->momentum[0];
      float momY     = simpar->momentum[1];
      float momZ     = simpar->momentum[2];
      float fposx    = simpar->fpos[0]*10.; //mm
      float fposy    = simpar->fpos[1]*10.; //mm
      float fposz    = simpar->fpos[2]*10.; //mm

      float momentum = sqrt(momX*momX+momY*momY+momZ*momZ);
      float sim_ang;
      if(momZ==0.     ){ sim_ang = 90.; }
      else             { sim_ang = atan( sqrt(momX*momX+momY*momY)/momZ )*180./PI ;}
      if(sim_ang < 0. ){ sim_ang = 180. + sim_ang; }
      if(sim_ang==180.){ sim_ang = 0.;}
      double ang_fold;
      if     (sim_ang>= 0.&&sim_ang< 90.) ang_fold = sim_ang;
      else if(sim_ang>=90.&&sim_ang<180.) ang_fold = 180.-sim_ang;
      else ang_fold = 180.;

      fposx -= C_B2MotherPosX + C_B2INGPosX;
      fposy -= C_B2MotherPosY + C_B2INGPosY;
      fposz -= C_B2MotherPosZ + C_B2INGPosZ;
      if(pdg==13){
        sim_muon = true;
        num_sim_mu++;
        mu_mom = momentum;
        mu_ang = sim_ang;
        if(fabs(fposx)<INGStopXY && fabs(fposy)<INGStopXY 
            && fposz>zposi(14,0,2,0) && fposz<zposi(14,0,INGStopPlnZ,0))
        {
          sim_stopING = true;
        }
        if(mu_mom>=0.4){sim_mu_400MeV=true;}
        if(mu_ang< 45.){sim_mu_45deg =true;}
        if(mu_ang< 30.){sim_mu_30deg =true;}
      }
      else if(pdg==211){ //charged pion
        if(momentum>=PION_MOM_TH[sim_mod]&&fabs(ang_fold)<PION_ANG_TH[sim_mod]){
          num_sim_pi++;
          sim_pion=true;
          if(momentum>pi_mom){
            pi_mom = momentum;
            pi_ang = sim_ang;
          }
        }
      }
      else if(pdg==2212){ //proton
        if(momentum>=PROTON_MOM_TH[sim_mod]&&fabs(ang_fold)<PROTON_ANG_TH[sim_mod]){
          num_sim_p ++;
          sim_proton=true;
          if(momentum>p_mom){
            p_mom = momentum;
            p_ang = sim_ang;
          }
        }
      }
    }
    if(num_sim_mu>0){
      if(num_sim_pi==0){
        if(num_sim_p==0){ inttype2=0; }
        else            { inttype2=1; }
      }
      else if(num_sim_pi==1){ inttype2=2; }
      else{ inttype2=3; }
    }

    if     (target== 7&&inttype2!=0){ continue; } //CC0pi0p
    else if(target== 8&&inttype2!=1){ continue; } //CC0piNp
    else if(target== 9&&inttype2!=2){ continue; } //CC1pi
    else if(target==10&&inttype2!=3){ continue; } //CCother
    else if(target==11&&inttype2!=0&&inttype2!=1){ continue; } //CC0pi


    // ========= For 8-bin Analysis ================
    bool truebin[8];
    bool detect_p_pi = (num_sim_pi>0) || (num_sim_p>0);
    bool detect_mu   =  sim_muon && mu_mom >= 0.4;
    bool cc_events   =  sim_muon &&  sim_inFV;
    bool sim_ooFV    =  sim_muon && !sim_inFV;
    bool nc_events   = !sim_muon;
    bool int_onCH    = (targetz==6);
    truebin[0] = cc_events && (!detect_mu || detect_p_pi);
    truebin[1] = cc_events && detect_mu && !detect_p_pi && (mu_ang>= 0. && mu_ang<  5.);
    truebin[2] = cc_events && detect_mu && !detect_p_pi && (mu_ang>= 5. && mu_ang< 10.);
    truebin[3] = cc_events && detect_mu && !detect_p_pi && (mu_ang>=10. && mu_ang< 15.);
    truebin[4] = cc_events && detect_mu && !detect_p_pi && (mu_ang>=15. && mu_ang< 20.);
    truebin[5] = cc_events && detect_mu && !detect_p_pi && (mu_ang>=20. && mu_ang< 25.);
    truebin[6] = cc_events && detect_mu && !detect_p_pi && (mu_ang>=25. && mu_ang< 30.);
    truebin[7] = cc_events && detect_mu && !detect_p_pi && (mu_ang>=30. && mu_ang<180.);
#ifdef DEBUG_B2PLOT
    // --- Checking the binning
    bool check_truebin = false;
    for(int ibin=0;ibin<8;ibin++){
      if(truebin[ibin]){check_truebin=true;break;}
    }
    if(!check_truebin&&cc_events){
      cout << "Found an event out of the true binning." << endl;
      exit(1);
    }
#endif

    for(int ibin=0;ibin<8;ibin++){
      if(truebin[ibin]){
        h[386+sim_mod]->Fill(ibin,norm_tot);
        if(int_onCH){ h[413+sim_mod]->Fill(ibin,norm_tot); }
        break;
      }
    }



    // ---------------------------------------------
    // Simulation : True Hit Information
    // ---------------------------------------------
    int nsimhit = evtsum->NSimHits();
    int nmuonhit[2][11];
    for(int ii=0;ii<2;ii++){
      for(int jj=0;jj<11;jj++){
        nmuonhit[ii][jj]=0;
      }
    }
    for(int isimhit=0;isimhit<nsimhit;isimhit++){
      SimHitSummary* simhit = evtsum->GetSimHit(isimhit);
      HitSummary   * hit    = evtsum->GetHit(isimhit);
      int simhitmod = hit->mod;
      int cmod = -1;
      if     (simhitmod==14){cmod=0;}
      else if(simhitmod==16){cmod=1;}
      else if(simhitmod==21){cmod=2;}
      else{continue;}

      if(abs(simhit->pdg)==13&&hit->mod==14&&hit->pln<11){
        nmuonhit[hit->view][hit->pln] = 1;
      }
      if(abs(simhit->pdg)==13){
        h[326+cmod]->Fill(hit->pe,norm_tot);
      }
      else{
        h[329+cmod]->Fill(hit->pe,norm_tot);
      }
    }
    int nmuoncnt = 0;
    for(int ii=0;ii<11;ii++){
      if(nmuonhit[0][ii]==1&&nmuonhit[1][ii]==1){nmuoncnt++;}
    }
    if(nmuoncnt>2){sim_hitING=true;}


    // ====================
    // Nu Flux Distribution
    if(sim_inFVarea){
      h [74+sim_mod]->Fill(nuE,flux_tot);
      if(t2krewfilename!="local"){
        h3[35+sim_mod]->Fill(0.,flux_tot);
      }
    }
    if(sim_inFV){
      h[302+sim_mod]->Fill(nuE,flux_tot);
      if(sim_mu_30deg&&sim_mu_400MeV){
        h[305+sim_mod]->Fill(nuE,flux_tot);
      }
    }

    // ======================
    // Total Cross Section
    h2[67+sim_mod]->Fill(nuE,totcrsne,norm*reweight);
    if(sim_inFV){
      if(sim_muon){
        h[296+sim_mod]->Fill(totcrsne,norm*reweight);
        if(sim_mu_30deg&&sim_mu_400MeV){
          h[299+sim_mod]->Fill(totcrsne,norm*reweight);
        }
      }
    }

    // ======================
    // All interaction in FV 
    if(sim_inFV&&sim_muon){
      h[ 77+sim_mod]->Fill(nuE,norm_tot);
      h[332+sim_mod]->Fill(mu_ang,norm_tot);
      if(
            (sim_mod==0&&targetz==26) //Fe
          ||(sim_mod==1&&targetz== 6) //C from CH
          ||(sim_mod==2&&targetz== 1) //H from H2O
          )        
      {
        h[314+sim_mod]->Fill(mu_ang,norm_tot);
      }
    }

    // ===========================
    // Interactions in Phase Space
    if(sim_inFV&&sim_muon&&sim_mu_400MeV&&sim_mu_30deg){
      h[287+sim_mod]->Fill(mu_ang,norm_tot);
      h3[38+sim_mod]->Fill(0.,norm_tot);
      if(
            (sim_mod==0&&targetz==26) //Fe
          ||(sim_mod==1&&targetz== 6) //C from CH
          ||(sim_mod==2&&targetz== 1) //H from H2O
          )        
      {
        h[308+sim_mod]->Fill(mu_ang,norm_tot);
      }
      else if(targetz!=-1){
        h[311+sim_mod]->Fill(mu_ang,norm_tot);
      }

      h[365+sim_mod]->Fill(inttype_id,norm_tot);
      h2[82+sim_mod]->Fill(inttype_id,inttype2,norm_tot);

      h2[85+sim_mod]->Fill(pi_mom,pi_ang,norm_tot);
      h[368+sim_mod]->Fill(pi_mom,norm_tot);
      if(!sim_pion){
        h2[88+sim_mod]->Fill(p_mom,p_ang,norm_tot);
        h[371+sim_mod]->Fill(p_mom,norm_tot);
      }
      h[374+sim_mod]->Fill(p_mom,norm_tot);
    }

    // =====================
    // CCinc interaction in FV
    if(sim_inFV&&sim_muon){
      h[80+sim_mod]->Fill(nuE   ,norm_tot);
      h[83+sim_mod]->Fill(mu_mom,norm_tot);
      h[86+sim_mod]->Fill(mu_ang,norm_tot);

      h2[6+sim_mod*NumSelec2]->Fill(mu_mom,mu_ang,norm_tot);

      if(sim_mu_400MeV){
        h[227+sim_mod]->Fill(nuE   ,norm_tot);
        h[230+sim_mod]->Fill(mu_mom,norm_tot);
        h[233+sim_mod]->Fill(mu_ang,norm_tot);
      }
      if(sim_mu_45deg){
        h[236+sim_mod]->Fill(nuE   ,norm_tot);
        h[239+sim_mod]->Fill(mu_mom,norm_tot);
        h[242+sim_mod]->Fill(mu_ang,norm_tot);
      }
      if(sim_mu_30deg){
        h[245+sim_mod]->Fill(nuE   ,norm_tot);
        h[248+sim_mod]->Fill(mu_mom,norm_tot);
        h[251+sim_mod]->Fill(mu_ang,norm_tot);
      }
    }

    // ========================================
    // CCinc interaction in FV & Hits in INGRID
    if(sim_inFV&&sim_muon&&sim_hitING){
      h[176+sim_mod]->Fill(nuE   ,norm_tot);
      h[179+sim_mod]->Fill(mu_mom,norm_tot);
      h[182+sim_mod]->Fill(mu_ang,norm_tot);
      h2[11+sim_mod*NumSelec2]->Fill(mu_mom,mu_ang,norm_tot);
    }

    // ========================================
    // CCinc interaction in FV & Stop in INGRID
    if(sim_inFV&&sim_muon&&sim_stopING){
      h[167+sim_mod]->Fill(nuE   ,norm_tot);
      h[170+sim_mod]->Fill(mu_mom,norm_tot);
      h[173+sim_mod]->Fill(mu_ang,norm_tot);
      h2[12+sim_mod*NumSelec2]->Fill(mu_mom,mu_ang,norm_tot);
    }




    // -------------------
    // Hits and 2D Recons
    // -------------------
    for(int icyc=4; icyc<13;icyc++){
      int b2mod[3] = {14,16,21};
      for(int imod=0;imod<3;imod++){
        if(icyc==13&&imod!=2){continue;}
        int nhit = evtsum->NModHits(b2mod[imod],icyc);

        vector<Hit> v_hit(0);
        for(int ihit=0;ihit<nhit;ihit++){ 
          HitSummary* hit = evtsum->GetModHit(ihit,b2mod[imod],icyc);
          int    view = hit->view;
          int    mod  = hit->mod;
          double pe   = hit->pe;
          double time = hit->time;
          if(pe<3.5){continue;}

          Hit c_hit;
          c_hit.clear();
          c_hit.mod  = mod;
          c_hit.view = view;
          c_hit.pe   = pe;
          c_hit.time = time;
          v_hit.push_back(c_hit);
        }
        fSortTime(v_hit);
        fHitTimeWindow(v_hit,50,5);
        fTimeClustering(v_hit,50);

        nhit = (int)v_hit.size();

        // =====================
        // Hit Timing Clustering
        for(int ihit=0;ihit<nhit;ihit++){
          double time = v_hit[ihit].time;
          double hittime;
          if(fdid==0){
            hittime = time - (581.*icyc+hittime_off[imod]);
            if(imod==2) hittime -= wg_hittime_corr;
          }
          else{
            hittime = time;
          }
          h[18+imod]->Fill(hittime,norm_tot);
          if(v_hit[ihit].used==1){
            h[21+imod]->Fill(hittime,norm_tot);
          }
          if(v_hit[ihit].id==1){
            h[24+imod]->Fill(hittime,norm_tot);
            h[27+imod]->Fill(v_hit[ihit].time-v_hit[ihit].tdc,norm_tot);
          }
        }

        for(int iview=0;iview<2;iview++){
          int nrecon = evtsum->NModTwoDimRecons(b2mod[imod],icyc,iview);
          for(int irecon=0;irecon<nrecon;irecon++){
            TwoDimReconSummary* recon = evtsum->GetModTwoDimRecon(irecon,b2mod[imod],icyc,iview);
            double clstime = recon->clstime;
            int    nrechit = recon->Nhits();
            float  intcpt  = recon->intcpt;
            float  slope   = recon->slope;

            // ==========================
            // Timing of Hits on 2D Track
            for(int irechit=0;irechit<nrechit;irechit++){
              HitSummary* rechit = recon->GetHit(irechit);
              double hittime = rechit->time - clstime;
              h[imod]->Fill(hittime,norm_tot);
            }

            // ===============================
            // 2D Track Matching (INGRID/PM Hits)
            if(imod==2){
              int ninghit = evtsum->NModHits(14,icyc);
              for(int iinghit=0;iinghit<ninghit;iinghit++){
                HitSummary* inghit = evtsum->GetModHit(iinghit,14,icyc);
                int ingview = inghit->view;
                int ingpln  = inghit->pln;
                int ingch   = inghit->ch;
                if(ingview!=iview){continue;}
                if(ingpln > 1    ){continue;}
                double zpos  = zposi (14,ingview,ingpln,ingch);
                double xypos = xyposi(14,ingview,ingpln,ingch);
                zpos += C_B2INGPosZ - C_B2WMPosZ;
                if(ingview==0){ xypos += C_B2INGPosY - C_B2WMPosY; }
                else          { xypos += C_B2INGPosX - C_B2WMPosX; }
                double ingdist = (intcpt + slope*zpos) - xypos;
                h[39]->Fill(ingdist,norm_tot);
              }
              int npmhit = evtsum->NModHits(16,icyc);
              for(int ipmhit=0;ipmhit<npmhit;ipmhit++){
                HitSummary* pmhit = evtsum->GetModHit(ipmhit,16,icyc);
                int pmview = pmhit->view;
                int pmpln  = pmhit->pln;
                int pmch   = pmhit->ch;
                if(pmview!=iview){continue;}
                if(pmpln < 16   ){continue;}
                double zpos  = zposi (16,pmview,pmpln,pmch);
                double xypos = xyposi(16,pmview,pmpln,pmch);
                zpos += C_B2CHPosZ - C_B2WMPosZ;
                if(pmview==0){ xypos += C_B2CHPosY - C_B2WMPosY; }
                else         { xypos += C_B2CHPosX - C_B2WMPosX; }
                double pmdist = (intcpt + slope*zpos) - xypos;
                h[40]->Fill(pmdist,norm_tot);
              }
            }
          }
        }
      }
    }//for(icyc)


    // -----------
    // Vertices
    // -----------
    int nvtx = evtsum->NThreeDimRecons();
    for(int ivtx=0;ivtx<nvtx;ivtx++){

      bool sel_vtx      = false;
      bool sel_ironpene = false;
      bool sel_beamtime = false;
      bool sel_veto     = false;
      bool sel_fvcut    = false;
      bool sel_trkangle = false;
      bool sel_ingside  = false;
      bool sel_ingdwst  = false;
      bool sel_onetrk   = false;
      bool sel_accept   = false;
      bool sim_vtxmod   = false;

      ThreeDimReconSummary *vtx = evtsum->GetThreeDimRecon(ivtx);
      int ntrk = vtx->Ntrack;
      if(ntrk< 0){ sel_vtx   =false; continue; }
      else       { sel_vtx   =true ;}
      if(ntrk==1){ sel_onetrk=true ;}

      int startmod0 = vtx->startmod[0];
      int imod = -1;
      if     (startmod0==14){imod=0;}
      else if(startmod0==16){imod=1;}
      else if(startmod0==21){imod=2;}
      else{continue;}

      if(sim_mod==-1||sim_mod==imod){sim_vtxmod=true;}

      double vertex_xy_1st[2] = {0.,0.};
      int    vertex_z_1st [2] = {0,0};
      double angle_1sttrk = -1.;
      int    ironpene = -1;
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

        if(startmod!=startmod0){continue;}

        double angle = fabs(atan(sqrt(pow(tan(thetax),2)+pow(tan(thetay),2))))*180./PI;
        if(itrk==0){ angle_1sttrk = angle; }
        //if(ivtx==0&&itrk==0){ angle_1sttrk = angle; }


        float vx  = startych;
        float vy  = startxch;
        float vz;
        float vz1 = zposi(startmod, 0, startxpln, 0); //side
        float vz2 = zposi(startmod, 1, startypln, 0); //top
        if(vz2<vz1){vz=vz2;}else{vz=vz1;}


        int startpln = -1;
        int endpln   = -1;
        if(itrk==0){

          if(stopmodx!=14||stopmody!=14){break;}

          //if(vz1<vz2){startpln=startxpln;}else{startpln=startypln;}
          startpln = startypln;

          // ================
          // Iron Penetration
          if(endxpln>endypln){endpln=endxpln;}else{endpln=endypln;}
          if(itrk==0){
            if(endpln>9){ ironpene=9; }else{ ironpene=endpln; }
            h[imod+3]->Fill(ironpene, norm_tot); 
          }


          // ================
          // Track Stop Pos
          h[107+imod]->Fill(endych);
          h[110+imod]->Fill(endxch);

          // ********************
          // Iron penetration cut.
          // ********************
          if     (imod==0&&endpln<IronPeneCut0){sel_ironpene=false;}
          else if(imod==1&&endpln<IronPeneCut1){sel_ironpene=false;}
          else if(imod==2&&endpln<IronPeneCut2){sel_ironpene=false;}
          else{sel_ironpene=true;}


          // *******************
          // INGRID Side Escape
          // *******************
          if(fabs(endxch)<INGStopXY && fabs(endych)<INGStopXY){ sel_ingside=true;  }
          else                                                { sel_ingside=false; }

          // *******************
          // INGRID Downstream 
          // *******************
          if(endpln<INGStopPlnZ){sel_ingdwst=true ;}
          else                  {sel_ingdwst=false;}

        }

        // ================
        // 3D Track Matching
        if(oneview==0){
          h[15+imod]->Fill(startxpln-startypln, norm_tot);
        }

        int nanahit = vtx->Nhits();
        double largestpe = 0.;
        double diff_exptime;
        for(int ianahit=0;ianahit<nanahit;ianahit++){
          HitSummary* anahit = vtx->GetHit(ianahit);
          int   anahitview = anahit->view;
          int   anahitmod  = anahit->mod;
          int   anahitpln  = anahit->pln;
          float anahittime = anahit->time;
          float anahitpe   = anahit->pe;
          float anahitcyc  = anahit->cyc;
          if(fdid==0){
            if(anahitmod==21) anahittime -= wg_hittime_corr;
          }

          int cmod = -1;
          if     (anahitmod==14){cmod=0;}
          else if(anahitmod==16){cmod=1;}
          else if(anahitmod==21){cmod=2;}
          else{continue;}

          if(anahitmod!=21&&largestpe<anahitpe){
            largestpe = anahitpe;
            if(fdid!=0){
              diff_exptime = anahittime;
            }else{
              diff_exptime = anahittime - (581.*anahitcyc + hittime_off[cmod]);
            }
          }

          // =====================
          // Hit Timing on 3D Track
          h[36+cmod]->Fill(anahittime, norm_tot);
        }

        // ****************
        // Beam Timing Cut
        if(itrk==0){
          if(fdid==0 && fabs(diff_exptime)<BeamTimingCut){
            sel_beamtime=true ;
          }
          else if(fdid!=0 && fabs(diff_exptime)<BeamTimingCut_MC){
            sel_beamtime=true;
          }
          else{sel_beamtime=false;}

          // ======================
          // Event Timing Diff from Exp Time
          h[53+imod]->Fill(diff_exptime, norm_tot);
        }
        sel_beamtime = sel_beamtime || beam_timing_mask;


        // ================
        // 2D Track Matching
        if(endxpln>1){
          //WM-ING
          if(diff_angx  !=-1000.){ h[ 6]->Fill(diff_angx  ,norm_tot);} 
          if(diff_posx  !=-1000.){ h[ 9]->Fill(diff_posx  ,norm_tot);} 
          if(diff_timex !=-1000.){ h[12]->Fill(diff_timex ,norm_tot);} 
          //PM-WM                                                
          if(diff_angx2 !=-1000.){ h[ 7]->Fill(diff_angx2 ,norm_tot);} 
          if(diff_posx2 !=-1000.){ h[10]->Fill(diff_posx2 ,norm_tot);} 
          if(diff_timex2!=-1000.){ h[13]->Fill(diff_timex2,norm_tot);} 
          //PM-ING                                               
          if(diff_angx3 !=-1000.){ h[ 8]->Fill(diff_angx3 ,norm_tot);} 
          if(diff_posx3 !=-1000.){ h[11]->Fill(diff_posx3 ,norm_tot);} 
          if(diff_timex3!=-1000.){ h[14]->Fill(diff_timex3,norm_tot);} 
        }                                                        
        if(endypln>1){                                           
          //WM-ING                                               
          if(diff_angy  !=-1000.){ h[ 6]->Fill(diff_angy  ,norm_tot);} 
          if(diff_posy  !=-1000.){ h[ 9]->Fill(diff_posy  ,norm_tot);} 
          if(diff_timey !=-1000.){ h[12]->Fill(diff_timey ,norm_tot);} 
          //PM-WM                                                
          if(diff_angy2 !=-1000.){ h[ 7]->Fill(diff_angy2 ,norm_tot);} 
          if(diff_posy2 !=-1000.){ h[10]->Fill(diff_posy2 ,norm_tot);} 
          if(diff_timey2!=-1000.){ h[13]->Fill(diff_timey2,norm_tot);} 
          //PM-ING                                               
          if(diff_angy3 !=-1000.){ h[ 8]->Fill(diff_angy3 ,norm_tot);} 
          if(diff_posy3 !=-1000.){ h[11]->Fill(diff_posy3 ,norm_tot);} 
          if(diff_timey3!=-1000.){ h[14]->Fill(diff_timey3,norm_tot);} 
        }


        // ================
        // Vertexing
        if(itrk==0){ 
          vertex_xy_1st[0] = startych; //x
          vertex_xy_1st[1] = startxch; //y
          vertex_z_1st [0] = startypln;
          vertex_z_1st [1] = startxpln;
        }
        else if(oneview==0){
          double diff_vtx_xy = sqrt(pow(startych-vertex_xy_1st[0],2)+pow(startxch-vertex_xy_1st[1],2));
          int    diff_vtx_z  = abs(startypln-vertex_z_1st[0]) + abs(startxpln-vertex_z_1st[1]);
          h[30+imod]->Fill(diff_vtx_z ,norm_tot);
          h[33+imod]->Fill(diff_vtx_xy,norm_tot);
        }

        if(itrk==0){

          // ================
          // Vertex Position 
          h[41+imod]->Fill(startpln ,norm_tot);
          h[44+imod]->Fill(startych ,norm_tot);
          h[47+imod]->Fill(startxch ,norm_tot);
          h2[imod*2+0]->Fill(vz1,vy ,norm_tot); //side
          h2[imod*2+1]->Fill(vz2,vx ,norm_tot); //top


          // *******************
          // Upstream Veto Cut
          // *******************
          if     (imod==0&&(startxpln<=VetoCutPln0||startypln<=VetoCutPln0)){ sel_veto=false;}
          else if(imod==1&&(startxpln<=VetoCutPln1||startypln<=VetoCutPln1)){ sel_veto=false;}
          else if(imod==2&&(startxpln<=VetoCutPln2||startypln<=VetoCutPln2)){ sel_veto=false;}
          else                                                              { sel_veto=true; }

          // *******************
          // Fiducial Volume Cut
          // *******************
          if     (imod==0&&fabs(startxch)<=FVcutY0&&fabs(startych)<=FVcutX0){sel_fvcut=true; }
          else if(imod==1&&fabs(startxch)<=FVcutY1&&fabs(startych)<=FVcutX1){sel_fvcut=true; }
          else if(imod==2&&fabs(startxch)<=FVcutY2&&fabs(startych)<=FVcutX2){sel_fvcut=true; }
          else                                                              {sel_fvcut=false;}


          // ***********************
          // Reconstructed Angle Cut
          // ***********************
          if(angle<30.){sel_trkangle=true; }
          else         {sel_trkangle=false;}

          // ***************************
          // Additional Acceptance Limit
          // ***************************
          double pos_pm_ing = zposi(14,0,2,0) + C_B2INGPosZ - (zposi(16,0,15,0)+C_B2CHPosZ);
          if     (startmod==16){sel_accept=true;}
          else if(startmod==21){
            double pos_ing  = pos_pm_ing + zposi(21,0,15,0) + ACCEPT_OFFSETZ;
            double x_at_ing = vx + tan(thetay)*(pos_ing-vz2);
            double y_at_ing = vy + tan(thetax)*(pos_ing-vz1);
            double offsetx = C_B2WMPosX - C_B2INGPosX - ACCEPT_OFFSETX;
            double offsety = C_B2WMPosY - C_B2INGPosY - ACCEPT_OFFSETY;
            x_at_ing += offsetx;
            y_at_ing += offsety;
            if(   fabs(x_at_ing)<C_INGScintiLength/2.
                &&fabs(y_at_ing)<C_INGScintiLength/2.
              )  {sel_accept=true; }
            else {sel_accept=false;}
          }
          else if(startmod==14){
            //double pos_ing  = zposi(14,0,2,0) + 1676.5; //1569.5;
            double pos_ing  = pos_pm_ing + zposi(14,0,8,0) + ACCEPT_OFFSETZ;
            double x_at_ing = vx + tan(thetay)*(pos_ing-vz2);
            double y_at_ing = vy + tan(thetax)*(pos_ing-vz1);
            double offsetx = C_B2CHPosX - C_B2INGPosX - ACCEPT_OFFSETX;
            double offsety = C_B2CHPosY - C_B2INGPosY - ACCEPT_OFFSETY;
            x_at_ing += offsetx;
            y_at_ing += offsety;
            if(   fabs(x_at_ing)<C_INGScintiLength/2.
                &&fabs(y_at_ing)<C_INGScintiLength/2.
              )  {sel_accept=true; }
            else {sel_accept=false;}
          }

          // ==========================================================
          // Vertex Position (After Iron Penetration & Beam Timing Cut)
          // & Track angle
          if(sel_ironpene&&sel_beamtime){
            h[56+imod]->Fill(startpln ,norm_tot);
            h[59+imod]->Fill(startych ,norm_tot);
            h[62+imod]->Fill(startxch ,norm_tot);
            if(sel_veto){
              h[65+imod]->Fill(startych ,norm_tot);
              h[68+imod]->Fill(startxch ,norm_tot);
              if(sel_fvcut){
                h[ 71+imod]->Fill(angle    ,norm_tot);
                h[188+imod]->Fill(startpln ,norm_tot);
                h[191+imod]->Fill(startych ,norm_tot);
                h[194+imod]->Fill(startxch ,norm_tot);
              }
            }
          }

          // ===================================
          // Track Stop Pos After Track Angle Cut
          if(   sel_vtx
              &&sel_ironpene
              &&sel_beamtime
              &&sel_veto
              &&sel_fvcut
              &&sel_trkangle)
          {
            h[113+imod]->Fill(endpln,norm_tot);
            h[116+imod]->Fill(endych,norm_tot);
            h[119+imod]->Fill(endxch,norm_tot);
            if(sel_ingside){
              h[185+imod]->Fill(endpln,norm_tot);
            }
          }

          // ==================
          // Reconstructed Angle
          h[50+imod]->Fill(angle,norm_tot);

          if(   sel_vtx
              &&sel_ironpene
              &&sel_beamtime
              &&sel_veto
              &&sel_fvcut)
          {
            h[158+imod]->Fill(angle,norm_tot);
            if(sel_trkangle){
              if(sel_ingside){
                h[161+imod]->Fill(angle,norm_tot);
                if(sel_ingdwst){
                  h[164+imod]->Fill(angle,norm_tot);
                }
              }
              if(sel_onetrk){
                h[206+imod]->Fill(angle,norm_tot);
              }
            }
          }


          if(   sel_vtx
              &&sel_ironpene
              &&sel_beamtime
              &&sel_veto
              &&sel_fvcut)
          {
            // ==================
            // Number of Tracks
            h[197+imod]->Fill(ntrk ,norm_tot);
            if(sel_trkangle){
              h[200+imod]->Fill(ntrk ,norm_tot);
              if(sel_ingside){
                if(sel_ingdwst){
                  h[203+imod]->Fill(ntrk ,norm_tot);
                }
              }
              if(sel_accept){
                h[284+imod]->Fill(ntrk ,norm_tot);
              }
            }

            // ===================
            // Particle Info
            int num_trkpdg[5] = {0,0,0,0,0};
            for(int ianahit=0;ianahit<nanahit;ianahit++){
              HitSummary* anahit = vtx->GetHit(ianahit);
              if(anahit->NSimHits()<1){continue;}
              SimHitSummary* anasimhit = anahit->GetSimHit(0);
              int anapdg = abs(anasimhit->pdg);
              int pdg_id;
              if     (anapdg==  13              ){pdg_id=0;} // muon
              else if(anapdg==  11||anapdg==  22){pdg_id=1;} // e/gamma
              else if(anapdg==2212||anapdg==2112){pdg_id=2;} // p/n
              else if(anapdg== 211||anapdg== 111){pdg_id=3;} // pion
              else                               {pdg_id=4;} // others
              h[209+imod]->Fill(pdg_id,norm_tot);
              num_trkpdg[pdg_id]++;
            }
            int trkpdg    = -1;
            int numtrkpdg = -1;
            for(int ii=0;ii<5;ii++){
              if(numtrkpdg<num_trkpdg[ii]){
                numtrkpdg = num_trkpdg[ii];
                trkpdg    = ii;
              }
            }
            h[212+imod]->Fill(trkpdg,norm_tot);
            if(sel_trkangle&&sel_onetrk){
              h[215+imod]->Fill(trkpdg ,norm_tot);
            }

            if(sel_trkangle&&sel_accept){
              h[290+imod]->Fill(trkpdg ,norm_tot);
              if(sel_onetrk){
                h[293+imod]->Fill(trkpdg ,norm_tot);
              }
            }

          }

          // ===============
          // True Vertex
          if(sel_vtx){
            h2[30+2*imod]->Fill(vertex[2]/10.,vertex[1]/10.,norm_tot);
            h2[31+2*imod]->Fill(vertex[2]/10.,vertex[0]/10.,norm_tot);
            if(   sel_ironpene
                &&sel_beamtime
                &&sel_veto
                &&sel_fvcut)
            {
              h2[36+2*imod]->Fill(vertex[2]/10.,vertex[1]/10.,norm_tot);
              h2[37+2*imod]->Fill(vertex[2]/10.,vertex[0]/10.,norm_tot);
            }
          }

          // =================
          // Resolution [Vertex]
          if(sim_muon){
            if(   sel_vtx
                &&sel_ironpene
                &&sel_beamtime)
            {
              int truemod = -1;
              if     (sim_mod==0){truemod=14;}
              else if(sim_mod==1){truemod=16;}
              else if(sim_mod==2){truemod=21;}
              int true_startpln = GetTrueStartPln(truemod,vertex[2]-CoT[sim_mod][2]);
              int reco_startpln = startpln;
              if     (sim_mod==0){true_startpln+=42;}
              else if(sim_mod==1){true_startpln+= 0;}
              else if(sim_mod==2){true_startpln+=18;}
              if     (imod   ==0){reco_startpln+=42;}
              else if(imod   ==1){reco_startpln+= 0;}
              else if(imod   ==2){reco_startpln+=18;}
              if(sim_inFVarea){
                h2[54]->Fill(true_startpln,reco_startpln,norm_tot);
              }
              if(sel_veto){
                if(imod==sim_mod){
                  h2[48+imod]->Fill(vertex[0]-CoT[imod][0],vx,norm_tot);
                  h2[51+imod]->Fill(vertex[1]-CoT[imod][1],vy,norm_tot);
                }
              }
            }
          }

          // =================
          // U & P Matrix [Angle]
          if(   sel_vtx
              &&sel_ironpene
              &&sel_beamtime
              &&sel_veto
              &&sel_fvcut)
          {
            if(   sim_muon
                &&sim_inFV
                &&sim_mu_400MeV)
            {
              if(sim_mod==imod){
                h2[42+imod]->Fill(mu_ang,angle,norm_tot);
                if(sel_accept){
                }
                if(sel_onetrk){
                  h2[45+imod]->Fill(mu_ang,angle,norm_tot);
                  if(sel_accept){
                  }
                }
              }
            }
          }

          // ===================
          // Light Yield
          if(   sel_vtx
              &&sel_ironpene
              &&sel_beamtime
              &&sel_veto
              &&sel_fvcut
              &&sel_accept
              &&sel_onetrk
            )
          {
            for(int ianahit=0;ianahit<nanahit;ianahit++){
              HitSummary* anahit = vtx->GetHit(ianahit);
              int   anahitview = anahit->view;
              int   anahitmod  = anahit->mod;
              int   anahitpln  = anahit->pln;
              float anahittime = anahit->time;
              float anahitpe   = anahit->pe;
              float anahitcyc  = anahit->cyc;
              int cmod = -1;
              if     (anahitmod==14){cmod=0;}
              else if(anahitmod==16){cmod=1;}
              else if(anahitmod==21){cmod=2;}
              else{continue;}
              h[323+cmod]->Fill(anahitpe,norm_tot);
            }
          }

        } //if(itrk==0)
      } //for(itrk)


      // **********************************************
      //   8-bin Analysis
      bool recobin[8];
      bool presel = sel_vtx && sel_ironpene && sel_beamtime && sel_veto && sel_fvcut && sel_accept;
      recobin[0] = presel && !sel_onetrk;
      recobin[1] = presel &&  sel_onetrk && (angle_1sttrk>= 0.&&angle_1sttrk<  5.); 
      recobin[2] = presel &&  sel_onetrk && (angle_1sttrk>= 5.&&angle_1sttrk< 10.); 
      recobin[3] = presel &&  sel_onetrk && (angle_1sttrk>=10.&&angle_1sttrk< 15.); 
      recobin[4] = presel &&  sel_onetrk && (angle_1sttrk>=15.&&angle_1sttrk< 20.); 
      recobin[5] = presel &&  sel_onetrk && (angle_1sttrk>=20.&&angle_1sttrk< 25.); 
      recobin[6] = presel &&  sel_onetrk && (angle_1sttrk>=25.&&angle_1sttrk< 30.); 
      recobin[7] = presel &&  sel_onetrk && (angle_1sttrk>=30.&&angle_1sttrk<180.); 
#ifdef DEBUG_B2PLOT
      bool check_recobin =false;
      for(int ibin=0;ibin<8;ibin++){
        if(recobin[ibin]){check_recobin=true;break;}
      }
      if(!check_recobin&&presel){
        cout << "Found an event out of reconstruction binning.." << endl;
        exit(1);
      }
#endif

      // ===============
      // Fill plots for 8-bin Analysis
      for(int ibin=0;ibin<8;ibin++){
        if(truebin[ibin]){
          if(presel                      ){ h[389+sim_mod]->Fill(ibin,norm_tot); }
          if(presel  &&sel_onetrk        ){ h[392+sim_mod]->Fill(ibin,norm_tot); }
          if(int_onCH&&presel            ){ h[416+sim_mod]->Fill(ibin,norm_tot); }
          if(int_onCH&&presel&&sel_onetrk){ h[419+sim_mod]->Fill(ibin,norm_tot); }
          break;
        }
      }
      for(int ibin=0;ibin<8;ibin++){
        if(recobin[ibin]){
          h[395+imod]->Fill(ibin,norm_tot);
          if(sel_onetrk            ){h[398+imod]->Fill(ibin,norm_tot);}
          if(sim_ooFV              ){h[401+imod]->Fill(ibin,norm_tot);}
          if(sim_ooFV &&sel_onetrk ){h[404+imod]->Fill(ibin,norm_tot);}
          if(nc_events             ){h[407+imod]->Fill(ibin,norm_tot);}
          if(nc_events&&sel_onetrk ){h[410+imod]->Fill(ibin,norm_tot);}
          break;
        }
      }
      for(int ibin=0;ibin<8;ibin++){
        for(int jbin=0;jbin<8;jbin++){
          if(recobin[jbin]&&truebin[ibin]){
            h2[100+sim_mod]->Fill(ibin,jbin,norm_tot);
            h2[106+sim_mod]->Fill(jbin,ibin,norm_tot);
            if(sel_onetrk){
              h2[103+sim_mod]->Fill(ibin,jbin,norm_tot);
              h2[109+sim_mod]->Fill(jbin,ibin,norm_tot);
            }
          }
        }
      }

      if(presel){
        h2[112+sim_mod]->Fill(vertex_z_1st[0],angle_1sttrk,norm_tot);
        h[428+imod]->Fill(ntrk,norm_tot);
        if(sel_onetrk){
          h[431+imod]->Fill(angle_1sttrk,norm_tot);
          if(angle_1sttrk<30){
            h[434+imod]   ->Fill(vertex_z_1st [0],norm_tot);
            h[437+imod]   ->Fill(vertex_xy_1st[0],norm_tot);
            h[440+imod]   ->Fill(vertex_xy_1st[1],norm_tot);
            h[443+sim_mod]->Fill(nuE     ,norm_tot);
            h[446+sim_mod]->Fill(mu_mom  ,norm_tot);
            h[449+sim_mod]->Fill(mu_ang  ,norm_tot);
            h[452+imod]   ->Fill(ironpene,norm_tot);
          }
        }
      }


      // ===============
      // Fill Plots for non-limited phase space
      if(sim_inFV&&sim_muon){ 
        if(sel_vtx){
          if(sel_ironpene){
            if(sel_beamtime){
              if(sel_veto){
                if(sel_fvcut){
                  h[89+sim_mod]->Fill(nuE   ,norm_tot);
                  h[92+sim_mod]->Fill(mu_mom,norm_tot);
                  h[95+sim_mod]->Fill(mu_ang,norm_tot);
                  if(sel_trkangle){
                    h[ 98+sim_mod]->Fill(nuE   ,norm_tot);
                    h[101+sim_mod]->Fill(mu_mom,norm_tot);
                    h[104+sim_mod]->Fill(mu_ang,norm_tot);
                    if(sel_onetrk){
                      //h3[10+imod*NumSelec3]->Fill(0.,norm_tot);
                      //h2[13+sim_mod*NumSelec2]->Fill(mu_mom,mu_ang,norm_tot);
                      h[218+sim_mod]->Fill(nuE   ,norm_tot);
                      h[221+sim_mod]->Fill(mu_mom,norm_tot);
                      h[224+sim_mod]->Fill(mu_ang,norm_tot);
                      if(sel_accept){
                        h[263+sim_mod]->Fill(nuE   ,norm_tot);
                        h[266+sim_mod]->Fill(mu_mom,norm_tot);
                        h[269+sim_mod]->Fill(mu_ang,norm_tot);
                        h2[64+sim_mod]->Fill(mu_mom,mu_ang,norm_tot);
                      }
                    }
                    if(sel_accept){
                      h[254+sim_mod]->Fill(nuE   ,norm_tot);
                      h[257+sim_mod]->Fill(mu_mom,norm_tot);
                      h[260+sim_mod]->Fill(mu_ang,norm_tot);
                      h2[61+sim_mod]->Fill(mu_mom,mu_ang,norm_tot);
                    }
                    if(sel_ingside){
                      //h3[8+imod*NumSelec3]->Fill(0.,norm_tot);
                      h[122+sim_mod]->Fill(nuE   ,norm_tot);
                      h[125+sim_mod]->Fill(mu_mom,norm_tot);
                      h[128+sim_mod]->Fill(mu_ang,norm_tot);
                      h2[9+sim_mod*NumSelec2]->Fill(mu_mom,mu_ang,norm_tot);
                      if(sel_ingdwst){ //INGRID stopping events
                        //h3[9+imod*NumSelec3]->Fill(0.,norm_tot);
                        h[131+sim_mod]->Fill(nuE   ,norm_tot);
                        h[134+sim_mod]->Fill(mu_mom,norm_tot);
                        h[137+sim_mod]->Fill(mu_ang,norm_tot);
                        h2[10+sim_mod*NumSelec2]->Fill(mu_mom,mu_ang,norm_tot);
                      }
                      else{ //Non-Stopping Event
                        h[149+sim_mod]->Fill(nuE   ,norm_tot);
                        h[152+sim_mod]->Fill(mu_mom,norm_tot);
                        h[155+sim_mod]->Fill(mu_ang,norm_tot);
                      }
                    }
                    else{ //Side Escaping Event
                      h[140+sim_mod]->Fill(nuE   ,norm_tot);
                      h[143+sim_mod]->Fill(mu_mom,norm_tot);
                      h[146+sim_mod]->Fill(mu_ang,norm_tot);
                    }
                  }
                }
              }
            }
          }
        }
      }

      bool signal   = sim_inFV && sim_muon && sim_mu_400MeV && sim_mu_30deg;
      bool bg_nc    = !sim_muon;
      bool bg_outfv = sim_muon && !sim_inFV;
      bool bg_outps = sim_muon &&  sim_inFV && ( !sim_mu_400MeV || !sim_mu_30deg );


      // ===============
      // Fill Number
      if(sel_vtx){
        h3[2+imod*NumSelec3]->Fill(0.,norm_tot);
        ///////////////////
        if(sel_ironpene){
          h3[3+imod*NumSelec3]->Fill(0.,norm_tot);
          ///////////////////
          if(sel_beamtime){
            h3[4+imod*NumSelec3]->Fill(0.,norm_tot);
            ///////////////////
            if(sel_veto){
              h3[5+imod*NumSelec3]->Fill(0.,norm_tot);
              ///////////////////
              if(sel_fvcut){
                h3[6+imod*NumSelec3]->Fill(0.,norm_tot);
                ///////////////////
                if(sel_accept){
                  h3[7+imod*NumSelec3]->Fill(0.,norm_tot);
                  ///////////////////
                  if(sel_onetrk){
                    h3[11+imod*NumSelec3]->Fill(0.,norm_tot);
                    ///////////////////
                    if(sel_trkangle){
                      h3[12+imod*NumSelec3]->Fill(0.,norm_tot);
                    }
                  }
                }
              }
            }
          }
        }
      }

      // ===============
      // Fill Plot
      if(sel_vtx){
        if(sel_ironpene){
          if(sel_beamtime){
            if(sel_veto){
              if(sel_fvcut){
                h2[7+imod*NumSelec2]->Fill(mu_mom,mu_ang,norm_tot);
                if(sel_trkangle){
                  h2[8+imod*NumSelec2]->Fill(mu_mom,mu_ang,norm_tot);
                  if(sel_accept){
                    h[272+imod]->Fill(angle_1sttrk,norm_tot);
                    if(signal  ){ 
                      h[335+imod]->Fill(nuE   ,norm_tot); //For detection eff.
                      h[338+imod]->Fill(mu_mom,norm_tot); //For detection eff.
                      h[341+imod]->Fill(mu_ang,norm_tot); //For detection eff.
                      h2[55+imod]->Fill(mu_ang,angle_1sttrk,norm_tot); //U-Matrix
                      h2[70+imod]->Fill(angle_1sttrk,mu_ang,norm_tot); //P-Matrix
                      if(targetz==6){
                        h[359+imod]->Fill(mu_ang,norm_tot); //For detection eff. on CH
                        h2[76+imod]->Fill(angle_1sttrk,mu_ang,norm_tot); //P-Matrix
                      }
                    } 
                    else if(bg_nc   ){ h[353+imod]->Fill(angle_1sttrk,norm_tot); } 
                    else if(bg_outfv){ h[317+imod]->Fill(angle_1sttrk,norm_tot); }
                    else if(bg_outps){ h[278+imod]->Fill(angle_1sttrk,norm_tot); }

                    ///////////////////
                    if(sel_onetrk){
                      h[275+imod]->Fill(angle_1sttrk,norm_tot);
                      if(signal){
                        h[344+imod]->Fill(nuE   ,norm_tot); //For detection eff.
                        h[347+imod]->Fill(mu_mom,norm_tot); //For detection eff.
                        h[350+imod]->Fill(mu_ang,norm_tot); //For detection eff.
                        h2[58+imod]->Fill(mu_ang,angle_1sttrk,norm_tot); //U-Matrix
                        h2[73+imod]->Fill(angle_1sttrk,mu_ang,norm_tot); //P-Matrix
                        if(targetz==6){
                          h[362+imod]->Fill(mu_ang,norm_tot); //For detection eff. on CH
                          h2[79+imod]->Fill(angle_1sttrk,mu_ang,norm_tot); //P-Matrix
                        }
                      }
                      else if(bg_nc   ){ h[356+imod]->Fill(angle_1sttrk,norm_tot); }
                      else if(bg_outfv){ h[320+imod]->Fill(angle_1sttrk,norm_tot); }
                      else if(bg_outps){ h[281+imod]->Fill(angle_1sttrk,norm_tot); }
                    }
                    else{ //Multi-track events
                      if(signal){
                        h2[94+sim_mod]->Fill(pi_mom,pi_ang,norm_tot);
                        h[377+sim_mod]->Fill(pi_mom,norm_tot);
                        if(!sim_pion){
                          h2[97+sim_mod]->Fill(p_mom,p_ang,norm_tot);
                          h[380+sim_mod]->Fill(pi_mom,norm_tot);
                        }
                        h[383+sim_mod]->Fill(pi_mom,norm_tot);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }

    }//for(ivtx)

    // =======================================
    // Beam Information
    int nbeamsum = evtsum->NBeamSummarys();
    if(nbeamsum>0){
      BeamInfoSummary* beamsum = evtsum->GetBeamSummary(0);
      double pot = beamsum->ct_np[4][0];
      h3[0]->Fill(0.,pot); //POT
      h3[1]->Fill(0.); //spill
    }


  }//for(ievt)

  if(t2krewfilename=="local"){
    for(int i=0;i<3;i++){
      TH1F *h_tmp = (TH1F*)readfile->Get(Form("h3_%d",i));
      if(FluxErr){
        double fluxerr = 1.;
        int nbins = h_tmp->GetNbinsX();
        for(int ibin=1;ibin<=nbins;ibin++){
          double nuE = h_tmp->GetXaxis()->GetBinLowEdge(ibin);
          int offset;
          if     (nutype==2){offset=20;}
          else if(nutype==1){offset= 0;}
          else{
            cout << "Unexpected neutrino type: " << nutype << endl;
            exit(1);
          }
          int id;
          if     (nuE< 3.){ id = (int)(nuE/0.2)      + offset; fluxerr = 1.+fluxwei[id];}
          else if(nuE< 4.){ id = (int)(nuE/1.0) + 15 + offset; fluxerr = 1.+fluxwei[id];}
          else if(nuE<10.){ id = (int)(nuE/2.0) + 16 + offset; fluxerr = 1.+fluxwei[id];}
          else if(nuE<30.){ id = (int)(nuE/20.) + 19 + offset; fluxerr = 1.+fluxwei[id];}

          double cont = h_tmp->GetBinContent(ibin);
          cont *= fluxerr;
          h_tmp->SetBinContent(ibin,cont);
        }
      }
      h3[35+i]->Fill(0.,h_tmp->Integral());
    }
  }

  if(evtsum) delete evtsum;

  outfile->cd();
  if(!NselOnly){
    for(int i=0; i<NumHist; i++){
      if(h[i]!=NULL) h[i]->Write();
    }
    for(int i=0; i<NumHist2; i++){
      if(h2[i]!=NULL) h2[i]->Write();
    }
    for(int i=0; i<NumHist3; i++){
      if(h3[i]!=NULL) h3[i]->Write();
    }
  }
  else{
    for(int i=386; i<NumHist; i++){
      if(h[i]!=NULL) h[i]->Write();
    }
    for(int i=100; i<NumHist2; i++){
      if(h2[i]!=NULL) h2[i]->Write();
    }
    for(int i=0; i<NumHist3; i++){
      if(h3[i]!=NULL) h3[i]->Write();
    }
  }

  outfile->Close();

  std::cout << "Done!!!" << std::endl;
  return 0;
}
