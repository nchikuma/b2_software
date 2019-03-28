#include "B2Plot.hxx"

//#define DRAW_COS

FileStat_t fs;
void PrintUsage(){
  std::cout << "./TrackEfficiency [options]" << std::endl;
  std::cout << "-i <filename.root>   : Input ROOT file, after 2d reconstruction." << std::endl;
  std::cout << "-o <filename.root>   : Outpu ROOT file." << std::endl;
  std::cout << "-m <fdid>            : Select fdid;  0-> data, 7->WM, 8->PM, 9->ING, ... " << std::endl;
  std::cout << "-b                   : To open all bad channels. Default: masked." << std::endl;
  std::cout << std::endl;
}

int main(int argc, char** argv){

  //------- Arguments -------
  int c = -1;
  std::string ifilename;
  std::string ofilename = "tmp.root";
  bool BadCh_Ana = true;
  bool irename = false;
  bool orename = false;
  int fdid = -1;
  int sim_mod = -1;

  while ((c = getopt(argc, argv, "i:o:m:bh")) != -1){
    switch(c){
      case 'i':
        ifilename = optarg;
        irename = true;
        break;
      case 'o':
        ofilename = optarg;
        orename = true;
        break;
      case 'm':
        fdid = atoi(optarg);
        if(fdid==7){sim_mod=2;}
        if(fdid==8){sim_mod=1;}
        if(fdid==9){sim_mod=0;}
        break;
      case 'b':
        BadCh_Ana = false;
        break;
      case 'h':
        PrintUsage();
        exit(0);
      default:
        PrintUsage();
        exit(0);
    }
  }


  if((fdid<0)||(!irename||!orename))
  {
    PrintUsage();
    exit(1);
  }


  if(gSystem->GetPathInfo(ifilename.c_str(),fs)){
    std::cout << "No such a file: "<< ifilename << std::endl;
    exit(1);
  }

  // ----------------------------------------------------
  // ------- Open root file and get tree information -------
  TFile* ifile = new TFile(ifilename.c_str(),"read");
  if(ifile->IsZombie()){
    std::cout << "Cannot open file : " << ifilename << std::endl;
    exit(1);
  }
  std::cout << "Open file : " << ifilename << std::endl;
  TTree   *tree  = (TTree*)ifile ->Get("tree");
  TBranch *evtbr = tree->GetBranch("fDefaultReco.");
  EventSummary* evtsum = new EventSummary();
  evtbr -> SetAddress(&evtsum);
  tree  -> SetBranchAddress("fDefaultReco.",&evtsum);

  int nevt = (int)tree -> GetEntries();
  std::cout << "Total # of events = " << nevt << std::endl; 

  // ------------------ Flux Tuning File ----------------
  char tunefile[300];
  string tunefilename = "~/b2_data/jnubeam/tunefile/tune_nd7_8_9.root";
  sprintf(tunefile,"%s",tunefilename.c_str());
  FluxTuning* fluxtune;
  if(fdid==7||fdid==8||fdid==9) fluxtune = new FluxTuning(fdid,tunefile);

  //------- Define histograms ------
  TFile *ofile = new TFile(ofilename.c_str(),"recreate");
  std::cout << "Output file name : " << ofilename << std::endl;
  std::string name, xtitle, modname, axisname;

  const int NMOD_TRKEFF   =  2; //
  TH1F *h1[NMOD_TRKEFF];
  TH1F *h2[NMOD_TRKEFF];
  for(int i=0;i<NMOD_TRKEFF;i++){
    if     (i==0){modname="WM"   ;} //WAGASCI
    else if(i==1){modname="PM"   ;} //Proton Module
#ifdef DRAW_COS
    axisname="Z";xtitle="Track angle from the Z axis [cos]";
    name  =  Form("h1_%s_axis%s",modname.c_str(),axisname.c_str());
    h1[i] =  new TH1F(name.c_str(),name.c_str(),100,-1.,1.);
    h1[i] -> GetXaxis()->SetTitle(xtitle.c_str());
    name  =  Form("h2_%s_axis%s",modname.c_str(),axisname.c_str());
    h2[i] =  new TH1F(name.c_str(),name.c_str(),100,-1.,1.);
    h2[i] -> GetXaxis()->SetTitle(xtitle.c_str());
#else
    axisname="Z";xtitle="Track angle from the Z axis";
    name  =  Form("h1_%s_axis%s",modname.c_str(),axisname.c_str());
    h1[i] =  new TH1F(name.c_str(),name.c_str(),180,0.,180.);
    h1[i] -> GetXaxis()->SetTitle(xtitle.c_str());
    name  =  Form("h2_%s_axis%s",modname.c_str(),axisname.c_str());
    h2[i] =  new TH1F(name.c_str(),name.c_str(),180,0.,180.);
    h2[i] -> GetXaxis()->SetTitle(xtitle.c_str());
#endif

  }



  //------- Event loop -------
  detdim   = new DetectorDimension();
  badch  = new INGRID_BadCh_mapping();
  badch -> set_BadCh(BadCh_Ana);
  //bool newcanv = true;
  int start_cyc,stop_cyc;
  if(fdid==0){start_cyc= 4;stop_cyc=11;}
  else       {start_cyc= 4;stop_cyc= 4;}

  //int n_count = 0;
  for(int ievt=0;ievt<nevt;ievt++){

    if(ievt%100==0){
      std::cout << "Event # is " << ievt << std::endl;
    }

    tree  -> GetEntry(ievt);	


    double norm_tot=1.,reweight=1.;
    double norm=1.,totcrsne=1.;
    double nuE=-1.;
    int    nutype=-1;

    int nsimver = evtsum->NSimVertexes();
    if(nsimver>0){
      SimVertexSummary* simver = evtsum->GetSimVertex(0);
      norm      = simver->norm;
      totcrsne  = simver->totcrsne;
      nuE       = simver->nuE;
      nutype    = simver->nutype;

      if(fdid==7||fdid==8||fdid==9){
        reweight = fluxtune->reweight(nuE,nutype);
      }
      if(sim_mod>=0){
        norm_tot = norm*totcrsne*1e-38*Avogadro*density[sim_mod]*thickness[sim_mod]*reweight;
      }
      else{
        float thick = 0.;
        if     (fdid==10){thick=thickness_wallbg;}
        else if(fdid==11){thick=thickness_ceiling;}
        else if(fdid==12){thick=thickness_floor;}
        else if(fdid==13){thick=thickness_pillar_r;}
        else if(fdid==14){thick=thickness_pillar_l;}
        norm_tot = norm*totcrsne*1e-38*Avogadro*density_wallbg*thick;
      }
    }

    for(int icyc=start_cyc;icyc<=stop_cyc;icyc++){
      for(int imod=0;imod<NMOD_TRKEFF;imod++){
        int c_mod = -1;
        int t_mod = -1;
        int max_diff_cyc = 0;
        //double THRES_DIST = -1.; 
        if     (imod==0){c_mod=21;t_mod=14;max_diff_cyc=1;} //WAGASCI
        else if(imod==1){c_mod=16;t_mod=14;max_diff_cyc=0;} //Proton Module
        double offset[3] = {0.,0.,0.};
        double fv_endz;
        if(imod==0){
          offset[0] = C_B2INGPosX-C_B2WMPosX;
          offset[1] = C_B2INGPosY-C_B2WMPosY;
          offset[2] = C_B2INGPosZ-C_B2WMPosZ;
          fv_endz   = 59.1;
        }
        else if(imod==1){
          offset[0] = C_B2INGPosX-C_B2CHPosX;
          offset[1] = C_B2INGPosY-C_B2CHPosY;
          offset[2] = C_B2INGPosZ-C_B2CHPosZ;
          fv_endz   = 306.;
        }
  
        int n2drec[2] = {0,0};
        for(int iview=0;iview<2;iview++){
          n2drec[iview] = evtsum->NModTwoDimRecons(t_mod,icyc,iview);
        }
        if(n2drec[0]!=1||n2drec[1]!=1){continue;}


        double slopex=-1.,slopey=-1.;
        double intcptx=-1.,intcpty=-1.;
        bool ingrid_recon = true;
        //int endpln_tmp = -1;
        for(int iview=0;iview<2;iview++){
          TwoDimReconSummary* reconsum = evtsum->GetModTwoDimRecon(0,t_mod,icyc,iview);
          int   startpln = reconsum->startpln;
          //int   endpln   = reconsum->endpln;
          double slope   = reconsum->slope;
          double intcpt  = reconsum->intcpt;
          if(iview==0){slopey=slope;intcpty=intcpt;}
          else        {slopex=slope;intcptx=intcpt;}

          if(startpln!=0){ingrid_recon=false;break;}
        } //iview
        if(!ingrid_recon){continue;}

        double cmod_x = (fv_endz-offset[2])*slopex+intcptx+offset[0];
        double cmod_y = (fv_endz-offset[2])*slopey+intcpty+offset[1];
        if(fabs(cmod_x)>350.||fabs(cmod_y)>350.){continue;}

        double norm = sqrt(1.+slopex*slopex+slopey*slopey);
        double cos  = 1./norm;
        double angle = acos(cos)*180./PI;

        //h1[imod]->Fill(cos,norm_tot);
        h1[imod]->Fill(angle,norm_tot);

        bool target_recon[2] = {false,false};
        for(int iview=0;iview<2;iview++){
          for(int diff_cyc=0;diff_cyc<=max_diff_cyc;diff_cyc++){
            int cyc = icyc+diff_cyc;
            int n2dreco = evtsum->NModTwoDimRecons(c_mod,cyc,iview);
            for(int irec=0;irec<n2dreco;irec++){
              TwoDimReconSummary* reconsum = evtsum->GetModTwoDimRecon(irec,c_mod,cyc,iview);
              double endxy  = reconsum->endxy;
              //double endpln = reconsum->endpln;
              double slope  = reconsum->slope;
              double endz   = reconsum->endz;
              double ang_diff,pos_diff;
              if(iview==0){
                ang_diff = fabs(atan(slope)-atan(slopey))*180./PI;
                pos_diff = 
                  fabs(slopey*(endz-offset[2])-(endxy-offset[1])+intcpty)
                  /sqrt(1.+slopey*slopey);
              }
              else{
                ang_diff = fabs(atan(slope)-atan(slopex))*180./PI;
                pos_diff =
                  fabs(slopex*(endz-offset[2])-(endxy-offset[0])+intcptx)
                  /sqrt(1.+slopex*slopex);
              }
              if(ang_diff<35&&pos_diff<150){target_recon[iview]=true;}
            } //for(irec)
          } //for(diff_cyc)
        }
        if(target_recon[0]&&target_recon[1]){
          //h2[imod]->Fill(cos,norm_tot);
          h2[imod]->Fill(angle,norm_tot);
        }
      } //for(imod)
    } //icyc
  } //ievt

  ofile->cd();
  for(int i=0;i<NMOD_TRKEFF;i++){
    h1[i]->Write();
    h2[i]->Write();
  }
  ofile->Write();
  ofile->Close();

  delete badch;
  delete detdim;
  return 0;
}
