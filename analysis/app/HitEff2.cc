#include "B2Plot.hxx"

FileStat_t fs;
void PrintUsage(){
  std::cout << "./HitEff2 [options]" << std::endl;
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
  TH1F *h1[NMOD_TRKEFF][2];
  TH1F *h2[NMOD_TRKEFF][2];
  TH1F *h3[NMOD_TRKEFF][2];
  TH1F *h4[NMOD_TRKEFF][2];

  for(int i=0;i<NMOD_TRKEFF;i++){
    if     (i==0){modname="WM"   ;} //WAGASCI
    else if(i==1){modname="PM"   ;} //Proton Module

    for(int grid=0;grid<2;grid++){
      if(i==1&&grid==1){continue;}

      name  =  Form("h1_%s_%01d",modname.c_str(),grid);
      h1[i][grid] =  new TH1F(name.c_str(),name.c_str(),1400,-700.,700.);
      h1[i][grid] -> GetXaxis()->SetTitle(xtitle.c_str());
      name  =  Form("h2_%s_%01d",modname.c_str(),grid);
      h2[i][grid] =  new TH1F(name.c_str(),name.c_str(),1400,-700.,700.);
      h2[i][grid] -> GetXaxis()->SetTitle(xtitle.c_str());

      name  =  Form("h3_%s_%01d",modname.c_str(),grid);
      h3[i][grid] =  new TH1F(name.c_str(),name.c_str(),1400,-700.,700.);
      h3[i][grid] -> GetXaxis()->SetTitle(xtitle.c_str());
      name  =  Form("h4_%s_%01d",modname.c_str(),grid);
      h4[i][grid] =  new TH1F(name.c_str(),name.c_str(),1400,-700.,700.);
      h4[i][grid] -> GetXaxis()->SetTitle(xtitle.c_str());
    }
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

    if(evtsum->NThreeDimRecons()!=1){continue;}
    if(evtsum->GetThreeDimRecon(0)->Ningtrack!=1){continue;}


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
  
        int n2drec0[2] = {0,0};
        int n2drec1[2] = {0,0};
        for(int iview=0;iview<2;iview++){
          n2drec0[iview] += evtsum->NModTwoDimRecons(t_mod,icyc,iview);
          for(int iicyc=0;iicyc<=max_diff_cyc;iicyc++){
            n2drec1[iview] += evtsum->NModTwoDimRecons(c_mod,icyc+iicyc,iview);
          }
        }
        if(n2drec0[0]!=1||n2drec0[1]!=1){continue;}
        if(n2drec1[0]!=1||n2drec1[1]!=1){continue;}


        double slopex=-1.,slopey=-1.;
        double intcptx=-1.,intcpty=-1.;
        bool ingrid_recon = true;
        //int endpln_tmp = -1;
        for(int iview=0;iview<2;iview++){
          TwoDimReconSummary* reconsum = evtsum->GetModTwoDimRecon(0,t_mod,icyc,iview);
          int   startpln = reconsum->startpln;
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


        int maxpln=-1,maxch=-1;
        double dist_max;
        if     (imod==0){maxpln= 8;maxch=80;dist_max=25.;}
        else if(imod==1){maxpln=18;maxch=32;dist_max=25.;}

        // Penetrated planes
        bool hitplane    [2][20][2]; //view,pln,grid
        bool hitplane_hit[2][20][2]; //view,pln,grid
        // Find channels close to the INGRID track
        bool channel    [2][18][80];
        bool channel_hit[2][18][80];
        for(int view=0;view<2;view++){
          for(int pln=0;pln<maxpln;pln++){
            for(int ch=0;ch<maxch;ch++){
              if(imod==1&&pln==0&&ch>=24){continue;}
              channel    [view][pln][ch] = false;
              channel_hit[view][pln][ch] = false;
            }
          }
          for(int pln=0;pln<20;pln++){
            for(int grid=0;grid<2;grid++){
              hitplane    [view][pln][grid] = false;
              hitplane_hit[view][pln][grid] = false;
            }
          }
        }


        double rec_iz=-1e+5,rec_fz=1e+5;
        double rec_iy=-1e+5,rec_fy=1e+5;
        double rec_ix=-1e+5,rec_fx=1e+5;
        for(int iicyc=0;iicyc<=max_diff_cyc;iicyc++){
          for(int view=0;view<2;view++){
            if(evtsum->NModTwoDimRecons(c_mod,icyc+iicyc,view)!=1){continue;}
            TwoDimReconSummary *recsum = evtsum->GetModTwoDimRecon(0,c_mod,icyc+iicyc,view);
            double iz  = recsum->startz;
            double fz  = recsum->endz;
            double ixy = recsum->startxy;
            double fxy = recsum->endxy;
            if(rec_iz<iz){rec_iz=iz;}
            if(rec_fz>fz){rec_fz=fz;}
            if(view==0){
              if(rec_iy<ixy){rec_iy=ixy;}
              if(rec_fy>fxy){rec_fy=fxy;}
            }else{
              if(rec_ix<ixy){rec_ix=ixy;}
              if(rec_fx>fxy){rec_fx=fxy;}
            }
          }
        }

        for(int iicyc=0;iicyc<=max_diff_cyc;iicyc++){
          for(int view=0;view<2;view++){
            if(evtsum->NModTwoDimRecons(c_mod,icyc+iicyc,view)!=1){continue;}
            for(int pln=0;pln<maxpln;pln++){
              for(int ch=0;ch<maxch;ch++){
                if(imod==1&&pln==0&&ch>=24){continue;}
                int grid=0;
                if(imod==0&&ch>=40){grid=1;}

                double tmpx,tmpy,tmpz;
                detdim->GetPosInMod(c_mod,pln,view,ch,&tmpx,&tmpy,&tmpz);
                if( (view==0
                    &&(tmpy-rec_iy)*(tmpy-rec_fy)<0
                    &&(tmpz-rec_iz)*(tmpz-rec_fz)<0)
                  ||(view==1
                    &&(tmpx-rec_ix)*(tmpx-rec_fx)<0
                    &&(tmpz-rec_iz)*(tmpz-rec_fz)<0)
                    )
                {
                  int recpln = pln;
                  if(grid==1){ recpln = ch%20; }
                  hitplane[view][recpln][grid] = true;
                }


                tmpx -= offset[0];
                tmpy -= offset[1];
                tmpz -= offset[2];

                double posx_ingtrk = intcptx+slopex*tmpz;
                double posy_ingtrk = intcpty+slopey*tmpz;

                double tmpxy;
                double slopexy,intcptxy;
                if(view==0){
                  if(fabs(posx_ingtrk+offset[0])>350.){continue;}
                  tmpxy    = tmpy;
                  slopexy  = slopey;
                  intcptxy = intcpty;
                }else{
                  if(fabs(posy_ingtrk+offset[1])>350.){continue;}
                  tmpxy    = tmpx;
                  slopexy  = slopex;
                  intcptxy = intcptx;
                }
                // Distance from INGRID track
                double dist = fabs(intcptxy+slopexy*tmpz-tmpxy)/sqrt(1.+slopexy*slopexy);
                if(dist<dist_max){
                  channel[view][pln][ch] = true;
                }
              }
            }
          }
        }

        double slopex2=-1.,slopey2=-1.;
        double intcptx2=-1.,intcpty2=-1.;
        for(int diff_cyc=0;diff_cyc<=max_diff_cyc;diff_cyc++){
          int cyc = icyc+diff_cyc;
          for(int view=0;view<2;view++){
            int n2dreco = evtsum->NModTwoDimRecons(c_mod,cyc,view);
            if(n2dreco!=1){continue;}
            for(int irec=0;irec<n2dreco;irec++){
              TwoDimReconSummary* reconsum = evtsum->GetModTwoDimRecon(irec,c_mod,cyc,view);
              double endxy  = reconsum->endxy;
              double slope  = reconsum->slope;
              double endz   = reconsum->endz;
              double slope2 = reconsum->slope;
              double intcpt2= reconsum->intcpt;
              if(view==0){slopey2=slope2;intcpty2=intcpt2;}
              else       {slopex2=slope2;intcptx2=intcpt2;}

              double ang_diff,pos_diff;
              if(view==0){
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
              int nrechit = reconsum->Nhits();
              for(int ihit=0;ihit<nrechit;ihit++){
                HitSummary* hitsum = reconsum->GetHit(ihit);
                int rechitpln = hitsum->pln;
                int rechitch  = hitsum->ch ;
                channel_hit[view][rechitpln][rechitch] = true;

                int grid = 0;
                if(imod==0&&rechitch>=40){grid=1;}
                int recpln = rechitpln;
                if(grid==1){recpln = rechitch%20;}
                hitplane_hit[view][recpln][grid] = true;

              } //for(ihit)
            }
          }
        } //for(diff_cyc)


        intcptx += offset[0] - slopex*offset[2];
        intcpty += offset[1] - slopey*offset[2];
        if     (imod==0){maxpln= 8;maxch=80;}
        else if(imod==1){maxpln=18;maxch=32;}
        for(int view=0;view<2;view++){
          for(int pln=0;pln<maxpln;pln++){
            for(int ch=0;ch<maxch;ch++){
              if(imod==1&&pln==0&&ch>=24){continue;}
              double posx=0.,posy=0.;
              double posx2=0.,posy2=0.;
              int grid = 0;
              if(imod==0&&ch>=40){grid=1;}
              get_hitpos_in_scinti(intcptx ,slopex ,intcpty ,slopey ,c_mod,view,pln,ch,posx ,posy );
              get_hitpos_in_scinti(intcptx2,slopex2,intcpty2,slopey2,c_mod,view,pln,ch,posx2,posy2);

              if(channel[view][pln][ch]){ 
                h3[imod][grid]->Fill(posx,norm_tot);
                if(channel_hit[view][pln][ch]){ 
                  h4[imod][grid]->Fill(posx,norm_tot);
                }
              }

              int recpln = pln;
              if(grid==1){recpln = ch%20;}
              if(hitplane[view][recpln][grid]){
                h1[imod][grid]->Fill(posx2,norm_tot);
                if(hitplane_hit[view][recpln][grid]){
                  h2[imod][grid]->Fill(posx2,norm_tot);
                }
              }
            }
          }
        }

      } //for(imod)
    } //icyc
  } //ievt

  ofile->cd();
  for(int i=0;i<NMOD_TRKEFF;i++){
    for(int grid=0;grid<2;grid++){
      if(i==1&&grid==1){continue;}
      h1[i][grid]->Write();
      h2[i][grid]->Write();
      h3[i][grid]->Write();
      h4[i][grid]->Write();
    }
  }
  ofile->Write();
  ofile->Close();

  delete badch;
  delete detdim;
  return 0;
}
