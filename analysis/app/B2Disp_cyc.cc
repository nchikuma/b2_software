#include "TArc.h"
#include "TBox.h"
#include "TPolyLine.h"
#include "TLine.h"

#include "ThreeDimRecon.hxx"

#define DEBUG_DISP
//#define DEBUG_DETDIM


#define OFFSET_CANV 1000.

bool debug      = false;
double PEpara   = 3.;//6.  ;
double PEth     = 0.1 ;
double LineWid  = 0.05; 
double LineWid2 = 1.5 ; 


void sci_ing(double x, double y, double x_len, double y_len, double deg);
void sci_par(double x, double y, double x_len, double y_len, double deg);
void sci_sci(double x, double y, double x1, double y1);
void sci_veto(double x, double y, double x_len, double y_len, double deg_veto);
void iron(double x, double y, double x_len, double y_len, double deg);

void DrawINGRID(int view, double x_center, double y_center, double z_center, double deg);
void DrawProtonModule(int view, double x_center, double y_center, double z_center);
void DrawWaterModule(int view, double x_center, double y_center, double z_center);
void DrawHits(int mod, int pln, int view, int ch, double pe,int color=2);
void DrawPoints(int mod,int view,double xy,double z, int color=1,int offset=1);

void drawx(int targetmod);
void drawy(int targetmod);

void tline(double iX,double iY,double fX,double fY);
void drawdotline(double iX,double iY,double fX,double fY);
void drawtext(const char* text, double x, double y,double font);

void PrintUsage(){
  std::cout << "-m Select mode:\n" 
    << "   0 : MC true info\n"
    << "   1 : 2D reconstruction\n"
    << "   2 : 3D recon and vertex" << std::endl;
  std::cout << "-i Input ROOT file." << std::endl;
  std::cout << "-t Target module:\n"
	      << "   5: OnAxisWaterModule\n" 
	      << "   7: WaterModule\n" 
	      << "   8: CHModule" << std::endl;
  std::cout << "-v : Activate simulation parameters." << std::endl;
  std::cout << "-e : Set a start event number (default=0)." << std::endl;
  std::cout << "-s : Set an interaction type." << std::endl;
  std::cout << "-k : Activate scanning mode."  << std::endl;
  std::cout << "-c : Cosmic trigger mode."  << std::endl;
  std::cout << "-b : Open all bad channels (Default: masked)" <<std::endl;
  std::cout << "-d : Debug mode" <<std::endl;
  std::cout << "-a : Check event by event" <<std::endl;
}


int main(int argc, char** argv){

  //------- Arguments -------
  int c = -1;
  std::string readfilename;
  int mode        = -1;
  int type        = -1;
  int target      = -1;
  int evt_start   = 0; 
  bool selecttype = false;
  bool SCANNING   = false;
  bool simulation = false;
  bool BadCh_Ana  = true;
  bool Cosmic     = false;
  bool allevt     = false;

  while ((c = getopt(argc, argv, "i:m:t:s:e:cabkvdh")) != -1){
    switch(c){
    case 'i':
      readfilename = optarg;
      break;
    case 'm':
      mode = atoi(optarg);
      break;
    case 't':
      target = atoi(optarg);
      break;
    case 's':
      type = atoi(optarg);
      selecttype = true;
      break;
    case 'e':
      evt_start = atoi(optarg);
      break;
    case 'k':
      SCANNING = true;
      gROOT->SetBatch(true);
      break;
    case 'c':
      Cosmic = true;
      break;
    case 'b':
      BadCh_Ana = true;
      break;
    case 'v':
      simulation = true;
      break;
    case 'd':
      debug = true;
      break;
    case 'a':
      allevt = true;
      break;
    case 'h':
      PrintUsage();
      exit(0);
    default:
      PrintUsage();
      exit(0);
    }
  }


  //------- Open root file and get tree information -------
  TROOT root("GUI","GUI");
  TApplication theApp("App",0,0);

  if(
      (argc<3)||
      (mode==-1)||
      (target!=5 && target!=7 && target!=8) )
  {
    PrintUsage();
    exit(1);
  }

  TFile* readfile = new TFile(readfilename.c_str(),"read");
  if(readfile->IsZombie()){
    std::cout << "Cannot open file : " << readfilename << std::endl;
    exit(1);
  }
  std::cout << "Open file : " << readfilename << std::endl;
  TTree* tree = (TTree*)readfile ->Get("tree");
  TBranch* evtbr = tree->GetBranch("fDefaultReco.");
  EventSummary* evtsum = new EventSummary();
  evtbr -> SetAddress(&evtsum);
  tree  -> SetBranchAddress("fDefaultReco.",&evtsum);

  int nevt = (int)tree -> GetEntries();
  std::cout << "Total # of events = " << nevt << std::endl; 


  //------- Set up -------
  double HistLeft  = 0.;
  double HistRight = 0.;
  int    HistBin   = 2050*100;
  if(target==7 || target==8){
    HistLeft  = -C_B2MotherSizeZ+300;
    HistRight =  C_B2MotherSizeZ-1800;
  }
  else if(target==5){
    HistLeft  = C_PMMotherPosZ   - C_PMSizeZ;
    HistRight = C_INGHMotherPosZ + C_INGSizeZ;
  }
  int start_cyc, stop_cyc;
  if     (Cosmic    ){start_cyc=14;stop_cyc=15;}
  else if(simulation){start_cyc= 4;stop_cyc= 4;} 
  else               {start_cyc= 4;stop_cyc=12;}

  badch =  new INGRID_BadCh_mapping();
  badch -> set_BadCh(BadCh_Ana);
  detdim  = new DetectorDimension();
  bool newcanv = true;
  TCanvas *c1,*c2;

  //------- Event loop -------
  for(int ievt = evt_start;ievt<nevt;ievt++){

    std::cout << "Event # is " << ievt << std::endl;

    //Information text in event display
    std::string textevt, texthits, texttrk, textingtrk, textinttype, textspill, textcyc;

    tree  -> GetEntry(ievt);	
    int nhits = evtsum->NHits();
    int spill = evtsum->nd280nspill;

    std::cout << "  Hits # is " << nhits << std::endl;

    //SimvertexSummary
    SimVertexSummary* simvertexsum;
    int inttype;
    double xnu,ynu,znu;
    if(simulation){
      simvertexsum  = evtsum->GetSimVertex(0);
      inttype   = simvertexsum->inttype;
      if(selecttype==true&&inttype!=type) continue;
      std::cout << "  Neutrino energy is "  << simvertexsum->nuE  << " [GeV]" << std::endl;
      std::cout << "  Interaction Type # is " << inttype << std::endl;
      xnu = simvertexsum->xnu*10.;
      ynu = simvertexsum->ynu*10.;
      znu = simvertexsum->znu*10.;
    }
    for(int icyc=start_cyc;icyc<=stop_cyc;icyc++){
      cout << "icyc="<<icyc <<endl;

      //------- Create canvas and histgram -------
      if (newcanv){
        gROOT->SetStyle("Plain");
        double canvas_norm = 0.8;
        double canvas_x = 600;
        double canvas_y = 800;
        if(debug){
          canvas_x = 1900; 
          canvas_y = 1000;
          PEpara   = 5.;
        }
        if(target==5){ canvas_x*=0.7; }
        c1 = new TCanvas("c1","c1",0,0,canvas_x*canvas_norm,canvas_y*canvas_norm);
        if(debug){
          c2 = new TCanvas("c2","c2",0,0,canvas_x*canvas_norm,canvas_y*canvas_norm);
        }
        newcanv = false;			
      }

      if(!allevt){
          //if( target==5&&
          //    (evtsum->NModTwoDimRecons(15,icyc,0)<1&&evtsum->NModTwoDimRecons(15,icyc,1)<1)
          //    //&&(evtsum->NModTwoDimRecons(3,icyc,0)<1&&evtsum->NModTwoDimRecons(3,icyc,1)<1)
          //      )
          //    //(evtsum->NModHits(15,icyc)<3)
          //    )
          //{continue;}
          //if( target==7&&
          //    //(evtsum->NModTwoDimRecons(14,icyc,0)<1&&evtsum->NModTwoDimRecons(14,icyc,1)<1)&&
          //    //(evtsum->NModTwoDimRecons(21,icyc,0)<1&&evtsum->NModTwoDimRecons(21,icyc,1)<1)
          //    //)
          //    (evtsum->NModHits(16,icyc)<1)
          //    )
          //{continue;}
        //if(evtsum->NThreeDimRecons()<1) {continue;}
        if(
            (evtsum->NModTwoDimRecons(14,icyc,0)<1||evtsum->NModTwoDimRecons(14,icyc,1)<1)
            &&
            (evtsum->NModTwoDimRecons(14,icyc+1,0)<1||evtsum->NModTwoDimRecons(14,icyc+1,1)<1)
            &&
            (evtsum->NModTwoDimRecons(14,icyc-1,0)<1||evtsum->NModTwoDimRecons(14,icyc-1,1)<1)
            ){continue;}
        if(evtsum->NThreeDimRecons()<1){continue;}
        if(evtsum->GetThreeDimRecon(0)->Ntrack<2){continue;}
        if(evtsum->GetThreeDimRecon(0)->startmod[0]!=21){ continue; }
        double thetax = evtsum->GetThreeDimRecon(0)->thetax[0]*PI/180.;
        double thetay = evtsum->GetThreeDimRecon(0)->thetay[0]*PI/180.;
        double angle = fabs(atan(sqrt(pow(tan(thetax),2)+pow(tan(thetay),2))))*180./PI;
        if(angle>35.||angle<20.) { continue; }
        //if(mode==2&&!allevt){
        //  if(target==7 && evtsum->NThreeDimRecons()<1) {continue;}
        //  //if(evtsum->NThreeDimRecons()<2){continue;}
        //  ThreeDimReconSummary *vtx  = evtsum->GetThreeDimRecon(0);
        //  bool sel_fv      = true;
        //  //bool sel_ingstop = true;
        //  if(vtx->startmod[0]!=21       ){ sel_fv=false; }
        //  //if(vtx->startmod[0]!=16       ){ sel_fv=false; }
        //  if(vtx->startxpln[0]<=4       ){ sel_fv=false; }
        //  if(vtx->startypln[0]<=4       ){ sel_fv=false; }
        //  //if(vtx->startxpln[0]<=0       ){ sel_fv=false; }
        //  //if(vtx->startypln[0]<=0       ){ sel_fv=false; }
        //  if(fabs(vtx->startxch[0])>350.){ sel_fv=false; }
        //  if(fabs(vtx->startych[0])>350.){ sel_fv=false; }
        //  //if(vtx->stopmodx[0] != 14     ){ sel_ingstop=false; }
        //  //if(vtx->stopmody[0] != 14     ){ sel_ingstop=false; }
        //  //if(fabs(vtx->endxch[0])>500.  ){ sel_ingstop=false; }
        //  //if(fabs(vtx->endych[0])>500.  ){ sel_ingstop=false; }
        //  //if(vtx->endxpln[0]>=9         ){ sel_ingstop=false; }
        //  //if(vtx->endypln[0]>=9         ){ sel_ingstop=false; }
      }

      bool sim_fv      = false;
      bool sim_ingstop = false;
      float iposx, iposy, iposz, fposx, fposy, fposz; 
      float mu_ang = 90.;
      int nsimp = evtsum->NSimParticles();
      for(int isimp=0;isimp<nsimp;isimp++){
        SimParticleSummary   *simp = evtsum->GetSimParticle(isimp);
        if(abs(simp->pdg)!=13){continue;}
        iposx    = simp->ipos[0]*10.; //mm
        iposy    = simp->ipos[1]*10.; //mm
        iposz    = simp->ipos[2]*10.; //mm
        fposx    = simp->fpos[0]*10.; //mm
        fposy    = simp->fpos[1]*10.; //mm
        fposz    = simp->fpos[2]*10.; //mm
        if(target==7){
          iposx -= C_B2MotherPosX + C_B2WMPosX;
          iposy -= C_B2MotherPosY + C_B2WMPosY;
          iposz -= C_B2MotherPosZ + C_B2WMPosZ;
        }else if(target==8){
          iposx -= C_B2MotherPosX + C_B2CHPosX;
          iposy -= C_B2MotherPosY + C_B2CHPosY;
          iposz -= C_B2MotherPosZ + C_B2CHPosZ;
        }
        fposx -= C_B2MotherPosX + C_B2INGPosX;
        fposy -= C_B2MotherPosY + C_B2INGPosY;
        fposz -= C_B2MotherPosZ + C_B2INGPosZ;
        float dirX  = simp->momentum[0];
        float dirY  = simp->momentum[1];
        float dirZ  = simp->momentum[2];
        if(dirZ==0.    ){ mu_ang = 90.; }
        else            { mu_ang = atan( sqrt(dirX*dirX+dirY*dirY)/dirZ )*180./PI;}
        if(mu_ang < 0. ){ mu_ang = 180. + mu_ang; }
        if(mu_ang==180.){ mu_ang = 0.;}

        bool tmp = true;
        if(fabs(iposx)>350.       ){ tmp = false; }
        if(fabs(iposy)>350.       ){ tmp = false; }
        if(target==7){
          if(iposz<zposi(21,1, 4,0)+1.5){ tmp = false; }
          if(iposz>zposi(21,0,15,0)-1.5){ tmp = false; }
        }else if(target==8){
          if(iposz<zposi(16,1, 1,0)+6.5){ tmp = false; }
          if(iposz>zposi(16,1,15,0)-6.5){ tmp = false; }
        }
        sim_fv = tmp;
        tmp = true;
        if(fabs(fposx)>500.       ){ tmp = false; }
        if(fabs(fposy)>500.       ){ tmp = false; }
        if(fposz<zposi(14,0,2,0)  ){ tmp = false; }
        if(fposz>zposi(14,0,9,0)  ){ tmp = false; }
        sim_ingstop = tmp;
      }

      //  //if(!sim_fv || sel_fv){continue;}
      //  //if(!sel_fv   ){continue;}
      //  if(!sim_fv   ){continue;}
      //  //if(mu_ang>20.){continue;}

      //  ////cout << "ipos={" << iposx << ", "<< iposy << ", "<< iposz << "}" <<endl;
      //  ////cout << "recon v={" 
      //  ////  << vtx->startych [0] << ", "<< vtx->startxch [0] << ", "
      //  ////  << vtx->startypln[0] << ", "<< vtx->startxpln[0] << "}" <<endl;
      //}


      TH1F *h,*v;
      if(!debug){
        h = new TH1F("","Side View",HistBin,HistLeft,HistRight);
        v = new TH1F("","Top View" ,HistBin,HistLeft,HistRight);
      }
      else{
        h = new TH1F("","Side View",HistBin,HistLeft,HistRight);
        v = new TH1F("","Top View" ,HistBin,HistLeft,HistRight);

      }
      if(!debug){
        h->SetMinimum(-C_INGSizeY-100.+200+OFFSET_CANV);
        h->SetMaximum( C_INGSizeY     +200+OFFSET_CANV);
      }
      else{
        h->SetMinimum(-C_INGSizeY-100.+200+OFFSET_CANV);
        h->SetMaximum( C_INGSizeY     +200+OFFSET_CANV);
      }
      h->GetXaxis()->SetLabelSize(0);
      h->GetYaxis()->SetLabelSize(0);
      h->SetStats(0);
      v->SetMinimum(-C_INGSizeX+OFFSET_CANV);
      v->SetMaximum(+C_INGSizeX+OFFSET_CANV);
      v->GetXaxis()->SetLabelSize(0);
      v->GetYaxis()->SetLabelSize(0);
      v->SetStats(0);
      gStyle->SetOptTitle(0);
      gStyle->SetPadTopMargin(0.01);
      gStyle->SetPadLeftMargin(0.01);
      gStyle->SetPadRightMargin(0.01);
      gStyle->SetPadBottomMargin(0.01);

      //------- Create pad for Top/Side view -------
      TPad *pad1;
      TPad *pad2;
      if(!debug){
        pad1 = new TPad("pad1","pad1",0.01,0.02,0.99,0.41);
        pad2 = new TPad("pad2","pad2",0.01,0.44,0.99,0.83);
      }else{
        pad1 = new TPad("pad1","pad1",0.01,0.02,0.99,0.83);
        pad2 = new TPad("pad2","pad2",0.01,0.02,0.99,0.83);
      }

      c1->cd();
      pad1->Draw();
      if(debug){c2->cd();}
      pad2->Draw();

      pad1->cd();
      h->Draw("AH");
      std::cout << "drawx: Side view " << std::endl;
      drawx(target); //Draw module in sideview
      pad2->cd();
      v->Draw("AH");
      std::cout << "drawy: Top view " << std::endl;
      drawy(target); //Draw module in topview

      c1->cd();
      if(debug){
        drawtext("sideview",0.1,0.84,0.04);
      }else{
        drawtext("sideview",0.1,0.41,0.04);
      }
      if(debug){c2->cd();}
      drawtext("topview" ,0.1,0.84,0.04);

      //if(simulation){
      //  if(target==7){
      //    xnu -= C_B2MotherPosX;
      //    ynu -= C_B2MotherPosY;
      //    znu -= C_B2MotherPosZ;
      //    pad1->cd();
      //    DrawPoints(21,0,ynu,znu,1,0);
      //    pad2->cd();
      //    DrawPoints(21,1,xnu,znu,1,0);
      //  }
      //  else if(target==8){
      //    xnu -= C_B2MotherPosX;
      //    ynu -= C_B2MotherPosY;
      //    znu -= C_B2MotherPosZ;
      //    pad1->cd();
      //    DrawPoints(16,0,ynu,znu,1,0);
      //    pad2->cd();
      //    DrawPoints(16,1,xnu,znu,1,0);
      //  }
      //  else if(target==5){
      //    pad1->cd();
      //    DrawPoints(15,0,ynu,znu,1,0);
      //    pad2->cd();
      //    DrawPoints(15,1,xnu,znu,1,0);
      //  }
      //}


      //------- Draw hit -------
      bool n2drecon = false;
      bool n3drecon = false;
      for(int imod=0;imod<5;imod++){
        int c_mod = -1;
        if     (imod==0){c_mod=3 ;} //INGRID mod3
        else if(imod==1){c_mod=15;} //INGRID Water Module
        else if(imod==2){c_mod=21;} //WAGASCI
        else if(imod==3){c_mod=16;} //Proton Module
        else if(imod==4){c_mod=14;} //B2 INGRID
        int nmodhits = evtsum->NModHits(c_mod,icyc);
        cout << "mod=" << c_mod <<" nhits=" << nmodhits << endl;
        if(target==7||target==8){
          if(imod<2){ continue; }
        }
        else if(target==5){
          if(imod>1){ continue; }
        }


        for(int ihits =0; ihits< nmodhits; ihits++){
          HitSummary* hitsum = evtsum->GetModHit(ihits,c_mod,icyc);

          int hitmod  = hitsum->mod;
          int hitview = hitsum->view;
          int hitpln  = hitsum->pln;
          int hitch   = hitsum->ch;
          double hitpe   = hitsum->pe;

#ifdef DEBUG_DISP
          std::cout << "Hits "
            << "\t mod:"  << hitmod
            << "\t pln:"  << hitpln
            << "\t view:" << hitview
            << "\t ch:"   << hitch
            << "\t pe:"   << hitpe
            << std::endl;
#endif
          if(hitview==0) pad1->cd();
          else           pad2->cd();
          DrawHits(hitmod,hitpln,hitview,hitch,hitpe);
        }

        if(mode==0&&simulation){
          int nsimpar = evtsum->NSimParticles();
          for(int isimpar=0;isimpar<nsimpar;isimpar++){
            SimParticleSummary * simpar = evtsum->GetSimParticle(isimpar);
            double ix,iy,iz;
            double fx,fy,fz;
            double dirx,diry,dirz;
            ix = simpar->ipos[0]*10.;
            iy = simpar->ipos[1]*10.;
            iz = simpar->ipos[2]*10.;
            fz = simpar->fpos[2]*10.;
            if(target==5){
              ix -= C_PMMotherPosX;
              iy -= C_PMMotherPosY;
              iz -= C_PMMotherPosZ;
              fz -= C_PMMotherPosZ;
            }
            else if(target==7||target==8){
              ix -= C_B2MotherPosX;
              iy -= C_B2MotherPosY;
              iz -= C_B2MotherPosZ;
              fz -= C_B2MotherPosZ;
            }
            dirx = simpar->momentum[0];
            diry = simpar->momentum[1];
            dirz = simpar->momentum[2];
            if(dirz!=0.){
              dirx/=dirz;
              diry/=dirz;
              fx = ix + dirx*(fz-iz);
              fy = iy + diry*(fz-iz);
            }
            else{
              fx = simpar->fpos[0];
              fy = simpar->fpos[1];
            }
            std::cout
              << " simpar:"<< isimpar
              << " pdg:" << simpar->pdg
              << " ix=" << ix
              << " iy=" << iy
              << " iz=" << iz
              << " fx=" << fx
              << " fy=" << fy
              << " fz=" << fz
              << std::endl;
            pad1->cd();
            iy += OFFSET_CANV;
            fy += OFFSET_CANV;
            drawdotline(iz,iy,fz,fy);
            pad2->cd();
            ix += OFFSET_CANV;
            fx += OFFSET_CANV;
            drawdotline(iz,ix,fz,fx);
          }
        }


        //2D reconstructed track
        if(mode==1){
          for(int iview=0;iview<2;iview++){
            int ntwodimrecon = evtsum->NModTwoDimRecons(c_mod,icyc,iview);
            //if(ntwodimrecon>0) n2drecon = true;
            n2drecon = true;

            cout << "  TwoDim Recon Track # is " << ntwodimrecon <<  endl;

            //------- Draw track -------
            for(int irecon = 0; irecon < ntwodimrecon ; irecon++){
              TwoDimReconSummary* reconsum = evtsum->GetModTwoDimRecon(irecon,c_mod,icyc,iview);

              int   view     = iview;
              int   hitmod   = reconsum->hitmod;
              int   startpln = reconsum->startpln;
              int   endpln   = reconsum->endpln;
              float startz   = reconsum->startz;
              float endz     = reconsum->endz;
              float startxy  = reconsum->startxy;
              float endxy    = reconsum->endxy;
              float angle    = reconsum->angle;
              float intcpt   = reconsum->intcpt;
              float slope    = reconsum->slope;
              int   clstime  = reconsum->clstime;
              std::cout << "===== track" << irecon+1 << ":"
                << " view=" << view
                << " hitmod=" << hitmod 
                << " angle=" << angle 
                << " intcpt=" << intcpt 
                << " slope=" << slope 
                << " clstime=" << clstime
                << std::endl
                << " startpln=" << startpln 
                << " startz=" << startz 
                << " startxy=" << startxy 
                << std::endl
                << " endpln=" << endpln 
                << " endz=" << endz 
                << " endxy=" << endxy 
                << std::endl;

              double iX, iY, iZ, fX, fY, fZ; 
              iX = startxy; 
              iY = startxy;
              iZ = startz;
              fX = endxy;
              fY = endxy;
              fZ = endz;

              if(view==0){     
                cout << "       iY-intcpt-slope*iZ = " << iY- intcpt - slope*iZ << endl;
                cout << "       fY-intcpt-slope*fZ = " << fY- intcpt - slope*fZ << endl;
              }
              else{
                cout << "       iX-intcpt-slope*iZ = " << iX- intcpt - slope*iZ << endl;
                cout << "       fX-intcpt-slope*fZ = " << fX- intcpt - slope*fZ << endl;
              }

              double OFFSET[3] = {0., 0., 0.};
              if(hitmod==MOD_B2_WM){
                OFFSET[0] = C_B2WMPosX;
                OFFSET[1] = C_B2WMPosY;
                OFFSET[2] = C_B2WMPosZ;
              }
              else if(hitmod==MOD_B2_CH||hitmod==MOD_PM){
                OFFSET[0] = C_B2CHPosX;
                OFFSET[1] = C_B2CHPosY;
                OFFSET[2] = C_B2CHPosZ;
              }
              else if(hitmod==MOD_B2_INGRID){
                OFFSET[0] = C_B2INGPosX;
                OFFSET[1] = C_B2INGPosY;
                OFFSET[2] = C_B2INGPosZ;
              }
              else if(hitmod==MOD_ONAXIS_WM){
                OFFSET[0] = C_PMMotherPosX;
                OFFSET[1] = C_PMMotherPosY;
                OFFSET[2] = C_PMMotherPosZ;
              }

              iX += OFFSET[0];
              iY += OFFSET[1];
              iZ += OFFSET[2]; 

              fX += OFFSET[0];
              fY += OFFSET[1];
              fZ += OFFSET[2];

              //std::cout << "      iX:"<<iX << " fX:"<<fX;
              //std::cout << "      iY:"<<iY << " fY:"<<fY;
              //std::cout << "      iZ:"<<iZ << " fZ:"<<fZ;
              //std::cout << std::endl;

              //Draw track for reconstruction
              if(view==0){  //Side view
                pad1->cd();
                tline(iZ,iY,fZ,fY);
              }
              if(view==1){  //Top view
                pad2->cd();
                tline(iZ,iX,fZ,fX);
              }


              int nrechit = reconsum->Nhits();
              //cout << "   --- Recon hits (" << nrechit <<") :";
              if(irecon==0){
                for(int irechit=0;irechit<nrechit;irechit++){
                  HitSummary* rechit = reconsum->GetHit(irechit);
                  int    hitmod  = rechit->mod;
                  int    hitview = rechit->view;
                  int    hitpln  = rechit->pln;
                  int    hitch   = rechit->ch;
                  double hitpe   = rechit->pe;
                  //cout
                  //  <<"{"<<hitmod  
                  //  <<" "<<hitview 
                  //  <<" "<<hitpln  
                  //  <<" "<<hitch   
                  //  <<" "<<hitpe   
                  //  <<"} ";
                  //DrawHits(hitmod,hitpln,hitview,hitch,hitpe,4);
                }
              }
              //cout << endl;

            }

          }
        } // if (mode==1)
      } //for(imod)

        //------- Vertex loop -------
      if(mode==2){
        int nvertex = evtsum->NThreeDimRecons();
        std::cout << "  Vertex # is " << nvertex << std::endl;

        //for(int ivertex=0; ivertex<nvertex; ivertex++){
        for(int ivertex=nvertex-1; ivertex>=0; ivertex--){
          ThreeDimReconSummary* anasum = evtsum->GetThreeDimRecon(ivertex);
          if(anasum->hitcyc!=icyc){continue;}
          int ntrack    = anasum->Ntrack;
          int ningtrack = anasum->Ningtrack;
          std::cout << "  Track # is " << ntrack 
            << "(ing:" << ningtrack << ")" << std::endl;
          if(ntrack>0) n3drecon = true;

          texttrk = Form("# of Track: %d",ntrack);

          //if(ntrack!=2) goto ENDEVENT;
          c1->cd();

          //------- Draw track -------
          vector<bool>::iterator it_ingtrk = anasum->ing_trk.begin();
          //for(int itrack = 0; itrack < ntrack ; itrack++){		
          for(int itrack = ntrack-1; itrack >=0 ; itrack--){		
            int   ing_startmod  = anasum->ing_startmod[itrack];
            int   startmod      = anasum->startmod    [itrack];
            int   stopmodx      = anasum->stopmodx    [itrack];
            int   stopmody      = anasum->stopmody    [itrack];
            int   ing_endmod    = anasum->ing_endmod  [itrack];
            int   startxpln     = anasum->startxpln   [itrack];
            int   startypln     = anasum->startypln   [itrack];
            float startxch      = anasum->startxch    [itrack];
            float startych      = anasum->startych    [itrack];
            int   ingendpln     = anasum->ing_endpln  [itrack];
            int   endxpln       = anasum->endxpln     [itrack];
            int   endypln       = anasum->endypln     [itrack];
            float endxch        = anasum->endxch      [itrack];
            float endych        = anasum->endych      [itrack];
            float thetax        = anasum->thetax      [itrack]*PI/180.;
            float thetay        = anasum->thetay      [itrack]*PI/180.;
            int   oneview       = anasum->oneview     [itrack];

#ifdef DEBUG_DISP
            std::cout << "  / track"      << itrack+1 
              << " ing_startmod:" << ing_startmod 
              << " startmod:"     << startmod
              << " stopmodx:"     << stopmodx
              << " stopmody:"     << stopmody
              << " startxpln:"    << startxpln 
              << " startypln:"    << startypln 
              << " startxch:"     << startxch 
              << " startych:"     << startych 
              << " ing_endmod:"   << ing_endmod
              << " endxpln:"      << endxpln
              << " endypln:"      << endypln
              << " endxch:"       << endxch
              << " endych:"       << endych
              << std::endl
              << "   ing_trk:"    << *it_ingtrk
              << " oneview:"      << oneview
              << " thetax:"       << anasum->thetax[itrack] 
              << " thetay:"       << anasum->thetay[itrack] 
              << std::endl;
#endif

            double iX, iY, iZ0, iZ1;
            double fX, fY, fZ0, fZ1;
            iZ0 = zposi(startmod, 0, startxpln, 0); //side
            iZ1 = zposi(startmod, 1, startypln, 0); //top
            fZ0 = zposi(stopmodx, 0, endxpln  , 0); //side
            fZ1 = zposi(stopmody, 1, endypln  , 0); //top
            iY = startxch; //side
            iX = startych; //top
            fY = endxch;   //side
            fX = endych;   //top

            //Draw start/stop points
            for(int view=0;view<2;view++){
              if(view==0){  //Side view
                pad1->cd();
                DrawPoints(startmod,view,iY,iZ0);
                DrawPoints(stopmodx,view,fY,fZ0);
              }
              if(view==1){  //Top view
                pad2->cd();
                DrawPoints(startmod,view,iX,iZ1);
                DrawPoints(stopmody,view,fX,fZ1);
              }
            }


#ifdef DEBUG_DISP
            std::cout
              << "  iX :" << iX  << "  fX :" << fX  
              << "  iY :" << iY  << "  fY :" << fY  
              << "  iZ0:" << iZ0 << "  fZ0:" << fZ0 
              << "  iZ1:" << iZ1 << "  fZ1:" << fZ1
              << std::endl;
#endif

            double OFFSET1[3] = {0.,0.,0.};
            double OFFSET2[3] = {0.,0.,0.};
            double OFFSET3[3] = {0.,0.,0.};

            if(startmod==MOD_B2_WM){
              OFFSET1[0] = C_B2WMPosX;
              OFFSET1[1] = C_B2WMPosY;
              OFFSET1[2] = C_B2WMPosZ;
            }
            else if(startmod==MOD_B2_CH||startmod==MOD_PM){
              OFFSET1[0] = C_B2CHPosX;
              OFFSET1[1] = C_B2CHPosY;
              OFFSET1[2] = C_B2CHPosZ;
            }
            else if(startmod==MOD_B2_INGRID){
              OFFSET1[0] = C_B2INGPosX;
              OFFSET1[1] = C_B2INGPosY;
              OFFSET1[2] = C_B2INGPosZ;
            }
            else if(startmod==MOD_ONAXIS_WM){
              OFFSET1[0] = C_PMMotherPosX;
              OFFSET1[1] = C_PMMotherPosY;
              OFFSET1[2] = C_PMMotherPosZ;
            }
            else if(startmod==MOD_INGRID_C){
              OFFSET1[0] = 0;
              OFFSET1[1] = 0;
              OFFSET1[2] = 0;
            }
            if(stopmodx==MOD_B2_WM){
              OFFSET2[0] = C_B2WMPosX;
              OFFSET2[1] = C_B2WMPosY;
              OFFSET2[2] = C_B2WMPosZ;
            }
            else if(stopmodx==MOD_B2_CH||stopmodx==MOD_PM){
              OFFSET2[0] = C_B2CHPosX;
              OFFSET2[1] = C_B2CHPosY;
              OFFSET2[2] = C_B2CHPosZ; //+100.;
            }
            else if(stopmodx==MOD_B2_INGRID){
              OFFSET2[0] = C_B2INGPosX;
              OFFSET2[1] = C_B2INGPosY;
              OFFSET2[2] = C_B2INGPosZ;
            }
            else if(stopmodx==MOD_ONAXIS_WM){
              OFFSET2[0] = C_PMMotherPosX;
              OFFSET2[1] = C_PMMotherPosY;
              OFFSET2[2] = C_PMMotherPosZ;
            }
            else if(stopmodx==MOD_INGRID_C){
              OFFSET2[0] = 0;
              OFFSET2[1] = 0;
              OFFSET2[2] = 0;
            }

            if(stopmody==MOD_B2_WM){
              OFFSET3[0] = C_B2WMPosX;
              OFFSET3[1] = C_B2WMPosY;
              OFFSET3[2] = C_B2WMPosZ;
            }
            else if(stopmody==MOD_B2_CH||stopmody==MOD_PM){
              OFFSET3[0] = C_B2CHPosX;
              OFFSET3[1] = C_B2CHPosY;
              OFFSET3[2] = C_B2CHPosZ; //+100.;
            }
            else if(stopmody==MOD_B2_INGRID){
              OFFSET3[0] = C_B2INGPosX;
              OFFSET3[1] = C_B2INGPosY;
              OFFSET3[2] = C_B2INGPosZ;
            }
            else if(stopmody==MOD_ONAXIS_WM){
              OFFSET3[0] = C_PMMotherPosX;
              OFFSET3[1] = C_PMMotherPosY;
              OFFSET3[2] = C_PMMotherPosZ;
            }
            else if(stopmody==MOD_INGRID_C){
              OFFSET3[0] = 0;
              OFFSET3[1] = 0;
              OFFSET3[2] = 0;
            }

            //X
            iX  += OFFSET1[0];
            fX  += OFFSET2[0];
            //Y
            iY  += OFFSET1[1];
            fY  += OFFSET3[1];
            //Z0
            iZ0 += OFFSET1[2];
            fZ0 += OFFSET2[2];
            //Z1
            iZ1 += OFFSET1[2];
            fZ1 += OFFSET3[2];

            //Calculate fX, fY
            fX  =  iX + tan(thetay)*(fZ1-iZ1); //top
            fY  =  iY + tan(thetax)*(fZ0-iZ0); //side

            bool htrk = (oneview==0||oneview==1);
            bool vtrk = (oneview==0||oneview==2);

            for(int VIEW=0; VIEW<2; VIEW++){
              if(VIEW==0&&htrk){
                pad1->cd();
                tline(iZ0,iY,fZ0,fY);
              }
              if(VIEW==1&&vtrk){
                pad2->cd();
                tline(iZ1,iX,fZ1,fX);
              }
            }

            int nanahit = anasum->NhitTs(itrack);
            for(int ianahit=0;ianahit<nanahit;ianahit++){
              HitSummary* rechit = anasum->GetHitTrk(ianahit,itrack);
              int    hitview = rechit->view;
              int    hitmod  = rechit->mod;
              int    hitpln  = rechit->pln;
              int    hitch   = rechit->ch;
              double hitpe   = rechit->pe;
              //cout
              //  <<"{"<<hitmod  
              //  <<" "<<hitview 
              //  <<" "<<hitpln  
              //  <<" "<<hitch   
              //  <<" "<<hitpe   
              //  <<"} ";
              int color = 2;
              if(ivertex==0&&itrack==0){color=4;}
              else{color=6;}
              if(hitview==0){pad1->cd();}
              else          {pad2->cd();}
              DrawHits(hitmod,hitpln,hitview,hitch,hitpe,color);
              //cout << endl;
            }

            it_ingtrk++;
          }
        }
      }// if (mode==2)



      c1->cd();
      textevt   = Form("Event # is %d",ievt);
      texthits  = Form("# of Hits is %d",nhits);
      textspill = Form("spill=%d",spill);
      textcyc   = Form("cyc=%d",icyc);
      drawtext(textevt     .c_str(),0.05,0.95 ,0.04 );
      drawtext(texthits    .c_str(),0.50,0.95 ,0.03 );
      drawtext(textspill   .c_str(),0.50,0.90 ,0.03 );
      drawtext(textcyc     .c_str(),0.50,0.85 ,0.03 );
      drawtext(readfilename.c_str(),0.03,0.005,0.015);
      if(debug){
        c2->cd();
        drawtext(textevt     .c_str(),0.05,0.95 ,0.04 );
        drawtext(texthits    .c_str(),0.50,0.95 ,0.03 );
        drawtext(textspill   .c_str(),0.50,0.90 ,0.03 );
        drawtext(textcyc     .c_str(),0.50,0.85 ,0.03 );
        drawtext(readfilename.c_str(),0.03,0.005,0.015);
      }

      if(simulation){
        if     (fabs(inttype)==1)  textinttype = "CCQE";
        else if(fabs(inttype)==2)  textinttype = "MEC2p2h";
        else if(fabs(inttype)>=11&&
            fabs(inttype)<=13)     textinttype = "CC1Pi";
        else if(fabs(inttype)==16) textinttype = "CC_coh";
        else if(fabs(inttype)==21||
            fabs(inttype)==26)     textinttype = "CC_dis";
        else if(fabs(inttype)>=31&&
            fabs(inttype)<=35)     textinttype = "NC_1Pi";
        else if(fabs(inttype)==36) textinttype = "NC_coh";
        else if(fabs(inttype)==41||
            fabs(inttype)==46)     textinttype = "NC_dis";
        else if(fabs(inttype)==52) textinttype = "NC_elastic";
        else		           textinttype = Form("IntType:%d",inttype);
        std::cout << textinttype << std::endl;
        c1->cd(0);
        drawtext(textinttype.c_str(),0.33,0.95,0.04);
      }

      c1->cd(0);
      c1->Update();
      if(debug){
        c2->cd(0);
        c2->Update();
      }
      if(
          (SCANNING)
          //|| (mode==1&&!n2drecon)
          //|| (mode==2&&!n3drecon)
          )
      {
        std::cout << "========== Next. ==========\n" << std::endl;
        if(SCANNING){
          std::string filename;
          filename = Form("png/target%d_event%d_cyc%d_mode%d.png",target,ievt,icyc,mode);
          c1->Print(filename.c_str());
          if(debug){ 
            filename = Form("pnt/target%d_event%d_cyc%d_mode%d_2.png",target,ievt,icyc,mode);
            c2->Print(filename.c_str());
          }
        }
      }
      else{
        printf("  Type \' n\' to move to next event.\n");
        printf("  Type \' s\' to save the event display.\n");
        printf("  Type \' q\' to quit.\n");
        printf("  Type \' e\' to get a event you want.\n");
        printf("  Type \' m\' to change mode.\n");
        printf("  Type \' b\' to go back the last cycle.\n");
        printf("  Type any other key to go to the next event.\n");

        while(1){
          char ans[8];
          cin >> ans;
          if(*ans=='n'){
            std::cout << ">> Next." << std::endl;
            break;
          }
          else if(*ans=='s'){
            //std::stringstream filename;
            //filename << "event" << ievt << textinttype<<"_"<<mode<<".eps";
            //c1->Print(filename.str().c_str());
            std::string filename;
            filename = Form("png/target%d_event%d_cyc%d_mode%d.png",target,ievt,icyc,mode);
            c1->Print(filename.c_str());
            if(debug){ 
              filename = Form("png/target%d_event%d_cyc%d_mode%d_2.png",target,ievt,icyc,mode);
              c2->Print(filename.c_str());
            }
            //break;
          }
          else if(*ans=='q') exit(0);
          else if(*ans=='e'){
            int eventnumber;
            std::cout << "  Type the number of event :";
            cin >> eventnumber;
            ievt = eventnumber - 1;
            icyc = 4;
            goto ENDEVENT;
          }
          else if(*ans=='b'){
            if(icyc==4){
              if(simulation){
                icyc = 4;
                ievt = ievt - 2;
              }
              else{
                icyc = 11;
                ievt = ievt -1;
              }
            }else{
              icyc = icyc - 2;
            }
            break;
          }

          else if(*ans=='m'){
            std::cout << "  Type the mode :";
            cin >> mode;
            icyc = icyc - 1;
            //ievt = ievt - 1;
            break;
          }

          printf("  Type \' n\' to move to next event.\n");
          printf("  Type \' s\' to save the event display.\n");
          printf("  Type \' q\' to quit.\n");
          printf("  Type \' e\' to get a event you want.\n");
          printf("  Type \' m\' to change mode.\n");
          printf("  Type any other key to go to the next event.\n");

        }
      }
ENDCYCLE:
      cout << "===== End of a cycle ===== " <<endl;
      c1->Clear();
      if(debug){c2->Clear();}
    } //icyc
ENDEVENT:
    cout << "===== End of an event ===== " <<endl;
    c1->Clear();
    if(debug){c2->Clear();}
  } //ievt

  delete detdim;
  delete badch;
  return 0;
}


void sci_ing(double x,double y,double x_len,double y_len,double deg_ing=0.){
  std::complex<double> pos_tmp[5] = {
    std::complex<double>(-x_len/2.,-y_len/2.),
    std::complex<double>( x_len/2.,-y_len/2.),
    std::complex<double>( x_len/2., y_len/2.),
    std::complex<double>(-x_len/2., y_len/2.),
    std::complex<double>(-x_len/2.,-y_len/2.)
  };
  Double_t x_tmp[5], y_tmp[5];
  for(int i=0; i<5; i++){
    pos_tmp[i] *= exp(std::complex<double>(0.,-deg_ing));
    x_tmp[i] = x + pos_tmp[i].real();
    y_tmp[i] = y + pos_tmp[i].imag();
  }

  TPolyLine *pline = new TPolyLine(5,x_tmp,y_tmp);
  pline->SetLineColor(kGreen);
  pline->SetFillStyle(0);
  pline->SetLineWidth(LineWid);
  pline->Draw("SAME");
};

void sci_par(double x,double y,double x_len,double y_len,double deg_par=0.){
  std::complex<double> pos_tmp[5] = {
    std::complex<double>(-x_len/2.,-y_len/2.),
    std::complex<double>( x_len/2.,-y_len/2.),
    std::complex<double>( x_len/2., y_len/2.),
    std::complex<double>(-x_len/2., y_len/2.),
    std::complex<double>(-x_len/2.,-y_len/2.)
  };
  Double_t x_tmp[5], y_tmp[5];
  for(int i=0; i<5; i++){
    pos_tmp[i] *= exp(std::complex<double>(0.,-deg_par));
    x_tmp[i] = x + pos_tmp[i].real();
    y_tmp[i] = y + pos_tmp[i].imag();
  }

  TPolyLine *pline = new TPolyLine(5,x_tmp,y_tmp);
  pline->SetLineColor(kYellow);
  pline->SetFillStyle(0);
  pline->SetLineWidth(LineWid);
  pline->Draw("SAME");
};

void sci_sci(double x,double y,double x1,double y1){
  TBox *b1=new TBox(x,y,x1,y1);
  b1->SetLineColor(kGreen+2);
  b1->SetFillStyle(0);
  b1->SetLineWidth(LineWid);
  b1->Draw("SAME");
};

void sci_veto(double x,double y,double x_len,double y_len,double deg_veto=0.){

  std::complex<double> pos_tmp[5] = {
    std::complex<double>(-x_len/2.,-y_len/2.),
    std::complex<double>( x_len/2.,-y_len/2.),
    std::complex<double>( x_len/2., y_len/2.),
    std::complex<double>(-x_len/2., y_len/2.),
    std::complex<double>(-x_len/2.,-y_len/2.)
  };

  Double_t x_tmp[5], y_tmp[5];
  for(int i=0;i<5;i++){
    pos_tmp[i] *= exp(std::complex<double>(0.,-deg_veto));
    x_tmp[i] = x + pos_tmp[i].real();
    y_tmp[i] = y + pos_tmp[i].imag();
  }

  TPolyLine *pline = new TPolyLine(5,x_tmp,y_tmp);
  pline->SetLineColor(kBlue);
  pline->SetLineWidth(LineWid);
  pline->Draw("SAME");

};

void iron(double x,double y,double x_len,double y_len,double deg_iron=0.){
  std::complex<double> pos_tmp[5] = {
    std::complex<double>(-x_len/2.,-y_len/2.),
    std::complex<double>( x_len/2.,-y_len/2.),
    std::complex<double>( x_len/2., y_len/2.),
    std::complex<double>(-x_len/2., y_len/2.),
    std::complex<double>(-x_len/2.,-y_len/2.)
  };
  Double_t x_tmp[5], y_tmp[5];

  for(int i=0; i<5; i++){
    pos_tmp[i] *= exp(std::complex<double>(0.,-deg_iron));
    x_tmp[i] = x + pos_tmp[i].real();
    y_tmp[i] = y + pos_tmp[i].imag();
  }

  TPolyLine *pline = new TPolyLine(5,x_tmp,y_tmp);
  pline->SetFillColor(17);
  pline->SetLineWidth(LineWid);
  pline->Draw("f SAME");
};

void watertank(double x,double y,double x1,double y1){
  TBox *b1 = new TBox(x,y,x1,y1);
  b1->SetFillColor(kBlack);
  b1->SetFillStyle(0);
  b1->SetLineWidth(LineWid);
  b1->Draw("SAME");
};

void DrawProtonModule(int view, double x_center, double y_center, double z_center){

#ifdef DEBUG_DISP
  std::cout << "DrawProtonModule" << std::endl;
#endif

  int pln,ch,ch_num;
  double x,y,z,xy,x_tmp,z_tmp,z_length,xy_length; 

  for(int pln=0; pln<22; pln++){
    if     (pln==0) ch_num = 24;
    else if(pln<18) ch_num = 32;
    else if(pln<22) ch_num = 17;

    if(pln<18){
      detdim->GetPosPM(pln,1-view,0,&x,&y,&z);
      x = x + x_center;
      y = y + y_center;
      z = z + z_center;

      if(view==0) xy = y;
      else        xy = x;

      sci_par(z,xy,C_PMScintiThick,C_PMScintiLength);
    }


    for(ch=0; ch<ch_num; ch++){
      detdim->GetPosPM(pln,view,ch,&x,&y,&z);
      x = x + x_center;
      y = y + y_center;
      z = z + z_center;

      if(view==0) xy = y;
      else        xy = x;

      if(pln<18){
        if(pln==0 || ch<8 || ch>=24) sci_ing(z,xy,C_INGScintiThick,C_INGScintiWidth);
        else                         sci_ing(z,xy,C_PMScintiThick ,C_PMScintiWidth);
      }
      else if((view==0 && (pln==18 || pln==20))||
	      (view==1 && (pln==19 || pln==21)))
        sci_veto(z,xy,C_INGScintiWidth,C_INGScintiThick);
    }

  }
}

void DrawINGRID(int view, double x_center, double y_center, double z_center, double deg=0.){

#ifdef DEBUG_DISP
  std::cout << "DrawINGRID" << std::endl;
#endif

  int pln,ch,ch_num;
  double x,y,z,xy,x_tmp,z_tmp,z_length,xy_length,rotate; 
  rotate = deg*PI/180.;

  for(int pln=0; pln<15; pln++){
    if(pln<11)      ch_num = 24;
    else if(pln<15) ch_num = 22;

    if(pln<11){
      detdim->GetPosING(8,pln,1-view,0,&x,&y,&z);
      x_tmp = cos(rotate)*x - sin(rotate)*z; 
      z_tmp = sin(rotate)*x + cos(rotate)*z; 

      x = x_tmp + x_center;
      y = y     + y_center;
      z = z_tmp + z_center;

      if(view==0) xy = y;
      else        xy = x;
      sci_par(z,xy,C_INGScintiThick,C_INGScintiLength,rotate);

      if(pln<9){
        detdim->GetPosING(8,pln,1-view,0,&x,&y,&z);
        x_tmp = cos(rotate)*x - sin(rotate)*z; 
        z_tmp = sin(rotate)*x + cos(rotate)*z; 

        x_tmp = x_tmp - sin(rotate)*(C_INGIronStart - C_INGPlnStart);
        z_tmp = z_tmp + cos(rotate)*(C_INGIronStart - C_INGPlnStart);

        x = x_tmp + x_center;
        y = y     + y_center;
        z = z_tmp + z_center;	  

        if(view==0) xy = y;
        else        xy = x;

        iron(z,xy,C_INGIronThick,C_INGIronXY,rotate);
      }
    }

    for(ch=0;ch<ch_num;ch++){
      detdim->GetPosING(8,pln,view,ch,&x,&y,&z);
      x_tmp = cos(rotate)*x - sin(rotate)*z; 
      z_tmp = sin(rotate)*x + cos(rotate)*z; 

      x = x_tmp + x_center;
      y = y     + y_center;
      z = z_tmp + z_center;

      if(view==0) xy = y;
      else        xy = x;
      z_length  = cos(rotate)*C_INGScintiThick/2. - sin(rotate)*C_INGScintiWidth/2.;
      xy_length = sin(rotate)*C_INGScintiThick/2. + cos(rotate)*C_INGScintiWidth/2.;

      if(pln<11){
        sci_ing(z,xy,C_INGScintiThick*0.99,C_INGScintiWidth,rotate);
      }
      else if((view==0&&(pln==13||pln==14))||
	      (view==1&&(pln==11||pln==12))){
        sci_veto(z,xy,C_INGScintiWidth,C_INGScintiThick,rotate);
      }
    }

  }
}

void DrawWaterModule(int mod, int view, double x_center, double y_center, double z_center){

#ifdef DEBUG_DISP
  std::cout << "DrawWaterModule" << std::endl;
#endif

  int  pln, ch, ch_num;
  double x, y, z, xy; 
  double z1, z2, xy1, xy2; 

  z1 = z_center-C_WMWaterTargetSizeZ/2.;
  z2 = z_center+C_WMWaterTargetSizeZ/2.;

  if(view==0){
    xy1 = y_center-C_WMWaterTargetSizeY/2.;
    xy2 = y_center+C_WMWaterTargetSizeY/2.;
  }
  else{
    xy1 = x_center-C_WMWaterTargetSizeX/2.;
    xy2 = x_center+C_WMWaterTargetSizeX/2.;
  }
  watertank(z1,xy1,z2,xy2);

  for(int pln=0; pln<8; pln++){

    detdim->GetPosWM(mod,pln,1-view,0,0,&x,&y,&z);
    x = x + x_center;
    y = y + y_center;
    z = z + z_center;	  

    if(view==0) xy = y;
    else        xy = x;

    sci_par(z,xy,C_WMScintiThick,C_WMScintiLength,0.);

    ch_num = 40;
    for(ch=0; ch<ch_num; ch++){

      detdim->GetPosWM(mod,pln,view,ch,0,&x,&y,&z);
      x = x + x_center;
      y = y + y_center;
      z = z + z_center;

      if(view==0) xy = y;
      else        xy = x;

      sci_ing(z,xy,C_WMScintiThick,C_WMScintiWidth,0.);

      if(ch<20){
	detdim->GetPosWM(mod,pln,view,ch+40,1,&x,&y,&z);
        x = x + x_center;
        y = y + y_center;
        z = z + z_center;

        if(view==0) xy = y;
        else        xy = x;
        sci_ing(z,xy,C_WMScintiWidth,C_WMScintiThick,0.);

	detdim->GetPosWM(mod,pln,view,ch+60,2,&x,&y,&z);
        x = x + x_center;
        y = y + y_center;
        z = z + z_center;

        if(view==0) xy = y;
        else        xy = x;
        sci_ing(z,xy,C_WMScintiWidth,C_WMScintiThick,0.);
      }
    }

  }

}

void DrawHits(int mod, int pln, int view, int ch, double pe, int color){
  if(badch->is_BadCh(mod,view,pln,ch)){ return; }

  double OFFSET[3] = {0.,0.,0.};

  //Offset from mod#3
  if(mod==MOD_INGRID_C){
    OFFSET[0] = 0.;
    OFFSET[1] = 0.;
    OFFSET[2] = 0.;
  }
  else if(mod==MOD_ONAXIS_WM){
    OFFSET[0] = C_PMMotherPosX;
    OFFSET[1] = C_PMMotherPosY;
    OFFSET[2] = C_PMMotherPosZ;
  }
  else if(mod==MOD_B2_WM){
    OFFSET[0] = C_B2WMPosX;
    OFFSET[1] = C_B2WMPosY;
    OFFSET[2] = C_B2WMPosZ;
  }
  else if(mod==MOD_PM||mod==MOD_B2_CH){
    OFFSET[0] = C_B2CHPosX;
    OFFSET[1] = C_B2CHPosY;
    OFFSET[2] = C_B2CHPosZ;
  }
  else if(mod==MOD_B2_INGRID){
    OFFSET[0] = C_B2INGPosX;
    OFFSET[1] = C_B2INGPosY;
    OFFSET[2] = C_B2INGPosZ;
  }

  OFFSET[1-view] += OFFSET_CANV;

  double x_center = OFFSET[0];
  double y_center = OFFSET[1];
  double z_center = OFFSET[2];


  double X = 0., Y = 0., Z = 0., R = 0., XY = 0.;
  detdim->GetPosInMod(mod,pln,view,ch,&X,&Y,&Z);

  X = X + x_center;
  Y = Y + y_center;
  Z = Z + z_center;

  if(view==0) XY = Y;
  else        XY = X;

  if(pe<PEth) R = 0.;
  else        R = sqrt(pe-PEth)*PEpara;

  //if(debug){ R = PEpara; }

  TArc *arc = new TArc(Z,XY,R);
  arc->SetFillColor(color);
  arc->SetLineColor(color);
  arc->SetLineWidth(LineWid);
  arc->Draw("SAME");
}

void DrawPoints(int mod,int view,double xy,double z, int color, int offset){
  double OFFSET[3] = {0.,0.,0.};

  //Offset from mod#3
  if(mod==MOD_INGRID_C){
    OFFSET[0] = 0.;
    OFFSET[1] = 0.;
    OFFSET[2] = 0.;
  }
  else if(mod==MOD_ONAXIS_WM){
    OFFSET[0] = C_PMMotherPosX;
    OFFSET[1] = C_PMMotherPosY;
    OFFSET[2] = C_PMMotherPosZ;
  }
  else if(mod==MOD_B2_WM){
    OFFSET[0] = C_B2WMPosX;
    OFFSET[1] = C_B2WMPosY;
    OFFSET[2] = C_B2WMPosZ;
  }
  else if(mod==MOD_PM||mod==MOD_B2_CH){
    OFFSET[0] = C_B2CHPosX;
    OFFSET[1] = C_B2CHPosY;
    OFFSET[2] = C_B2CHPosZ;
  }
  else if(mod==MOD_B2_INGRID){
    OFFSET[0] = C_B2INGPosX;
    OFFSET[1] = C_B2INGPosY;
    OFFSET[2] = C_B2INGPosZ;
  }

  OFFSET[1-view] += OFFSET_CANV;

  double x_center = OFFSET[0];
  double y_center = OFFSET[1];
  double z_center = OFFSET[2];


  if(offset==1){
    if(view==0) xy += y_center;
    if(view==1) xy += x_center;
    z += z_center;
  }
  else{
    xy += OFFSET_CANV;
  }

  double R = 50.;

  TArc *arc = new TArc(z,xy,R);
  arc->SetFillColor(color);
  arc->SetFillStyle(3002);
  arc->SetLineColor(color);
  arc->SetLineWidth(LineWid);
  arc->Draw("SAME");
}

//Sideview
void drawx(int targetmod){
  int modules[3];
  int pln, ch, ch_num;
  double x, y, z;
  double offset_x, offset_y, offset_z;

  if(targetmod==7 || targetmod==8){
    modules[0] = MOD_B2_WM;
    modules[1] = MOD_B2_CH;
    modules[2] = MOD_B2_INGRID;
  }
  else if(targetmod==5){
    modules[0] = MOD_ONAXIS_WM;
    modules[1] = MOD_INGRID_C;
    modules[2] = -1;
  }

  for(int mod=0; mod<3; mod++){
    if(modules[mod]==MOD_INGRID_C){
      offset_x = C_INGHMotherPosX;
      offset_y = C_INGHMotherPosY + OFFSET_CANV;
      offset_z = C_INGHMotherPosZ;
      DrawINGRID(0,offset_x,offset_y,offset_z);
    }
    if(modules[mod]==MOD_ONAXIS_WM){
      offset_x = C_PMMotherPosX;
      offset_y = C_PMMotherPosY + OFFSET_CANV;
      offset_z = C_PMMotherPosZ;
      DrawWaterModule(MOD_ONAXIS_WM,0,offset_x,offset_y,offset_z);
    }
    if(modules[mod]==MOD_B2_WM){
      offset_x = C_B2WMPosX;
      offset_y = C_B2WMPosY + OFFSET_CANV;
      offset_z = C_B2WMPosZ;
      DrawWaterModule(MOD_B2_WM,0,offset_x,offset_y,offset_z);
    }
    if(modules[mod]==MOD_B2_CH){
      offset_x = C_B2CHPosX;
      offset_y = C_B2CHPosY + OFFSET_CANV;
      offset_z = C_B2CHPosZ;
      DrawProtonModule(0,offset_x,offset_y,offset_z);
    }
    if(modules[mod]==MOD_B2_INGRID){
      offset_x = C_B2INGPosX;
      offset_y = C_B2INGPosY + OFFSET_CANV;
      offset_z = C_B2INGPosZ;
      DrawINGRID(0,offset_x,offset_y,offset_z);
    }
  }
};

//Topview
void drawy(int targetmod){
  int modules[3];
  int pln, ch, ch_num;
  double x, y, z;
  double offset_x, offset_y, offset_z;

  if(targetmod==7 || targetmod==8){
    modules[0] = MOD_B2_WM;
    modules[1] = MOD_B2_CH;
    modules[2] = MOD_B2_INGRID;
  }
  else if(targetmod==5){
    modules[0] = MOD_ONAXIS_WM;
    modules[1] = MOD_INGRID_C;
    modules[2] = -1;
  }

  for(int mod=0; mod<3; mod++){
    if(modules[mod]==MOD_INGRID_C){
      offset_x = C_INGHMotherPosX + OFFSET_CANV;
      offset_y = C_INGHMotherPosY;
      offset_z = C_INGHMotherPosZ;
      DrawINGRID(1,offset_x,offset_y,offset_z);
    }
    if(modules[mod]==MOD_ONAXIS_WM){
      offset_x = C_PMMotherPosX + OFFSET_CANV;
      offset_y = C_PMMotherPosY;
      offset_z = C_PMMotherPosZ;
      DrawWaterModule(MOD_ONAXIS_WM,1,offset_x,offset_y,offset_z);
    }
    if(modules[mod]==MOD_B2_WM){
      offset_x = C_B2WMPosX + OFFSET_CANV;
      offset_y = C_B2WMPosY;
      offset_z = C_B2WMPosZ;
      DrawWaterModule(MOD_B2_WM,1,offset_x,offset_y,offset_z);
    }
    if(modules[mod]==MOD_B2_CH){
      offset_x = C_B2CHPosX + OFFSET_CANV;
      offset_y = C_B2CHPosY;
      offset_z = C_B2CHPosZ;
      DrawProtonModule(1,offset_x,offset_y,offset_z);
    }
    if(modules[mod]==MOD_B2_INGRID){
      offset_x = C_B2INGPosX + OFFSET_CANV;
      offset_y = C_B2INGPosY;
      offset_z = C_B2INGPosZ;
      DrawINGRID(1,offset_x,offset_y,offset_z);
    }
  }
};

void tline(double iX,double iY,double fX,double fY){
  iY += OFFSET_CANV;
  fY += OFFSET_CANV;
  TLine *l1 = new TLine(iX,iY,fX,fY);
  l1->SetLineWidth(LineWid2);
  l1->Draw("SAME");
};

void drawdotline(double iX,double iY,double fX,double fY){
  TLine *l1 = new TLine(iX,iY,fX,fY);
  l1->SetLineStyle(2);
  l1->SetLineColor(kBlue);
  l1->SetLineWidth(LineWid2);
  l1->Draw("SAME");
};

void drawtext(const char* text,double x, double y,double font){
  TText *t1 = new TText(x,y,text);
  t1->SetTextSize(font);  
  t1->Draw("SAME");
}
