#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TSpectrum.h"
#include "TTree.h"
#include "TArc.h"
#include "TBox.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TText.h"
#include "TGaxis.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <stdio.h>

#include "FluxTuning.h"

//#include "DetectorDimension.hh"
#include "EVENTSUMMARY.h"
//#include "TwoDimReconSummary.h"
#include "Const.hh"

const int NumHist  = 24;
std::stringstream title, xtitle, ytitle;

void SetHistBins(int i, int* nbin, double* min, double* max);
void SetHistBins2(int i, int* nbin1, double* min1, double* max1, int* nbin2, double* min2, double* max2);
void PrintHistList();
void TitleSetting(int i, int nbin, double max, double min);
void TitleSetting2(int i);
void chToCoord(int mod, int pln, int view, int ch, double *x, double *y, double *z);

const int ccothertype[9] = {11,12,13,15,16,21,22,23,26};
const int nctype[15] = {31,32,33,34,36,38,39,41,42,43,44,45,46,51,52};

int main(int argc, char** argv){

	const int fdid = 2;
	char tuningfile[300]; 
	char tuningfile_posi[300]; 
	char tuningfile_nega[300]; 
	//sprintf(tuningfile_posi,"/gpfs/fs03/t2k/beam/work/nchikuma/data/jnubeam/tunefile/tune_data13.root");
	//sprintf(tuningfile_nega,"/gpfs/fs03/t2k/beam/work/nchikuma/data/jnubeam/tunefile/tune_data13.root");
	//sprintf(tuningfile_posi,"/gpfs/fs03/t2k/beam/work/nchikuma/data/jnubeam/tunefile/test_posi250_fluxtune.root");
	//sprintf(tuningfile_posi,"/gpfs/fs03/t2k/beam/work/nchikuma/data/jnubeam/tunefile/test_nega250_fluxtune.root");

	//*****arguments*******************//
	int c = -1;
	std::stringstream readfilename, outputfilename;
	int horn = -1;
	int module = -1;
	while ((c = getopt(argc, argv, "i:o:t:m:")) != -1){
		switch(c){
			case 'i':
				readfilename << optarg;
				break;
			case 'o':
				outputfilename << optarg;
				break;
			case 't':
				horn = atoi(optarg);
				break;
			case 'm':
				module = atoi(optarg);
				break;
		}
	}



	//  ==== open root file and get tree ==== //
	TROOT root("GUI","GUI");
	TApplication theApp("App",0,0);
	gROOT->SetStyle("Plain");

	if(argc<4){
		std::cout << "-i Input ROOT file." << std::endl;
		std::cout << "-o Output ROOT file." << std::endl;
		std::cout << "-t Horn mode. <1 or 2>" << std::endl;
		std::cout << "-m Module. <7:WaterModule, 8:CHModule, 9:B2BG,"
		          << " 10:B2ING>" << std::endl;
	
		exit(1);
	}
	
	//if(horn==1)      sprintf(tuningfile,tuningfile_posi);
	//else if(horn==2) sprintf(tuningfile,tuningfile_nega);
	//else{
	//	std::cout << "Select tuning mode:\n"
	//		  << "  1 : fdid=2, +250kA\n"
	//		  << "  2 : fdid=2, -250kA\n"
	//		  << std::endl;
	//	exit(1);
	//}
	//FluxTuning* fluxtune = new FluxTuning(fdid,tuningfile);


	PrintHistList();

	TH1F *h[NumHist];
	for(int i=0;i<NumHist;i++){
	       	h[i] = NULL;
		int nbin;
		double min, max;
		SetHistBins(i,&nbin,&min,&max);
		TitleSetting(i,nbin,max,min);

		h[i] = new TH1F(title.str().c_str(),title.str().c_str(),nbin,min,max);
		h[i]->GetXaxis()->SetTitle(xtitle.str().c_str());
		h[i]->GetYaxis()->SetTitle(ytitle.str().c_str());
	}

	TFile* readfile = new TFile(readfilename.str().c_str(),"read");
	if(readfile->IsZombie()){
	  std::cout << "Cannot open file : " << readfilename.str().c_str() << std::endl;
	  exit(1);
	}
	std::cout << "Open file : " << readfilename.str().c_str() << std::endl;
	TTree* tree = (TTree*)readfile ->Get("tree");
	TBranch* evtbr = tree->GetBranch("fDefaultReco.");
	Int_t nevt = (int)tree -> GetEntries();
	std::cout << "Total # of events = " << nevt << std::endl;
	
	//  ==== open output root file ==== //
	TFile* outfile = new TFile(outputfilename.str().c_str(),"recreate");


	// ===== Constants ========== //
	Float_t density_h2o = 1.00;      //water        //g/cm3
	Float_t density_ch  = 1.03;      //scintillator //g/cm3
	Float_t density_fe  = 7.874;     //iron         //g/cm3
	Float_t thickness_tank = C_WMWaterTargetSizeZ/10.;                          //Water Module //cm
	Float_t area_tank      = C_WMWaterTargetSizeX/10.*C_WMWaterTargetSizeY/10.; //Water Module //cm^2
	Float_t fiducial_size[3]  = {95.,95.,thickness_tank};
	Float_t CenterOfTarget[3] = {0.,0.,0.};
	Float_t CoT_WM[3] = {C_B2MotherPosX/10.+C_B2WMPosX/10.,C_B2MotherPosY/10.+C_B2WMPosY/10.,C_B2MotherPosZ/10.+C_B2WMPosZ/10.};
	Float_t CoT_CH[3] = {C_B2MotherPosX/10.+C_B2CHPosX/10.,C_B2MotherPosY/10.+C_B2CHPosY/10.,C_B2MotherPosZ/10.+C_B2CHPosZ/10.};

	Float_t density   = 0.;
	Float_t area      = 0.;
	Float_t thickness = 0.;
	
	density  = density_fe;
	area      = fiducial_size[0]*fiducial_size[1];
	thickness = fiducial_size[2]; 


	//  ***** event loop start ************  ///
	EventSummary*       evtsum         = new EventSummary();
	SimVertexSummary*   simvertexsum   = new SimVertexSummary();
	SimParticleSummary* simparticlesum = new SimParticleSummary();
	HitSummary*         hitsum         = new HitSummary();	
	SimHitSummary*      simhitsum      = new SimHitSummary();	
	for(int ievt = 0;ievt<nevt;ievt++){

	  if(ievt%1000==0) std::cout << ">> event " << ievt << std::endl;
	  //std::cout << ">> event " << ievt << std::endl;
  
	  // ===== variables ========== //
  	  Int_t nhits, nsimparticles, ntwodrecon, nthreedrecon, ingtrknum, trknum;
	  // --- Normalization ---
  	  Float_t Norm_tot, Flux_tot, rewight;
	  // --- SimVertex ---
	  Float_t norm, totcrsne, nuE, inttype, nutype, vertex[3];
	  // --- SimParticle ---
  	  Float_t simpar_pdg, iX, iY, iZ, fX, fY, fZ, mom[3];
	  Float_t momentum, angle, range, angle_muon, momentum_muon, momentum_proton, momentum_pion;
	  Float_t muon_stop[3], range_proton, range_pion;
	  bool    muon_simpar=false, proton_simpar=false, pion_simpar=false;
	  // --- Hit,SimHit ---
  	  Int_t   hitmod, hitview, hitpln, hitch, hitpdg;
	  Float_t hitpe, hitedep;
	  Int_t   NumModHit_muon[6][2]; //0:WM, 1:CH, 2:B2ING

	  
	  // ==== get eventsummary ==== //
          evtbr -> SetAddress(&evtsum);
	  tree  -> SetBranchAddress("fDefaultReco.",&evtsum);
	  tree  -> GetEntry(ievt);
	  	
	  // ==== get simvertexsummar ==== //
	  simvertexsum = evtsum->GetSimVertex(0);
	  norm     = simvertexsum->norm;
	  totcrsne = simvertexsum->totcrsne;
	  nuE      = simvertexsum->nuE;
	  inttype  = simvertexsum->inttype;
	  nutype   = simvertexsum->nutype;
	  vertex[0]= simvertexsum->xnu;
	  vertex[1]= simvertexsum->ynu;
	  vertex[2]= simvertexsum->znu;



	  // === Normalization factors ======= //
	  //reweight = fluxtune->reweight(nuE,nutype);
	  //Norm_tot = norm*totcrsne*1e-38*6.02*1e+23*density*thickness*reweight;
	  //Flux_tot = norm/area*reweight;
	  Norm_tot = norm*totcrsne*1e-38*6.02*1e+23*density*thickness;	
	                           // DEBUG rewight for data13/15 must be tuned later.
	  Flux_tot = norm/area;



	  // ==== get hitsummary ==== //
	  nhits = evtsum->NHits();
	  for(int j=0;j<6;j++){ for(int i=0;i<2;i++)  NumModHit_muon[j][i] = 0;}
	  for(int ihits=0;ihits<nhits;ihits++){
	    hitsum  = evtsum->GetHit(ihits);		

	    hitmod  = hitsum->mod;
	    hitview = hitsum->view;
	    hitpln  = hitsum->pln;
	    hitch   = hitsum->ch;
	    hitpe   = hitsum->pe;

	    simhitsum = hitsum->GetSimHit(0);		
	    hitedep   = simhitsum->edeposit;
	    hitpdg    = simhitsum->pdg;

	    if(fabs(hitpdg)==13 && (hitmod==21)) NumModHit_muon[0][hitview]++;
	    if(fabs(hitpdg)==13 && (hitmod==22)) NumModHit_muon[1][hitview]++;
	    if(fabs(hitpdg)==13 && (hitmod==14)) NumModHit_muon[2][hitview]++;
	    if(fabs(hitpdg)==13 && (hitmod==24)) NumModHit_muon[3][hitview]++;
	    if(fabs(hitpdg)==13 && (hitmod==25)) NumModHit_muon[4][hitview]++;
	    if(fabs(hitpdg)==13 && (hitmod==26)) NumModHit_muon[5][hitview]++;
	  }
	  
	  // ==== get simparticles ==== //
	  range_proton = 0.;
	  range_pion   = 0.;
	  nsimparticles = evtsum->NSimParticles();
	  for(int isimparticles = 0;isimparticles<nsimparticles;isimparticles++){
	     simparticlesum = evtsum->GetSimParticle(isimparticles);

	     simpar_pdg = simparticlesum->pdg;
	     iX         = simparticlesum->ipos[0];
	     iY         = simparticlesum->ipos[1];
	     iZ         = simparticlesum->ipos[2];
	     fX         = simparticlesum->fpos[0];
	     fY         = simparticlesum->fpos[1];
	     fZ         = simparticlesum->fpos[2];
	     mom[0]     = simparticlesum->momentum[0];
	     mom[1]     = simparticlesum->momentum[1];
	     mom[2]     = simparticlesum->momentum[2];

	     momentum = sqrt(mom[0]*mom[0]+mom[1]*mom[1]+mom[2]*mom[2]);
	     angle    = acos(mom[2]/momentum)*180./PI;
	     range    = sqrt((fX-iX)*(fX-iX)+(fY-iY)*(fY-iY)+(fZ-iZ)*(fZ-iZ));
	     if(abs(simpar_pdg)==13){
		     momentum_muon  = momentum;
		     angle_muon     = angle;
		     muon_stop[0]   = fX;
		     muon_stop[1]   = fY;
		     muon_stop[2]   = fZ;
		     muon_simpar    = true;
	     }
	     if(abs(simpar_pdg)==2212&&range>range_proton){
		     momentum_proton = momentum;
		     range_proton    = range;
		     proton_simpar   = true;
	     }
	     if(abs(simpar_pdg)==211&&range>range_pion){
		     momentum_pion = momentum;
		     range_pion    = range;
		     pion_simpar   = true;
	     }
	  }
	  
	  //If muon stops in FV if H2O or CH
	  bool StopInFV[2][3] = {{false,false,false},{false,false,false}};
	  bool muonstopWM = false;
	  bool muonstopCH = false;
	  if(muon_simpar){
	     for(int i=0;i<3;i++){
	             if(fabs(muon_stop[i]-CoT_WM[i])<fiducial_size[i]/2.) StopInFV[0][i]=true;
	             if(fabs(muon_stop[i]-CoT_CH[i])<fiducial_size[i]/2.) StopInFV[1][i]=true;
	     }
	     if(StopInFV[0][0]&&StopInFV[0][1]&&StopInFV[0][2]) muonstopWM = true;
	     if(StopInFV[1][0]&&StopInFV[1][1]&&StopInFV[1][2]) muonstopCH = true;
	  }

	  // =========== Fill Histogram =============== //

	  bool draw_wmmuhit    = false;
	  bool draw_chmuhit    = false;
	  bool draw_inghit     = false;
	  bool draw_ingmuhit = false;
	  if(NumModHit_muon[0][0]>=3 && NumModHit_muon[0][1]>=3) draw_wmmuhit    = true;
	  if(NumModHit_muon[1][0]>=3 && NumModHit_muon[1][1]>=3) draw_chmuhit    = true;
	  if(NumModHit_muon[2][0]>=3 && NumModHit_muon[2][1]>=3) draw_ingmuhit = true;
	  
	  if     (module==B2ING) draw_inghit = draw_ingmuhit;

	  for(int p2=0;p2<4;p2++){
		  bool draw1 = false;
		  if     (p2==0){if(draw_wmmuhit)                draw1 = true;}
		  else if(p2==1){if(draw_chmuhit)                draw1 = true;}
		  else if(p2==2){if(draw_wmmuhit && muonstopWM)  draw1 = true;}
		  else if(p2==3){if(draw_chmuhit && muonstopCH)  draw1 = true;}
		  
		  if(draw1){	
			  if(muon_simpar){
				  h[p2]  ->Fill(nuE          ,Norm_tot); 
				  h[p2+4]->Fill(momentum_muon,Norm_tot); 
				  h[p2+8]->Fill(angle_muon   ,Norm_tot);
				  if(draw_inghit){ 
					  h[p2+12]->Fill(nuE          ,Norm_tot);
					  h[p2+16]->Fill(momentum_muon,Norm_tot);
					  h[p2+20]->Fill(angle_muon   ,Norm_tot);
				  }
			  }
			  
		  }
	  }

	} // Event Loop end
	if(evtsum) delete evtsum;

	for(int i=0;i<NumHist;i++){
		if(h[i]!=NULL) h[i]->Write();
	}

	outfile->Close();
	if(outfile) delete outfile;


	std::cout << "Done!!!" << std::endl;
	return 0;
}



void TitleSetting(int i, int nbin, double max, double min){

   title.str("");
   xtitle.str("");
   ytitle.str("");

   if(i>=0&&i<=3){
	title << "AllInteraction_nuE";
     	if     (i==0 )title<<"_BGWM";
     	else if(i==1 )title<<"_BGCH";
     	else if(i==2 )title<<"_BGWMstop";
     	else if(i==3 )title<<"_BGCHstop";
     	xtitle << "Neutrino energy [GeV]";
	ytitle << "Events [/10^{21}POT/"<<(int)((max-min)/nbin*1000) <<"MeV]";
   }
   else if(i>=4&&i<=7){
	title << "AllInteraction_muP";
     	if     (i==4 )title<<"_BGWM";
     	else if(i==5 )title<<"_BGCH";
     	else if(i==6 )title<<"_BGWMstop";
     	else if(i==7 )title<<"_BGCHstop";
     	xtitle << "Muon momentum [GeV/c]";
	ytitle << "Events [/10^{21}POT/"<<(int)((max-min)/nbin*1000) <<"MeV]";
   }
   else if(i>=8&&i<=11){
	title << "AllInteraction_muAng";
     	if     (i==8 )title<<"_BGWM";
     	else if(i==9 )title<<"_BGCH";
     	else if(i==10)title<<"_BGWMstop";
     	else if(i==11)title<<"_BGCHstop";
     	xtitle << "Muon angle to beam [deg]";
	ytitle << "Events [/10^{21}POT/"<<(int)((max-min)/nbin) <<"deg]";
   }
   else if(i>=12&&i<=15){
	title << "IngMuon3Hits_nuE";
     	if     (i==12)title<<"_BGWM";
     	else if(i==13)title<<"_BGCH";
     	else if(i==14)title<<"_BGWMstop";
     	else if(i=15)title<<"_BGCHstop";
     	xtitle << "Neutrino energy [GeV]";
	ytitle << "Events [/10^{21}POT/"<<(int)((max-min)/nbin*1000) <<"MeV]";
   }
   else if(i>=16&&i<=19){
	title << "IngMuon3Hits_muP";
     	if     (i==16)title<<"_BGWM";
     	else if(i==17)title<<"_BGCH";
     	else if(i==18)title<<"_BGWMstop";
     	else if(i==19)title<<"_BGCHstop";
     	xtitle << "Muon momentum [GeV/c]";
	ytitle << "Events [/10^{21}POT/"<<(int)((max-min)/nbin*1000) <<"MeV]";
   }
   else if(i>=20&&i<=23){
	title << "IngMuon3Hits_muAng";
     	if     (i==20)title<<"_BGWM";
     	else if(i==21)title<<"_BGCH";
     	else if(i==22)title<<"_BGWMstop";
     	else if(i==23)title<<"_BGCHstop";
     	xtitle << "Muon angle to beam [deg]";
	ytitle << "Events [/10^{21}POT/"<<(int)((max-min)/nbin) <<"deg]";
   }

}

void PrintHistList(){
	std::cout << "Histograms:" 
			  << "\n  == 3 hits/more in WM, 3 hits/more in CH, stop in WM, stop in CH== "
			  << "\n   0  - 3    : All Interactions, neutrino energy" 
			  << "\n   4  - 7    : All Interactions, muon momentum" 
			  << "\n   8  - 11   : All Interactions, muon angle" 
			  << "\n   12 - 15   : Requiring 3 hits or more of Muon in INGRID, neutrino energy" 
			  << "\n   16 - 19   : requiring 3 hits or more of muon in ingrid, muon momentum" 
			  << "\n   20 - 23   : requiring 3 hits or more of muon in ingrid, muon angle" 
			  << std::endl;

}

void SetHistBins(int i, int* nbin, double* min, double* max){

	if (i>=0  &&i<=3  ) { *nbin=200; *min=0.; *max=10. ;} //GeV       //Flux & nuE for each Int
	else if (i>=4  &&i<=7  ) { *nbin=200; *min=0.; *max=5.  ;} //GeV       //Momentum
	else if (i>=8  &&i<=11 ) { *nbin=180; *min=0.; *max=180.;} //deg       //angle
	else if (i>=12 &&i<=15 ) { *nbin=200; *min=0.; *max=10. ;} //GeV       //Flux & nuE for each Int
	else if (i>=16 &&i<=19 ) { *nbin=200; *min=0.; *max=5.  ;} //GeV       //Momentum
	else if (i>=20 &&i<=23 ) { *nbin=180; *min=0.; *max=180.;} //deg       //angle



	else { *nbin=100; *min=0.; *max=50.;}
}
