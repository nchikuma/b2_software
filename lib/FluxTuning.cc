#include "FluxTuning.h"

FluxTuning::FluxTuning(int fdid, char* tunefile){
	FDID = fdid;
	fTFRatio = new TFile(tunefile);

	tune();
	std::cout << "#### FluxTuning ####"  << std::endl;
}

FluxTuning::~FluxTuning()
{
}

double FluxTuning::reweight(float Enu, int flav){

  double wratio;
  int intnue=-1;

  //ununiform binning
  for(int ibin=0;ibin<NBIN;ibin++){
    if(Enu<tune_xbins[ibin]){
	intnue = ibin;
	break;
    }
    if(ibin==NBIN-1){
	intnue = ibin;
	break;
    }
  }
  wratio = rwrationu[FDID][intnue][flav];
  return wratio;

}

/////////////////////////
void FluxTuning::tune(){
	
  TGaxis::SetMaxDigits(4);
  char histname[300];
  char histname2[300];
  for(int iflav=1;iflav<=4;iflav++){
    if     (iflav == Cnumu ) sprintf(Flavor,"numu" );
    else if(iflav == Cnumub) sprintf(Flavor,"numub");
    else if(iflav == Cnue  ) sprintf(Flavor,"nue"  ); 
    else if(iflav == Cnueb ) sprintf(Flavor,"nueb" );
    else{
	    std::cout << "Wrong flavor id: " << iflav << std::endl;
	    exit(1);
    }

    sprintf(histname,"nd%d_tune_%s",FDID,Flavor);
    sprintf(histname2,"nd%d_nom_%s",FDID,Flavor);
    trationu[FDID][iflav-1]  = (TH1D*)fTFRatio->Get(histname);
    tratio2nu[FDID][iflav-1] = (TH1D*)fTFRatio->Get(histname2);

    for(int j=0;j<NBIN;j++){
      	  tune_xbins[j]=trationu[FDID][iflav-1]->GetBinCenter(j+1) 
          	           + trationu[FDID][iflav-1]->GetBinWidth(j+1)/2.;
            
      	  if(tratio2nu[FDID][iflav-1]->GetBinContent(j+1)!=0)
            	  rwrationu[FDID][j][iflav-1] = (trationu[FDID][iflav-1]->GetBinContent(j+1))
          	                                /(tratio2nu[FDID][iflav-1]->GetBinContent(j+1));
      	  else
            	  rwrationu[FDID][j][iflav-1]=1;
    }
  }
}

ClassImp(FluxTuning)
