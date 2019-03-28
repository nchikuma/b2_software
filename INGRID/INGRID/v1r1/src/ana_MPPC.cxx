#ifndef _ANA_MPPC_C
#define _ANA_MPPC_C

#include "ana_MPPC.hxx"

ana_MPPC::ana_MPPC(){
  without_sigma=2.5;
}


Bool_t ana_MPPC::analysis_pedestal(TH1F *noise_hist){
  TF1 *gaus = new TF1("gaus","gaus",HIST_MIN,HIST_MAX);
  //get maximum bin X position
  fBinWidth         =  noise_hist  -> GetBinWidth    (150);
  fBinMax           =  noise_hist  -> GetMaximumBin  ();
  fBinEdge          =  noise_hist  -> GetBinLowEdge  (1);
  fXMax             =  fBinEdge + fBinWidth * fBinMax;
  pedestal_peak_pos =  fXMax;
  //fit around pedestal_peak_pos and get sigma and peak number
  gaus->SetParLimits(0,0,100000000);
  gaus->SetParLimits(1,0,300);
  gaus->SetParLimits(2,0,10);
  noise_hist->Fit("gaus","qn","",(pedestal_peak_pos-FITRANGE),(pedestal_peak_pos+FITRANGE));
  pedestal_peak_pos    = gaus->GetParameter(1);
  pedestal_peak_height = gaus->GetParameter(0);
  pedestal_peak_sigma  = gaus->GetParameter(2);
  gaus  -> Clear();
  delete gaus;
  return true; 
}

Bool_t ana_MPPC::analysis_again(TH1F *noise_hist, bool draw = false){
  TF1*        dgaus = new TF1("dgaus","gaus(0)+gaus(3)");
  TSpectrum*  peak1 = new TSpectrum();
  //Int_t      numpeak= 2;
  number_of_entries = noise_hist->GetEntries();
  //double       par[6];


  bool flag=false;
  //peak1->Search(noise_hist,2,"goff",0.01);
  char dop[300];
  sprintf(dop,"goff");
  if(draw){memset(dop,'\0',sizeof(dop));}

  peak1->Search(noise_hist,2,dop,0.01);
  if(peak1->GetNPeaks()<2){
    peak1->Search(noise_hist,3,dop,0.01);
    if(peak1->GetNPeaks()<2){
      peak1->Search(noise_hist,2,dop,0.005);
      if(peak1->GetNPeaks()<2){
	peak1->Search(noise_hist,3,dop,0.005);
	if( peak1->GetNPeaks()<=1 ){
	  onepe_peak_pos = *(peak1->GetPositionX()) + 7;
	  flag = true;
	}
      }
    }
  }

  pedestal_peak_pos = *(peak1->GetPositionX());
  if(!flag){
    onepe_peak_pos    = *(peak1->GetPositionX()+1);
    if(  onepe_peak_pos - pedestal_peak_pos > 25 &&
	 peak1->GetNPeaks()>=3 ){
      cout << "ereplace" << endl; 
     onepe_peak_pos =  *(peak1->GetPositionX()+2);
    }
  }
  if( draw )
    cout << "ped:" << pedestal_peak_pos << "\t" << onepe_peak_pos << endl;

  if(pedestal_peak_pos>onepe_peak_pos){
    Double_t aaa = pedestal_peak_pos;
    pedestal_peak_pos = onepe_peak_pos;
    onepe_peak_pos=aaa;
  }

  dgaus -> SetParameter(0, 1000);
  dgaus -> SetParameter(1, pedestal_peak_pos );
  dgaus -> SetParameter(2, 2 );
  dgaus -> SetParameter(3, 100);
  dgaus -> SetParameter(4, onepe_peak_pos );
  dgaus -> SetParameter(5, 2 );
  dgaus -> SetParLimits( 2, 1, 5 );//sigma
  dgaus -> SetParLimits( 5, 1, 5 );
  dgaus -> SetParLimits( 5, pedestal_peak_pos + 2, onepe_peak_pos + 10  );
  
  sprintf(dop,"qn0");
  if(draw){memset(dop,'\0',sizeof(dop));}

  noise_hist->Fit("dgaus",dop,"",pedestal_peak_pos-FITRANGE,onepe_peak_pos+FITRANGE);
  float temp = dgaus->GetParameter(4) - dgaus->GetParameter(1); 
  if( temp < 3 || temp > 28 ||
      fabs( dgaus->GetParameter(1) - pedestal_peak_pos )>3 ||
      fabs( dgaus->GetParameter(4) - onepe_peak_pos    )>3 
      ){
    for(int ipar=0; ipar<6; ipar++){
      dgaus -> SetParLimits( 0, 10, 1000 );//heigt
      dgaus -> SetParLimits( 3, 10, 100 );
      dgaus -> SetParLimits( 2, 1, 5 );//sigma
      dgaus -> SetParLimits( 5, 1, 5 );
      dgaus -> SetParLimits( 1, 0, 300 );//position
      dgaus -> SetParLimits( 4, 0, 300 );
 
    }
    noise_hist->Fit("dgaus",dop,"",pedestal_peak_pos-FITRANGE,onepe_peak_pos+FITRANGE);
    temp = dgaus->GetParameter(4) - dgaus->GetParameter(1); 
  
    if( temp < 3 || temp > 30 ){
      dgaus -> SetParLimits( 1, pedestal_peak_pos-3, pedestal_peak_pos+3 );//position
      dgaus -> SetParLimits( 4, onepe_peak_pos-3, onepe_peak_pos+3 );

      noise_hist->Fit("dgaus",dop,"",pedestal_peak_pos-4,onepe_peak_pos+4);
 
    }
  
  }

  pedestal_peak_pos     = dgaus->GetParameter(1);
  pedestal_peak_height  = dgaus->GetParameter(0);
  pedestal_peak_sigma   = dgaus->GetParameter(2);
  onepe_peak_pos        = dgaus->GetParameter(4);
  onepe_peak_height     = dgaus->GetParameter(3);
  onepe_peak_sigma      = dgaus->GetParameter(5);
  fChisquare            = dgaus->GetChisquare();
  fNDF                  = dgaus->GetNDF();
  gain = onepe_peak_pos - pedestal_peak_pos;
  if(gain<0)gain = -1.0 * gain;

  dgaus -> Clear();
  peak1 -> Clear();

  delete dgaus;
  delete peak1;

  if(pedestal_peak_height<onepe_peak_height){
    return false;
  }
  if(pedestal_peak_sigma>10||onepe_peak_sigma>10){
    return false;
  }
  //if(fChisquare/fNDF>15){
  //return false;
  //}

  return true;


}


Bool_t ana_MPPC::analysis_old_version(TH1F *noise_hist, int mod, bool draw = false){
  int h_min,h_max,h_fitrange;
  if(mod==15){ h_min=HIST_MIN_WM; h_max=HIST_MAX_WM; h_fitrange=FITRANGE_WM; }
  else       { h_min=HIST_MIN;    h_max=HIST_MAX;    h_fitrange=FITRANGE;    }

  TF1*         gaus = new TF1("gaus","gaus",h_min,h_max);
  TF1*        dgaus = new TF1("dgaus","gaus(0)+gaus(3)",h_min,h_max);
  TSpectrum*  peak1 = new TSpectrum();
  number_of_entries = noise_hist->GetEntries();
  double       par[6];

  bool flag=false;
  char dop2[300];
  sprintf(dop2,"goff");
  if(draw){memset(dop2,'\0',sizeof(dop2));}

  if(mod!=15){ peak1->Search(noise_hist,2  ,dop2,0.01        );}
  else       { peak1->Search(noise_hist,1.4,dop2,100./100000.);}

  if(peak1->GetNPeaks()<2){
    if(mod!=15){ peak1->Search(noise_hist,3  ,dop2,0.01        );}
    else       { peak1->Search(noise_hist,1.5,dop2,100./100000.);}

    if(peak1->GetNPeaks()<2){
      return false;
      //if(mod!=15){ peak1->Search(noise_hist,2,dop2,0.005); }
      //if(peak1->GetNPeaks()<2){
      //  if(mod!=15){ peak1->Search(noise_hist,3,dop2,0.005); }
      //  if( peak1->GetNPeaks()==1 ){
      //    //std::cout << "###can not see 2 peaks by TSpectrum###" << std::endl;
      //    onepe_peak_pos = *(peak1->GetPositionX()) + 7;
      //    flag = true;
      //  }
      //  else if( peak1->GetNPeaks()==0 ){
      //    //std::cout << "###can not see 1 peak by TSpectrum###" << std::endl;
      //    return false;
      //  }
      //}
    }
  }

  pedestal_peak_pos = *(peak1->GetPositionX());
  if(!flag){ onepe_peak_pos = *(peak1->GetPositionX()+1); }
  
  if(pedestal_peak_pos>onepe_peak_pos){
    Double_t aaa = pedestal_peak_pos;
    pedestal_peak_pos = onepe_peak_pos;
    onepe_peak_pos=aaa;
  }
   
  char dop[300];
  sprintf(dop,"qn");
  if(draw){memset(dop,'\0',sizeof(dop));}
 
  noise_hist->Fit("gaus",dop,"",
      (pedestal_peak_pos-h_fitrange),
      (pedestal_peak_pos+h_fitrange));
  gaus->GetParameters(par);

  noise_hist->Fit("gaus",dop,"",
      (onepe_peak_pos-h_fitrange),
      (onepe_peak_pos+h_fitrange));
  gaus ->GetParameters(par+3);
  gaus ->Clear();

  // if the 1pe pos stays in the range
  if( onepe_peak_pos-h_fitrange > par[4] || onepe_peak_pos+h_fitrange < par[4] ){
    if(mod==15){par[3] = par[0] * 0.01;} //adjust the parameter for 1pe height
    else       {par[3] = par[0] * 0.1 ;}
    par[4] = par[1] + 8;                 //for center pos 
    par[5] = par[2];                     //for sigma
  }

  dgaus -> SetParameters(par);
  dgaus -> SetParLimits(0, par[0]*0.1, par[0]*10  );
  dgaus -> SetParLimits(1, par[1]-10 , par[1]+10  );
  dgaus -> SetParLimits(2, par[2]*0.5, par[2]*1.5 );
  dgaus -> SetParLimits(3, par[3]*0.1, par[3]*10  );
  dgaus -> SetParLimits(4, par[4]-10 , par[4]+10  );
  dgaus -> SetParLimits(5, par[5]*0.5, par[5]*1.5 );


  noise_hist->Fit("dgaus",dop,"",pedestal_peak_pos-h_fitrange,onepe_peak_pos+h_fitrange);
  
  float temp = dgaus->GetParameter(4) - dgaus->GetParameter(1); 
  if( temp < 3 || temp > 22 ){
    //std::cout << "reset limit and fit again" << std::endl;
    for(int ipar=0; ipar<6; ipar++){
      dgaus -> SetParLimits( 0, 0, 10000000 );//heigt
      dgaus -> SetParLimits( 3, 0, 10000000 );
      dgaus -> SetParLimits( 2, 0, 100 );//sigma
      dgaus -> SetParLimits( 5, 0, 100 );
      dgaus -> SetParLimits( 1, 0, 300 );//position
      dgaus -> SetParLimits( 4, 0, 300 );
    }
    noise_hist->Fit("dgaus",dop,"",pedestal_peak_pos-h_fitrange,onepe_peak_pos+h_fitrange);
  }

  pedestal_peak_height  = dgaus->GetParameter(0);
  pedestal_peak_pos     = dgaus->GetParameter(1);
  pedestal_peak_sigma   = dgaus->GetParameter(2);
  onepe_peak_height     = dgaus->GetParameter(3);
  onepe_peak_pos        = dgaus->GetParameter(4);
  onepe_peak_sigma      = dgaus->GetParameter(5);
  fChisquare            = dgaus->GetChisquare();
  fNDF                  = dgaus->GetNDF();

  gain = onepe_peak_pos - pedestal_peak_pos;
  if(gain<0)gain = -1.0 * gain;

  dgaus -> Clear();
  peak1 -> Clear();
  delete gaus;
  delete dgaus;
  delete peak1;

  if(pedestal_peak_height<onepe_peak_height     ){ return false; }
  if(pedestal_peak_sigma>10||onepe_peak_sigma>10){ return false; }

  return true;
}

void ana_MPPC::analysis_no_mppc(TH1F *pedestal_hist){
  TF1 *gaus = new TF1("gaus","gaus",HIST_MIN,HIST_MAX);
  TSpectrum *peak1 = new TSpectrum();
  Double_t par[3];
  peak1->Search(pedestal_hist,4,"goff",0.05);
  pedestal_peak_pos = *(peak1->GetPositionX());
  //debug peak serch
  if(pedestal_peak_pos>MAXPEDESTAL||pedestal_peak_pos<MINPEDESTAL||onepe_peak_pos<MINONEPE||onepe_peak_pos>MAXONEPE){
    peak1->Search(pedestal_hist,3,"goff",0.05);
    pedestal_peak_pos = *(peak1->GetPositionX());
    onepe_peak_pos = *(peak1->GetPositionX()+1);
  }
  pedestal_hist->Fit("gaus","qn","",(pedestal_peak_pos-FITRANGE*2),(pedestal_peak_pos+FITRANGE*2));
  gaus->GetParameters(par);
  pedestal_peak_pos = gaus->GetParameter(1);
  pedestal_peak_height = gaus->GetParameter(0);
  pedestal_peak_sigma = gaus->GetParameter(2);
  onepe_peak_pos = 0;
  onepe_peak_height = 0;
  onepe_peak_sigma = 0;
  gain=0;
}



Double_t ana_MPPC::get_gain(){
  return gain;
}
Double_t ana_MPPC::get_pedestal(){
  return pedestal_peak_pos;
}
Double_t ana_MPPC::get_pedestal_sigma(){
  return pedestal_peak_sigma;
}
Double_t ana_MPPC::get_pedestal_height(){
  return pedestal_peak_height;
}
Double_t ana_MPPC::get_onepe(){
  return onepe_peak_pos;
}

Double_t ana_MPPC::get_noise(){
  TF1 *gaus = new TF1("gaus","gaus",0,1000);
  gaus -> SetParameters ( pedestal_peak_height,
			  pedestal_peak_pos,
			  pedestal_peak_sigma );
  number_of_pedestal_events 
    = gaus->Integral ( pedestal_peak_pos - 5.0 * pedestal_peak_sigma,
		       pedestal_peak_pos + 5.0 * pedestal_peak_sigma );

  if ( number_of_pedestal_events<10 || number_of_pedestal_events > number_of_entries )
    return -444;
  noise   = 1.0 * log( 1.0 * number_of_entries / number_of_pedestal_events ) / GATE * 1000;
  gaus    ->Clear();
  delete gaus;
  return noise;

}

Double_t ana_MPPC::get_noise(Int_t entry){
  TF1 *gaus = new TF1("gaus","gaus",0,1000);
  gaus->SetParameters(pedestal_peak_height,pedestal_peak_pos,pedestal_peak_sigma);
  number_of_pedestal_events = gaus->Integral(pedestal_peak_pos-5.0*pedestal_peak_sigma,pedestal_peak_pos+5.0*pedestal_peak_sigma);

  if(number_of_pedestal_events<10||number_of_pedestal_events>entry)return -444;
  noise= 1.0*log(1.0*entry/number_of_pedestal_events)/GATE*1000;
  mean_pe = 1.0*log(1.0*entry/number_of_pedestal_events);
  gaus  ->Clear();
  delete gaus;
  return noise;
}

Double_t ana_MPPC::get_crosstalk_and_afterpulse(){
  TF1 *gaus = new TF1("gaus","gaus",0,1000);
  gaus->SetParameters(pedestal_peak_height,pedestal_peak_pos,pedestal_peak_sigma);
  number_of_pedestal_events = gaus->Integral(pedestal_peak_pos-5.0*pedestal_peak_sigma,pedestal_peak_pos+5.0*pedestal_peak_sigma);
  gaus->SetParameters(onepe_peak_height,onepe_peak_pos,onepe_peak_sigma);
  number_of_onepe_events = gaus->Integral(onepe_peak_pos-5.0*onepe_peak_sigma,onepe_peak_pos+5.0*onepe_peak_sigma);
  if(number_of_pedestal_events<10||number_of_pedestal_events>number_of_entries)return -444;
  crosstalk_and_afterpulse= 1.0-number_of_onepe_events/(number_of_pedestal_events*log(1.0*number_of_entries/number_of_pedestal_events));
  gaus -> Clear();
  delete gaus;
  return crosstalk_and_afterpulse;
}


Double_t ana_MPPC::get_crosstalk_and_afterpulse(Int_t entry){

  TF1 *gaus = new TF1("gaus","gaus",0,1000);
  gaus->SetParameters(pedestal_peak_height,pedestal_peak_pos,pedestal_peak_sigma);
  number_of_pedestal_events = gaus->Integral(pedestal_peak_pos-5.0*pedestal_peak_sigma,pedestal_peak_pos+5.0*pedestal_peak_sigma);
  gaus->SetParameters(onepe_peak_height,onepe_peak_pos,onepe_peak_sigma);
  number_of_onepe_events = gaus->Integral(onepe_peak_pos-5.0*onepe_peak_sigma,onepe_peak_pos+5.0*onepe_peak_sigma);
  if(number_of_pedestal_events<10||number_of_pedestal_events>entry)return -444;
  crosstalk_and_afterpulse= 1.0-number_of_onepe_events/(number_of_pedestal_events*log(1.0*entry/number_of_pedestal_events));
  gaus -> Clear();
  delete gaus;
  return crosstalk_and_afterpulse;
}

Double_t ana_MPPC::get_mean_pe(){
  return mean_pe;
}

Bool_t ana_MPPC::analysis(TH1F *noise_hist){
  TF1 *gaus = new TF1("gaus","gaus",HIST_MIN,HIST_MAX);
  TF1 *dgaus = new TF1("dgaus","gaus(0)+gaus(3)",HIST_MIN,HIST_MAX);
  Double_t par[6];

  fBinWidth=noise_hist->GetBinWidth(150);
  fBinMax=noise_hist->GetMaximumBin();
  fBinEdge=noise_hist->GetBinLowEdge(1);
  fXMax=fBinEdge+fBinWidth*fBinMax;
  pedestal_peak_pos=fXMax;//define maximum X position as pedestal peak

  //fit around pedestal_peak_pos and get sigma and peak number
  gaus->SetParLimits(0,0,100000000);
  gaus->SetParLimits(1,0,300);
  gaus->SetParLimits(2,0,10);
  noise_hist->Fit("gaus","qn","",(pedestal_peak_pos-FITRANGE),(pedestal_peak_pos+FITRANGE));
  pedestal_peak_pos = gaus->GetParameter(1);
  pedestal_peak_height = gaus->GetParameter(0);
  pedestal_peak_sigma = gaus->GetParameter(2);
  gaus->GetParameters(par);

  //define new histgram without pedestal peak
  Int_t fMinX_noise_hist_wo_pedestal=static_cast<Int_t>(pedestal_peak_pos+pedestal_peak_sigma*without_sigma);
 
  TH1F *noise_hist_wo_pedestal; 
  noise_hist_wo_pedestal = new TH1F("noise_hist_wo_pedestal","noise_hist_wo_pedestal",HIST_MAX-fMinX_noise_hist_wo_pedestal,fMinX_noise_hist_wo_pedestal,HIST_MAX);

  for(Int_t bin=0;bin<HIST_MAX-fMinX_noise_hist_wo_pedestal;bin++){
    Int_t ftempbincontent=noise_hist->GetBinContent(bin+fMinX_noise_hist_wo_pedestal-HIST_MIN);
    for(Int_t nevent=0;nevent<ftempbincontent;nevent++){
      noise_hist_wo_pedestal->Fill(fMinX_noise_hist_wo_pedestal+bin-1);
    }
  }

  //get maximum bin X position 
  fBinWidth=noise_hist_wo_pedestal->GetBinWidth(HIST_MAX-1);
  fBinMax=noise_hist_wo_pedestal->GetMaximumBin();
  fBinEdge=noise_hist_wo_pedestal->GetBinLowEdge(1);
  fXMax=fBinEdge+fBinWidth*fBinMax;
  onepe_peak_pos=fXMax;//define maximum X position as onepe peak
  gaus->SetParLimits(0,0,100000000);
  gaus->SetParLimits(1,0,300);
  gaus->SetParLimits(2,0,10);
  noise_hist_wo_pedestal->Fit("gaus","qn","",(onepe_peak_pos-FITRANGE),(onepe_peak_pos+FITRANGE));
  onepe_peak_pos = gaus->GetParameter(1);
  onepe_peak_height = gaus->GetParameter(0);
  onepe_peak_sigma = gaus->GetParameter(2);
  gaus->GetParameters(par+3);

  if(onepe_peak_sigma<0||onepe_peak_height<0)return false;
  //fit by double gaussiun
  dgaus->SetParameters(par);


  noise_hist->Fit("dgaus","qn","",(pedestal_peak_pos-FITRANGE),(onepe_peak_pos+FITRANGE));
  pedestal_peak_pos = dgaus->GetParameter(1);
  pedestal_peak_height = dgaus->GetParameter(0);
  pedestal_peak_sigma = dgaus->GetParameter(2);
  onepe_peak_pos = dgaus->GetParameter(4);
  onepe_peak_height = dgaus->GetParameter(3);
  onepe_peak_sigma = dgaus->GetParameter(5);
  fChisquare=dgaus->GetChisquare();
  fNDF=dgaus->GetNDF();
  gain=onepe_peak_pos-pedestal_peak_pos;

  delete noise_hist_wo_pedestal;
  delete gaus;
  delete dgaus;


  if(pedestal_peak_height<onepe_peak_pos){
    return false;
  }

  return true;
}


Bool_t ana_MPPC::analysis_logain(TH1F *noise_hist){
  Int_t ped_bin   = ( noise_hist -> GetMaximumBin() );
  logain_pedestal =   noise_hist -> GetBinLowEdge( ped_bin );
  return true; 
}


#endif
