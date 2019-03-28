#include "DetectorResponse.hh"
#include <CLHEP/Random/Randomize.h>
#include <algorithm>
#include <TRandom.h>
#include "G4SystemOfUnits.hh"

using namespace std;

DetectorResponse::DetectorResponse(DetectorDimension* fdim)
{

  //=====================================
  //=======     PE correction     =======
  //=====================================
  //Definition of "pe_ratio"
  for(int imod=0; imod<24; imod++){
    for(int iview=0; iview<2; iview++){
      for(int ipln=0; ipln<18; ipln++){
        for(int ich=0; ich<80; ich++){
          peratio[imod][iview][ipln][ich]  = 1.;
        }
      }
    }
    for(int itype=0;itype<2;itype++){
      attleng[imod][itype]  = C_FiberAttLeng;
      reflec [imod][itype]  = 0.;
    }
  }

  //PE correction
  cout<<"***** Set PE correction value *****"<<endl;
  string pedata_n    =Form("%s/pefiles/t2krun9/pemm_data.txt",MY_PATH);
  string pedata_mc_n =Form("%s/pefiles/t2krun9/pemm_mc.txt"  ,MY_PATH);
  ifstream pedata   (pedata_n   .c_str());
  ifstream pedata_mc(pedata_mc_n.c_str());
  if(!pedata||!pedata_mc){
    cout << "Error: Cannot open pe file (data)" << endl;
    cout << pedata_n    << endl;
    cout << pedata_mc_n << endl;
    exit(1);
  }
  int imod, iview, ipln, ich;
  int imod2, iview2, ipln2, ich2;
  double pe, pe2;
  while(pedata >> imod >> iview >> ipln >> ich >> pe){
    if(pedata_mc >> imod2 >> iview2 >> ipln2 >> ich2 >> pe2){
      if(imod==imod2&&iview==iview2&&ipln==ipln2&&ich==ich2){
        if(pe2!=0){
          peratio[imod][iview][ipln][ich] = pe/pe2;
        }
        else{
          cout << "Error: Light yield in MC is zero!" << endl;
          cout 
            << "mod" << imod << " view" << iview << " pln" << ipln << " ch" << ich << " peratio" 
            <<  peratio[imod][iview][ipln][ich] << endl;
          exit(1);
        }
      }
      else{
        cout << "Error: Mod/View/Pln/Ch is NOT consistent b/w data and MC." << endl;
        cout 
          << "mod" << imod << " view" << iview << " pln" << ipln << " ch" << ich
          << endl;
        exit(1);
      }
    }
    else{
      cout << "Error: No corresponding MC info." << endl;
      exit(1);
    }
  }
  pedata   .close();
  pedata_mc.close();

  //Fiber attenuation length & reflection paramter
  cout <<"****** Set Attenuation length & Reflection parameter *****" << endl;
  string attfile_n = Form("%s/pefiles/t2krun9/fiber_att.txt",MY_PATH);
  ifstream attfile(attfile_n.c_str());
  if(!attfile){
    cout << "Error: No attenuation length file." << endl;
    cout << attfile_n << endl;
    exit(1);
  }
  int att_mod, att_type;
  double attlen,reflection;
  while(attfile >> att_mod >> att_type >> attlen >> reflection){
    attleng[att_mod][att_type] = attlen/10.;
    reflec [att_mod][att_type] = reflection;
  }
  attfile.close();


  //for(int imod=0; imod<24; imod++){
  //  for(int iview=0; iview<2; iview++){
  //    for(int ipln=0; ipln<18; ipln++){
  //      for(int ich=0; ich<80; ich++){
  //        cout 
  //          << " mod:"  << imod
  //          << " view:" << iview
  //          << " pln:"  << ipln
  //          << " ch:"   << ich
  //          << " ratio:"<< peratio[imod][iview][ipln][ich]
  //          << endl;
  //      }
  //    }
  //  }
  //}


  detdim = fdim;
    
}

DetectorResponse::~DetectorResponse()
{
}

void DetectorResponse::ApplyScintiResponse(G4double* edep, G4double* length, G4Track* aTrack)
{
  //Quenching
  BirksSaturation(edep, length, aTrack);
  
  return;
}

void DetectorResponse::ApplyScintiResponse2(G4double* edep, G4Track* aTrack)
{
  //Quenching
  BirksSaturation2(edep, aTrack);
  
  return;
}

void DetectorResponse::BirksSaturation(G4double* edep, G4double* length, G4Track* aTrack)
{
  
  //G4double              kineticE = aTrack->GetKineticEnergy();
  G4ParticleDefinition* particle = aTrack->GetDefinition();
  //G4Material*           material = aTrack->GetMaterial();


  if(particle->GetPDGCharge()==0){
    return;
  }

  G4double dedx = (*edep)/(*length); //Same method as INGRIDWaterModule

  (*edep) = (*edep)/(1. + C_Corr_Birks*dedx);

  return;
}

void DetectorResponse::BirksSaturation2(G4double* edep, G4Track* aTrack)
{
  
  G4double              kineticE = aTrack->GetKineticEnergy();
  G4ParticleDefinition* particle = aTrack->GetDefinition();
  G4Material*           material = aTrack->GetMaterial();


  if(particle->GetPDGCharge()==0){
    return;
  }

  G4double dedx = emcal.GetDEDX(kineticE, particle, material)/(MeV/cm);

  (*edep) = (*edep)/(1. + C_Corr_Birks*dedx);

  return;
}

void DetectorResponse::ApplyLightCollection(G4double* edep, G4int mod, G4int pln, G4int view, G4int ch, G4ThreeVector hitpos){

  distToFiber = 0.;

  G4double X_tmp, Y_tmp, Z_tmp;
  X_tmp = hitpos[0]/cm;
  Y_tmp = hitpos[1]/cm;
  Z_tmp = hitpos[2]/cm;

  //H-INGRID, V-INGRID, B2-INGRID
  if(mod<NUMINGMOD){
    detdim->GetPosING(3,pln,view,ch,&fiberpos[0],&fiberpos[1],&fiberpos[2]);
    if(view==TopView){
      distToFiber = fabs(X_tmp - fiberpos[0]/10.);
    }
    else if(view==SideView){
      distToFiber = fabs(Y_tmp - fiberpos[1]/10.);
    }
    *edep *= exp(-1.*distToFiber/C_ScintiAttLeng)*C_INGRIDFactor;
  }
  //Proton Module (On-axis/B2)
  else if(mod==MOD_PM || mod==MOD_B2_CH){
    detdim->GetPosPM(pln,view,ch,&fiberpos[0],&fiberpos[1],&fiberpos[2]);

    if(view==TopView){
      distToFiber = fabs(X_tmp - fiberpos[0]/10.);
    }
    else if(view==SideView){
      distToFiber = fabs(Y_tmp - fiberpos[1]/10.);
    }
    
    if(pln!=0&&ch>=8&&ch<24){
      *edep *= exp(-1.*distToFiber/C_ScintiAttLeng)*C_ScibarFactor;
    }
    else{
      *edep *= exp(-1.*distToFiber/C_ScintiAttLeng)*C_INGRIDFactor;
    }
  }
  //INGRIDWaterModule, WAGASCI water module
  else if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
    int grid, gridch;
    detdim->GetWMGridCh(pln,view,ch,&grid,&gridch);
    detdim->GetPosWM(mod,pln,view,ch,grid,&fiberpos[0],&fiberpos[1],&fiberpos[2]);
    if(view==TopView){
      if(grid==0){
	distToFiber = fabs(X_tmp - (fiberpos[0]-C_WMScintiHoleShift1)/10.);
      }
      else{
	distToFiber = fabs(Z_tmp - (fiberpos[2]+C_WMScintiHoleShift1)/10.);
      }
    }
    else{
      if(grid==0){
	distToFiber = fabs(Y_tmp - (fiberpos[1]-C_WMScintiHoleShift1)/10.);
      }
      else{
	distToFiber = fabs(Z_tmp - (fiberpos[2]-C_WMScintiHoleShift1)/10.);
      }
    }
    *edep *= exp(-1.*distToFiber/C_WMScintiAttLeng)*C_WMFactor;
  }
  
  *edep *= peratio[mod][view][pln][ch];
  
  return;
}
//void DetectorResponse::ApplyFiberResponse(G4double* edep, G4double* time, G4int mod, G4int view, G4ThreeVector hitpos)
void DetectorResponse::ApplyFiberResponse(G4double* edep, G4double* time, G4int mod, G4int view, G4int pln, G4int ch, G4ThreeVector hitpos)
{

  distToMPPC  = 0.;
  distToMPPC2 = 0.; //concerning reflection

  double fiber_attlen = C_FiberAttLeng;
  double fiber_reflec = 0.;
  G4double X_tmp, Y_tmp, Z_tmp;
  X_tmp =  hitpos[0]/cm;
  Y_tmp =  hitpos[1]/cm;
  Z_tmp =  hitpos[2]/cm;

  //H-INGRID, V-INGRID, B2-INGRID
  if(mod<NUMINGMOD){
    if(view==TopView){
      distToMPPC  = fabs(C_INGScintiLength/10./2. - Y_tmp);
      distToMPPC2 = fabs(C_INGScintiLength/10.    + Y_tmp);
    }
    else{
      distToMPPC  = fabs(C_INGScintiLength/10./2. + X_tmp);
      distToMPPC2 = fabs(C_INGScintiLength/10.    - X_tmp);
    }
    fiber_attlen = attleng[mod][0];
    fiber_reflec = reflec [mod][0];
  }
  //Proton Module (On-axis/B2)
  else if(mod==MOD_PM || mod==MOD_B2_CH){
    if(view==TopView){
      distToMPPC  = fabs(C_PMScintiLength/10./2. - Y_tmp);
      distToMPPC2 = fabs(C_PMScintiLength/10.    + Y_tmp);
    }
    else{
      distToMPPC  = fabs(C_PMScintiLength/10./2. + X_tmp);
      distToMPPC2 = fabs(C_PMScintiLength/10.    - X_tmp);
    }
    if(pln==0||ch<8||ch>=24){ fiber_attlen = attleng[mod][1];fiber_reflec = reflec[mod][1];}
    else                    { fiber_attlen = attleng[mod][0];fiber_reflec = reflec[mod][0];}
  }
  //INGRIDWaterModule, WAGASCI watermodule
  else if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
    if(view==TopView){
      distToMPPC  = fabs(C_WMScintiLength/10./2. - Y_tmp) + C_WMFiberMargin/10.;
      distToMPPC2 = fabs(C_WMScintiLength/10.    + Y_tmp) + C_WMFiberMargin/10.;
    }
    else{
      distToMPPC  = fabs(C_WMScintiLength/10./2. + X_tmp) + C_WMFiberMargin/10.;
      distToMPPC2 = fabs(C_WMScintiLength/10.    - X_tmp) + C_WMFiberMargin/10.;
    }
    fiber_attlen = attleng[mod][0];
    fiber_reflec = reflec [mod][0];
  }
  
  //Attenuation
  //*edep *= exp(-1.*distToMPPC/C_FiberAttLeng);
  *edep *= exp(-1.*distToMPPC/fiber_attlen) + fiber_reflec*exp(-1.*distToMPPC2/fiber_attlen);

  //Delay in fiber
  *time += C_TransTimeInFiber*distToMPPC;

  return;
}

void DetectorResponse::ApplyFiberResponseV(G4double* edep, G4double* time, G4int pln, G4ThreeVector hitpos)
{

  distToMPPC = 0.;
  
  //Attenuation
  *edep *= exp(-1.*distToMPPC/C_FiberAttLeng);
  
  //Delay in fiber
  *time += C_TransTimeInFiber*distToMPPC;
  
  return;
}

void DetectorResponse::ApplyMPPCResponse(G4double edep, G4double* pe, G4int mod=0)
{
  G4double npe;
  
  //Energy to PE
  if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
    npe = edep*(C_WMMeV2PE);
  }
  else{
    npe = edep*(C_MeV2PE);
  }

  //MPPC linearity
  if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
    npe = C_MPPCPixel * (1. - exp( C_WMEff_PDE * npe / C_MPPCPixel ));
  }
  else{
    npe = C_MPPCPixel * (1. - exp( C_Eff_PDE * npe / C_MPPCPixel ));
  }
  
  //Poisson statistics & 1 pe resolution
  npe = (int) CLHEP::RandPoisson::shoot(npe); //Common setting for INGRID and WAGASCI

  ////Modified cross-talk and after pulse
  if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
    npe += (int) CLHEP::RandPoisson::shoot(npe*C_WMCrossAfterRate);
  }
  else{
    npe += (int) CLHEP::RandPoisson::shoot(npe*C_CrossAfterRate);
  }

  if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
    npe = gRandom->Gaus(npe, npe*C_WMPixelGainVari);
  }
  else{
    npe = gRandom->Gaus(npe, npe*C_PixelGainVari);
  }

  *pe = npe;


  return;
}

void DetectorResponse::ApplyADCResponse(G4double *pe, G4double *lope, G4int* adc, G4int* loadc, G4int mod=0)
{

  G4double adc_tmp, loadc_tmp,Q,loQ;
  
  //PE to ADC
  if(mod==MOD_B2_WM){
    adc_tmp   = C_WMPedestal + (*pe)*C_WMGain;
    loadc_tmp = C_Pedestal   + (*pe)*C_LowGain*C_LowGain_corr;
  }
  else if(mod==MOD_ONAXIS_WM){
    adc_tmp   = C_WMPedestal + (*pe)*C_WMGain;
    loadc_tmp = C_Pedestal   + (*pe)*C_LowGain*C_LowGain_corr;
  }
  else{
    adc_tmp   = C_Pedestal + (*pe)*C_Gain;
    loadc_tmp = C_Pedestal + (*pe)*C_LowGain*C_LowGain_corr;
  }  

  if(mod!=MOD_B2_WM){

    //Electronics noise
    adc_tmp   = gRandom->Gaus(adc_tmp  ,C_ElecNoise);
    loadc_tmp = gRandom->Gaus(loadc_tmp,C_LowElecNoise);

    //ADC to Charge
    Q   = (adc_tmp)  /C_ADCtoCharge;
    loQ = (loadc_tmp)/C_LowADCtoCharge;

    //Non linearlity of high gain ADC
    if(Q<C_NonLinADCTh[0]){
      *adc=C_NonLinADC1[0]*Q+C_NonLinADC2[0];
    }
    else if(Q<C_NonLinADCTh[1]){
      *adc=C_NonLinADC1[1]*Q+C_NonLinADC2[1];
    }
    else if(Q<C_NonLinADCTh[2]){
      *adc=C_NonLinADC1[2]*Q+C_NonLinADC2[2];
    }
    else if(Q<C_NonLinADCTh[3]){
      *adc=C_NonLinADC1[3]*Q+C_NonLinADC2[3];
    }
    else{
      *adc=C_ADCSaturation;
    }

    //Non linearlity of low gain ADC
    if(Q<C_LowNonLinADCTh[0]){
      *adc=C_LowNonLinADC1[0]*Q+C_LowNonLinADC2[0];
    }
    else if(Q<C_LowNonLinADCTh[1]){
      *adc=C_LowNonLinADC1[1]*Q+C_LowNonLinADC2[1];
    }
    else if(Q<C_LowNonLinADCTh[2]){
      *adc=C_LowNonLinADC1[2]*Q+C_LowNonLinADC2[2];
    }
    else if(Q<C_LowNonLinADCTh[3]){
      *adc=C_LowNonLinADC1[3]*Q+C_LowNonLinADC2[3];
    }
    else{
      *adc=C_LowADCSaturation;
    }

    //ADC to PE
    if(mod==MOD_ONAXIS_WM){
      *pe   = (float)((*adc)   - C_WMPedestal)/C_WMGain;
      *lope = (float)((*loadc) - C_Pedestal)/C_LowGain;
    }
    else{
      *pe   = (float)((*adc)   - C_Pedestal)/C_Gain;
      *lope = (float)((*loadc) - C_Pedestal)/C_LowGain;
    }
  }
}

void DetectorResponse::ApplyTDCResponse(G4double time, G4int* tdc)
{
  *tdc = 0;
}
