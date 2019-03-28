#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RandGauss.h"
#include "Randomize.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"


//#define SANDMU0
#define SANDMU1
//#define COSMIC0
#define COSMIC1

#ifdef COSMIC1
const double distributionX[2][100][2] =
{
  {
    {-1.00, 0.132},
    {-0.98, 0.029},
    {-0.96, 0.026},
    {-0.94, 0.027},
    {-0.92, 0.028},
    {-0.90, 0.032},
    {-0.88, 0.034},
    {-0.86, 0.045},
    {-0.84, 0.044},
    {-0.82, 0.039},
    {-0.80, 0.070},
    {-0.78, 0.043},
    {-0.76, 0.054},
    {-0.74, 0.072},
    {-0.72, 0.067},
    {-0.70, 0.052},
    {-0.68, 0.068},
    {-0.66, 0.102},
    {-0.64, 0.090},
    {-0.62, 0.103},
    {-0.60, 0.105},
    {-0.58, 0.090},
    {-0.56, 0.174},
    {-0.54, 0.095},
    {-0.52, 0.118},
    {-0.50, 0.142},
    {-0.48, 0.178},
    {-0.46, 0.144},
    {-0.44, 0.152},
    {-0.42, 0.186},
    {-0.40, 0.214},
    {-0.38, 0.168},
    {-0.36, 0.237},
    {-0.34, 0.208},
    {-0.32, 0.227},
    {-0.30, 0.236},
    {-0.28, 0.304},
    {-0.26, 0.285},
    {-0.24, 0.315},
    {-0.22, 0.350},
    {-0.20, 0.351},
    {-0.18, 0.383},
    {-0.16, 0.391},
    {-0.14, 0.401},
    {-0.12, 0.404},
    {-0.10, 0.392},
    {-0.08, 0.368},
    {-0.06, 0.336},
    {-0.04, 0.337},
    {-0.02, 0.697},
    {0.00, 1.000},
    {0.02, 0.315},
    {0.04, 0.319},
    {0.06, 0.343},
    {0.08, 0.372},
    {0.10, 0.390},
    {0.12, 0.409},
    {0.14, 0.406},
    {0.16, 0.418},
    {0.18, 0.398},
    {0.20, 0.422},
    {0.22, 0.387},
    {0.24, 0.360},
    {0.26, 0.389},
    {0.28, 0.318},
    {0.30, 0.320},
    {0.32, 0.300},
    {0.34, 0.347},
    {0.36, 0.245},
    {0.38, 0.331},
    {0.40, 0.296},
    {0.42, 0.239},
    {0.44, 0.244},
    {0.46, 0.278},
    {0.48, 0.250},
    {0.50, 0.210},
    {0.52, 0.173},
    {0.54, 0.291},
    {0.56, 0.168},
    {0.58, 0.174},
    {0.60, 0.189},
    {0.62, 0.155},
    {0.64, 0.201},
    {0.66, 0.138},
    {0.68, 0.112},
    {0.70, 0.240},
    {0.72, 0.137},
    {0.74, 0.108},
    {0.76, 0.083},
    {0.78, 0.117},
    {0.80, 0.075},
    {0.82, 0.074},
    {0.84, 0.075},
    {0.86, 0.060},
    {0.88, 0.056},
    {0.90, 0.050},
    {0.92, 0.048},
    {0.94, 0.044},
    {0.96, 0.046},
    {0.98, 0.262}
  },
  {
    {-1.00, 0.013},
    {-0.98, 0.002},
    {-0.96, 0.002},
    {-0.94, 0.001},
    {-0.92, 0.002},
    {-0.90, 0.002},
    {-0.88, 0.002},
    {-0.86, 0.003},
    {-0.84, 0.004},
    {-0.82, 0.004},
    {-0.80, 0.006},
    {-0.78, 0.006},
    {-0.76, 0.009},
    {-0.74, 0.010},
    {-0.72, 0.012},
    {-0.70, 0.008},
    {-0.68, 0.009},
    {-0.66, 0.024},
    {-0.64, 0.029},
    {-0.62, 0.027},
    {-0.60, 0.016},
    {-0.58, 0.016},
    {-0.56, 0.027},
    {-0.54, 0.020},
    {-0.52, 0.024},
    {-0.50, 0.038},
    {-0.48, 0.039},
    {-0.46, 0.043},
    {-0.44, 0.035},
    {-0.42, 0.046},
    {-0.40, 0.051},
    {-0.38, 0.052},
    {-0.36, 0.058},
    {-0.34, 0.064},
    {-0.32, 0.077},
    {-0.30, 0.073},
    {-0.28, 0.085},
    {-0.26, 0.108},
    {-0.24, 0.111},
    {-0.22, 0.137},
    {-0.20, 0.149},
    {-0.18, 0.169},
    {-0.16, 0.171},
    {-0.14, 0.188},
    {-0.12, 0.198},
    {-0.10, 0.208},
    {-0.08, 0.230},
    {-0.06, 0.248},
    {-0.04, 0.289},
    {-0.02, 0.581},
    {0.00, 1.000},
    {0.02, 0.326},
    {0.04, 0.296},
    {0.06, 0.297},
    {0.08, 0.281},
    {0.10, 0.281},
    {0.12, 0.263},
    {0.14, 0.246},
    {0.16, 0.233},
    {0.18, 0.206},
    {0.20, 0.180},
    {0.22, 0.144},
    {0.24, 0.137},
    {0.26, 0.103},
    {0.28, 0.094},
    {0.30, 0.096},
    {0.32, 0.070},
    {0.34, 0.066},
    {0.36, 0.058},
    {0.38, 0.057},
    {0.40, 0.049},
    {0.42, 0.038},
    {0.44, 0.055},
    {0.46, 0.048},
    {0.48, 0.047},
    {0.50, 0.028},
    {0.52, 0.021},
    {0.54, 0.035},
    {0.56, 0.018},
    {0.58, 0.023},
    {0.60, 0.040},
    {0.62, 0.041},
    {0.64, 0.035},
    {0.66, 0.014},
    {0.68, 0.014},
    {0.70, 0.022},
    {0.72, 0.020},
    {0.74, 0.024},
    {0.76, 0.015},
    {0.78, 0.015},
    {0.80, 0.009},
    {0.82, 0.008},
    {0.84, 0.008},
    {0.86, 0.006},
    {0.88, 0.005},
    {0.90, 0.005},
    {0.92, 0.005},
    {0.94, 0.004},
    {0.96, 0.004},
    {0.98, 0.024}
  }
};

const double distributionY[2][100][2] =
{
  {
    {-1.00, 0.247},
    {-0.98, 0.385},
    {-0.96, 0.380},
    {-0.94, 0.334},
    {-0.92, 0.329},
    {-0.90, 0.293},
    {-0.88, 0.284},
    {-0.86, 0.278},
    {-0.84, 0.217},
    {-0.82, 0.146},
    {-0.80, 0.222},
    {-0.78, 0.121},
    {-0.76, 0.167},
    {-0.74, 0.140},
    {-0.72, 0.127},
    {-0.70, 0.085},
    {-0.68, 0.100},
    {-0.66, 0.243},
    {-0.64, 0.144},
    {-0.62, 0.147},
    {-0.60, 0.112},
    {-0.58, 0.099},
    {-0.56, 0.149},
    {-0.54, 0.081},
    {-0.52, 0.122},
    {-0.50, 0.126},
    {-0.48, 0.105},
    {-0.46, 0.089},
    {-0.44, 0.069},
    {-0.42, 0.087},
    {-0.40, 0.086},
    {-0.38, 0.060},
    {-0.36, 0.067},
    {-0.34, 0.059},
    {-0.32, 0.048},
    {-0.30, 0.041},
    {-0.28, 0.041},
    {-0.26, 0.038},
    {-0.24, 0.031},
    {-0.22, 0.029},
    {-0.20, 0.023},
    {-0.18, 0.021},
    {-0.16, 0.016},
    {-0.14, 0.014},
    {-0.12, 0.013},
    {-0.10, 0.012},
    {-0.08, 0.011},
    {-0.06, 0.012},
    {-0.04, 0.012},
    {-0.02, 0.064},
    {0.00, 0.104},
    {0.02, 0.018},
    {0.04, 0.017},
    {0.06, 0.017},
    {0.08, 0.016},
    {0.10, 0.017},
    {0.12, 0.017},
    {0.14, 0.019},
    {0.16, 0.024},
    {0.18, 0.025},
    {0.20, 0.029},
    {0.22, 0.030},
    {0.24, 0.037},
    {0.26, 0.041},
    {0.28, 0.041},
    {0.30, 0.051},
    {0.32, 0.062},
    {0.34, 0.065},
    {0.36, 0.056},
    {0.38, 0.086},
    {0.40, 0.092},
    {0.42, 0.073},
    {0.44, 0.086},
    {0.46, 0.110},
    {0.48, 0.135},
    {0.50, 0.136},
    {0.52, 0.090},
    {0.54, 0.182},
    {0.56, 0.126},
    {0.58, 0.142},
    {0.60, 0.201},
    {0.62, 0.202},
    {0.64, 0.365},
    {0.66, 0.195},
    {0.68, 0.163},
    {0.70, 0.332},
    {0.72, 0.272},
    {0.74, 0.384},
    {0.76, 0.264},
    {0.78, 0.400},
    {0.80, 0.266},
    {0.82, 0.351},
    {0.84, 0.422},
    {0.86, 0.458},
    {0.88, 0.480},
    {0.90, 0.552},
    {0.92, 0.604},
    {0.94, 0.755},
    {0.96, 1.000},
    {0.98, 0.957}
  },
  {
    {-1.00, 0.016},
    {-0.98, 0.006},
    {-0.96, 0.004},
    {-0.94, 0.003},
    {-0.92, 0.007},
    {-0.90, 0.013},
    {-0.88, 0.052},
    {-0.86, 0.039},
    {-0.84, 0.023},
    {-0.82, 0.017},
    {-0.80, 0.067},
    {-0.78, 0.030},
    {-0.76, 0.079},
    {-0.74, 0.017},
    {-0.72, 0.020},
    {-0.70, 0.020},
    {-0.68, 0.030},
    {-0.66, 0.288},
    {-0.64, 0.068},
    {-0.62, 0.093},
    {-0.60, 0.095},
    {-0.58, 0.114},
    {-0.56, 0.133},
    {-0.54, 0.089},
    {-0.52, 0.277},
    {-0.50, 0.192},
    {-0.48, 0.138},
    {-0.46, 0.125},
    {-0.44, 0.104},
    {-0.42, 0.150},
    {-0.40, 0.166},
    {-0.38, 0.121},
    {-0.36, 0.124},
    {-0.34, 0.151},
    {-0.32, 0.117},
    {-0.30, 0.097},
    {-0.28, 0.099},
    {-0.26, 0.101},
    {-0.24, 0.071},
    {-0.22, 0.075},
    {-0.20, 0.060},
    {-0.18, 0.050},
    {-0.16, 0.032},
    {-0.14, 0.027},
    {-0.12, 0.019},
    {-0.10, 0.017},
    {-0.08, 0.014},
    {-0.06, 0.013},
    {-0.04, 0.017},
    {-0.02, 0.048},
    {0.00, 0.076},
    {0.02, 0.027},
    {0.04, 0.020},
    {0.06, 0.018},
    {0.08, 0.020},
    {0.10, 0.023},
    {0.12, 0.029},
    {0.14, 0.035},
    {0.16, 0.049},
    {0.18, 0.056},
    {0.20, 0.069},
    {0.22, 0.067},
    {0.24, 0.094},
    {0.26, 0.099},
    {0.28, 0.092},
    {0.30, 0.119},
    {0.32, 0.173},
    {0.34, 0.141},
    {0.36, 0.134},
    {0.38, 0.207},
    {0.40, 0.199},
    {0.42, 0.133},
    {0.44, 0.149},
    {0.46, 0.183},
    {0.48, 0.297},
    {0.50, 0.426},
    {0.52, 0.108},
    {0.54, 0.151},
    {0.56, 0.135},
    {0.58, 0.106},
    {0.60, 0.133},
    {0.62, 0.160},
    {0.64, 0.750},
    {0.66, 0.195},
    {0.68, 0.116},
    {0.70, 0.145},
    {0.72, 0.141},
    {0.74, 0.461},
    {0.76, 0.171},
    {0.78, 0.306},
    {0.80, 0.116},
    {0.82, 0.144},
    {0.84, 0.170},
    {0.86, 0.206},
    {0.88, 0.140},
    {0.90, 0.212},
    {0.92, 0.383},
    {0.94, 0.660},
    {0.96, 1.000},
    {0.98, 0.592}
  }
};
#endif

PrimaryGeneratorAction::
  PrimaryGeneratorAction(Neut *neut0,RunAction* rac, EventAction* evt,int nd,int flavor0)
:runaction(rac)
{
  eventaction      = evt;
  neut_file        = neut0;
  module_mode      = nd;
  neutrino_flavor  = flavor0;
  particleGun      = new G4ParticleGun(1);
  particleTable    = G4ParticleTable::GetParticleTable();

  runaction->NotEntry = 0;
}

//
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  if(particleGun!=NULL){
    delete particleGun; particleGun=NULL;
  }
}

//
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
#ifdef DEBUG_PRINGEN
  G4cout << " ============ Primary Generator : GeneratePrimaries ===========  " << G4endl;
#endif

  SecondaryVector Secondary;
  Neut   *neut        = neut_file;
  int    fdid         = 0;
  int    mode         = 0;
  float  pos[3]       = {-1000.,-1000.,-1000.};
  int    ID           = -1;
  float  direction[3] = {-1000.,-1000.,-1000.};
  int    vertex_flag  = 0;
  int    flayer       = -1 ;
  float  dz = -999;
  double prob         = -1.;//Probability to hit scintillator for Proton Module
  int    scitype      = -1 ;//Scintillator type

  //======================================================
  //=======     Particle-Gun: Sand-muon mode     =======
  //======================================================
  if(neutrino_flavor==11){

    G4ParticleDefinition* particle;
    particle = particleTable->FindParticle(13);
    particleGun->SetParticleDefinition(particle);
    G4float mass  = particle->GetPDGMass();

    G4float mumom = 1*GeV;
    G4float energy = sqrt(mumom*mumom+mass*mass) - mass;
    particleGun->SetParticleEnergy(energy);
    float pi_calc = 4.0*atan(1.0);
    float thetamax_calc = 70./180.*pi; //PDG2017
    float rand_theta_calc;
    float theta_calc, cos2theta_calc;
    float costheta_calc;
    float phi_calc;
    float xpos_calc;
    float ypos_calc;
    float zpos_calc;
    double cos_calc;
    float px;
    float py;
    float pz;

#ifdef SANDMU0
    while(1){
      costheta_calc = G4UniformRand();
      phi_calc = 2*pi_calc*G4UniformRand();
      cos2theta_calc = costheta_calc*costheta_calc;
      rand_theta_calc = G4UniformRand();
      if(cos2theta_calc>rand_theta_calc){
        theta_calc = acos(costheta_calc);
        break;
      }
    }
    px = sin(phi_calc)*sin(theta_calc);
    py = -cos(theta_calc);
    pz = cos(phi_calc)*sin(theta_calc);

    //Position setting
    xpos_calc = (G4UniformRand()-0.5)*650 + B2WM_OFFSETX;
    ypos_calc = 200. + B2WM_OFFSETY;
    zpos_calc = (G4UniformRand()-0.5)*650 + B2WM_OFFSETZ;

#else
#ifdef SANDMU1

    //Oriented to the WAGASCI
    double const_a = 17.1519;
    theta_calc = G4UniformRand()*2.*PI;
    while(true){
      cos_calc = G4UniformRand()*0.65+0.35; //0.35 to 1
      double prob0 = G4UniformRand()*1.;
      if(prob0<exp(-1.*const_a*(1.-cos_calc))){break;}
    }
    px = sqrt(1.-cos_calc*cos_calc)*cos(theta_calc);
    py = sqrt(1.-cos_calc*cos_calc)*sin(theta_calc);
    pz = cos_calc;

    while(true)
    {
      xpos_calc = (G4UniformRand()-0.5)*400.;// + B2WM_OFFSETX;
      ypos_calc = (G4UniformRand()-0.5)*400.;// + B2WM_OFFSETY;
      zpos_calc =                       -60.;// + B2WM_OFFSETZ;

      double tmpx = xpos_calc+B2WM_OFFSETX-B2ING_OFFSETX;
      double tmpy = ypos_calc+B2WM_OFFSETY-B2ING_OFFSETY;
      double tmpz = zpos_calc+B2WM_OFFSETZ-B2ING_OFFSETZ;

      double prod  = (xpos_calc*px + ypos_calc*py + zpos_calc*pz);
      double norm1 = xpos_calc*xpos_calc+ypos_calc*ypos_calc+zpos_calc*zpos_calc;
      double norm2 = px*px + py*py + pz*pz;
      double cos2  = prod*prod/norm1/norm2;
      if(cos2==0.){continue;}
      double dist2 = norm1*(1.-cos2)/cos2;
      if(dist2>50.*50.){continue;}

      prod  = (tmpx*px + tmpy*py + tmpz*pz);
      norm1 = tmpx*tmpx+tmpy*tmpy+tmpz*tmpz;
      norm2 = px*px + py*py + pz*pz;
      cos2  = prod*prod/norm1/norm2;
      if(cos2==0.){continue;}
      dist2 = norm1*(1.-cos2)/cos2;
      if(dist2<50.*50.){break;}
    }

    int imod = -1;
    while(1){
      imod = (int) (G4UniformRand()*2.);
      if(imod==0||imod==1){break;}
    }
    if(imod==0){
      xpos_calc += B2WM_OFFSETX;
      ypos_calc += B2WM_OFFSETY;
      zpos_calc += B2WM_OFFSETZ;
    }
    else if(imod==1){
      xpos_calc += PM_OFFSETX;
      ypos_calc += PM_OFFSETY;
      zpos_calc += PM_OFFSETZ;
    }
#endif
#endif

    particleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
    G4float xpos = xpos_calc*cm;
    G4float ypos = ypos_calc*cm;
    G4float zpos = zpos_calc*cm;
    particleGun->SetParticlePosition(G4ThreeVector(xpos, ypos, zpos));
    particleGun->SetParticleTime(0.0*ns);

    SimVertexSummary* simvertex = new SimVertexSummary();
    simvertex -> Clear();
    simvertex -> nutype   = neutrino_flavor;
    simvertex -> inttype  = -100;
    simvertex -> nuE      = mumom/GeV;
    simvertex -> xnu      = xpos/cm;
    simvertex -> ynu      = ypos/cm;
    simvertex -> znu      = zpos/cm;
    simvertex -> gmomx.push_back(px);
    simvertex -> gmomy.push_back(py);
    simvertex -> gmomz.push_back(pz);
    //simvertex -> mod      = ID;
    simvertex -> norm     =1.;
    simvertex -> totcrsne		= 1.;

    runaction  -> GetEvtSum() -> AddSimVertex( simvertex );    
    particleGun->GeneratePrimaryVertex(anEvent);
    return;
  }

  //====================================================
  //=======     Particle-Gun: Cosmic-muon mode     =======
  //====================================================
  if(neutrino_flavor==12){

    G4ParticleDefinition* particle;
    particle = particleTable->FindParticle(13);
    particleGun->SetParticleDefinition(particle);
    G4float mass  = particle->GetPDGMass();
    G4float mumom = 2*GeV;
    G4float energy = sqrt(mumom*mumom+mass*mass) - mass;
    particleGun->SetParticleEnergy(energy);
    float pi_calc = 4.0*atan(1.0);
    float rand_theta_calc;
    float theta_calc, cos2theta_calc;
    float costheta_calc;
    float phi_calc;
    float xpos_calc;
    float ypos_calc;
    float zpos_calc;
    double cos_calc;
    float px;
    float py;
    float pz;


#ifdef COSMIC0
    costheta_calc = 2*(G4UniformRand()-0.5);
    phi_calc = 2*pi_calc*G4UniformRand();
    theta_calc=acos(costheta_calc);

    G4double pos_select = G4UniformRand();
    if(pos_select>=0. && pos_select<=0.75){ //Wall
      xpos_calc = (G4UniformRand()-0.5)*50. + B2WM_OFFSETX;
      ypos_calc = (G4UniformRand()-0.5)*50. + B2WM_OFFSETY;
      zpos_calc = (G4UniformRand()-0.5)*10. + (B2WM_OFFSETZ - 150.);
    }
    else if(pos_select>0.75 && pos_select<=0.89){ //Left pillar
      px = -cos(theta_calc);
      py = sin(phi_calc)*sin(theta_calc);
      pz = cos(phi_calc)*sin(theta_calc);
      xpos_calc = -373.5;
      ypos_calc = (G4UniformRand()-0.5)*100 + B2WM_OFFSETY;
      zpos_calc = (G4UniformRand()-0.5)*100 + B2WM_OFFSETZ - 50.;
    }
    else if(pos_select>0.89 && pos_select<=0.99){ //Right pillar
      px = cos(theta_calc);
      py = sin(phi_calc)*sin(theta_calc);
      pz = cos(phi_calc)*sin(theta_calc);
      xpos_calc = -773.5;
      ypos_calc = (G4UniformRand()-0.5)*100 + B2WM_OFFSETY;
      zpos_calc = (G4UniformRand()-0.5)*100 + B2WM_OFFSETZ - 50.;
    }
    else if(pos_select>0.99 && pos_select<=1.){ //Ceiling
      px = sin(phi_calc)*sin(theta_calc);
      py = -cos(theta_calc);
      pz = cos(phi_calc)*sin(theta_calc);
      xpos_calc = (G4UniformRand()-0.5)*50 + B2WM_OFFSETX + 50.;
      ypos_calc = 387. + B2WM_OFFSETY;
      zpos_calc = (G4UniformRand()-0.5)*50 + B2WM_OFFSETZ - 150.;
    }
#else
#ifdef COSMIC1
    //Oriented to the WAGASCI
    double a_calc,b_calc;
    double norm_ab;
    int  on_off_axis;
    int  imod;
    while(true){
      imod = (int)(G4UniformRand()*5.);
      if(
             imod==0  //WAGASCI
          || imod==1  //INGRIDWM
          || imod==2  //INGmod3
          || imod==3  //PM
          || imod==4  //B2ING
        )
      {
        if     (imod==0){on_off_axis==0;} //B2
        else if(imod==1){on_off_axis==1;} //SS
        else if(imod==2){on_off_axis==1;} //SS
        else if(imod==3){on_off_axis==0;} //B2
        else if(imod==4){on_off_axis==0;} //B2
        break;
      }
    }
    while(true){
      double r_calc = G4UniformRand()*1.0; //0. to 1.
      double t_calc = G4UniformRand()*2.*PI;
      a_calc = r_calc*cos(t_calc);
      b_calc = r_calc*sin(t_calc);
      double a_prob = G4UniformRand()*1.;
      double b_prob = G4UniformRand()*1.;
      int  a_cos_id = (a_calc+1.)/0.02;
      int  b_cos_id = (b_calc+1.)/0.02;
      if(a_cos_id>=100) a_cos_id=99;
      if(a_cos_id<   0) a_cos_id= 0;
      if(b_cos_id>=100) b_cos_id=99;
      if(b_cos_id<   0) b_cos_id= 0;
      norm_ab = a_calc*a_calc+b_calc*b_calc;
      if( a_prob<distributionX[on_off_axis][a_cos_id][1]&&
          b_prob<distributionY[on_off_axis][b_cos_id][1]&&
          norm_ab<=1.
          )
      {break;}
    }
    px = a_calc;
    py = b_calc;
    pz = sqrt(1.-norm_ab);

    while(true){
      double thickness =  46.;
      double width     = 100.;
      if     (imod==0){ thickness = 46.; width = 100.;}
      else if(imod==1){ thickness = 46.; width = 100.;}
      else if(imod==2){ thickness =107.; width = 120.;}
      else if(imod==3){ thickness = 80.; width = 120.;}
      else if(imod==4){ thickness =107.; width = 120.;}
      double tmpx = (G4UniformRand()-0.5)*width;
      double tmpz = (G4UniformRand()-0.5)*thickness;
      double tmpy;
      if(py==0.){tmpy=(G4UniformRand()-0.5)*100.;}
      else      {tmpy=-(tmpx*px+tmpz*pz)/py;}
      if(fabs(tmpy)<=width/2.){
        xpos_calc = tmpx - px*300.;
        ypos_calc = tmpy - py*300.;
        zpos_calc = tmpz - pz*300.;
        break;
      }
    }
    if(imod==0){
      xpos_calc += B2WM_OFFSETX;
      ypos_calc += B2WM_OFFSETY;
      zpos_calc += B2WM_OFFSETZ;
    }
    else if(imod==1){
      xpos_calc += PM_OFFSETX;
      ypos_calc += PM_OFFSETY;
      zpos_calc += PM_OFFSETZ;
    }
    else if(imod==2){
      xpos_calc += 0;
      ypos_calc += 0;
      zpos_calc += 0;
    }
    else if(imod==3){
      xpos_calc += B2CH_OFFSETX;
      ypos_calc += B2CH_OFFSETY;
      zpos_calc += B2CH_OFFSETZ;
    }
    else if(imod==4){
      xpos_calc += B2ING_OFFSETX;
      ypos_calc += B2ING_OFFSETY;
      zpos_calc += B2ING_OFFSETZ;
    }

#endif
#endif

    particleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
    G4float xpos = xpos_calc*cm;
    G4float ypos = ypos_calc*cm;
    G4float zpos = zpos_calc*cm;
    particleGun->SetParticlePosition(G4ThreeVector(xpos, ypos, zpos));
    particleGun->SetParticleTime(0.0*ns);

    SimVertexSummary* simvertex = new SimVertexSummary();
    simvertex -> Clear();
    simvertex -> nutype   = neutrino_flavor;
    simvertex -> inttype  = -100;
    simvertex -> nuE      = mumom/GeV;
    simvertex -> xnu      = xpos/cm;
    simvertex -> ynu      = ypos/cm;
    simvertex -> znu      = zpos/cm;
    simvertex -> gmomx.push_back(px);
    simvertex -> gmomy.push_back(py);
    simvertex -> gmomz.push_back(pz);
    //simvertex -> mod      = ID;
    simvertex -> norm			= 1.;
    simvertex -> totcrsne		= 1.;
    runaction  -> GetEvtSum() -> AddSimVertex( simvertex );

    particleGun->GeneratePrimaryVertex(anEvent);
    return;
  }

  //============================================================
  //=======     Particle-Gun: Muon/Pion/Proton-Injection             =======
  //============================================================
  if(neutrino_flavor==13||neutrino_flavor==14||neutrino_flavor==15){

    G4ParticleDefinition* particle;
    if     (neutrino_flavor==13){ particle = particleTable->FindParticle(13); }
    else if(neutrino_flavor==14){ particle = particleTable->FindParticle(211); }
    else if(neutrino_flavor==15){ particle = particleTable->FindParticle(2212); }
    particleGun->SetParticleDefinition(particle);
    G4float mass  = particle->GetPDGMass();

    G4float mumom = G4UniformRand()*2*GeV+0.1;
    G4float energy = sqrt(mumom*mumom+mass*mass) - mass;
    particleGun->SetParticleEnergy(energy);

    //Angle setting
    float theta_calc;
    float phi_calc;
    float xpos_calc;
    float ypos_calc;
    float zpos_calc;
    float px;
    float py;
    float pz;

    theta_calc = G4UniformRand()*PI;
    phi_calc   = G4UniformRand()*2.*PI;

    px = sin(theta_calc)*sin(phi_calc);
    py = sin(theta_calc)*cos(phi_calc);
    pz = cos(theta_calc);

    //Position setting
    double thickness, startz;
    double offset[3];
    while(true){
      int imod = (int)(G4UniformRand()*3.);
      if(imod==0){ 
        thickness =  20.5167; 
        startz    = -14.6067;
        offset[0] = B2WM_OFFSETX;
        offset[1] = B2WM_OFFSETY;
        offset[2] = B2WM_OFFSETZ;
        break;
      }
      else if(imod==1){ 
        thickness =  63.1; 
        startz    = -32.5;
        offset[0] = B2CH_OFFSETX;
        offset[1] = B2CH_OFFSETY;
        offset[2] = B2CH_OFFSETZ;
        break;
      }
      else if(imod==2){
        thickness =  73.9; 
        startz    = -41.8;
        offset[0] = B2ING_OFFSETX;
        offset[1] = B2ING_OFFSETY;
        offset[2] = B2ING_OFFSETZ;
        break;
      }
    }
    xpos_calc = (G4UniformRand()-0.5)*70.;
    ypos_calc = (G4UniformRand()-0.5)*70.;
    zpos_calc = G4UniformRand()*thickness + startz;
    xpos_calc += offset[0];
    ypos_calc += offset[1];
    zpos_calc += offset[2];

    particleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
    G4float xpos = xpos_calc*cm;
    G4float ypos = ypos_calc*cm;
    G4float zpos = zpos_calc*cm;
    particleGun->SetParticlePosition(G4ThreeVector(xpos, ypos, zpos));
    particleGun->SetParticleTime(0.0*ns);

    SimVertexSummary* simvertex = new SimVertexSummary();
    simvertex -> Clear();
    simvertex -> nutype   = neutrino_flavor;
    simvertex -> inttype  = -100;
    simvertex -> nuE      = mumom/GeV;
    simvertex -> xnu      = xpos/cm;
    simvertex -> ynu      = ypos/cm;
    simvertex -> znu      = zpos/cm;
    simvertex -> gmomx.push_back(px);
    simvertex -> gmomy.push_back(py);
    simvertex -> gmomz.push_back(pz);
    //simvertex -> mod      = ID;
    simvertex -> norm     = 1.;
    simvertex -> totcrsne = 1.;

    //Side angle
    if(direction[2]>=0){
      double alpha = acos( direction[2] / sqrt( pow(direction[1],2) + pow(direction[2],2) ) );
      simvertex -> thetaX = (-alpha)*180./pi;
    }
    else if(direction[2]<0){
      double alpha = acos( direction[2] / sqrt( pow(direction[1],2) + pow(direction[2],2) ) );
      simvertex -> thetaX = (pi-alpha)*180./pi;
    }

    //Top angle
    if(direction[2]>=0){
      double Theta = acos( direction[2] / sqrt( pow(direction[0],2) + pow(direction[2],2) ) );
      if(direction[0]>=0){
        simvertex -> thetaY = (Theta)*180./pi;
      }
      else if(direction[0]<0){
        simvertex -> thetaY = (-Theta)*180./pi;
      }
    }
    else if(direction[2]<0){
      double Theta = acos( direction[2] / sqrt( pow(direction[0],2) + pow(direction[2],2) ) );
      if(direction[0]>=0){
        simvertex -> thetaY = (-(pi-Theta))*180./pi;
      }
      else if(direction[0]<0){
        simvertex -> thetaY = (pi-Theta)*180./pi;
      }
    }

    //Angle along with X-axis
    double AngleX = 180./pi*(acos( direction[0] / sqrt( pow(direction[0],2) + pow(direction[1],2) + pow(direction[2],2) ) )); 
    simvertex -> angleX = AngleX;
    if(AngleX>90.){
      simvertex->angleX = 180.-AngleX;
    }

    //Angle along with Y-axis
    double AngleY = 180./pi*(acos( -direction[1] / sqrt( pow(direction[0],2) + pow(direction[1],2) + pow(direction[2],2) ) )); 
    simvertex -> angleY = AngleY;
    if(AngleY>90.){
      simvertex->angleY = 180.-AngleY;
    }

    //Angle along with Z-axis
    double AngleZ = 180./pi*(acos( direction[2] / sqrt( pow(direction[0],2) + pow(direction[1],2) + pow(direction[2],2) ) )); 
    simvertex -> angleZ = AngleZ;
    if(AngleZ>90.){
      simvertex->angleZ = 180.-AngleZ;
    }


    runaction  -> GetEvtSum() -> AddSimVertex( simvertex );    
    particleGun->GeneratePrimaryVertex(anEvent);
    return;
  }

  // Start loop of neut file
  while(1){

    mode = neut->NtupleReadEvent(Secondary,direction);
#ifdef DEBUG_PRINGEN
    G4cout << "mode: " << mode << G4endl;
#endif
    if(mode==-1111){
      G4cout <<"Aboart Run (mode =" << mode << G4endl;
      G4RunManager* runManager = G4RunManager::GetRunManager();
      eventaction->SetWriteFlag(-1); 
      runManager->AbortRun(1);

      return;
    }

    //Check neutrino flavor
    int neutrino_flavor_tmp = (int)(((neut->Vector).Neutrino.ProductionMode)/10);

    if(neutrino_flavor_tmp != neutrino_flavor){
#ifdef DEBUG_PRIMGEN
      G4cout << " === This neutrino id : " << neutrino_flavor_tmp
        << " ,not selected now === " << G4endl;
#endif
      continue;
    }

    //Define neutrino interaction vertex
    fdid = (neut->Vector).Neutrino.FDID;

    //X-Y vertex
    pos[0] = (neut->Vector).Neutrino.x;
    pos[1] = (neut->Vector).Neutrino.y;

#ifdef DEBUG_PRIMGEN
    G4cout << " ======= PrimaryGenerator: Vertex =======" << G4endl;
    G4cout << "fdid: " << fdid   << G4endl;
    G4cout << "X: " << pos[0] << G4endl;
    G4cout << "Y: " << pos[1] << G4endl;
#endif

    //======================================
    //==== Z-Vertex for INGRID Signal MC ===
    //====             and B2 INGRID MC  ===
    //======================================
    //if( fdid==3 || fdid==4 || (fdid==7&&(module_mode==B2ING))) 
    if( fdid==3 || fdid==4 || (fdid==9&&module_mode==B2ING) )
    {

      //Select vertex in Fe or scinti 
      //double ratio = total_mass_sci / (total_mass_fe + total_mass_sci);

      //if(ratio > G4UniformRand()){
      //  vertex_flag = 1; // flag = 1 -> vertex in scinti
      //}
      //else{
      //  vertex_flag = 0;
      //}
      vertex_flag=0; //iron only

      //INGRID mod#3
      if( fdid==3){
        if(vertex_flag==0){
          flayer = (int)(9*G4UniformRand());       // 0 < 9*rand < 9
          while(flayer<0||flayer>8){ flayer = (int)(9*G4UniformRand()); }
          pos[2] = width_fe*(G4UniformRand()-0.5); // -6.5/2. -- 6.5/2.cm   
          pos[2] = pos[2] + iron_start03 + GAP03*flayer;    
        }
        else if(vertex_flag==1){
          flayer = (int)(11*G4UniformRand());         // 0 < 11*rand < 11
          while(flayer<0||flayer>10){ flayer = (int)(11*G4UniformRand()); }

          pos[2] = width_sci*2*(G4UniformRand()-0.5); // -1 -- 1cm
          //For a pair of X-Y layer
          pos[2] = pos[2] + scinti_start03 + GAP03*flayer + width_sci*0.5;
        }
      }
      //INGRID modules (except for mod#3)
      else{
        if(vertex_flag==0){//iron
          flayer = (int)(9*G4UniformRand());       // 0 < 9*rand < 9
          while(flayer<0||flayer>8){ flayer = (int)(9*G4UniformRand()); }
          pos[2] = width_fe*(G4UniformRand()-0.5); // -6.5/2. -- 6.5/2.cm
          pos[2] = pos[2] + iron_start + GAP*flayer;    
        }
        else if(vertex_flag==1){//scinti
          flayer = (int)(11*G4UniformRand());       // 0 < 9*rand < 9
          while(flayer<0||flayer>10){ flayer = (int)(11*G4UniformRand()); }
          pos[2] = width_sci*2*(G4UniformRand()-0.5); // -1. -- 1.cm 
          //For a pair of X-Y layer
          pos[2] = pos[2] + scinti_start + GAP*flayer + width_sci*0.5;
        }
      }
    }

    //========================================================
    //=== Z-Vertex for ProtonModule (On-axis/B2) Signal MC ===
    //========================================================
    //else if( (fdid==2&&module_mode==PROTON) || ((fdid==7||fdid==8) && module_mode==B2CH))
    else if( 
        ( fdid==2 && module_mode==PROTON) || 
        ( fdid==8 && module_mode==B2CH  ) )
    {

      //double frontratio =  total_mass_front_pm / (total_mass_front_pm + total_mass_sci_pm);
      double frontratio =  total_mass_front_pm / total_mass_sci_pm;

      //Vertex_flag: 0->front vetor
      //             1-> middle part (Scibar+Scibar)
      //             2-> X-middle,Y-side (INGRID+Scibar)
      //             3-> X-side,Y-middle (INGRID+Scibar)
      //             4-> side region (INGRID+INGRID)
      // Normalized by thickness of Scibar scintillator
      // ==> Normalized by INGRIDScinti x2 + SciBarScinti x2 x17pln
      if( frontratio > (G4UniformRand()) ){
        //vertex_flag=0; prob=sciing_region/scibar_region;
        vertex_flag=0; prob=1.;
      }
      else if(fabs(pos[0])<=20&&fabs(pos[1])<=20){
        vertex_flag=1; prob=1.;
      }
      else if(fabs(pos[0])<=20){
        vertex_flag=2; prob=sciing_region/scibar_region;
      }      
      else if(fabs(pos[1])<=20){
        vertex_flag=3; prob=sciing_region/scibar_region;
      }      
      else{
        vertex_flag=4; prob=ingrid_region/scibar_region;
      }

      if(prob<(G4UniformRand())){
        continue;
      }

      //Pos: [cm]
      if(vertex_flag==0){
        flayer = (int)(2*(G4UniformRand()));         // 0 < 2*rand < 2
        while(flayer<0||flayer>1){ flayer = (int)(2*(G4UniformRand()));}
        pos[2] = ingrid_width*(G4UniformRand()-0.5); //-1/2. -- 1/2.cm
        pos[2] = pos[2] + Pscinti_start + flayer*distance_pln;    
      }
      else if(vertex_flag==1) {
        flayer = (int)(34*(G4UniformRand()));        // 0 < 34*rand < 34
        while(flayer<0||flayer>33){ flayer = (int)(34*(G4UniformRand()));}
        pos[2] = scibar_width*(G4UniformRand()-0.5); // -1.3/2. -- 1.3/2.cm
        pos[2] = pos[2] + Pscinti_start + distance_first + distance_pln*(flayer+1);
      }
      else if(vertex_flag==2){
        scitype= (int)((G4UniformRand())/scibar_width*(ingrid_width+scibar_width));
        while(scitype<0||scitype>1){
          scitype= (int)((G4UniformRand())/scibar_width*(ingrid_width+scibar_width));
        }
        flayer = (int)(17*(G4UniformRand())); // 0 < 17*rand < 17
        while(flayer<0||flayer>16){ flayer = (int)(17*(G4UniformRand()));}
        if(scitype==0){
          pos[2] = scibar_width*(G4UniformRand()-0.5); 
          pos[2] = pos[2] + Pscinti_start + distance_first + distance_pln*(flayer*2+2);
        }
        else{
          pos[2] = ingrid_width*(G4UniformRand()-0.5);
          pos[2] = pos[2] + Pscinti_start + distance_first + distance_pln*(flayer*2+1);
        }
      }
      else if(vertex_flag==3){
        scitype= (int)((G4UniformRand())/scibar_width*(ingrid_width+scibar_width));
        while(scitype<0||scitype>1){
          scitype= (int)((G4UniformRand())/scibar_width*(ingrid_width+scibar_width));
        }
        flayer = (int)(17*(G4UniformRand())); // 0 < 17*rand < 17
        while(flayer<0||flayer>16){ flayer = (int)(17*(G4UniformRand()));}
        if(scitype==0){
          pos[2] = scibar_width*(G4UniformRand()-0.5);
          pos[2] = pos[2] + Pscinti_start + distance_first + distance_pln*(flayer*2+1);
        }
        else{
          pos[2] = ingrid_width*(G4UniformRand()-0.5);
          pos[2] = pos[2] + Pscinti_start + distance_first + distance_pln*(flayer*2+2);
        }
      }
      else if(vertex_flag==4){
        flayer = (int)(34*(G4UniformRand())); // 0 < 34*rand < 34
        while(flayer<0||flayer>33){ flayer = (int)(34*(G4UniformRand()));}
        pos[2] = ingrid_width*(G4UniformRand()-0.5);
        pos[2] = pos[2] + Pscinti_start + distance_first + distance_pln*(flayer+1);
      }
    }

    //===================================================
    //==== Z-Vertex for On-Axis WaterModule Signal MC ===
    //===================================================
    else if(fdid==2 && module_mode==INGWaterModule)
    {
      pos[2]=watertank_z*(G4UniformRand()-0.5); 
    }
    else if(fdid==3 && module_mode==INGWaterModuleBG)
    {
      pos[2]=altank_z*(G4UniformRand()-0.5); // -60.6/2. --- 60.6/2.cm
    }

    //===================================================
    //==== Z-Vertex for B2 Modules Signal MC ============
    //===================================================
    else if(fdid==7 && module_mode==B2Water)
    {
      pos[2]=watertank_z*(G4UniformRand()-0.5);
    }

    //===============================================
    //==== Z-Vertex for B2 Background MC ============
    //===============================================
    if((fdid==1||fdid==7) && module_mode==B2Wall){
      if (fabs(pos[1])<980){ //Vertex_x  > 0 along Neut X axis
        G4double lx = fabs(pos[0] + 322.2);    //lx = distance from hall center along x axis 

        if(lx<HallRadiusMin){
          dz = sqrt((HallRadiusMin+500)*(HallRadiusMin+500) - lx*lx) - sqrt(HallRadiusMin*HallRadiusMin - lx*lx);
          pos[2] = 170 - sqrt(HallRadiusMin*HallRadiusMin - lx*lx) - G4UniformRand()*dz;
          //G4cout << "pos[2] (up) : " << pos[2] << G4endl;
        }
        if((HallRadiusMin <= lx) && (lx < HallRadiusMin +500)){	
          dz = sqrt( (HallRadiusMin+500)*(HallRadiusMin+500) - lx*lx );
          pos[2] = 170 - G4UniformRand()*dz;
          //G4cout << "pos[2] (down) : " << pos[2] << G4endl;
        }
      }
      else continue;
    }

    //----- Right-pillar BG ----- (Kin 2017/11/08)
    if(fdid==6 &&  module_mode == B2Right_Pillar){
      pos[2] = (G4UniformRand())*400.0 -200;  // cm   
      // set for dz = -200 ~ 200,  is length of B2Pillar along z axis throughly, later hamidashi cut
    }                                               

    //----- Left-pillar BG ----- (Kin 2017/11/08)
    if(fdid==6 &&  module_mode == B2Left_Pillar){
      if( (250.5<pos[0]) && (pos[0]<350.5) && (pos[1]>-112.4) && (pos[1]<187.6)  ){  // Vertex_x  > 0 along Neut X axis
        pos[2] = (G4UniformRand())*713.3 -513.3;              // cm
        // set for dz = -513.3 ~ 200,  is length of B2Pillar along z axis throughly, later hamidashi cut
        // dz = -513.3 means left mini pillar edge.
      }
      else continue;    
    }                                                

    //----- Ceiling BG ----- (Kin 2017/11/08)
    if(fdid==5 && module_mode==B2Ceiling ){
      pos[2] = 200.0 - (G4UniformRand())*(B2WM_OFFSETZ - 170 + HallRadiusMin);   
      // set for dz (reference point is B2INGRID ) 
      //       = 200 ~ to halledge( -(B2INGRIDz + hallradius_min - 170cm) ),
      //       ( 170cm is Z_offset from hall center in Geant coordinate  later hamidashi space cut )
    }

    //----- Floor BG ----- (Kin 2017/11/08)
    if(fdid==5 && module_mode==B2Floor){
      pos[2] = 200.0 - (G4UniformRand())*(B2WM_OFFSETZ + HallRadiusMin - 170);   
      // set for dz (reference point is B2INGRID ) 
      //       = 200 ~ to halledge( -(B2INGRIDz + hallradius_min - 170cm) ),
      //       ( 170cm is Z_offset from hall center in Geant coordinate  later hamidashi space cut )
    }




#ifdef DEBUG_PRIMGEN
    G4cout << "Z:" << pos[2] << G4endl;
#endif

    //Define module ID
    ID = -1;
    if(fdid==7)
    {    
      if(module_mode==B2Water){
        if( fabs(pos[0]) <= watertank_x/2. &&
            fabs(pos[1]) <= watertank_y/2. ) {
          ID = 21;
          goto NEXTSTEP;
        }
      }
    }
    else if(
        (fdid==2 && module_mode == PROTON) ||
        (fdid==8 && module_mode == B2CH  ) )
    {
      if( fabs(pos[0]) <= scintileng_pm/2. &&
          fabs(pos[1]) <= scintileng_pm/2. ) {
        ID = 16;
        goto NEXTSTEP;
      }
    }
    else if(fdid==9 && module_mode == B2ING){
      if( fabs(pos[0]) <= scintileng_ing/2. &&
          fabs(pos[1]) <= scintileng_ing/2. ) 
      {
        ID = 14;
        goto NEXTSTEP;
      }
    }
    else if(fdid==1 && module_mode==B2Wall){
      if(fabs(pos[1])<=980) {
	ID = 31;
	goto NEXTSTEP;
      }
    }
    else if(fdid==3 && module_mode == HORIZONTAL) {
      for( int m=0;m<7;m++ ) {
        if( fabs(pos[0]+INGRID_MODSpace*(3-m)) <= scintileng_ing/2. &&
            fabs(pos[1]) <= scintileng_ing/2. ) {
          ID = m;
          goto NEXTSTEP;
        }
      }
    }
    else if(fdid==4 && module_mode == VERTICAL) {
      for( int m=0;m<7;m++ ) {
        if( fabs(pos[0]) <= scintileng_ing/2. &&
            fabs(pos[1]+INGRID_MODSpace*(3-m)) <= scintileng_ing/2. ) {
          ID = m+7;
          goto NEXTSTEP;
        }
      }
    }
    else if(fdid==2 && module_mode == INGWaterModule) {
      if( fabs(pos[0]) <= watertank_x/2. &&
          fabs(pos[1]) <= watertank_y/2. ) {
        ID = 15;
        goto NEXTSTEP;
      }
    }
    else if(fdid==5 || fdid==6){
      if(module_mode==B2Ceiling){
	//G4double lx = fabs(pos[0] + C_B2CeilingPosX);
	//G4double lz = fabs(pos[2] + C_B2CeilingPosZ);
	G4double lx = fabs(pos[0] - 322.2);
	G4double lz = fabs(pos[2] + 170.);
	//if(lx*lx + lz*lz < 1000*1000 && fabs(pos[1]) <= 50)
	if(fabs(pos[1]) <= 50)
        {
	  ID = 32;
	  goto NEXTSTEP;
	}
      }
      else if(module_mode==B2Floor){
	//G4double lx = fabs(pos[0] + C_B2CeilingPosX);
	G4double lx = fabs(pos[0] - 322.2);
	//if( (lx*lx + pos[2]*pos[2] < HallRadiusMin*HallRadiusMin ) && (-C_HallDirtHeight < pos[1] + C_B2FloorPosY) && (pos[1] + C_B2FloorPosY < -594.343) )
	if(pos[1]>=-300 && pos[1]<=0)
        {	
	  ID = 33;
	  goto NEXTSTEP;
	}
      }
      else if(module_mode==B2Right_Pillar){
	if( fabs(pos[0]) <= 50 && fabs(pos[1]) <= 200 ) {
	  ID = 34;
	  goto NEXTSTEP;
	}
      }
      else if(module_mode==B2Left_Pillar){
	//if( ((pos[2]>-200)&&(pos[2]<200)) || ((pos[2]>-513.3)&&(pos[2]<-413.3)))
	if(fabs(pos[1]) <= 200) 
        {
	  ID = 35;
	  goto NEXTSTEP;
	}
      }
    }



    //Count events which have vertex out of modules
#ifdef DEBUG_PRIMGEN
    G4cout << "##### Interaction vertex is out of modules #####" << G4endl;
    G4cout << "##### Skip this event                      #####" << G4endl;
#endif
    runaction->NotEntry++; 

  }//End while loop


  //--------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------

NEXTSTEP:

  //Offset for surface of each module [cm]
  ////V-INGRID
  if(ID > 6 && ID < 14){
    pos[2] = pos[2] + INGRID_VMODOFFSET;
  }

  ////On-axis Proton Module
  //else if(ID == 16)  pos[2] = pos[2] + PM_OFFSET;

  ////INGRIDWaterModule
  //else if(ID == 20)  pos[2] = pos[2] + WM_OFFSET;
  else if(ID == 15)  pos[2] = pos[2] + WM_OFFSET;

  ////BG for INGRIDWaterModule
  else if(ID == 18)  pos[2] = pos[2] + WMBG_OFFSET;

  ////WAGASCI water module
  else if(ID == 21){
    pos[0] = pos[0] + B2WM_OFFSETX;
    pos[1] = pos[1] + B2WM_OFFSETY;
    pos[2] = pos[2] + B2WM_OFFSETZ;
  }

  ////B2 Proton Module
  else if(ID == 16){
    pos[0] = pos[0] + B2CH_OFFSETX;
    pos[1] = pos[1] + B2CH_OFFSETY;
    pos[2] = pos[2] + B2CH_OFFSETZ;
  }

  ////B2-INGRID
  else if(ID == 14){
    pos[0] = pos[0] + B2ING_OFFSETX;
    pos[1] = pos[1] + B2ING_OFFSETY;
    pos[2] = pos[2] + B2ING_OFFSETZ;
  }

  //Upstream-Wall BG
  else if(ID==31){
    pos[0] = pos[0];
    pos[1] = pos[1] + 17.38;
    pos[2] = pos[2];
  }

  //Ceiling BG
  else if(ID==32){
    pos[0] = pos[0] + B2WM_OFFSETX + 250.;
    pos[1] = pos[1] - 157.3;
    //pos[2] = pos[2] + B2WM_OFFSETZ - 200;
    pos[2] = pos[2] + B2WM_OFFSETZ;
  }

  //Floor BG
  else if(ID==33){
    pos[0] = pos[0] + B2WM_OFFSETX + 100.;
    //pos[1] = pos[1] - 797.17;
    pos[1] = pos[1] + B2WM_OFFSETY - 100.;
    pos[2] = pos[2] + B2WM_OFFSETZ;
  }

  //Right pillar BG
  else if(ID==34){
    pos[0] = pos[0] - 823.5;
    pos[1] = pos[1] - 400.;
    pos[2] = pos[2] + 233.5;
  }

  //Left pillar BG
  else if(ID==35){
    pos[0] = pos[0] - 624.;
    pos[1] = pos[1] - 400.;
    pos[2] = pos[2] + 233.5;
  }


  // Input Neut file info to output ROOT class
  neut->ID = ID;
  for(int i=0;i<3;i++) (runaction->vertex)[i] = pos[i];

  SimVertexSummary* simvertex = new SimVertexSummary();
  simvertex -> Clear();
  simvertex -> nutype   = neutrino_flavor;
  simvertex -> inttype  = (neut->Vector).Primary.Mode;
  simvertex -> nuE      = (neut->Vector).Neutrino.Energy;
  simvertex -> xnu      = pos[0];
  simvertex -> ynu      = pos[1];
  simvertex -> znu      = pos[2];
  simvertex -> mod      = ID;
  simvertex -> norm     =  (neut->Vector).Neutrino.Norm;
  simvertex -> totcrsne	= (neut->Vector).neutcrs.Totcrsne;

  //for Al density
  if(ID==18){
    simvertex -> norm			= (neut->Vector).Neutrino.Norm * 2.7;
  }

  runaction  -> GetEvtSum() -> AddSimVertex( simvertex );

  G4cout.precision( 3 );

#ifdef DEBUG_PRIMGEN
  G4cout << "\n=== Neutrino Information from Jnubeam ===" << G4endl;
  G4cout << "Norm: " <<  (neut->Vector).Neutrino.Norm << G4endl;
  G4cout << "Totcrsne: " <<  (neut->Vector).neutcrs.Totcrsne << G4endl;
  G4cout << "ParentID: " << (neut->Vector).Neutrino.ParentID;
  G4cout << "  Neut Production Mode: " << (neut->Vector).Neutrino.ProductionMode;
  G4cout << "  Neutrino.FDID: " << (neut->Vector).Neutrino.FDID << G4endl;
  G4cout << "Neut interaction Mode: " << (neut->Vector).Primary.Mode << G4endl;
  G4cout << "Energy[GeV]: " << (neut->Vector).Neutrino.Energy;
  G4cout << "  Direction: {" << direction[0] << "," << direction[1] << "," << direction[2] << "}" << G4endl;
  G4cout << "Vertex(cm): {" << pos[0] << ", "<< pos[1] << ", "<< pos[2] << "}";
  G4cout << "  Module: " << ID << "\n\n";
#endif


#ifndef SUPPRESS_FOR_EXT_BG
  // Input Neut info for T2KReWeight to SK__h1 class
  runaction -> numnu = (neut->Vector).Primary.NumParticle;
  runaction -> mode  = (neut->Vector).Primary.Mode;
  for ( int i = 0; i<50; i++ ) {
    runaction -> ipnu[i] = (neut->Vector).Primary.ParticleID[i];
    runaction -> pnu[i]  = (neut->Vector).Primary.AbsMomentum[i];
    for ( int j = 0 ; j < 3 ; j++ ){
      runaction -> dirnu[i][j] = (neut->Vector).Primary.Momentum[i][j] / (neut->Vector).Primary.AbsMomentum[i];
    }
  }

  runaction -> Crsx   = (neut->Vector).Crs.Crsx;
  runaction -> Crsy   = (neut->Vector).Crs.Crsy;
  runaction -> Crsz   = (neut->Vector).Crs.Crsz;
  runaction -> Crsphi = (neut->Vector).Crs.Crsphi;

  runaction -> Nvert = (neut->Vector).Fsi.Nvert;
  for (int ivert=0; ivert<150; ivert++) {
    runaction -> Iflgvert[ivert] = (neut->Vector).Fsi.Iflgvert[ivert];
    for (int j=0; j<3; j++)
      runaction -> Posvert[ivert][j] = (neut->Vector).Fsi.Posvert[ivert][j];
  }

  runaction -> Nvcvert = (neut->Vector).Fsi.Nvcvert;
  for (int ip=0; ip<900; ip++) {
    runaction -> Abspvert[ip]  = (neut->Vector).Fsi.Abspvert[ip];
    runaction -> Abstpvert[ip] = (neut->Vector).Fsi.Abstpvert[ip];
    runaction -> Ipvert[ip]    = (neut->Vector).Fsi.Ipvert[ip];
    runaction -> Iverti[ip]    = (neut->Vector).Fsi.Iverti[ip];
    runaction -> Ivertf[ip]    = (neut->Vector).Fsi.Ivertf[ip];
    for (int j=0; j<3; j++)
      runaction -> Dirvert[ip][j] = (neut->Vector).Fsi.Dirvert[ip][j];
  }
  runaction -> Fsiprob = (neut->Vector).Fsi.Fsiprob;
  runaction -> Numbndn = (neut->Vector).target_info.Numbndn;
  runaction -> Numbndp = (neut->Vector).target_info.Numbndp;
  runaction -> Numfrep = (neut->Vector).target_info.Numfrep;
  runaction -> Numatom = (neut->Vector).target_info.Numatom;
  runaction -> Ibound  = (neut->Vector).Fsi.Ibound;
  runaction -> Npvc    = (neut->Vector).Secondary.NumParticle;
  for (int i=0; i<100; i++) {
    runaction -> Ipvc[i]    = (neut->Vector).Secondary.ParticleID[i];
    runaction -> Ichvc[i]   = (neut->Vector).Secondary.TrackingFlag[i];
    runaction -> Iorgvc[i]  = (neut->Vector).Secondary.ParentID[i];
    runaction -> Iflvc[i]   = (neut->Vector).Secondary.InteractionCode[i];
    runaction -> Abspvc[i]  = (neut->Vector).Secondary.AbsMomentum[i];
    for (int j=0; j<3; j++)
      runaction -> Pvc[i][j]     = (neut->Vector).Secondary.Momentum[i][j];
  }
#endif

  for(int ipart=0; ipart<Secondary.NumParticle; ipart++) {
    if( Secondary.TrackingFlag[ipart]==1 ) {
#ifdef DEBUG_PRIMGEN
      G4cout << "Particle:" << (neut->Vector).Secondary.ParticleID[ipart];
      G4cout << "  Index: " << ipart;
      G4cout << "  Parent Index: " << (neut->Vector).Secondary.ParentID[ipart] -1 << "\n";
      G4cout << "Tracking Flag: " << (neut->Vector).Secondary.TrackingFlag[ipart];
      G4cout << "  Interaction code: " << (neut->Vector).Secondary.InteractionCode[ipart] << "\n";
      G4cout << " Momentum[MeV/c]:";
      for (int k=0;k<3;k++)   G4cout << (neut->Vector).Secondary.Momentum[ipart][k]*MeV << " ";
      G4cout << "\n";
#endif

      // Set ParticleGun with information above ------------------------------------------


      G4ParticleDefinition* particle;
      particle = particleTable->FindParticle(Secondary.ParticleID[ipart]);

      double nvec[3];
      for(int ixyz=0; ixyz<3; ixyz++)
        nvec[ixyz] = Secondary.Momentum[ipart][ixyz]/ Secondary.AbsMomentum[ipart];
      G4ThreeVector dir(nvec[0], nvec[1], nvec[2]);

      G4double mass = particle->GetPDGMass();
      G4double mom = Secondary.AbsMomentum[ipart]*MeV;
      G4double energy = sqrt(mass*mass+mom*mom) - mass;

      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(pos[0]*cm,pos[1]*cm,pos[2]*cm));
      particleGun->SetParticleMomentumDirection(dir);
      particleGun->SetParticleEnergy(energy);
      particleGun->SetParticleTime(0.0*ns);
      particleGun->GeneratePrimaryVertex(anEvent);

#ifdef DEBUG_PRINGEN
      G4cout << "ipart: " << ipart << "\n";
      G4cout << "PID:" << (neut->Vector).Secondary.ParticleID[ipart] << "\n";
      G4cout << "Tracking Flag: " << (neut->Vector).Secondary.TrackingFlag[ipart] << "\n";
      G4cout << "  Kinetic Energy[MeV]: " << energy << G4endl;;
      G4cout << "  Momentum:";
      for (int k=0;k<3;k++)   G4cout << (neut->Vector).Secondary.Momentum[ipart][k]*MeV << " ";
      G4cout << " [MeV]\n";
#endif

    } // end of if
  } // end of for loop
}
