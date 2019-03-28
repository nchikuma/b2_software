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

////#define DEBUG_PRIMGEN

PrimaryGeneratorAction::
  //PrimaryGeneratorAction(Neut *neut0,RunAction* rac, EventAction* evt,int nd,int flavor0)
  PrimaryGeneratorAction(RunAction* rac, EventAction* evt,int nd,int flavor0)
:runaction(rac)
{
  eventaction      = evt;
  //neut_file        = neut0;
  genie_tree       = runaction->GetGenieTree();
  module_mode      = nd;
  neutrino_flavor  = flavor0;
  particleGun      = new G4ParticleGun(1);
  particleTable    = G4ParticleTable::GetParticleTable();
  runaction->NotEntry = 0;

  iEvt   = -1;
  NumEvt = genie_tree->GetEntries();
  SetBranchGENIE(genie_tree,genie_para);

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
  //Neut   *neut        = neut_file;
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


  // Start loop of neut file
  while(1){
    iEvt++;

    if(iEvt%100==0){
      cout << "event : " << iEvt << endl;
    }

    if(iEvt>=NumEvt){
      G4cout <<"Aboart Run (mode =" << mode << G4endl;
      G4RunManager* runManager = G4RunManager::GetRunManager();
      eventaction->SetWriteFlag(-1); 
      runManager->AbortRun(1);
      return;
    }


    genie_tree->GetEntry(iEvt);

    //Check neutrino flavor
    //int neutrino_flavor_tmp = (int)(((neut->Vector).Neutrino.ProductionMode)/10);

    //Define neutrino interaction vertex
    //fdid = (neut->Vector).Neutrino.FDID;
    fdid = module_mode;

    //X-Y vertex
    //pos[0] = (neut->Vector).Neutrino.x;
    //pos[1] = (neut->Vector).Neutrino.y;
    pos[0] = genie_para.VtxX*100.;
    pos[1] = genie_para.VtxY*100.;

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
    if( fdid==3 || fdid==4 || (fdid==9&&module_mode==B2ING) )
    {
      vertex_flag=0; //iron only
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

    //========================================================
    //=== Z-Vertex for ProtonModule (On-axis/B2) Signal MC ===
    //========================================================
    //else if( (fdid==2&&module_mode==PROTON) || ((fdid==7||fdid==8) && module_mode==B2CH))
    else if( 
        ( fdid==2 && module_mode==PROTON) || 
        ( fdid==8 && module_mode==B2CH  ) )
    {
      double frontratio =  total_mass_front_pm / total_mass_sci_pm;
      if     ( frontratio > (G4UniformRand())   ){ vertex_flag=0; prob=1.; }
      else if(fabs(pos[0])<=20&&fabs(pos[1])<=20){ vertex_flag=1; prob=1.; }
      else if(fabs(pos[0])<=20                  ){ vertex_flag=2; prob=sciing_region/scibar_region; }      
      else if(fabs(pos[1])<=20                  ){ vertex_flag=3; prob=sciing_region/scibar_region; }      
      else                                       { vertex_flag=4; prob=ingrid_region/scibar_region; }

      if(prob<(G4UniformRand())){ continue; }

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
    else if(fdid==2 && module_mode==INGWaterModule  ){ pos[2]=watertank_z*(G4UniformRand()-0.5); }
    else if(fdid==3 && module_mode==INGWaterModuleBG){ pos[2]=altank_z*(G4UniformRand()-0.5);    }

    //===================================================
    //==== Z-Vertex for B2 Modules Signal MC ============
    //===================================================
    else if(fdid==7 && module_mode==B2Water) { pos[2]=watertank_z*(G4UniformRand()-0.5); }

    //===============================================
    //==== Z-Vertex for B2 Background MC ============
    //===============================================
    if((fdid==1||fdid==7) && module_mode==B2Wall){
      if (fabs(pos[1])<980){ 
        G4double lx = fabs(pos[0] + 322.2);   
        if(lx<HallRadiusMin){
          dz = sqrt((HallRadiusMin+500)*(HallRadiusMin+500) - lx*lx) - sqrt(HallRadiusMin*HallRadiusMin - lx*lx);
          pos[2] = 170 - sqrt(HallRadiusMin*HallRadiusMin - lx*lx) - G4UniformRand()*dz;
        }
        if((HallRadiusMin <= lx) && (lx < HallRadiusMin +500)){	
          dz = sqrt( (HallRadiusMin+500)*(HallRadiusMin+500) - lx*lx );
          pos[2] = 170 - G4UniformRand()*dz;
        }
      }
      else continue;
    }

    //----- Right-pillar BG -----
    if(fdid==6 &&  module_mode == B2Right_Pillar){ pos[2] = (G4UniformRand())*400.0 -200; }
    //----- Left-pillar BG -----
    if(fdid==6 &&  module_mode == B2Left_Pillar){
      if( (250.5<pos[0]) && (pos[0]<350.5) && (pos[1]>-112.4) && (pos[1]<187.6)  ){ 
        pos[2] = (G4UniformRand())*713.3 -513.3;
      }
      else continue;    
    }                                                
    //----- Ceiling BG -----
    if(fdid==5 && module_mode==B2Ceiling ){ pos[2] = 200.0 - (G4UniformRand())*(B2WM_OFFSETZ - 170 + HallRadiusMin); }
    //----- Floor BG -----
    if(fdid==5 && module_mode==B2Floor){ pos[2] = 200.0 - (G4UniformRand())*(B2WM_OFFSETZ + HallRadiusMin - 170);   }



#ifdef DEBUG_PRIMGEN
    G4cout << "Z:" << pos[2] << G4endl;
#endif

    //Define module ID
    ID = -1;
    if(fdid==7)
    {    
      if(module_mode==B2Water){
        if( fabs(pos[0]) <= watertank_x/2. &&
            fabs(pos[1]) <= watertank_y/2. ) { ID = 21; goto NEXTSTEP; }
      }
    }
    else if(
        (fdid==2 && module_mode == PROTON) ||
        (fdid==8 && module_mode == B2CH  ) )
    {
      if( fabs(pos[0]) <= scintileng_pm/2. &&
          fabs(pos[1]) <= scintileng_pm/2. ) { ID = 16; goto NEXTSTEP; }
    }
    else if(fdid==9 && module_mode == B2ING){
      if( fabs(pos[0]) <= scintileng_ing/2. &&
          fabs(pos[1]) <= scintileng_ing/2. ) { ID = 14; goto NEXTSTEP; }
    }
    else if(fdid==1 && module_mode==B2Wall){
      if(pos[0]>400.){continue;} 
      if(pos[1]>400.||pos[1]<-800.){continue;} 
      if(fabs(pos[1])<=980) { ID = 31; goto NEXTSTEP; }
    }
    else if(fdid==3 && module_mode == HORIZONTAL) {
      for( int m=0;m<7;m++ ) {
        if( fabs(pos[0]+INGRID_MODSpace*(3-m)) <= scintileng_ing/2. &&
            fabs(pos[1]) <= scintileng_ing/2. ) { ID = m; goto NEXTSTEP; }
      }
    }
    else if(fdid==4 && module_mode == VERTICAL) {
      for( int m=0;m<7;m++ ) {
        if( fabs(pos[0]) <= scintileng_ing/2. &&
            fabs(pos[1]+INGRID_MODSpace*(3-m)) <= scintileng_ing/2. ) { ID = m+7; goto NEXTSTEP; }
      }
    }
    else if(fdid==2 && module_mode == INGWaterModule) {
      if( fabs(pos[0]) <= watertank_x/2. &&
          fabs(pos[1]) <= watertank_y/2. ) { ID = 15; goto NEXTSTEP; }
    }
    else if(fdid==5 || fdid==6){
      if(module_mode==B2Ceiling){
	G4double lx = fabs(pos[0] - 322.2);
	G4double lz = fabs(pos[2] + 170.);
	if(fabs(pos[1]) <= 50) { ID = 32; goto NEXTSTEP; }
      }
      else if(module_mode==B2Floor){
	G4double lx = fabs(pos[0] - 322.2);
	if(pos[1]>=-300 && pos[1]<=0) {	ID = 33; goto NEXTSTEP; }
      }
      else if(module_mode==B2Right_Pillar){
	if( fabs(pos[0]) <= 50 && fabs(pos[1]) <= 200 ) { ID = 34; goto NEXTSTEP; }
      }
      else if(module_mode==B2Left_Pillar){
	if(fabs(pos[1]) <= 200) { ID = 35; goto NEXTSTEP; }
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
  if(ID > 6 && ID < 14){ pos[2] = pos[2] + INGRID_VMODOFFSET; }
  ////INGRIDWaterModule
  else if(ID == 15){ pos[2] = pos[2] + WM_OFFSET; }
  ////BG for INGRIDWaterModule
  else if(ID == 18){ pos[2] = pos[2] + WMBG_OFFSET; }
  ////WAGASCI water module
  else if(ID == 21){ pos[0] = pos[0] + B2WM_OFFSETX; pos[1] = pos[1] + B2WM_OFFSETY; pos[2] = pos[2] + B2WM_OFFSETZ; }
  ////B2 Proton Module
  else if(ID == 16){ pos[0] = pos[0] + B2CH_OFFSETX; pos[1] = pos[1] + B2CH_OFFSETY; pos[2] = pos[2] + B2CH_OFFSETZ; }
  ////B2-INGRID
  else if(ID == 14){ pos[0] = pos[0] + B2ING_OFFSETX; pos[1] = pos[1] + B2ING_OFFSETY; pos[2] = pos[2] + B2ING_OFFSETZ; }
  //Upstream-Wall BG
  else if(ID==31){ pos[0] = pos[0]; pos[1] = pos[1] + 17.38; pos[2] = pos[2]; }
  //Ceiling BG
  else if(ID==32){ pos[0] = pos[0] + B2WM_OFFSETX + 250.; pos[1] = pos[1] - 157.3; pos[2] = pos[2] + B2WM_OFFSETZ; }
  //Floor BG
  else if(ID==33){ pos[0] = pos[0] + B2WM_OFFSETX + 100.; pos[1] = pos[1] + B2WM_OFFSETY - 100.; pos[2] = pos[2] + B2WM_OFFSETZ; }
  //Right pillar BG
  else if(ID==34){ pos[0] = pos[0] - 823.5; pos[1] = pos[1] - 400.; pos[2] = pos[2] + 233.5; }
  //Left pillar BG
  else if(ID==35){ pos[0] = pos[0] - 624.; pos[1] = pos[1] - 400.; pos[2] = pos[2] + 233.5; }


  // Input Neut file info to output ROOT class
  //neut->ID = ID;
  for(int i=0;i<3;i++) (runaction->vertex)[i] = pos[i];

  SimVertexSummary* simvertex = new SimVertexSummary();
  simvertex -> Clear();
  simvertex -> nutype   = neutrino_flavor;
  //simvertex -> inttype  = (neut->Vector).Primary.Mode;
  if     (genie_para.CodeNeut!=0                     ){ simvertex -> inttype = genie_para.CodeNeut;  }
  else if(genie_para.Neutrino== 14&&genie_para.IsMec){ simvertex -> inttype =  2;                    }
  else if(genie_para.Neutrino==-14&&genie_para.IsMec){ simvertex -> inttype = -2;                    }
  //simvertex -> nuE      = (neut->Vector).Neutrino.Energy;
  simvertex -> nuE      = genie_para.Ev;
  simvertex -> xnu      = pos[0];
  simvertex -> ynu      = pos[1];
  simvertex -> znu      = pos[2];
  simvertex -> mod      = ID;
  //simvertex -> norm     = (neut->Vector).Neutrino.Norm;
  simvertex -> norm     = 1.;
  //simvertex -> totcrsne	= (neut->Vector).neutcrs.Totcrsne;
  simvertex -> totcrsne	= 1.;


  runaction  -> GetEvtSum() -> AddSimVertex( simvertex );

  G4cout.precision( 3 );



//#ifndef SUPPRESS_FOR_EXT_BG
//  // Input Neut info for T2KReWeight to SK__h1 class
//  runaction -> numnu = (neut->Vector).Primary.NumParticle;
//  runaction -> mode  = (neut->Vector).Primary.Mode;
//  for ( int i = 0; i<50; i++ ) {
//    runaction -> ipnu[i] = (neut->Vector).Primary.ParticleID[i];
//    runaction -> pnu[i]  = (neut->Vector).Primary.AbsMomentum[i];
//    for ( int j = 0 ; j < 3 ; j++ ){
//      runaction -> dirnu[i][j] = (neut->Vector).Primary.Momentum[i][j] / (neut->Vector).Primary.AbsMomentum[i];
//    }
//  }
//
//  runaction -> Crsx   = (neut->Vector).Crs.Crsx;
//  runaction -> Crsy   = (neut->Vector).Crs.Crsy;
//  runaction -> Crsz   = (neut->Vector).Crs.Crsz;
//  runaction -> Crsphi = (neut->Vector).Crs.Crsphi;
//
//  runaction -> Nvert = (neut->Vector).Fsi.Nvert;
//  for (int ivert=0; ivert<150; ivert++) {
//    runaction -> Iflgvert[ivert] = (neut->Vector).Fsi.Iflgvert[ivert];
//    for (int j=0; j<3; j++)
//      runaction -> Posvert[ivert][j] = (neut->Vector).Fsi.Posvert[ivert][j];
//  }
//
//  runaction -> Nvcvert = (neut->Vector).Fsi.Nvcvert;
//  for (int ip=0; ip<900; ip++) {
//    runaction -> Abspvert[ip]  = (neut->Vector).Fsi.Abspvert[ip];
//    runaction -> Abstpvert[ip] = (neut->Vector).Fsi.Abstpvert[ip];
//    runaction -> Ipvert[ip]    = (neut->Vector).Fsi.Ipvert[ip];
//    runaction -> Iverti[ip]    = (neut->Vector).Fsi.Iverti[ip];
//    runaction -> Ivertf[ip]    = (neut->Vector).Fsi.Ivertf[ip];
//    for (int j=0; j<3; j++)
//      runaction -> Dirvert[ip][j] = (neut->Vector).Fsi.Dirvert[ip][j];
//  }
//  runaction -> Fsiprob = (neut->Vector).Fsi.Fsiprob;
//  runaction -> Numbndn = (neut->Vector).target_info.Numbndn;
//  runaction -> Numbndp = (neut->Vector).target_info.Numbndp;
//  runaction -> Numfrep = (neut->Vector).target_info.Numfrep;
//  runaction -> Numatom = (neut->Vector).target_info.Numatom;
//  runaction -> Ibound  = (neut->Vector).Fsi.Ibound;
//  runaction -> Npvc    = (neut->Vector).Secondary.NumParticle;
//  for (int i=0; i<100; i++) {
//    runaction -> Ipvc[i]    = (neut->Vector).Secondary.ParticleID[i];
//    runaction -> Ichvc[i]   = (neut->Vector).Secondary.TrackingFlag[i];
//    runaction -> Iorgvc[i]  = (neut->Vector).Secondary.ParentID[i];
//    runaction -> Iflvc[i]   = (neut->Vector).Secondary.InteractionCode[i];
//    runaction -> Abspvc[i]  = (neut->Vector).Secondary.AbsMomentum[i];
//    for (int j=0; j<3; j++)
//      runaction -> Pvc[i][j]     = (neut->Vector).Secondary.Momentum[i][j];
//  }
//#endif


  if(genie_para.IsCC){
      G4ParticleDefinition* particle;
      if     (genie_para.Neutrino== 14){ particle = particleTable->FindParticle( 13);} 
      else if(genie_para.Neutrino==-14){ particle = particleTable->FindParticle(-13);} 
      else{
        cout << "ERROR : No muon was found in CC interaaction" << endl;
        exit(1);
      }
      double nvec[3];
      nvec[0] = genie_para.Pxl;
      nvec[1] = genie_para.Pyl;
      nvec[2] = genie_para.Pzl;
      G4ThreeVector dir(nvec[0], nvec[1], nvec[2]);
      double nvec_norm=0.;
      for(int i=0;i<3;i++){nvec_norm+=nvec[i]*nvec[i];}

      G4double mass   = particle->GetPDGMass();
      G4double mom    = sqrt(nvec_norm)*GeV;
      G4double energy = sqrt(mass*mass+mom*mom) - mass;

      particleGun->SetParticleDefinition(particle);
      particleGun->SetParticlePosition(G4ThreeVector(pos[0]*cm,pos[1]*cm,pos[2]*cm));
      particleGun->SetParticleMomentumDirection(dir);
      particleGun->SetParticleEnergy(energy);
      particleGun->SetParticleTime(0.0*ns);
      particleGun->GeneratePrimaryVertex(anEvent);

  }

  for(int ipart=0; ipart<genie_para.Nf; ipart++) {

    // Set ParticleGun with information above ------------------------------------------
    G4ParticleDefinition* particle;
    particle = particleTable->FindParticle(genie_para.Pdgf[ipart]);

    double nvec[3];
    nvec[0] = genie_para.Pxf[ipart];
    nvec[1] = genie_para.Pyf[ipart];
    nvec[2] = genie_para.Pzf[ipart];
    G4ThreeVector dir(nvec[0], nvec[1], nvec[2]);
    double nvec_norm=0.;
    for(int i=0;i<3;i++){nvec_norm+=nvec[i]*nvec[i];}

    G4double mass   = particle->GetPDGMass();
    G4double mom    = sqrt(nvec_norm)*GeV;
    G4double energy = sqrt(mass*mass+mom*mom) - mass;


    particleGun->SetParticleDefinition(particle);
    particleGun->SetParticlePosition(G4ThreeVector(pos[0]*cm,pos[1]*cm,pos[2]*cm));
    particleGun->SetParticleMomentumDirection(dir);
    particleGun->SetParticleEnergy(energy);
    particleGun->SetParticleTime(0.0*ns);
    particleGun->GeneratePrimaryVertex(anEvent);

  } // end of for loop
}
