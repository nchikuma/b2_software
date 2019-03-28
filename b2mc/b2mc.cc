#include "G4RunManager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
//#include "QGSP.hh"
#include "QGSP_BERT.hh"
#include "QGSP_BIC.hh"
//#include "QGSP_BERT_HP.hh"
//#include "G4NeutronHPData.hh"
#include "DetectorConstruction.hh"
#include "DetectorDimension.hh"
//#include "ND280mPhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "Neut.hh"
#include "Const.hh"
#include "G4VisExecutive.hh"
#include "G4HepRep.hh"
#include "G4HepRepFile.hh"
#include <string>

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

// ====================================================================
//     main
// ====================================================================
int main(int argc, char** argv) 
{
  char neutfile[500];
  char output[500];
  char cmd[500];

  int nd = 0;
  int batch = 0;
  int flav = 0;
  // 1:numu, 2:numubar, 3:nue, 4:nuebar
  int target = 0;
  // 1:ch, 2:water,3:fe

  int nevent = 0;
  int NumberOfEvent = 0;

  int c = -1;
  while ((c = getopt(argc, argv, "ho:t:e:i:m:b:f:")) != -1) {
    switch(c){
    case 'o':
      snprintf(output,sizeof(output),"%s",optarg);
      break;
    case 'i':
      snprintf(neutfile,sizeof(neutfile),"%s",optarg);
      break;
    case 'm':
      nd = atoi(optarg);
      break;
    case 'f':
      flav = atoi(optarg);
      break;
    case 't':
      target = atoi(optarg);
      break;
    case 'e':
      nevent = atoi(optarg);
      break;
    case 'b':
      batch = 1;
      sprintf(cmd,"%s",optarg);
      break;
    case 'h':
      std::cerr << "o : output root file name" << std::endl;
      std::cerr << "i : input neut file" << std::endl;
      std::cerr << "m : module" <<std::endl;
      std::cerr << "    2 ProtonModule" << std::endl;
      std::cerr << "    3 INGRID Horizontal" << std::endl;
      std::cerr << "    4 INGRID Vertical" << std::endl;
      std::cerr << "    5 On-Axis WaterModule" << std::endl;
      std::cerr << "    6 On-Axis WaterModule BG" << std::endl;
      std::cerr << "    7 B2 WaterModule" << std::endl;
      std::cerr << "    8 B2 CHModule" << std::endl;
      std::cerr << "    9 B2 INGRID" << std::endl;
      std::cerr << "   10 B2Wall"    << std::endl;
      std::cerr << "f : flavor" <<std::endl;
      std::cerr << "   1 numu" << std::endl;
      std::cerr << "   2 anti-numu" << std::endl;
      std::cerr << "   3 nue" << std::endl;
      std::cerr << "   4 anti-nue" << std::endl;
      std::cerr << "b : batch command" << std::endl;
      exit(1);
    }
  }

  if( nd==0 ) {
    G4cout << "Select horizontal or vertical or proton module" << G4endl;
    G4cout << "   2 ProtonModule" << G4endl;
    G4cout << "   3 INGRID Horizontal" << G4endl;
    G4cout << "   4 INGRID Vertical" << G4endl;
    G4cout << "   5 On-Axis WaterModule" << G4endl;
    G4cout << "   6 On-Axis WaterModule BG" << G4endl;
    G4cout << "   7 B2 WaterModule" << G4endl;
    G4cout << "   8 B2 CHModule" << G4endl;
    G4cout << "   9 B2 BG study" << G4endl;
    exit(1);
  }
  if( flav==0 ) {
    G4cout << "Select neutrino flavor" << G4endl;
    G4cout << "   1 numu" << G4endl;
    G4cout << "   2 anti-numu" << G4endl;
    G4cout << "   3 nue" << G4endl;
    G4cout << "   4 anti-nue" << G4endl;
    G4cout << "  11 sand-muon" << G4endl;
    G4cout << "  12 cosmic-muon" << G4endl;
    G4cout << "  13 muon-injection" << G4endl;
    exit(1);
  }

  DetectorDimension *detdim = new DetectorDimension();

  // run manager
  G4cout << "Open runManager..." << G4endl;
  G4RunManager* runManager= new G4RunManager;

  // Neut initialization
  Neut *neut = new Neut;
  if(flav==11 || flav==12 || flav==13 || flav==14 || flav==15){
    if(flav==11) G4cout << "Sand-muon mode (WAGASCI)" << G4endl;
    if(flav==12) G4cout << "Cosmic-muon mode" << G4endl;
    if(flav==13) G4cout << "Muon test mode" << G4endl;
    if(flav==14) G4cout << "Pion test mode" << G4endl;
    if(flav==15) G4cout << "Proton test mode" << G4endl;
    NumberOfEvent = nevent;
    if(NumberOfEvent==0){
      G4cout << "<<< Error: Input the number of events >>>" << G4endl;
      exit(1);
    }
    G4cout << "NumberOfEvent: " << NumberOfEvent << G4endl;
  }
  else{
    G4cout << "Open Neut..." << G4endl;
    NumberOfEvent = neut->NtupleReadInit(neutfile);
    G4cout << "NumberOfEvent :" << NumberOfEvent << G4endl;
  }

  // detector setup
  //runManager-> SetUserInitialization(new DetectorConstruction(nd,shift_INGpos));
  runManager-> SetUserInitialization(new DetectorConstruction(nd,detdim));
  G4cout << "Detector Init OK" << G4endl;

  // particles and physics processes
  //runManager-> SetUserInitialization(new ND280mPhysicsList);
  //runManager-> SetUserInitialization(new QGSP_BIC);
  runManager-> SetUserInitialization(new QGSP_BERT);
  G4cout << "PhysicsList Init OK" << G4endl;

  G4cout << "output : " << output << endl;
  RunAction * rac = new RunAction(output);
  G4cout << "RunAction init OK" << G4endl;

  EventAction * evt = new EventAction(rac,detdim);
  G4cout << "EventAction init OK" << G4endl;

  TrackingAction * tra = new TrackingAction(rac, evt);
  runManager->SetUserAction(tra);
  G4cout << "TrackingAction init OK" << G4endl;

  runManager-> SetUserAction(new PrimaryGeneratorAction(neut, rac, evt, nd, flav));
  G4cout << "PrimaryGenerator init OK" << G4endl;

  // user action classes... (optional)
  runManager-> SetUserAction(rac);
  runManager-> SetUserAction(evt);

#ifdef G4VIS_USE
  // initialize visualization package
  G4VisManager* visManager= new G4VisExecutive;
  visManager-> Initialize();
  G4cout <<"visualization init OK" << G4endl;
#endif

  // Initialize G4 kernel
  runManager-> Initialize();

  // get the pointer to the UI manager and set verbosities
  G4UImanager* UI= G4UImanager::GetUIpointer();

  G4cout << "*******************************************************" << G4endl;

  // Define (G)UI terminal for interactive mode  
  if(batch==1)
    { 
      // G4UIterminal is a (dumb) terminal.
      G4UIsession * session = 0;
      G4String command = "/control/execute ";
      G4String macro = cmd;
      session = new G4UIterminal(new G4UItcsh);
      UI->ApplyCommand(command+macro);
      session->SessionStart();
      delete session;
    }
  else 
    {
      runManager->BeamOn(NumberOfEvent);
    }

  //terminating...
  #ifdef G4VIS_USE
    delete visManager;
  #endif

  delete runManager;  //G4cout << G4endl;
  return 0;

}
