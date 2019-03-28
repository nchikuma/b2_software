#include "TrackingAction.hh"
#include  "DetectorConstruction.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include  "G4VProcess.hh"
#include <iostream>
#include <vector>
#include "G4SystemOfUnits.hh"


TrackingAction::TrackingAction(RunAction *rac,EventAction *evt)
  :runaction(rac)
{
  eventaction = evt;
  simpart = new SimParticleSummary();
}

TrackingAction::~TrackingAction()
{
  if(simpart){
    delete simpart;
  }
}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
#ifdef DEBUG_TRKACT
  G4cout << "\n------ In PreUserTrackingAction ------" << G4endl;
#endif

  int id = aTrack->GetParentID();

  if(id==0){
    runaction->GetEvtSum()->GetSimVertex(0)->targetz = aTrack->GetLogicalVolumeAtVertex()->GetMaterial()->GetElement(0)->GetZ();
  }

  for(int i=0; i<3; i++){
    posi[i] = (aTrack->GetPosition()/cm)[i];
    momi[i] = (float)(aTrack->GetMomentum())[i];

    simpart -> ipos[i]     = (float)(aTrack->GetPosition()/cm)[i];
    simpart -> momentum[i] = (float)(aTrack->GetMomentum()/GeV)[i];
    simpart -> dir[i]      = (float)(aTrack->GetMomentumDirection())[i];
  }

  initE = aTrack->GetTotalEnergy();
  
  simpart -> ipos[3]     = aTrack->GetGlobalTime()*ns;
  simpart -> momentum[3] = aTrack->GetKineticEnergy()/GeV;
  simpart -> iposflag = 1;
    
  if(id==0){
#ifdef DEBUG_TRKACT
    G4cout << "init  Pos(cm) " 
	   << (aTrack->GetPosition()/cm)[0] << " " 
	   << (aTrack->GetPosition()/cm)[1] << " "
	   << (aTrack->GetPosition()/cm)[2] << G4endl; 
    G4cout << "init  Momentum(GeV) " 
	   << (aTrack->GetMomentum()/GeV)[0] << " " 
	   << (aTrack->GetMomentum()/GeV)[1] << " "
	   << (aTrack->GetMomentum()/GeV)[2] << " "
	   << (aTrack->GetKineticEnergy()/GeV) << G4endl; 
#endif
  }
}


void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  int parentID = aTrack->GetParentID();
  int trackID  = aTrack->GetTrackID();
  int pid      = aTrack->GetDefinition()->GetPDGEncoding();  

  G4double trackL = aTrack->GetTrackLength();

  simpart -> parentid = parentID;
  simpart -> trackid  = trackID;
  simpart -> pdg      = pid;
  simpart -> length   = trackL/cm;

  for(int i=0; i<3; i++){
    simpart -> fpos[i] = (float)(aTrack->GetPosition()/cm)[i];
  }

  simpart -> fposflag = 1;

  if(parentID==0){// parentID = 0: first point of track
	   
    float posf[3];
    float momf[3];

    G4cout.precision( 4 );
	   
    for(int i=0; i<3; i++){
      posf[i] = (aTrack-> GetPosition()/cm)[i];
      momf[i] = (float)(aTrack->GetMomentum()/MeV)[i];

      simpart -> fpos[i] = (float)(aTrack->GetPosition()/cm)[i];
    }

#ifdef DEBUG_TRKACT
    G4cout << "Track:" << trackID;
    G4cout << "  " << aTrack->GetDefinition()->GetParticleName()  
	   << "(" << aTrack->GetDefinition()->GetPDGEncoding() << ")\n";  
    G4cout << " Ei:" << initE
	   << " Momi:{" << momi[0] << "," << momi[1] << "," << momi[2] << "}"
	   << " Start:{"<< posi[0] << "," << posi[1] << "," << posi[2] << "}"
	   << " Stop:{"<< posf[0] << "," << posf[1] << "," << posf[2] << "}"
	   << " Length:" << trackL/cm << G4endl;
#endif

    runaction -> GetEvtSum() -> AddSimParticle (simpart);
	   
  }
   
}
