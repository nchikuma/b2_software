#include "G4VPhysicalVolume.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "VetoSD.hh"

VetoSD::VetoSD(const G4String& name,DetectorDimension *fdim)
  : G4VSensitiveDetector(name)
{
  collectionName.insert("vetoHitCollection");
  detresp = new DetectorResponse(fdim);
}


VetoSD::~VetoSD()
{
    if(detresp!=NULL) delete detresp;
}


void VetoSD::Initialize(G4HCofThisEvent* HCTE)
{

#ifdef DEBUG_VETOHITSD
  G4cout << "Initialization of vetoSD" << G4endl;
#endif

  vetoHitCollection = 
    new VetoHitsCollection(SensitiveDetectorName,collectionName[0]);
  TotalvetoDep = 0.;
 
  static int HCID = -1;
  if(HCID<0) HCID = GetCollectionID(0); 
  HCTE->AddHitsCollection(HCID,vetoHitCollection);
}



G4bool VetoSD::ProcessHits(G4Step* aStep, 
				G4TouchableHistory* ROhist)
{
  G4Track* track = aStep->GetTrack();
  G4double edep  = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return false;
  TotalvetoDep += edep;
 
  const G4VTouchable* Touchable = aStep->GetPreStepPoint()->GetTouchable();
  G4int detID                   = Touchable->GetVolume(0)->GetCopyNo();
  G4int trackID                 = track->GetTrackID();
  G4int PDG                     = track->GetDefinition()->GetPDGEncoding();
  G4ThreeVector hitPos_True     = aStep->GetPreStepPoint()->GetPosition();
  G4double hittime              = aStep->GetPreStepPoint()->GetGlobalTime();
  G4ThreeVector hitPos_Center   = Touchable->GetTranslation(0);

  // True X,Y pos, and Z-center of veto-bar 
  G4ThreeVector hitPos;
  hitPos[0] = hitPos_True[0]; // X
  hitPos[1] = hitPos_True[1]; // Y
  hitPos[2] = hitPos_Center[2]; // Z


#ifdef DEBUG_VETOHITSD
  G4cout << " ========= VetoSD ========== " << G4endl;
  G4cout << " detID:" << detID 
	 << " edep:" << edep   
       	 << " trackID:" << trackID
	 << " PDG:" << PDG
	 << " hitPos:{"<<hitPos[0]<<","<<hitPos[1]<<","<< hitPos[2]<<"}"
	 << " hittime:" << hittime
      	 << G4endl;
#endif
  G4double edep_q = edep;
  detresp->ApplyScintiResponse2(&edep_q,track);

  VetoHit* aHit 
    = new VetoHit(detID,PDG,trackID,edep,edep_q,hitPos,hittime);
  
  VetoHit* bHit;
  
  for(int k=0;k<vetoHitCollection->entries();k++){
      	  bHit = (*vetoHitCollection)[k];
      	  if(bHit->CompareID(*aHit)){
	    	  bHit->AddEdep(edep,edep_q);
		  if(bHit->isFaster(*aHit)) { 
			  bHit->SetTime( aHit->GetTime() );
			  bHit->SetParticle( aHit->GetParticle() );
		  }
	    	  return true;
      	  }
  }
  
  vetoHitCollection->insert( aHit );
  
  return true;
  
}

void VetoSD::EndOfEvent(G4HCofThisEvent* HCTE)
{
#ifdef DEBUG_VETOHITSD	
	G4cout << "VetoSD::EndOfEvent(G4HCofThisEvent* HCTE)" << G4endl;
#endif

	VetoHit *cHit;

	G4double edep_tmp;
	G4double time_tmp;
	G4ThreeVector posinmod;
	G4int pln;
	G4int adc;
	G4int loadc;
	G4double pe;
	G4double lope;

	// apply ingrid response
	for(G4int k=0;k<vetoHitCollection->entries();k++) {
		cHit = (*vetoHitCollection)[k];

		edep_tmp = cHit->GetEdepQ();
		time_tmp = cHit->GetTime();
		posinmod = cHit->GetPosInMod();
		pln = cHit->GetPln();

		//apply fiber attenuation
		detresp->ApplyFiberResponseV(&edep_tmp,&time_tmp,pln,posinmod);

		//convert edep -> p.e. & cross & after pulse
		detresp->ApplyMPPCResponse(edep_tmp,&pe,0);

		//apply ADC responce
		detresp->ApplyADCResponse(&pe,&lope,&adc,&loadc,0);

		//fill variable to hitcollection
		cHit->SetPE(pe);
		cHit->SetLOPE(lope);
		cHit->SetDelayTime(time_tmp);
	}
}

void VetoSD::DrawAll()
{
  for(G4int k=0; k < vetoHitCollection->entries(); k++)
   (*vetoHitCollection)[k]->Draw(); 
}

void VetoSD::PrintAll()
{
   for(G4int k=0; k < vetoHitCollection->entries(); k++)
     (*vetoHitCollection)[k]->Print(); 
}

