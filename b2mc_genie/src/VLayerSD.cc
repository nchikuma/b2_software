#include "G4VPhysicalVolume.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "VLayerSD.hh"

VLayerSD::VLayerSD(const G4String& name,DetectorDimension* fdim)
  : G4VSensitiveDetector(name)
{
  collectionName.insert("vlayerHitCollection");
  detdim  = fdim;
  detresp = new DetectorResponse(detdim);
}


VLayerSD::~VLayerSD()
{
    if(detresp!=NULL) delete detresp;
}


void VLayerSD::Initialize(G4HCofThisEvent* HCTE)
{

#ifdef DEBUG_VHITSD
  G4cout << "Initialization of vlayerSD" << G4endl;
#endif

  vlayerHitCollection = 
    new VLayerHitsCollection(SensitiveDetectorName,collectionName[0]);
  TotalvlayerDep = 0.;
  
  static int HCID = -1;
  if(HCID<0) HCID = GetCollectionID(0); 
  HCTE->AddHitsCollection(HCID,vlayerHitCollection);
}



G4bool VLayerSD::ProcessHits(G4Step* aStep, 
				G4TouchableHistory* ROhist)
{
  G4Track* track  = aStep->GetTrack();
  G4double edep   = aStep->GetTotalEnergyDeposit();
  //G4double length = aStep->GetStepLength();
  G4double length = aStep->GetStepLength()/cm;
  if(edep==0.) return false;
  TotalvlayerDep += edep;

  const G4VTouchable* Touchable = aStep->GetPreStepPoint()->GetTouchable();
  G4int detID                   = Touchable->GetVolume(0)->GetCopyNo();
  G4int trackID                 = track->GetTrackID();
  G4int PDG                     = track->GetDefinition()->GetPDGEncoding();
  G4ThreeVector hitPos          = aStep->GetPreStepPoint()->GetPosition();
  G4double hittime              = aStep->GetPreStepPoint()->GetGlobalTime();

#ifdef DEBUG_VHITSD
  G4cout << " ========= VlayerSD ========== " << G4endl;
  G4cout << " detID:" << detID 
	 << " edep:" << edep   
       	 << " length:" << length 
       	 << " trackID:" << trackID
	 << " PDG:" << PDG
	 << " hitPos:{"<<hitPos[0]<<","<<hitPos[1]<<","<< hitPos[2]<<"}"
	 << " hittime:" << hittime
      	 << G4endl;
#endif
  
  G4double edep_q = edep;
  detresp->ApplyScintiResponse(&edep_q,&length,track);

  VLayerHit* aHit 
    = new VLayerHit(detID,PDG,trackID,edep,edep_q,hitPos,hittime,detdim);
    
  VLayerHit* bHit;
 
  for(int k=0;k<vlayerHitCollection->entries();k++){
    bHit = (*vlayerHitCollection)[k];

    if(bHit->CompareID(*aHit)){
         bHit->AddEdep(edep,edep_q);

	 if(bHit->isFaster(*aHit)) { 
		 bHit->SetTime( aHit->GetTime() );
	 }
	 if(bHit->LargerEdep(*aHit)) { 
		 bHit->SetParticle(aHit->GetParticle()); 
	 }
   	 return true;
    }
  }

  vlayerHitCollection->insert( aHit );

  return true;
}

void VLayerSD::EndOfEvent(G4HCofThisEvent* HCTE)
{
#ifdef DEBUG_VHITSD
  G4cout << " ========= VlayerSD::EndOfEvent ========== " << G4endl;
#endif

	VLayerHit *cHit;

	G4double edep_tmp;
	G4double time_tmp;
	G4ThreeVector posinmod;
	G4int mod;
	G4int view;
	G4int pln;
	G4int ch;
	G4int adc;
	G4int loadc;
	G4double pe;
	G4double lope;

	// apply detector response
	for(G4int k=0;k<vlayerHitCollection->entries();k++) {
		cHit = (*vlayerHitCollection)[k];
		edep_tmp = cHit->GetEdepQ();
		time_tmp = cHit->GetTime();
		posinmod = cHit->GetPosInMod();
		mod      = cHit->GetMod();
		view     = cHit->GetView();
		pln      = cHit->GetPln();
		ch       = cHit->GetCh();
#ifdef DEBUG_VHITSD
		G4cout << " edep_tmp:"  << edep_tmp
		       << " time_tmp:"  << time_tmp
		       << " mod:"  << mod
		       << " view:"  << view
		       << " pln:"  << pln
		       << " ch:"  << ch
		       << " pos:{" << posinmod[0] << "," << posinmod[1] << "," << posinmod[2] << "}"
		       << " adc:"  << adc
		       << " loadc:"  << loadc
		       << " pe:"  << pe
		       << " lope:"  << lope
		       << G4endl;
#endif

		//apply light collection
		detresp->ApplyLightCollection(&edep_tmp,mod,pln,view,ch,posinmod);

		//apply fiber attenuation
		detresp->ApplyFiberResponse(&edep_tmp,&time_tmp,mod,view,pln,ch,posinmod);

		//convert edep -> p.e. & cross & after pulse
		detresp->ApplyMPPCResponse(edep_tmp,&pe,mod);

		//apply ADC responce
		detresp->ApplyADCResponse(&pe,&lope,&adc,&loadc,mod);

		//fill variable to hitcollection
		cHit->SetPE(pe);
		cHit->SetLOPE(lope);
		cHit->SetDelayTime(time_tmp);
	}
}

void VLayerSD::DrawAll()
{
  for(G4int k=0; k < vlayerHitCollection->entries(); k++)
   (*vlayerHitCollection)[k]->Draw(); 
}

void VLayerSD::PrintAll()
{
   for(G4int k=0; k < vlayerHitCollection->entries(); k++)
     (*vlayerHitCollection)[k]->Print(); 
}

