#include "G4VPhysicalVolume.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "HLayerSD.hh"

HLayerSD::HLayerSD(const G4String& name,DetectorDimension* fdim)
  : G4VSensitiveDetector(name)
{
  collectionName.insert("hlayerHitCollection");
  detdim  = fdim;
  detresp = new DetectorResponse(detdim);
}


HLayerSD::~HLayerSD()
{
    if(detresp!=NULL) delete detresp;
}


void HLayerSD::Initialize(G4HCofThisEvent* HCTE)
{

#ifdef DEBUG_HHITSD
  G4cout << " ========= HlayerSD::Initialize ======== " << G4endl;
#endif

  hlayerHitCollection = 
    new HLayerHitsCollection(SensitiveDetectorName,collectionName[0]);
  TotalhlayerDep = 0.;
  
  static int HCID = -1;
  if(HCID<0) HCID = GetCollectionID(0); 
  HCTE->AddHitsCollection(HCID,hlayerHitCollection);
}



G4bool HLayerSD::ProcessHits(G4Step* aStep, 
				G4TouchableHistory* ROhist)
{
  G4Track* track  = aStep->GetTrack();
  G4double edep   = aStep->GetTotalEnergyDeposit();
  //G4double length = aStep->GetStepLength();
  G4double length = aStep->GetStepLength()/cm;
  if(edep==0.) return false;
  TotalhlayerDep += edep;

  const G4VTouchable* Touchable = aStep->GetPreStepPoint()->GetTouchable();
  G4int detID                   = Touchable->GetVolume(0)->GetCopyNo();
  G4int trackID                 = track->GetTrackID();
  G4int PDG                     = track->GetDefinition()->GetPDGEncoding();
  G4ThreeVector hitPos          = aStep->GetPreStepPoint()->GetPosition();
  G4double hittime              = aStep->GetPreStepPoint()->GetGlobalTime();
 

#ifdef DEBUG_HHITSD
  G4cout << " ========= HlayerSD::ProcessHits ========== " << G4endl;
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

  HLayerHit* aHit 
    = new HLayerHit(detID,PDG,trackID,edep,edep_q,hitPos,hittime,detdim);
  
  HLayerHit* bHit;
  
  for(int k=0;k<hlayerHitCollection->entries();k++){
	  bHit = (*hlayerHitCollection)[k];
	  
	  if(bHit->CompareID(*aHit)){
		  bHit->AddEdep(edep,edep_q);
		  
		  if(bHit->isFaster(*aHit)) { 
			  bHit->SetTime(aHit->GetTime()); 
		  }
		  if(bHit->LargerEdep(*aHit)) { 
			  bHit->SetParticle(aHit->GetParticle()); 
		  }
		  return true;
	  }
  }
  hlayerHitCollection->insert( aHit );
  
  return true;
}

void HLayerSD::EndOfEvent(G4HCofThisEvent* HCTE)
{

	HLayerHit *cHit;
	
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
	
#ifdef DEBUG_HHITSD
	      	G4cout << " ========= HlayerSD::EndOfEvent ========== " << G4endl;
	      	G4cout << " # of hitcollection " << hlayerHitCollection->entries() << G4endl;
#endif
	// apply ingrid response
	for(G4int k=0;k<hlayerHitCollection->entries();k++) {
		cHit = (*hlayerHitCollection)[k];		
      		edep_tmp = cHit->GetEdepQ();
      		time_tmp = cHit->GetTime();
      		posinmod = cHit->GetPosInMod();
      		mod      = cHit->GetMod();
      		view     = cHit->GetView();
      		pln      = cHit->GetPln();
      		ch       = cHit->GetCh();
#ifdef DEBUG_HHITSD
		G4cout << " edep_tmp:"  << edep_tmp
		       << " time_tmp:"  << time_tmp
		       << " mod:"  << mod
		       << " view:"  << view
		       << " pln:"  << pln
		       << " ch:"  << ch
		       << " adc:"  << adc
		       << " loadc:"  << loadc
		       << " pe:"  << pe
		       << " lope:"  << lope
		       << G4endl;
#endif

      		//apply light collection
      		detresp->ApplyLightCollection(&edep_tmp,mod,pln,view,ch,posinmod);
		
      		//apply fiber attenuation
      		//detresp->ApplyFiberResponse(&edep_tmp,&time_tmp,mod,view,posinmod);
		detresp->ApplyFiberResponse(&edep_tmp,&time_tmp,mod,view,pln,ch,posinmod);
		
      		//convert edep -> p.e. & cross-talk & after-pulse
      		detresp->ApplyMPPCResponse(edep_tmp,&pe,mod);
		
      		//apply ADC responce
		//detresp->ApplyADCResponse(&pe,&lope,&adc,&loadc);
		detresp->ApplyADCResponse(&pe,&lope,&adc,&loadc,mod);
		
      		//fill variable to hitcollection
      		cHit->SetPE(pe);
      		cHit->SetLOPE(lope);
      		cHit->SetDelayTime(time_tmp);
	}


}

void HLayerSD::DrawAll()
{
  for(G4int k=0; k < hlayerHitCollection->entries(); k++)
   (*hlayerHitCollection)[k]->Draw(); 
}

void HLayerSD::PrintAll()
{
   for(G4int k=0; k < hlayerHitCollection->entries(); k++)
     (*hlayerHitCollection)[k]->Print(); 
}

