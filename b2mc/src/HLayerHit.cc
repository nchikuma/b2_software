#include "HLayerHit.hh"

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"

#include "DetectorDimension.hh"


// allocator
G4Allocator<HLayerHit> HLayerHitAllocator;

HLayerHit::HLayerHit(G4int dID0, G4int P0, G4int trkid, G4double e, G4double eq, G4ThreeVector pos, G4double t,DetectorDimension *fdim) 
{
  detID    = dID0;
  trackID  = trkid;
  edep     = e;
  edep_q   = eq;
  Particle = P0;
  position = pos;
  time     = t;
  detdim   = fdim;

  mod = dID0/NORMMOD;
  pln = (dID0-NORMMOD*mod)/NORMPLN;
  ch  = (dID0-NORMMOD*mod)%NORMPLN;

  view = SideView;

#ifdef DEBUG_HHIT
  G4cout << " ========== HLayerHit::HLayerHit ==============" << G4endl;
  G4cout <<  " detID:" << dID0 << " mod:" << mod << " pln:" << pln << " view:" << view <<  " ch:" << ch << G4endl;
#endif

  if(mod<7){
	  ModOffset[0] = (C_INGHMotherPosX + C_INGStart + C_INGSpace*mod)*mm;
	  ModOffset[1] =  C_INGHMotherPosY*mm;
	  ModOffset[2] =  C_INGHMotherPosZ*mm;
  }
  else if(mod<14){
	  ModOffset[0] =  C_INGVMotherPosX*mm;
	  ModOffset[1] = (C_INGVMotherPosY + C_INGStart + C_INGSpace*(mod-7))*mm;
	  ModOffset[2] =  C_INGVMotherPosZ*mm;
  }
  //else if(mod==15||mod==16){
  else if(mod==15){
	  ModOffset[0] = C_PMMotherPosX*mm;
	  ModOffset[1] = C_PMMotherPosY*mm;
	  ModOffset[2] = C_PMMotherPosZ*mm;
  }
  else if(mod==21){
	  ModOffset[0] = (C_B2MotherPosX + C_B2WMPosX)*mm;
	  ModOffset[1] = (C_B2MotherPosY + C_B2WMPosY)*mm;
	  ModOffset[2] = (C_B2MotherPosZ + C_B2WMPosZ)*mm;
  }
  //else if(mod==22){
  else if(mod==16){
	  ModOffset[0] = (C_B2MotherPosX + C_B2CHPosX)*mm;
	  ModOffset[1] = (C_B2MotherPosY + C_B2CHPosY)*mm;
	  ModOffset[2] = (C_B2MotherPosZ + C_B2CHPosZ)*mm;
  }
  else if(mod==14){
	  ModOffset[0] = (C_B2MotherPosX + C_B2INGPosX)*mm;
	  ModOffset[1] = (C_B2MotherPosY + C_B2INGPosY)*mm;
	  ModOffset[2] = (C_B2MotherPosZ + C_B2INGPosZ)*mm;
  }

  posinmod[0] = pos[0] - ModOffset[0];
  posinmod[1] = pos[1] - ModOffset[1];
  posinmod[2] = pos[2] - ModOffset[2];

  gridcell_id_x1=-1;
  gridcell_id_x2=-1;
  gridcell_id_y1=-1;
  gridcell_id_y2=-1;

  if(mod==15 || mod==21){
      detdim->GetWMGridCellID(mod,pln,view,ch,posinmod[0]/mm,posinmod[1]/mm,posinmod[2]/mm,
		              &gridcell_id_x1,&gridcell_id_x2,&gridcell_id_y1,&gridcell_id_y2);
  }   


}

HLayerHit::HLayerHit(G4int dID0, G4double e,G4int P0) 
{
  detID = dID0;
  edep = e;
  Particle = P0;
}

HLayerHit::HLayerHit(G4int dID0, G4double e) 
{
  detID = dID0;
  edep = e;
}


HLayerHit::~HLayerHit() 
{
}

HLayerHit::HLayerHit(const HLayerHit &right)
  : G4VHit()
{
  detID = right.detID;
  edep       = right.edep;
  Particle = right.Particle;

  for(int i=0;i<3;i++) position[i] = right.position[i];
  eventID = right.eventID;

}

const HLayerHit& HLayerHit::operator=(const HLayerHit &right)
{
  detID    = right.detID;
  edep     = right.edep;
  Particle = right.Particle; 
  
  for(int i=0;i<3;i++) position[i] = right.position[i];
  eventID = right.eventID;

  return *this;
}

G4int HLayerHit::operator==(const HLayerHit &right) const
{
  return (this==&right) ? 1 : 0;
}  

G4int HLayerHit::CompareID(const HLayerHit right) 
{
  return (detID==right.detID) ? 1 : 0;
}

G4int HLayerHit::CompareP(const HLayerHit right) 
{
  return (Particle==right.Particle) ? 1 : 0;
}

G4int HLayerHit::isFaster(const HLayerHit right)
{
  return (time>right.time) ? 1 : 0;
}

G4int HLayerHit::LargerEdep(const HLayerHit right)
{
  //return (Particle==11||right.Particle!=11) ? 1 : 0;
  if(Particle==11) return 1;
  else if(right.Particle==11) return 0;
  else return (edep<right.edep) ? 1 : 0;
}

void HLayerHit::Draw()
{
    double size = edep + 7.;
    if( edep>9 ) size = 16.;
  #ifdef VIS_USE
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager)
    {
      G4ThreeVector pos; 
      for(int i=0;i<3;i++) pos[i] = position[i]*mm;
      G4Circle circle(pos);
      circle.SetScreenSize(size); 
      circle.SetFillStyle(G4Circle::filled);
      G4Colour colour(0.,1.,1.3); //cyan
      G4VisAttributes attribs(colour);
      circle.SetVisAttributes(attribs);
      pVVisManager->Draw(circle);
    }
  #endif
}

void HLayerHit::Print()
{
  G4cout.precision(4);
  
  G4cout << " Mod:" << mod 
	 << " detID:" << detID
	 << " Pln:" << pln 
	 << " Ch:" << ch 
	 << " Time:" << time
	 << " Edep:" << edep
	 << " p.e.:" << pe
	 << " PID:" << Particle
	 << " Trk:" << trackID
	 << " pos:{" << position[0]/cm << ", "<< position[1]/cm << ", "<< position[2]/cm<< "(cm)}"
	 << G4endl;
}
