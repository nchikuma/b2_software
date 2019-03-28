#ifndef VetoSD_H
#define VetoSD_H
 
#include "G4VSensitiveDetector.hh"

#include "VetoHit.hh"
#include "DetectorResponse.hh"
#include "Const.hh"
#include "DetectorDimension.hh"

class G4HCofThisEvent;
class G4Step;

class VetoSD : public G4VSensitiveDetector {
  G4THitsCollection<VetoHit>* vetoHitCollection;
  DetectorResponse *detresp;
  G4double TotalvetoDep;
  G4int vetoID[300];
  G4double DEP[300];
  G4int fHIT;
 
public:
  VetoSD(const G4String& name,DetectorDimension*);
  virtual ~VetoSD();

  // virtual methods
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  virtual void Initialize(G4HCofThisEvent* HCTE);
  virtual void EndOfEvent(G4HCofThisEvent* HCTE);

  virtual void DrawAll();
  virtual void PrintAll(); 

  
  inline G4int GetNhit()
  { return vetoHitCollection->entries(); }
  inline G4double GetTotalDep()
  { return TotalvetoDep; }
  


};

#endif
