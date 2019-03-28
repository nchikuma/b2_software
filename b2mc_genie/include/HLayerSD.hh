#ifndef HLayerSD_H
#define HLayerSD_H
 
#include "G4VSensitiveDetector.hh"

#include "HLayerHit.hh"
#include "DetectorResponse.hh"
#include "DetectorDimension.hh"
#include "Const.hh"

class G4HCofThisEvent;
class G4Step;

class HLayerSD : public G4VSensitiveDetector {
  G4THitsCollection<HLayerHit>* hlayerHitCollection;
  DetectorResponse *detresp;
  G4double TotalhlayerDep;

  public:
    HLayerSD(const G4String& name,DetectorDimension*);
    virtual ~HLayerSD();

    virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
    virtual void Initialize(G4HCofThisEvent* HCTE);
    virtual void EndOfEvent(G4HCofThisEvent* HCTE);

    virtual void DrawAll();
    virtual void PrintAll(); 


    inline G4int GetNhit()
    { return hlayerHitCollection->entries(); }
    inline G4double GetTotalDep()
    { return TotalhlayerDep; }
  private:
    DetectorDimension *detdim;
};

#endif
