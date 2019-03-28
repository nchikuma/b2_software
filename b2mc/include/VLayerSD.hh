#ifndef VLayerSD_H
#define VLayerSD_H
 
#include "G4VSensitiveDetector.hh"

#include "VLayerHit.hh"
#include "DetectorResponse.hh"
#include "DetectorDimension.hh"
#include "Const.hh"

class G4HCofThisEvent;
class G4Step;

class VLayerSD : public G4VSensitiveDetector {
  G4THitsCollection<VLayerHit>* vlayerHitCollection;
  G4double TotalvlayerDep;
  DetectorResponse *detresp;
 
  public:
    VLayerSD(const G4String& name, DetectorDimension*);
    virtual ~VLayerSD();

    virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
    virtual void Initialize(G4HCofThisEvent* HCTE);
    virtual void EndOfEvent(G4HCofThisEvent* HCTE);

    virtual void DrawAll();
    virtual void PrintAll(); 

     
    inline G4int GetNhit()
    { return vlayerHitCollection->entries(); }
    inline G4double GetTotalDep()
    { return TotalvlayerDep; }
  private:
    DetectorDimension *detdim;

};

#endif
