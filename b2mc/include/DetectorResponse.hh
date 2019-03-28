#ifndef DETECTOR_RESPONSE
#define DETECTOR_RESPONSE

#include "G4ThreeVector.hh"
#include "G4EmCalculator.hh"
#include "G4Track.hh"

#include "DetectorDimension.hh"
#include "Const.hh"

class DetectorResponse {
public:
  DetectorResponse(DetectorDimension*);
  ~DetectorResponse();


  //void ApplyScintiResponse(G4double* edep, G4Track* aTrack);
  void ApplyScintiResponse(G4double* edep, G4double* length, G4Track* aTrack);
  void ApplyScintiResponse2(G4double* edep, G4Track* aTrack); //For VetoSD.cc
  void ApplyLightCollection(G4double* edep, G4int mod, G4int pln, G4int view, G4int ch, G4ThreeVector hitpos);
  //void ApplyFiberResponse(G4double* edep, G4double* time, G4int mod, G4int view, G4ThreeVector hitpos);
  void ApplyFiberResponse(G4double* edep, G4double* time, G4int mod, G4int view, G4int pln, G4int ch, G4ThreeVector hitpos);
  void ApplyFiberResponseV(G4double* edep, G4double* time, G4int pln, G4ThreeVector hitpos);
  void ApplyMPPCResponse(G4double edep, G4double* pe, G4int mod);
  void ApplyADCResponse(G4double* pe, G4double* lope, G4int* adc, G4int* loadc, G4int mod);
  void ApplyTDCResponse(G4double time, G4int* tdc);

  double peratio [24][2][18][80]; //[mod][view][pln][ch]
  double attleng [24][2]; //[mod][scinti_type]
  double reflec  [24][2]; //[mod][scinti_type]

  G4double distToFiber, distToMPPC,distToMPPC2;
  G4ThreeVector fiberpos;

private:
  DetectorDimension *detdim;
  G4EmCalculator emcal;

  void BirksSaturation(G4double* edeposit, G4double* length, G4Track* aTrack);
  void BirksSaturation2(G4double* edeposit, G4Track* aTrack); //For VetoSD.cc

};

#endif
