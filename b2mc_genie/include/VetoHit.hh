#ifndef VetoHit_h
#define VetoHit_h 1
 
#include "G4ThreeVector.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

#include "Const.hh"

class VetoHit : public G4VHit 
{
public:
  VetoHit(G4int, G4int, G4int, G4double, G4double, G4ThreeVector, G4double);
  VetoHit(G4int, G4double, G4int);
  VetoHit(G4int, G4double);

  ~VetoHit();
  VetoHit(const VetoHit &right);
  const VetoHit& operator=(const VetoHit &right);
  G4int operator==(const VetoHit &right) const;
  G4int CompareID(const VetoHit right) ;
  G4int CompareP(const VetoHit right);
  G4int isFaster(const VetoHit right);

  // new/delete operators
  inline void *operator new(size_t);
  inline void operator delete(void  *aHit);

  void Draw();
  void Print();
  
private:
    G4int detID;
    G4int trackID;
    G4double edep;
    G4double edep_q;
    G4double pe;
    G4double lope;
    G4int Particle; 
    G4ThreeVector position;
    G4ThreeVector posinmod;
    G4ThreeVector ModOffset;
    G4double time;
    G4double delay_time;
    G4int eventID;

    G4int mod;
    G4int view;
    G4int pln;
    G4int ch;

  //set/get functions 
public:
    inline void AddEdep(G4double e, G4double eq) { edep += e; edep_q += eq; }
    inline void SetEventID(G4int event_ID){ eventID = event_ID;}
    inline void SetPE(G4double p) { pe = p; }
    inline void SetLOPE(G4double p) { lope = p; }
    inline void SetParticle(G4int p) { Particle = p; }
    inline void SetTime(G4double t) { time = t; }
    inline void SetDelayTime(G4double dt) { delay_time = dt; }
    
    inline G4int GetDetID() { return detID; }
    inline G4int GetTrackID() { return trackID; }
    inline G4int GetParticle() { return Particle; }
    inline G4double GetEdep() { return edep; }
    inline G4double GetEdepQ() { return edep_q; }
    inline G4double GetTime() { return time; }
    inline G4double GetDelayTime() { return delay_time; }
    inline G4ThreeVector GetPosition() { return position; }
    inline G4ThreeVector GetPosInMod() { return posinmod; }
    inline G4int GetEventID() { return eventID; }
    inline G4int GetMod() { return mod; }
    inline G4int GetView() { return view; }
    inline G4int GetPln() { return pln; }
    inline G4int GetCh() { return ch; }
    inline G4double GetPE() { return pe; }
    inline G4double GetLOPE() { return lope; }

};


// ====================================================================
// inline functions
// ====================================================================

// externally instanciated.
typedef G4THitsCollection<VetoHit> VetoHitsCollection;

extern G4Allocator<VetoHit> VetoHitAllocator; 

inline void* VetoHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) VetoHitAllocator.MallocSingle();
  return aHit;
}

inline void VetoHit::operator delete(void *aHit)
{
  VetoHitAllocator.FreeSingle((VetoHit*) aHit);
}

#endif
