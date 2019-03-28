#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "RunAction.hh"
#include "Const.hh"
#include "DetectorDimension.hh"
#include <vector>

class G4Event;

class EventAction : public G4UserEventAction
{
public:

  EventAction(RunAction*); 
  EventAction(RunAction*,DetectorDimension*); 
  ~EventAction();

public:
  void BeginOfEventAction(const G4Event* anEvent);
  void EndOfEventAction(const G4Event* anEvent);

  inline void SetTrackID(int i) { TrackID = i; }
  int GetTrackID() { return TrackID; }
  inline void SetWriteFlag(int j) { Flag = j; }

private:
  RunAction* runaction;
  int TrackID;
  int Flag;
    
  SimHitSummary *simhitsum;
  HitSummary *hitsum;
  DetectorDimension *detdim;

};

#endif
