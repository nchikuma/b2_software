#ifndef INGRID_TRACKING_ACTION_H
#define INGRID_TRACKING_ACTION_H 1

#include "G4UserTrackingAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"

#include "SimParticleSummary.h"

class TrackingAction : public G4UserTrackingAction {

public:
  TrackingAction(RunAction*,EventAction*);
  virtual ~TrackingAction();
   
  virtual void PreUserTrackingAction(const G4Track*);
  virtual void PostUserTrackingAction(const G4Track*);

private:
  RunAction* runaction;
  EventAction* eventaction;
  SimParticleSummary *simpart;
  int Tracking_Flag;

    float posi[3];
    float momi[3];

    double initE;
};

#endif
