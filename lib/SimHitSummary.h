#ifndef SIMHITSUMMARY_H
#define SIMHITSUMMARY_H
#include <iostream>
#include "TObject.h"


//......................................................................

class SimHitSummary : public TObject
{
public:
    SimHitSummary();
    virtual ~SimHitSummary();
    
    void Clear   (Option_t* option="");
    void Print();

public:
    
    float edeposit;             // Energy deposit (MeV)
    int   trackid;              // ID of GEANT4 track responsible for hit
    int   pdg;                  // PDG particle ID of GEANT4 track responsible for hit

private:

    ClassDef(SimHitSummary, 2) //  Hit Sim Summary
};

#endif // SIMHITSUMMARY_H
////////////////////////////////////////////////////////////////////////
