#include "SimParticleSummary.h"

//......................................................................

SimParticleSummary::SimParticleSummary() 
{ 
    trackid = -1;
    parentid = -1;
    pdg = 0;
    iposflag = -1;
    fposflag = -1;
    for (unsigned int i=0; i<4; ++i) {
        momentum[i] = -1.e5;
        ipos[i] = -1.e5;
        fpos[i] = -1.e5;
    }
    for (unsigned int i=0; i<3; ++i) {
        dir[i] = -1.e5;
    }
    //length  = -1.e5;
}

//......................................................................

SimParticleSummary::~SimParticleSummary() 
{ 
}

//......................................................................

void SimParticleSummary::Clear(Option_t* option)
{    
}

//......................................................................

void SimParticleSummary::Print()
{
    ;
}

//......................................................................

ClassImp(SimParticleSummary)

////////////////////////////////////////////////////////////////////////

