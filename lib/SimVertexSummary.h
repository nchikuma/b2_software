#ifndef SIMVERTEXSUMMARY_H
#define SIMVERTEXSUMMARY_H
#include <iostream>
#include <vector>
#include "TObject.h"

using namespace std;

//......................................................................

class SimVertexSummary : public TObject
{
public:
    SimVertexSummary();
    virtual ~SimVertexSummary();
    
    void Clear   (Option_t* option="");
    void Print(Option_t*) const;
    
public:
// Flux information
    int nutype;              // Neutrino type (generator-specific numbering
                             // scheme)
    float numomentum[3];     // Neutrino 3-momentum (GeV/c)
    float nuE;               // Neutrino Energy (GeV)
    float norm;							 // Weighting factor of Jnubeam				
    int ng;                  // # of parents (# of generations)
    vector<int> gpid;        // particle ID of each ancestor particle, gpid[0]=primary proton
    vector<int> gmec;        // particle production mechanism
    vector<float> gposx;     // interaction vertex point of each ancestor
    vector<float> gposy;     // interaction vertex point of each ancestor
    vector<float> gposz;     // interaction vertex point of each ancestor
    vector<float> gmomx;     // directional momentum
    vector<float> gmomy;     // directional momentum
    vector<float> gmomz;     // directional momentum
    vector<float> gcosbm;    // cosine of angle between ancestor and beam

// Neutrino interaction information
    int targeta;             // Atomic weight of target nucleus (TO BE FIXED)
    int targetz;             // Atomic number of target nucleus
    int targettype;          // Neutrino target type (generator-specific)
    float pfsurf;            // Fermi surface momentum (GeV/c)
    float vnuclini;          // Nuclear potential for the target initial state (GeV)
    int inttype;             // Neutrino interaction type (generator-specific
    int mod;                 // Interaction vertex module
    float xnu;               // interaction vertex x (from jnubeam)
    float ynu;               // interaction vertex y (from jnubeam)
    float znu;               // interaction vertex z (from jnubeam)
    float totcrsne;					 // Total cross section [10^-38 cm^2]

    double dirX;
    double dirY;
    double dirZ;

    double thetaX;
    double thetaY;
    double Phi_calc;

    double angleX;
    double angleY;
    double angleZ;

private:

    ClassDef(SimVertexSummary, 5) //  Simulation (generator + detector mc) 
                                    //  interaction vertex Summary
};

#endif // SIMVERTEXSUMMARY_H
////////////////////////////////////////////////////////////////////////
