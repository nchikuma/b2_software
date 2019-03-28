#include "SimVertexSummary.h"

//......................................................................

SimVertexSummary::SimVertexSummary() 
{ 
  nutype= -1;
  nuE   =  0;
  norm = -1.e5;
  for(int i=0;i<3;i++) numomentum[i] = -1.e5;

  mod   = -1;
  targeta = -1;
  targetz = -1;
  targettype = -1;
  vnuclini = -1.e5;
  pfsurf = -1.e5;
  inttype = -1;
  totcrsne = -1.e5;
  xnu = -1.e5;
  ynu = -1.e5;
  znu = -1.e5;

  dirX = 0;
  dirY = 0;
  dirZ = 0;

  thetaX = 0;
  thetaY = 0;
  Phi_calc = 0;

  angleX = 0;
  angleY = 0;
  angleZ = 0;
}

//......................................................................

SimVertexSummary::~SimVertexSummary() 
{ 
}

//......................................................................

void SimVertexSummary::Clear(Option_t* option)
{    
  nutype= -1;
  nuE   =  0;
  norm = -1.e5;
  for(int i=0;i<3;i++) numomentum[i] = -1.e5;

  mod   = -1;
  targeta = -1;
  targetz = -1;
  targettype = -1;
  vnuclini = -1.e5;
  pfsurf = -1.e5;
  inttype = -1;
  totcrsne = -1.e5;
  xnu = -1.e5;
  ynu = -1.e5;
  znu = -1.e5;

  dirX = 0;
  dirY = 0;
  dirZ = 0;

  thetaX = 0;
  thetaY = 0;
  Phi_calc = 0;

  angleX = 0;
  angleY = 0;
  angleZ = 0;

  ng = -1;
  gpid.clear();
  gmec.clear();
  gcosbm.clear();
  gposx.clear();
  gposy.clear();
  gposz.clear();
  gmomx.clear();
  gmomy.clear();
  gmomz.clear();
}

//......................................................................

void SimVertexSummary::Print(Option_t* option) const
{
}

//......................................................................

ClassImp(SimVertexSummary)

////////////////////////////////////////////////////////////////////////
