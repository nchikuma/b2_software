#ifndef _INC_FLUX_TUNING
#define _INC_FLUX_TUNING

#define NMOD 10
#define NBIN 200
#define NFLAV 4

#define NCENTER 4

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdio>
#include <vector>
#include <algorithm>
using namespace std;

#include "TFile.h"
#include "TH1D.h"
#include "TGaxis.h"

#include "TRefArray.h"
#include "TRef.h"

const int Cnumu =1;
const int Cnumub=2;
const int Cnue  =3;
const int Cnueb =4;

//const int FDID_PM = 2; //Proton Module
//const int FDID_INGh = 3; //INGRID horizontal
//const int FDID_INGv = 4; //INGRID vertical
//const int FDID_Wall = 7; //Wall neutrino study

////////////////////////
class FluxTuning : public TObject
{
private:
	TH1D*  trationu[NMOD][NFLAV];
	TH1D*  tratio2nu[NMOD][NFLAV];
	double rwrationu[NMOD][NBIN][NFLAV];
	double tune_xbins[NBIN];
	double tuneratio;

	TFile* fTFRatio;
	int FDID;
	char Flavor[10];

	void tune ();
	
public:
	FluxTuning(int fdid, char* tunefile);
	virtual ~FluxTuning();

	double reweight(float Enu, int flav);

private:
	ClassDef(FluxTuning, 1);

};


#endif
