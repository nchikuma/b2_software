#ifndef __HEADER__
#define __HEADER__


//##### Standard C++ lib. ######
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include <math.h>
#include <stdlib.h>
using namespace std;
//##### Root Library ###########
#include <TROOT.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TEventList.h>
#include <TBranch.h>
#include <TSystem.h>
//##### INGRID Library #########
#include "EVENTSUMMARY.h"
#include "../bsd/BeamData/BeamData.h"

string data_calib_dir  = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_calib" ;
string data_dst_dir    = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_dst"   ;
string data_dst_dir2   = "/home/t2k/nchikuma/b2_data2/data_dst"   ;
string data_cosmic_dir = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/data_cosmic";
string data_bsd_dir    = "/home/t2k/nchikuma/b2_data/bsd";


const static float CTtimeBase     = 1.08e-6; //[sec]base time of CT
const static int   GapbwBunch     = 581 ;    //[nsec]
const static int   TDCOffset      = 300;     //[nsec]
const static int   bunch1st_cyc   = 4;
const static int   fIntTime       = 480;     //[nsec]
const static int   fRstTime       = 100;     //[nsec]
const static int   fExpBeamTime   = 200;     //[nsec]
const static int   beamontime     = 100;     //ontime


#endif
