#ifndef _DetectorDimension_H
#define _DetectorDimension_H

#include"Const.hh"

#include<iostream>
#include<sstream>
#include<fstream>
#include "TFile.h"
#include "TTree.h"
using namespace std;

class DetectorDimension{
private:
  TFile* f;
  TTree* t;
  double position_xy[2][8][80];
  double position_z[2][8][80];
  double wposition_xy[2][8][80];
  double wposition_z[2][8][80];
public:
  DetectorDimension();

  ~DetectorDimension(){};

  //Position in World
  void GetPos(int mod, int pln, int view, int ch, double *x, double *y, double *z);
  void GetPosInMod(int mod, int pln, int view, int ch, double *x, double *y, double *z);

  //Position in each mother volume
  void GetPosING(int mod, int pln, int view, int ch, double *x, double *y, double *z);
  void GetPosPM(int pln, int view, int ch, double *x, double *y, double *z);
  void GetPosWM(int mod, int pln, int view, int ch, int grid, double *x, double *y, double *z);

  // WaterModule channel rearrangement
  bool GetWMGridCh(int pln, int view, int ch, int *grid, int *gridch);
  bool GetWMGridCellID(int mod, int pln, int view, int ch, double posx, double posy, double posz, //Add "mod" | Kin, 2018/03/31
		       int* gridcell_id_x1, int* gridcell_id_x2, int* gridcell_id_y1, int* gridcell_id_y2);
  void GetScintiID(int mod, int view, int pln, int gridcell_id, int* cross_n, int* ch);

  int INGNumVetoPln;

  //For TwoDimRecon
  bool  GetReconPlnCh(int mod, int view, int pln, int ch, int axis,
		      int* reconpln, int* reconch);
  bool  GetRawPlnChGrid(int mod, int view, int reconpln, int reconch, int axis,
			int* pln, int* gridch, int* grid);
  int   GetChMax(int mod, int view, int pln, int axis);
  int   GetPlnMax(int mod, int view, int pln, int axis);
  bool  GetPos_TwoDimRecon(int mod, int view, int reconpln, int reconch, int axis,
			   double* posxy, double* posz);
  float GetScintiWidth(int mod, int view, int reconpln, int reconch, int axis);
  float GetScintiThick(int mod, int view, int reconpln, int reconch, int axis);



};
#endif
