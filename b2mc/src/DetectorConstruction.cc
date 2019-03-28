#include <stdio.h>
#include <iostream>
#include <vector>

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4TwoVector.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"

DetectorConstruction::DetectorConstruction(int MODE)
{
  mode       = MODE;
}

DetectorConstruction::DetectorConstruction(int MODE, DetectorDimension* fdim)
{
  detdim  = fdim;
  mode    = MODE;
}

DetectorConstruction::~DetectorConstruction()
{
}


G4VPhysicalVolume* DetectorConstruction::Construct()
{

  //Initialise materials
  DefineMaterial();

  //Initialise detector parameters
  DefineStructures();

  int start_mod,stop_mod;
  if(mode==HORIZONTAL)           { start_mod=0; stop_mod=7; }
  else if(mode==VERTICAL)        { start_mod=7; stop_mod=14; }
  else if(mode==PROTON)          { start_mod=0; stop_mod=17; }
  else if(mode==INGWaterModule||
	  mode==INGWaterModuleBG){ start_mod=0; stop_mod=23; }
  else if(mode==B2Water||
	  mode==B2CH||
	  mode==B2ING)           { start_mod=0; stop_mod=36; }
  else if(mode==B2Wall        ||
	  mode==B2Ceiling     ||
	  mode==B2Floor       ||
	  mode==B2Right_Pillar||
	  mode==B2Left_Pillar)   { start_mod=0; stop_mod=36; }
  else                           { start_mod=0; stop_mod=36; }
  
  // ===========================================================
  // ========== Rotation Matrix  ===============================
  // ===========================================================
  //Rotation matrix for veto planes
  zrot  = new G4RotationMatrix(G4ThreeVector(0,0,1.),90.*degree);
  zrot2 = new G4RotationMatrix(G4ThreeVector(0,0,1.),-90.*degree);
  //rotations in reversal order
  xrot  = new G4RotationMatrix(G4ThreeVector(1,0,0), 90.*degree);
  xrot2 = new G4RotationMatrix(G4ThreeVector(1,0,0),-90.*degree);
  yrot  = new G4RotationMatrix(G4ThreeVector(0,1,0), (90.+sINGdeg)*degree);
  yrot2 = new G4RotationMatrix(G4ThreeVector(0,1,0),-(90.+sINGdeg)*degree);
  //for Y grid layer WaterModule
  rotgridY_v = new G4RotationMatrix(G4ThreeVector(0,0,1),180.*degree);
  rotgridY_v->rotateY(-90*degree);
  rotgridY_v->rotateX(-90*degree);
  rotgridY_h = new G4RotationMatrix(G4ThreeVector(1,0,0),-90.*degree);
  rotgridY_h->rotateY(-90*degree);
  //for Xgrid layer WaterModule
  rotgridX_v = new G4RotationMatrix(G4ThreeVector(0,1,0),-90.*degree);
  rotgridX_v->rotateX(-90*degree);
  rotgridX_h = new G4RotationMatrix(G4ThreeVector(0,0,1),180.*degree);
  rotgridX_h->rotateX(-90*degree);
  rotgridX_h->rotateY(-90*degree);


  // ========================================================================
  // ============= Mother Volumes ===========================================
  // ========================================================================
  //World volume
  experimentalHall_box = new G4Box("Hall",WorldSizeX,WorldSizeY,WorldSizeZ);
  worldLV              = new G4LogicalVolume(experimentalHall_box, Air, "world_log",0,0,0);
  worldPV              = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),worldLV,"world",0,false,0,check_MotherLV);

  //Hall Dirt volume
  HallDirtSLD    = new G4Tubs("halldirt", HallDirtRadiusMin, HallDirtRadiusMax,
			      HallDirtHeight, HallDirtSPhi, HallDirtDPhi);
  HallDirtLV     = new G4LogicalVolume(HallDirtSLD, Concrete, "HallDirt",0,0,0);
  HallDirtPV     = new G4PVPlacement(xrot,G4ThreeVector(HallDirtPosX,HallDirtPosY,HallDirtPosZ),
				     HallDirtLV,"HallDirt",worldLV,false,0,check_MotherLV);
  //HallDirtVisAtt = new G4VisAttributes(G4Color(1.,1.,1.));
  //HallDirtVisAtt -> SetForceWireframe(true);
  //HallDirtLV     -> SetVisAttributes(HallDirtVisAtt);

  //Horizontal volume (for the 7 horizontal modules)
  horizontalHall_box  = new G4Box("HorizontalHall",INGHMotherSizeX, INGHMotherSizeY, INGHMotherSizeZ);
  horizontalLV        = new G4LogicalVolume(horizontalHall_box, Air, "ingrid_h_log",0,0,0);
  horizontalPV        = new G4PVPlacement(0,G4ThreeVector(INGHMotherPosX,INGHMotherPosY,INGHMotherPosZ),
					  horizontalLV,"ingrid_h",worldLV,false,0,check_MotherLV);
  //Vertical volume (for the 7 vertical modules)
  verticalHall_box  = new G4Box("VerticalHall",INGVMotherSizeX, INGVMotherSizeY, INGVMotherSizeZ);
  verticalLV        = new G4LogicalVolume(verticalHall_box, Air, "ingrid_v_log",0,0,0);
  verticalPV        = new G4PVPlacement(0,G4ThreeVector(INGVMotherPosX,INGVMotherPosY,INGVMotherPosZ),
					verticalLV,"ingrid_v",worldLV,false,0,check_MotherLV);
  //Proton Moduel volume
  protonHall_box = new G4Box("ProtonHall",PMMotherSizeX, PMMotherSizeY, PMMotherSizeZ);
  ProtonLV       = new G4LogicalVolume(protonHall_box, Air, "pm_log",0,0,0);
  ProtonPV       = new G4PVPlacement(0,G4ThreeVector(PMMotherPosX,PMMotherPosY,PMMotherPosZ),
				     ProtonLV,"pm",worldLV,false,0,check_MotherLV);

  //B2 modules volume
  B2Hall_box = new G4Box("B2Hall",B2MotherSizeX, B2MotherSizeY, B2MotherSizeZ);
  B2HallLV   = new G4LogicalVolume(B2Hall_box, Air, "b2hall_log",0,0,0);
  B2HallPV   = new G4PVPlacement(0,G4ThreeVector(B2MotherPosX,B2MotherPosY,B2MotherPosZ),
				 B2HallLV,"b2hall",worldLV,false,0,check_MotherLV);

  // ========================================================================
  // ============= MODULE ===================================================
  // ========================================================================
  INGmodule_box  = new G4Box("INGModule" ,INGSizeX, INGSizeY, INGSizeZ);
  PMmodule_box   = new G4Box("PModule"   ,PMSizeX, PMSizeY, PMSizeZ);
  WMmodule_box   = new G4Box("WMModule"  ,WMSizeX, WMSizeY, WMSizeZ);
  //moduleLV[NUMMAX_MODULE]; //DEBUG

  for (int imod=start_mod;imod<stop_mod;imod++){
    bool notplace = false;
    for(int iarray=0;iarray<(int)sizeof(NotPlaceMod);iarray++){ //DEBUG
      if(NotPlaceMod[iarray]==imod){notplace = true; break;}
    }
    if(notplace) continue;

    //Set 7 INDRID horizontal modules:
    if( imod>=0 && imod<=MOD_INGRID_H ) {
      moduleLV[imod] = new G4LogicalVolume(INGmodule_box, Air, "module_log_ing_h");
      sprintf(name,"module%d",imod);
      posX = INGStart + INGSpace*imod;
      new G4PVPlacement(0,G4ThreeVector(posX,0.,0.),moduleLV[imod],name,horizontalLV,false,imod,check_moduleLV);
    }
    //Set 7 INGRID vertical modules
    else if( imod>MOD_INGRID_H && imod<=MOD_INGRID_V ) {
      moduleLV[imod] = new G4LogicalVolume(INGmodule_box, Air, "module_log_ing_v");
      sprintf(name,"module%d",imod);
      posY = INGStart + INGSpace*(imod-7);    
      new G4PVPlacement(0,G4ThreeVector(0.,posY,0.),moduleLV[imod],name,verticalLV,false,imod,check_moduleLV);     
    }
    //Set proton module (Temporary commented out)
    else if(imod==MOD_PM){
      sprintf(name,"module%d",imod);
      moduleLV[imod] = new G4LogicalVolume(PMmodule_box, Air, "module_log_pm");
      new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),moduleLV[imod],name,ProtonLV,false,imod,check_moduleLV);
    }

    //Set B2 modules ---------------------------------------------
    // WAGASCI water module
    else if(imod==MOD_B2_WM){
      sprintf(name,"module%d",imod);
      moduleLV[imod] = new G4LogicalVolume(WMmodule_box, Air, "module_log_b2wm");
      new G4PVPlacement(0,G4ThreeVector(B2WMPosX,B2WMPosY,B2WMPosZ),moduleLV[imod],name,B2HallLV,false,
			imod,check_moduleLV);
    }
    // Proton module @B2
    else if(imod==MOD_B2_CH){
      sprintf(name,"module%d",imod);
      moduleLV[imod] = new G4LogicalVolume(PMmodule_box, Air, "module__log_b2ch");
      new G4PVPlacement(0,G4ThreeVector(B2CHPosX,B2CHPosY,B2CHPosZ),moduleLV[imod],name,B2HallLV,false,
			imod,check_moduleLV);
    }
    // Downstream INGRID
    else if(imod==MOD_B2_INGRID){
      sprintf(name,"module%d",imod);
      moduleLV[imod] = new G4LogicalVolume(INGmodule_box, Air, "module_log_b2ing_d1");
      new G4PVPlacement(0,G4ThreeVector(B2INGPosX,B2INGPosY,B2INGPosZ),moduleLV[imod],name,B2HallLV,false,
			imod,check_moduleLV);
    }
  }

  // === IRON BLOCK =========================================  
  iron_block = new G4Box("Iron",INGIronXY/2.,INGIronXY/2.,INGIronThick/2.);
  ironLV     = new G4LogicalVolume(iron_block,Fe,"ironLV");
  //ironVisAtt = new G4VisAttributes(G4Color(0.7,0.,0.7)); // magenta
  //ironVisAtt -> SetForceSolid(true);
  //ironLV     -> SetVisAttributes(ironVisAtt);
  // === Water BLOCK =========================================  
  Water_box     = new G4Box("Water_box",WMWaterTargetSizeX/2.,WMWaterTargetSizeY/2.,WMWaterTargetSizeZ/2.);
  WaterLV       = new G4LogicalVolume(Water_box,Water,"WaterLV");
  //WaterVisAtt   = new G4VisAttributes(G4Colour(0.,0.,1.));
  //WaterLV       -> SetVisAttributes(WaterVisAtt);
  B2WaterLV     = new G4LogicalVolume(Water_box,Water,"B2WaterLV");
  //B2WaterVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.));
  //B2WaterLV     -> SetVisAttributes(B2WaterVisAtt);
  // === Air BLOCK with teh same size as Water Block ============  
  Air_box     = new G4Box("Air_box",WMWaterTargetSizeX/2.,WMWaterTargetSizeY/2.,WMWaterTargetSizeZ/2.);
  B2AirLV     = new G4LogicalVolume(Air_box,Air,"B2AirLV");
  // CH BLOCK fitting in each cell in WaterModule  ==============================  
  CH_box     = new G4Box("CH_box",WMCHCellTargetSizeX/2.,WMCHCellTargetSizeY/2.,WMCHCellTargetSizeZ/2.);
  B2CHLV     = new G4LogicalVolume(CH_box,Scinti,"B2CHLV");
  //B2CHVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.));
  //B2CHLV     -> SetVisAttributes(B2CHVisAtt);
  // Al BLOCK ============================================
  Al_box   = new G4Box("Al_box",WMAlTankSizeX/2.,WMAlTankSizeY/2.,WMAlTankSizeZ/2.);
  AlLV     = new G4LogicalVolume(Al_box,SUS304,"AlLV");
  //AlVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.));
  //AlLV     -> SetVisAttributes(AlVisAtt);

  // ========================================================================
  // ============= SCINTILLATORS FOR TRACKING PLANES ========================
  // ============= INGRID and ProtonModule           ========================
  // ========================================================================
  // INGRID 
  std::vector<G4TwoVector> vdim;
  std::vector<G4TwoVector> hdim;
  //SciBar type for Proton module
  std::vector<G4TwoVector> svdim;
  std::vector<G4TwoVector> shdim;
  //INGRID type for Proton module
  std::vector<G4TwoVector> ivdim;
  std::vector<G4TwoVector> ihdim;
  std::vector<G4ExtrudedSolid::ZSection> zsec;
  for(int iver=0; iver<8; iver++){
    vdim.push_back(  G4TwoVector( INGScintiVertexX[iver]  , INGScintiVertexY[iver]  ) );
    hdim.push_back(  G4TwoVector( INGScintiVertexY[7-iver], INGScintiVertexX[7-iver]) );
    svdim.push_back( G4TwoVector( PMScintiVertexX[iver]   , PMScintiVertexY[iver]   ) );
    shdim.push_back( G4TwoVector( PMScintiVertexY[7-iver] , PMScintiVertexX[7-iver] ) );
    ivdim.push_back( G4TwoVector( INGScintiVertexX[iver]  , INGScintiVertexY[iver]  ) );
    ihdim.push_back( G4TwoVector( INGScintiVertexY[7-iver], INGScintiVertexX[7-iver]) );
  }
  zsec.push_back( G4ExtrudedSolid::ZSection(-INGScintiLength/2., G4TwoVector(0*mm,0*mm), 1) );
  zsec.push_back( G4ExtrudedSolid::ZSection( INGScintiLength/2., G4TwoVector(0*mm,0*mm), 1) );

  //INGRID type for INGRID
  //vertical scintillator
  vscint_tmp       = new G4ExtrudedSolid("vscint_tmp", vdim, zsec);
  vsci_hole        = new G4EllipticalTube("vsci_hole",INGScintiHoleDia_a/2.,INGScintiHoleDia_b/2.,INGScintiLength/2.);
  vscint_int       = new G4SubtractionSolid("vscint_int",  vscint_tmp, vsci_hole);
  vscint_intLV     = new G4LogicalVolume(vscint_int,Scinti,"vscint_intLV");
  //vscint_intVisAtt = new G4VisAttributes(G4Color(0.,1.0,0.));
  //vscint_intVisAtt -> SetForceSolid(true);
  //vscint_intLV     -> SetVisAttributes(vscint_intVisAtt);
  //horizontal scintillator
  hscint_tmp       = new G4ExtrudedSolid("hscint_tmp", hdim, zsec);
  hsci_hole        = new G4EllipticalTube("hsci_hole",INGScintiHoleDia_b/2.,INGScintiHoleDia_a/2.,INGScintiLength/2.);
  hscint_int       = new G4SubtractionSolid("hscint_int",  hscint_tmp, hsci_hole);
  hscint_intLV     = new G4LogicalVolume(hscint_int,Scinti,"hscint_intLV");
  //hscint_intVisAtt = new G4VisAttributes(G4Color(0.,1.,0.));
  //hscint_intVisAtt -> SetForceSolid(true);
  //hscint_intLV     -> SetVisAttributes(hscint_intVisAtt);
  //SciBar type for proton module
  //vertical scintillator
  vscint2_tmp       = new G4ExtrudedSolid("vscint2_tmp", svdim, zsec);
  vsci_hole2        = new G4EllipticalTube("vsci_hole2",PMScintiHoleDia/2.,PMScintiHoleDia/2.,PMScintiLength/2.);
  vscint2_int       = new G4SubtractionSolid("vscint2_int", vscint2_tmp, vsci_hole2);
  vscint2_intLV     = new G4LogicalVolume(vscint2_int,Scinti,"vscint2_intLV");
  //vscint2_intVisAtt = new G4VisAttributes(G4Color(0.,1.0,0.));
  //vscint2_intVisAtt -> SetForceSolid(true);
  //vscint2_intLV     -> SetVisAttributes(vscint2_intVisAtt);
  //horizontal scintillator
  hscint2_tmp       = new G4ExtrudedSolid("hscint2_tmp", shdim, zsec);
  hsci_hole2        = new G4EllipticalTube("hsci_hole2",PMScintiHoleDia/2.,PMScintiHoleDia/2.,PMScintiLength/2.);
  hscint2_int       = new G4SubtractionSolid("hscint2_int",  hscint2_tmp, hsci_hole2);
  hscint2_intLV     = new G4LogicalVolume(hscint2_int,Scinti,"hscint2_intLV");
  //hscint2_intVisAtt = new G4VisAttributes(G4Color(0.,1.,0.));
  //hscint2_intVisAtt -> SetForceSolid(true);
  //hscint2_intLV     -> SetVisAttributes(hscint2_intVisAtt);
  //INGRID type for proton module
  //vertical scintillator
  vscint3_tmp       = new G4ExtrudedSolid("vscint3_tmp", ivdim, zsec);
  vsci_hole3        = new G4EllipticalTube("vsci_hole3",INGScintiHoleDia_a/2.,INGScintiHoleDia_b/2.,INGScintiLength/2.);
  vscint3_int       = new G4SubtractionSolid("vscint3_int",  vscint3_tmp, vsci_hole3);
  vscint3_intLV     = new G4LogicalVolume(vscint3_int,Scinti,"vscint3_intLV");
  //vscint3_intVisAtt = new G4VisAttributes(G4Color(0.,1.0,0.));
  //vscint3_intVisAtt -> SetForceSolid(true);
  //vscint3_intLV     -> SetVisAttributes(vscint3_intVisAtt);
  //horizontal scintillato
  hscint3_tmp       = new G4ExtrudedSolid("hscint3_tmp", ihdim, zsec);
  hsci_hole3        = new G4EllipticalTube("hsci_hole3",INGScintiHoleDia_b/2.,INGScintiHoleDia_a/2.,INGScintiLength/2.);
  hscint3_int       = new G4SubtractionSolid("hscint3_int",  hscint3_tmp, hsci_hole3);
  hscint3_intLV     = new G4LogicalVolume(hscint3_int,Scinti,"hscint3_intLV");
  //hscint3_intVisAtt = new G4VisAttributes(G4Color(0.,1.,0.));
  //hscint3_intVisAtt -> SetForceSolid(true);
  //hscint3_intLV     -> SetVisAttributes(hscint3_intVisAtt);


  // ========================================================================
  // ========== SCINTILLATORS FOR WaterModule ===============================
  // ========================================================================
  //Scintillator dimension for WaterModule
  std::vector<G4TwoVector> wvdim;
  std::vector<G4TwoVector> whdim;
  std::vector<G4TwoVector> wvdim2;
  std::vector<G4TwoVector> whdim2;
  std::vector<G4ExtrudedSolid::ZSection> wzsec;
  for(int iver=0; iver<8; iver++){
    wvdim .push_back( G4TwoVector( WMScintiVertexX[iver]  , WMScintiVertexY[iver]  ) );
    whdim .push_back( G4TwoVector( WMScintiVertexY[7-iver], WMScintiVertexX[7-iver]) );
    wvdim2.push_back( G4TwoVector( WMScintiVertexX2[iver] , WMScintiVertexY[iver]  ) );
    whdim2.push_back( G4TwoVector( WMScintiVertexY[7-iver], WMScintiVertexX2[7-iver]));
  }
  wzsec.push_back( G4ExtrudedSolid::ZSection(-WMScintiLength/2., G4TwoVector(0*mm,0*mm), 1) );
  wzsec.push_back( G4ExtrudedSolid::ZSection( WMScintiLength/2., G4TwoVector(0*mm,0*mm), 1) );

  //vertical scintillator
  vwscint_tmp    = new G4ExtrudedSolid("vwscint_tmp", wvdim, wzsec);
  wmhole_v       = new G4Box("wmhole_v",WMScintiHoleWidth/2.,WMScintiHoleThick/2.,WMScintiHoleLength/2.);
  vwscint_int    = new G4SubtractionSolid("vwscint_int", vwscint_tmp, wmhole_v, 0,
					  G4ThreeVector(-WMScintiHoleShift1,-WMScintiHoleShift2,0*mm));
  vwscint_intLV     = new G4LogicalVolume(vwscint_int,Scinti,"vwscint_intLV");
  //vwscint_intVisAtt = new G4VisAttributes(G4Color(0.,1.,0.));
  //vwscint_intVisAtt -> SetForceSolid(true);
  //vwscint_intLV     -> SetVisAttributes(vwscint_intVisAtt);
  //horizontal scintillator
  hwscint_tmp       = new G4ExtrudedSolid("hwscint_tmp", whdim, wzsec);
  wmhole_h          = new G4Box("wmhole_h",WMScintiHoleThick/2.,WMScintiHoleWidth/2.,WMScintiHoleLength/2.);
  hwscint_int       = new G4SubtractionSolid("hwscint_int", hwscint_tmp, wmhole_h, 0,
					     G4ThreeVector(WMScintiHoleShift2,-WMScintiHoleShift1,0*mm));
  hwscint_intLV     = new G4LogicalVolume(hwscint_int,Scinti,"hwscint_intLV");
  //hwscint_intVisAtt = new G4VisAttributes(G4Color(0.,1.,0.));
  //hwscint_intVisAtt -> SetForceSolid(true);
  //hwscint_intLV     -> SetVisAttributes(hwscint_intVisAtt);


  //vertical grid scintillator
  vwgridscint_tmp  = new G4ExtrudedSolid("vwgridscint_tmp", wvdim2, wzsec);
  vwgridscint_int  = new G4SubtractionSolid("vwgridscint_int", vwgridscint_tmp, wmhole_v, 0,
					    G4ThreeVector(-WMScintiHoleShift1,-WMScintiHoleShift2,0*mm));
  if(mode==INGWaterModule){
    wmslit_v         = new G4Box("wmslit_v",WMScintiSlitWidth/2.,WMScintiSlitThick/2.,WMScintiSlitLength2/2.); 
  }
  //if(mode==B2Water){
  else{
    wmslit_v         = new G4Box("wmslit_v",WMScintiSlitWidth/2.,WMScintiSlitThick/2.,WMScintiSlitLength/2.); 
  }
  for(int islit=0;islit<WMNumGridCh;islit++){
    vwgridscint_tmp2 = new G4SubtractionSolid("", vwgridscint_int, wmslit_v, 0, 
					      G4ThreeVector(WMScintiSlitShift,0*mm,WMGridChStart+islit*WMScintiSlitStep) );
    vwgridscint_int  = vwgridscint_tmp2;
  }
  vwgridscint_intLV     = new G4LogicalVolume(vwgridscint_int,Scinti,"vwgridscint_intLV");
  //vwgridscint_intVisAtt = new G4VisAttributes(G4Color(0.,1.,0.));
  //vwgridscint_intVisAtt -> SetForceSolid(true);
  //vwgridscint_intLV     -> SetVisAttributes(vwgridscint_intVisAtt);

  //horizontal grid scintillator
  hwgridscint_tmp = new G4ExtrudedSolid("hwgridscint_tmp", whdim2, wzsec);
  hwgridscint_int = new G4SubtractionSolid("hwgridscint_int", hwgridscint_tmp, wmhole_h, 0, 
					   G4ThreeVector(WMScintiHoleShift2,-WMScintiHoleShift1,0*mm));
  if(mode==INGWaterModule){
    wmslit_h        = new G4Box("wmslit_h",WMScintiSlitThick/2.,WMScintiSlitWidth/2.,WMScintiSlitLength2/2.);
  }
  //if(mode==B2Water){
  else{
    wmslit_h        = new G4Box("wmslit_h",WMScintiSlitThick/2.,WMScintiSlitWidth/2.,WMScintiSlitLength/2.);
  }
  for(int islit=0;islit<WMNumGridCh;islit++){
    hwgridscint_tmp2 = new G4SubtractionSolid("", hwgridscint_int, wmslit_h, 0, 
					      G4ThreeVector(0*mm,WMScintiSlitShift,WMGridChStart+islit*WMScintiSlitStep));
    hwgridscint_int = hwgridscint_tmp2;
  }
  hwgridscint_intLV     = new G4LogicalVolume(hwgridscint_int,Scinti,"hwgridscint_intLV");
  //hwgridscint_intVisAtt = new G4VisAttributes(G4Color(0.,1.,0.));
  //hwgridscint_intVisAtt -> SetForceSolid(true);
  //hwgridscint_intLV     -> SetVisAttributes(hwgridscint_intVisAtt);


  // ========================================================================
  // ========== SCINTILLATORS FOR VETO PLANES ===============================
  // ========================================================================
  //Long Veto plane for INGRID
  Lveto_box = new G4Box("Lveto_box",INGLVetoLength/2.,INGLVetoThick/2.,INGLVetoWidth/2.);
  LvetoLV   = new G4LogicalVolume(Lveto_box,Scinti,"LvetoLV");
  //LvetoLV   -> SetVisAttributes(G4VisAttributes::Invisible);
  //Short Veto plane for INGRID
  Sveto_box = new G4Box("Sveto_box",INGSVetoLength/2.,INGSVetoThick/2.,INGSVetoWidth/2.);
  SvetoLV   = new G4LogicalVolume(Sveto_box,Scinti,"SvetoLV");
  //SvetoLV   -> SetVisAttributes(G4VisAttributes::Invisible);
  //Proton module Long Veto plane
  PLveto_box = new G4Box("PLveto_box",PMLVetoLength/2.,PMLVetoThick/2.,PMLVetoWidth/2.);
  PLvetoLV   = new G4LogicalVolume(PLveto_box,Scinti,"PLvetoLV");
  //PLvetoLV   -> SetVisAttributes(G4VisAttributes::Invisible);
  //Proton module Short Veto plane
  PSveto_box = new G4Box("PSveto_box",PMSVetoLength/2.,PMSVetoThick/2.,PMSVetoWidth/2.);
  PSvetoLV   = new G4LogicalVolume(PSveto_box,Scinti,"PSvetoLV");
  //PSvetoLV   -> SetVisAttributes(G4VisAttributes::Invisible);


  // ========================================================================
  // ========== B2 Background sources =======================================
  // ========================================================================
  // ========================================================================
  //Ceiling
  B2ceiling_tube = new G4Tubs("B2Ceiling_box",B2CeilingRadiusMin, B2CeilingRadiusMax, B2CeilingSizeX, B2CeilingSizeY, B2CeilingSizeZ);
  B2Vbox = new G4Box("Vbox",B2VboxSizeX, B2VboxSizeY, B2VboxSizeZ);
  B2Ceiling_tube = new G4SubtractionSolid("B2Ceiling", B2ceiling_tube, B2Vbox, 0, G4ThreeVector(B2TransPosX,B2TransPosY,B2TransPosZ));
  B2CeilingLV = new G4LogicalVolume(B2Ceiling_tube, Concrete, "b2ceiling_log");
  B2CeilingPV = new G4PVPlacement(xrot,G4ThreeVector(B2CeilingPosX,B2CeilingPosY,B2CeilingPosZ),B2CeilingLV,"b2ceiling",worldLV,false,0,true);
  //G4VisAttributes* B2CeilingVisArt = new G4VisAttributes(G4Colour(0., 1.0, 1.0));
  //B2CeilingVisArt->SetForceSolid(true);
  //B2CeilingLV->SetVisAttributes(B2CeilingVisArt);
  //Floor
  B2Floor_tube = new G4Tubs("B2Floor",B2FloorRadiusMin, B2FloorRadiusMax, B2FloorSizeX, B2FloorSizeY, B2FloorSizeZ);
  B2FloorLV    = new G4LogicalVolume(B2Floor_tube, Concrete, "b2floor_log");
  B2FloorPV    =   new G4PVPlacement(xrot,G4ThreeVector(B2FloorPosX,B2FloorPosY,B2FloorPosZ),B2FloorLV,"b2floor",worldLV,false,0,true);
  //G4VisAttributes* B2FloorVisArt = new G4VisAttributes(G4Colour(0.7, 1., 0.4));
  //B2FloorVisArt->SetForceSolid(true);
  //B2FloorLV->SetVisAttributes(B2FloorVisArt);
  //Right pillar
  B2Right_Pillar_box = new G4Box("B2Right_Pillar",B2Right_PillarSizeX,B2Right_PillarSizeY,B2Right_PillarSizeZ);
  B2Right_PillarLV = new G4LogicalVolume(B2Right_Pillar_box, Concrete, "b2right_pillarlog",0,0,0);
  B2Right_PillarPV = new G4PVPlacement(0,G4ThreeVector(B2Right_PillarPosX,B2Right_PillarPosY,B2Right_PillarPosZ),B2Right_PillarLV,"b2right_pillar",worldLV,false,0,true);
  //G4VisAttributes* B2Right_PillarVisArt = new G4VisAttributes(G4Colour(1., 0., 1.));
  //B2Right_PillarVisArt->SetForceSolid(true);
  //B2Right_PillarLV->SetVisAttributes(B2Right_PillarVisArt);
  //Left pillar
  B2Left_Pillar_box = new G4Box("B2Left_Pillar",B2Left_PillarSizeX,B2Left_PillarSizeY,B2Left_PillarSizeZ);
  B2Left_PillarLV = new G4LogicalVolume(B2Left_Pillar_box, Concrete, "b2left_pillarlog",0,0,0);
  B2Left_PillarPV = new G4PVPlacement(0,G4ThreeVector(B2Left_PillarPosX,B2Left_PillarPosY,B2Left_PillarPosZ),B2Left_PillarLV,"b2left_pillar",worldLV,false,0,true);
  //G4VisAttributes* B2Left_PillarVisArt = new G4VisAttributes(G4Colour(1., 0., 1.));
  //B2Left_PillarVisArt->SetForceSolid(true);
  //B2Left_PillarLV->SetVisAttributes(B2Left_PillarVisArt);
  //Small left pillar
  B2Small_Pillar_box = new G4Box("B2Small_Pillar",B2Small_PillarSizeX,B2Small_PillarSizeY,B2Small_PillarSizeZ);
  B2Small_PillarLV = new G4LogicalVolume(B2Small_Pillar_box, Concrete, "b2small_pillarlog",0,0,0);
  B2Small_PillarPV = new G4PVPlacement(0,G4ThreeVector(B2Small_PillarPosX,B2Small_PillarPosY,B2Small_PillarPosZ),B2Small_PillarLV,"b2small_pillar",worldLV,false,0,true);
  //G4VisAttributes* B2Small_PillarVisArt = new G4VisAttributes(G4Colour(1., 0., 1.));
  //B2Small_PillarVisArt->SetForceSolid(true);
  //B2Small_PillarLV->SetVisAttributes(B2Small_PillarVisArt);




  // ====================================================================================================================
  // =================================  POSITIONNING OF ALL THE ELEMENTS ================================================
  // ====================================================================================================================

  for(int imod=start_mod;imod<stop_mod;imod++){
    bool notplace = false;
    for(int iarray=0;iarray<(int)sizeof(NotPlaceMod);iarray++){ //DEBUG
      if(NotPlaceMod[iarray]==imod){notplace = true; break;}
    }
    if(notplace) continue;


    // ============ INGRID =================
    if(imod<=MOD_INGRID_H)      this->PlaceINGRID(moduleLV[imod],imod,imod);
    else if(imod<=MOD_INGRID_V) this->PlaceINGRID(moduleLV[imod],imod,imod);

    // ============ ProtonModule =================
    //Temporary commented out
    //else if(imod==MOD_PM&&mode!=INGWaterModule) this->PlaceProtonModule(moduleLV[imod],imod);

    // ============ On-axis Water Module =================
    //else if(imod==MOD_ONAXIS_WM&&mode!=PROTON) this->PlaceWaterModule(moduleLV[MOD_PM],WaterLV,imod);
    //else if(imod==MOD_ONAXIS_WM) this->PlaceWaterModule(moduleLV[MOD_PM],WaterLV,imod);

    //else if(imod==MOD_ONAXIS_WM && mode==INGWaterModuleBG){
    //  new G4PVPlacement(0,G4ThreeVector(0,0,0),AlLV   ,"alflame",moduleLV[MOD_PM],false,0,0); //al frame
    //  new G4PVPlacement(0,G4ThreeVector(0,0,0),WaterLV,"watertarget",AlLV,false,0,0); //water target
    //}

    // ============      B2 modules      =================
    // WAGASCI water module
    else if(imod==MOD_B2_WM){
      new G4PVPlacement(0,G4ThreeVector(0,0,0),AlLV ,"alflame",moduleLV[imod],false,0,check_WM); //al frame
      //this->PlaceWaterModule(moduleLV[imod],B2WaterLV,imod);
      this->PlaceWaterModule(AlLV,B2WaterLV,imod);
    }

    // Proton module @B2
    else if(imod==MOD_PM) this->PlaceProtonModule(moduleLV[MOD_B2_CH],imod);

    // B2-INGRID
    else if(imod==MOD_B2_INGRID) this->PlaceINGRID(moduleLV[imod],8,imod);

  }//for imod

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  G4String hlayerSDname = "/hlayerSD";
  ahlayerSD = new HLayerSD( hlayerSDname ,detdim);
  SDman->AddNewDetector( ahlayerSD );
  hscint_intLV      -> SetSensitiveDetector( ahlayerSD );//INGRID
  hscint2_intLV     -> SetSensitiveDetector( ahlayerSD );//SciBar type in PM
  hscint3_intLV     -> SetSensitiveDetector( ahlayerSD );//INGRID type in PM 
  hwscint_intLV     -> SetSensitiveDetector( ahlayerSD ); //wagasci scinti h
  hwgridscint_intLV -> SetSensitiveDetector( ahlayerSD ); //wagasci scinti h

  G4String vlayerSDname = "/vlayerSD";
  avlayerSD = new VLayerSD( vlayerSDname , detdim);
  SDman->AddNewDetector( avlayerSD );
  vscint_intLV      -> SetSensitiveDetector( avlayerSD );//INGRID
  vscint2_intLV     -> SetSensitiveDetector( avlayerSD );//SciBar type
  vscint3_intLV     -> SetSensitiveDetector( avlayerSD );//INGRID type
  vwscint_intLV     -> SetSensitiveDetector( avlayerSD );//wagasci scinti v
  vwgridscint_intLV -> SetSensitiveDetector( avlayerSD );//wagasci scinti v

  G4String vetoSDname = "/vetoSD";
  avetoSD = new VetoSD( vetoSDname , detdim);
  SDman->AddNewDetector( avetoSD );

  SvetoLV  -> SetSensitiveDetector ( avetoSD ); //Veto for INGRID
  LvetoLV  -> SetSensitiveDetector ( avetoSD ); //Veto for INGRID
  PSvetoLV -> SetSensitiveDetector( avetoSD );//Veto for Proton Module
  PLvetoLV -> SetSensitiveDetector( avetoSD );//Veto for Proton Module

  return worldPV;

}




//___________________________________________________________________________________________________________

void DetectorConstruction::DefineMaterial()
{
  G4double a_atm;  // atomic mass
  G4double z_atm;  // atomic number
  G4double density;
  G4String name_atm, symbol;
  G4int nel;

  a_atm = 14.01*g/mole;
  G4Element* elN  = new G4Element(name_atm="Nitrogen",  symbol="N", z_atm=7.,  a_atm);
  a_atm = 16.00*g/mole;
  G4Element* elO  = new G4Element(name_atm="Oxigen",    symbol="O", z_atm=8.,  a_atm);
  a_atm = 1.01*g/mole;
  G4Element* elH  = new G4Element(name_atm="Hydrogen",  symbol="H", z_atm=1.,  a_atm);
  a_atm = 12.01*g/mole;
  G4Element* elC  = new G4Element(name_atm="Carbon",    symbol="C", z_atm=6.,  a_atm);
  a_atm = 28.1*g/mole;
  G4Element* elSi = new G4Element(name_atm="Silicon",   symbol="Si", z_atm=14., a_atm);
  a_atm = 40.1*g/mole;
  G4Element* elCa = new G4Element(name_atm="Calusium",  symbol="Ca", z_atm=20., a_atm);
  a_atm = 23.0*g/mole;
  G4Element* elNa = new G4Element(name_atm="Sodium",    symbol="Na", z_atm=11., a_atm);
  a_atm = 55.8*g/mole;
  G4Element* elFe = new G4Element(name_atm="Iron",      symbol="Fe", z_atm=26., a_atm);
  a_atm = 27.0*g/mole;
  G4Element* elAl = new G4Element(name_atm="Aluminium", symbol="Al", z_atm=13., a_atm);
  a_atm = 58.69*g/mole;
  G4Element* elNi = new G4Element(name_atm="Nickel",    symbol="Ni", z_atm=28., a_atm);
  a_atm = 51.99*g/mole;
  G4Element* elCr = new G4Element(name_atm="Chromium",  symbol="Cr", z_atm=24., a_atm);

  //Air
  density = 1.29*mg/cm3;
  Air = new G4Material(name_atm="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

  //Iron
  a_atm = 55.845*g/mole;
  density = 7.86*g/cm3;
  Fe = new G4Material(name_atm="Iron", z_atm=26., a_atm, density);

  //Water
  density = 1.000*g/cm3;
  Water = new G4Material(name_atm="Water",density,nel=2);
  Water->AddElement(elH,2);
  Water->AddElement(elO,1);

  //Al
  a_atm = 26.98*g/mole;
  density = 2.7*g/cm3;
  Al = new G4Material(name_atm="Aluminum", z_atm=13., a_atm, density);

  //Scintillator
  density = 1.032*g/cm3;
  Scinti = new G4Material(name_atm="Scintillator", density, nel=2);
  Scinti->AddElement(elC, 9);
  Scinti->AddElement(elH, 10);

  //Concrete
  density = 2.2*g/cm3;
  Concrete = new G4Material(name_atm="Concrete", density, nel=6);
  Concrete->AddElement(elO,  .53);
  Concrete->AddElement(elSi, .335);
  Concrete->AddElement(elCa, 0.06);
  Concrete->AddElement(elNa, 0.015);
  Concrete->AddElement(elFe, 0.02);
  Concrete->AddElement(elAl, 0.04);

  //SUS304 for Water tank
  density = 7.93*g/cm3;
  SUS304 = new G4Material(name_atm="SUS304", density, nel=3);
  SUS304->AddElement(elFe, 0.72);
  SUS304->AddElement(elCr, 0.19);
  SUS304->AddElement(elNi, 0.09);
}

//___________________________________________________________________________________________________________
void DetectorConstruction::PlaceINGRID(G4LogicalVolume* motherLV, int ingmod, int mod){

  // 11 planes of scintillator per module-----------------------------------------
  for(int ipln=0;ipln<INGNumPln;ipln++){
    for(int ich=0;ich<INGNumCh;ich++){ 

      //vertical scintillator
      detdim->GetPosING(ingmod,ipln,1,ich,&x,&y,&z);
      posX = x*mm; 
      posY = y*mm;
      posZ = z*mm;
      sprintf(name,"vlayer[%d][%d][%d]",mod,ipln,ich);
#ifdef DEBUG_DETCONST
      G4cout << name << " {" << posX << "," << posY << "," << posZ << "}" << G4endl;
#endif
      new G4PVPlacement(xrot,G4ThreeVector(posX,0.,posZ),vscint_intLV,name,
			motherLV,false,ich+ipln*NORMPLN+mod*NORMMOD,check_ING);

      //horizontal scintillator
      detdim->GetPosING(ingmod,ipln,0,ich,&x,&y,&z);
      posX = x*mm; 
      posY = y*mm;
      posZ = z*mm;
      sprintf(name,"hlayer[%d][%d][%d]",mod,ipln,ich);
#ifdef DEBUG_DETCONST
      G4cout << name << " {" << posX << "," << posY << "," << posZ << "}" << G4endl;
#endif
      new G4PVPlacement(yrot,G4ThreeVector(0.,posY,posZ),hscint_intLV,name,
			motherLV,false,ich+ipln*NORMPLN+mod*NORMMOD,check_ING); 
    }
    // 9 iron-blocs per module -----------------------------------------------------
    if(ipln<INGNumPln-2){
      detdim->GetPosING(ingmod,ipln,0,0,&x,&y,&z);
      posZ = z*mm;
      sprintf(name,"iron[%d][%d]",mod,ipln);

      posZ = posZ + (INGIronStart-INGPlnStart);
      if(mod==3) posZ = posZ - 1.*mm;
      new G4PVPlacement(0,G4ThreeVector(0.,0.,posZ),ironLV,name,motherLV,false,0,check_ING);   
    }
  }

  //4 veto planes in the first modules (mod#0,mod#7)------------------
  //3 veto planes on the other 12 modules-------------------------
  if (ingmod>=0 && ingmod<14){
    // 22 veto-bars per veto-plane
    for(int ich=0;ich<INGNumVetoCh;ich++){
      posZ = INGVetoStartZ + ich*INGSVetoWidth;

      char vetoname[4][22];
      sprintf(vetoname[0],"veto[%d][0][%d]",mod,ich);
      sprintf(vetoname[1],"veto[%d][1][%d]",mod,ich);
      sprintf(vetoname[2],"veto[%d][2][%d]",mod,ich);
      sprintf(vetoname[3],"veto[%d][3][%d]",mod,ich);

      if(ingmod>=0&&ingmod<14) new G4PVPlacement(0    ,G4ThreeVector(INGVetoPos1X,INGVetoPos1Y,posZ),
						 LvetoLV,vetoname[3],motherLV,false,NORMMOD*mod+ich+3*NORMPLN,check_ING); // UP
      if(ingmod<=7)            new G4PVPlacement(0    ,G4ThreeVector(INGVetoPos2X,INGVetoPos2Y,posZ),
						 SvetoLV,vetoname[2],motherLV,false,NORMMOD*mod+ich+2*NORMPLN,check_ING); // DOWN
      if(ingmod>=0&&ingmod<14) new G4PVPlacement(zrot ,G4ThreeVector(INGVetoPos3X,INGVetoPos3Y,posZ),
						 LvetoLV,vetoname[1],motherLV,false,NORMMOD*mod+ich+1*NORMPLN,check_ING); // LEFT
      if(ingmod==0||ingmod>=7) new G4PVPlacement(zrot ,G4ThreeVector(INGVetoPos4X,INGVetoPos4Y,posZ),
						 LvetoLV,vetoname[0],motherLV,false,NORMMOD*mod+ich+0*NORMPLN,check_ING); // RIGTH
    } 
  }
}

//___________________________________________________________________________________________________________
void DetectorConstruction::PlaceProtonModule(G4LogicalVolume* motherLV, int mod){

  // 18 tracking planes --------------------------------------------------------------------------
  for(int ipln=0;ipln<PMNumPln;ipln++){
    for(int ich=0;ich<PMNumCh;ich++){
      //First plane similar to an INGRID module plane arrangement
      if(ipln==0 && ich<PMNumCh1 ){

        //vertical scintillator
        detdim->GetPosPM(ipln,1,ich,&x,&y,&z);
        posX = x*mm;
        posY = y*mm;
        posZ = z*mm;
        sprintf(name,"vlayer[%d][%d][%d]",mod,ipln,ich);
        new G4PVPlacement(xrot,G4ThreeVector(posX,0.,posZ),vscint3_intLV,name,motherLV,false,
			  ich+ipln*NORMPLN+mod*NORMMOD,check_PM);

        //horizontal scintillator
        detdim->GetPosPM(ipln,0,ich,&x,&y,&z);
        posX = x*mm;
        posY = y*mm;
        posZ = z*mm;
        sprintf(name,"hlayer[%d][%d][%d]",mod,ipln,ich);
        new G4PVPlacement(yrot,G4ThreeVector(0.,posY,posZ),hscint3_intLV,name,motherLV,false,
			  ich+ipln*NORMPLN+mod*NORMMOD,check_PM);
      }
      else if(ipln>0){

        //vertical scintillator
        detdim->GetPosPM(ipln,1,ich,&x,&y,&z);
        posX = x*mm;
        posY = y*mm;
        posZ = z*mm;
        sprintf(name,"vlayer[%d][%d][%d]",mod,ipln,ich);
        if(ich>=PMNumCh_side && ich<PMNumCh_side+PMNumCh_mid)
          new G4PVPlacement(xrot,G4ThreeVector(posX,0.,posZ),vscint2_intLV,name,motherLV,false,
			    ich+ipln*NORMPLN+mod*NORMMOD,check_PM);
        else
          new G4PVPlacement(xrot,G4ThreeVector(posX,0.,posZ),vscint3_intLV,name,motherLV,false,
			    ich+ipln*NORMPLN+mod*NORMMOD,check_PM);

        //horizontal scintillator
        detdim->GetPosPM(ipln,0,ich,&x,&y,&z);
        posX = x*mm;
        posY = y*mm;
        posZ = z*mm;
        sprintf(name,"hlayer[%d][%d][%d]",mod,ipln,ich);
        if(ich>=PMNumCh_side && ich<PMNumCh_side+PMNumCh_mid)
          new G4PVPlacement(yrot,G4ThreeVector(0,posY,posZ),hscint2_intLV,name,motherLV,false,
			    ich+ipln*NORMPLN+mod*NORMMOD,check_PM);
        else
          new G4PVPlacement(yrot,G4ThreeVector(0,posY,posZ),hscint3_intLV,name,motherLV,false,
			    ich+ipln*NORMPLN+mod*NORMMOD,check_PM);
      }
    }
  }
  //Veto Planes --------------------------------------------------------------------------------
  for(int ich=0;ich<PMNumVetoCh;ich++){
    posZ = PMVetoStartZ + ich*PMSVetoWidth;

    char vetoname[4][22];
    sprintf(vetoname[0],"veto[%d][0][%d]",mod,ich);
    sprintf(vetoname[1],"veto[%d][1][%d]",mod,ich);
    sprintf(vetoname[2],"veto[%d][2][%d]",mod,ich);
    sprintf(vetoname[3],"veto[%d][3][%d]",mod,ich);

    new G4PVPlacement(0    ,G4ThreeVector(PMVetoPos1X,PMVetoPos1Y,posZ),PLvetoLV,vetoname[3],motherLV,false,  // UP
		      ich+2*NORMPLN+mod*NORMMOD,check_PM);
    new G4PVPlacement(0    ,G4ThreeVector(PMVetoPos2X,PMVetoPos2Y,posZ),PSvetoLV,vetoname[2],motherLV,false,  // DOWN
		      ich+0*NORMPLN+mod*NORMMOD,check_PM);
    new G4PVPlacement(zrot ,G4ThreeVector(PMVetoPos3X,PMVetoPos3Y,posZ),PLvetoLV,vetoname[1],motherLV,false,  // LEFT
		      ich+1*NORMPLN+mod*NORMMOD,check_PM);
    new G4PVPlacement(zrot ,G4ThreeVector(PMVetoPos4X,PMVetoPos4Y,posZ),PLvetoLV,vetoname[0],motherLV,false,  // RIGTH
		      ich+3*NORMPLN+mod*NORMMOD,check_PM);
  }   

}

//___________________________________________________________________________________________________________
void DetectorConstruction::PlaceWaterModule(G4LogicalVolume* motherLV, G4LogicalVolume* targetLV, int mod){

  sprintf(name,"waterLV[%d]",mod);
  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),targetLV,name,motherLV,false,0,check_WM); 

  for(int ipln=0;ipln<WMNumPln;ipln++){
    for(int ich=0;ich<WMNumXYLayerCh;ich++){ 
      //vertical layer
      if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
	detdim->GetPosWM(mod,ipln,1,ich,0, &x, &y, &z);
      }
      posX = x*mm; 
      posY = y*mm; 
      posZ = z*mm;
      sprintf(name,"vlayer[%d][%d][%d]",mod,ipln,ich);
      if(check_WM){
        cout << name << " " << posX << " " << posY << " " << posZ << endl;
      }
      new G4PVPlacement(xrot2,G4ThreeVector(posX,posY,posZ),vwscint_intLV,name,targetLV,false,
			ich+ipln*NORMPLN+mod*NORMMOD,check_WM);
      //horizontal layer
      if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
	detdim->GetPosWM(mod,ipln,0,ich,0, &x, &y, &z);
      }
      posX = x*mm; 
      posY = y*mm; 
      posZ = z*mm;
      sprintf(name,"hlayer[%d][%d][%d]",mod,ipln,ich); 
      if(check_WM){
        cout << name << " " << posX << " " << posY << " " << posZ << endl;
      }
      new G4PVPlacement(yrot2,G4ThreeVector(posX,posY,posZ),hwscint_intLV,name,targetLV,false,
			ich+ipln*NORMPLN+mod*NORMMOD,check_WM);

      if(ich<20){
        //vertical grid layer
        if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
	  detdim->GetPosWM(mod,ipln,1,ich+40,1, &x, &y, &z);
	}
        posX = x*mm; 
        posY = y*mm; 
        posZ = z*mm;
        sprintf(name,"vlayer1[%d][%d][%d]",mod,ipln,ich+40);
        if(check_WM){
          cout << name << " " << posX << " " << posY << " " << posZ << endl;
        }
        new G4PVPlacement(rotgridX_v,G4ThreeVector(posX,posY,posZ),vwgridscint_intLV,name,targetLV,false,
			  ich+40+ipln*NORMPLN+mod*NORMMOD,check_WM);
        if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
	  detdim->GetPosWM(mod,ipln,1,ich+60,2, &x, &y, &z);
	}
        posX = x*mm; 
        posY = y*mm; 
        posZ = z*mm;
        sprintf(name,"vlayer2[%d][%d][%d]",mod,ipln,ich+60);
        if(check_WM){
          cout << name << " " << posX << " " << posY << " " << posZ << endl;
        }
        new G4PVPlacement(rotgridY_v,G4ThreeVector(posX,posY,posZ),vwgridscint_intLV,name,targetLV,false,
			  ich+60+ipln*NORMPLN+mod*NORMMOD,check_WM);
        //horizontal grid layer
        if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
	  detdim->GetPosWM(mod,ipln,0,ich+40,1, &x, &y, &z);
	}
        posX = x*mm; 
        posY = y*mm; 
        posZ = z*mm;
        sprintf(name,"hlayer1[%d][%d][%d]",mod,ipln,ich+40); 
        if(check_WM){
          cout << name << " " << posX << " " << posY << " " << posZ << endl;
        }
        new G4PVPlacement(rotgridX_h,G4ThreeVector(posX,posY,posZ),hwgridscint_intLV,name,targetLV,false,
			  ich+40+ipln*NORMPLN+mod*NORMMOD,check_WM);
        if(mod==MOD_ONAXIS_WM||mod==MOD_B2_WM){
	  detdim->GetPosWM(mod,ipln,0,ich+60,2, &x, &y, &z);
	}
        posX = x*mm; 
        posY = y*mm; 
        posZ = z*mm;
        sprintf(name,"hlayer2[%d][%d][%d]",mod,ipln,ich+60); 
        if(check_WM){
          cout << name << " " << posX << " " << posY << " " << posZ << endl;
        }
        new G4PVPlacement(rotgridY_h,G4ThreeVector(posX,posY,posZ),hwgridscint_intLV,name,targetLV,false,
			  ich+60+ipln*NORMPLN+mod*NORMMOD,check_WM);
      }
    }
  }


}

//___________________________________________________________________________________________________________
void DetectorConstruction::PlaceCHModule(G4LogicalVolume* motherLV, G4LogicalVolume* targetLV, int mod){

  this->PlaceWaterModule(motherLV,targetLV,mod);

  int ich;
  for(int ipln=0;ipln<WMNumPln;ipln++){
    for(int ilayer=1;ilayer<=2;ilayer++){
      for(int igridch_h=0;igridch_h<WMNumGridCh;igridch_h++){
        for(int igridch_v=0;igridch_v<WMNumGridCh;igridch_v++){
          if(igridch_h==0||igridch_v==0) continue; 

          ich = igridch_h + (ilayer+1)*WMNumGridCh;

          detdim->GetPosWM(MOD_B2_WM,ipln,0,ich,ilayer, &x, &y, &z);
          posY = y*mm - WMScintiSlitStep/2.;

          ich = igridch_v + (ilayer+1)*WMNumGridCh;
          detdim->GetPosWM(MOD_B2_WM,ipln,1,ich,ilayer, &x, &y, &z);
          posX = x*mm - WMScintiSlitStep/2.;
          posZ = z*mm;

          sprintf(name,"chtarget[%d][%d][%d][%d][%d]",mod,ipln,ilayer,igridch_h,igridch_v);
          new G4PVPlacement(0,G4ThreeVector(posX,posY,posZ),B2CHLV,name,targetLV,false,0,check_CHM);
        }
      }
    }
  }

}


//___________________________________________________________________________________________________________

void DetectorConstruction::DefineStructures()
{

  WorldSizeX = C_WorldSizeX*mm; 
  WorldSizeY = C_WorldSizeY*mm;
  WorldSizeZ = C_WorldSizeZ*mm;

  // =============================================
  // =========      ND Hall           ============
  // =============================================
  HallDirtPosX     = C_HallDirtPosX     *mm;
  HallDirtPosY     = C_HallDirtPosY     *mm;
  HallDirtPosZ     = C_HallDirtPosZ     *mm;
  PillarPosX       = C_PillarPosX       *mm;
  PillarPosY       = C_PillarPosY       *mm;
  PillarPosZ       = C_PillarPosZ       *mm;
  B2FloorPosY      = C_B2FloorPosY      *mm;
  HallDirtRadiusMin= C_HallDirtRadiusMin*mm; 
  HallDirtRadiusMax= C_HallDirtRadiusMax*mm;
  HallDirtHeight   = C_HallDirtHeight   *mm;
  HallDirtSPhi     = C_HallDirtSPhi     ; 
  HallDirtDPhi     = C_HallDirtDPhi     ;
  PillarX          = C_PillarX          *mm; 
  PillarY          = C_PillarY          *mm;
  PillarZ          = C_PillarZ          *mm;

  // =============================================
  // =========    Mother Volumes      ============
  // =============================================
  INGHMotherPosX = C_INGHMotherPosX         *mm;
  INGHMotherPosY = C_INGHMotherPosY         *mm;
  INGHMotherPosZ = C_INGHMotherPosZ         *mm;
  INGVMotherPosX = C_INGVMotherPosX         *mm;
  INGVMotherPosY = C_INGVMotherPosY         *mm;
  INGVMotherPosZ = C_INGVMotherPosZ         *mm;
  PMMotherPosX   = C_PMMotherPosX           *mm;
  PMMotherPosY   = C_PMMotherPosY           *mm;
  PMMotherPosZ   = C_PMMotherPosZ           *mm;
  B2MotherPosX   = C_B2MotherPosX           *mm;
  B2MotherPosY   = C_B2MotherPosY           *mm;
  B2MotherPosZ   = C_B2MotherPosZ           *mm;
  INGHMotherSizeX= C_INGHMotherSizeX        *mm; 
  INGHMotherSizeY= C_INGHMotherSizeY        *mm;
  INGHMotherSizeZ= C_INGHMotherSizeZ        *mm;
  INGVMotherSizeX= C_INGVMotherSizeX        *mm; 
  INGVMotherSizeY= C_INGVMotherSizeY        *mm;
  INGVMotherSizeZ= C_INGVMotherSizeZ        *mm;
  PMMotherSizeX  = C_PMMotherSizeX          *mm; 
  PMMotherSizeY  = C_PMMotherSizeY          *mm;
  PMMotherSizeZ  = C_PMMotherSizeZ          *mm;
  B2MotherSizeX  = C_B2MotherSizeX          *mm; 
  B2MotherSizeY  = C_B2MotherSizeY          *mm;
  B2MotherSizeZ  = C_B2MotherSizeZ          *mm;

  // ==========================================================================
  // ================== Modules Volume  =======================================
  // ==========================================================================
  INGSpace   = C_INGSpace                   *mm;
  INGStart   = C_INGStart                   *mm;
  B2WMPosX   = C_B2WMPosX                   *mm;
  B2WMPosY   = C_B2WMPosY                   *mm;
  B2WMPosZ   = C_B2WMPosZ                   *mm;
  B2CHPosX   = C_B2CHPosX                   *mm;
  B2CHPosY   = C_B2CHPosY                   *mm;
  B2CHPosZ   = C_B2CHPosZ                   *mm;
  B2INGPosX  = C_B2INGPosX                  *mm;
  B2INGPosY  = C_B2INGPosY                  *mm;
  B2INGPosZ  = C_B2INGPosZ                  *mm;
  INGSizeX   = C_INGSizeX                   *mm;
  INGSizeY   = C_INGSizeY                   *mm;
  INGSizeZ   = C_INGSizeZ                   *mm;
  PMSizeX    = C_PMSizeX                    *mm;
  PMSizeY    = C_PMSizeY                    *mm;
  PMSizeZ    = C_PMSizeZ                    *mm;
  WMSizeX    = C_WMSizeX                    *mm;
  WMSizeY    = C_WMSizeY                    *mm;
  WMSizeZ    = C_WMSizeZ                    *mm;

  // =============================================
  // =========          INGRID        ============
  // =============================================
  INGNumPln          = C_INGNumPln          ;
  INGNumCh           = C_INGNumCh           ;
  INGNumVetoCh       = C_INGNumVetoCh       ;  
  INGPlnDist         = C_INGPlnDist         *mm;
  INGPlnStart        = C_INGPlnStart        *mm;
  INGChStart         = C_INGChStart         *mm;
  INGIronStart       = C_INGIronStart       *mm;
  INGIronThick       = C_INGIronThick       *mm;
  INGIronXY          = C_INGIronXY          *mm;
  INGScintiLength    = C_INGScintiLength    *mm;
  INGScintiWidth     = C_INGScintiWidth     *mm;
  INGScintiThick     = C_INGScintiThick     *mm;
  INGScintiHoleDia_a = C_INGScintiHoleDia_a *mm;
  INGScintiHoleDia_b = C_INGScintiHoleDia_b *mm;
  INGScintiVertexX[0]   = C_INGScintiVertexX[0]   *mm;
  INGScintiVertexX[1]   = C_INGScintiVertexX[1]   *mm;
  INGScintiVertexX[2]   = C_INGScintiVertexX[2]   *mm;
  INGScintiVertexX[3]   = C_INGScintiVertexX[3]   *mm;
  INGScintiVertexX[4]   = C_INGScintiVertexX[4]   *mm;
  INGScintiVertexX[5]   = C_INGScintiVertexX[5]   *mm;
  INGScintiVertexX[6]   = C_INGScintiVertexX[6]   *mm;
  INGScintiVertexX[7]   = C_INGScintiVertexX[7]   *mm;
  INGScintiVertexY[0]   = C_INGScintiVertexY[0]   *mm;
  INGScintiVertexY[1]   = C_INGScintiVertexY[1]   *mm;
  INGScintiVertexY[2]   = C_INGScintiVertexY[2]   *mm;
  INGScintiVertexY[3]   = C_INGScintiVertexY[3]   *mm;
  INGScintiVertexY[4]   = C_INGScintiVertexY[4]   *mm;
  INGScintiVertexY[5]   = C_INGScintiVertexY[5]   *mm;
  INGScintiVertexY[6]   = C_INGScintiVertexY[6]   *mm;
  INGScintiVertexY[7]   = C_INGScintiVertexY[7]   *mm;
  INGLVetoLength     = C_INGLVetoLength     *mm;
  INGLVetoWidth      = C_INGLVetoWidth      *mm;
  INGLVetoThick      = C_INGLVetoThick      *mm;
  INGVetoStartZ      = C_INGVetoStartZ      *mm;
  INGSVetoLength     = C_INGSVetoLength     *mm;
  INGSVetoWidth      = C_INGSVetoWidth      *mm;
  INGSVetoThick      = C_INGSVetoThick      *mm;
  INGVetoPos1X       = C_INGVetoPos1X       *mm;
  INGVetoPos1Y       = C_INGVetoPos1Y       *mm;
  INGVetoPos2X       = C_INGVetoPos2X       *mm;
  INGVetoPos2Y       = C_INGVetoPos2Y       *mm;
  INGVetoPos3X       = C_INGVetoPos3X       *mm;
  INGVetoPos3Y       = C_INGVetoPos3Y       *mm;
  INGVetoPos4X       = C_INGVetoPos4X       *mm;
  INGVetoPos4Y       = C_INGVetoPos4Y       *mm;


  // =============================================
  // =========      Proton Module     ============
  // =============================================
  PMNumPln           = C_PMNumPln           ;
  PMNumVetoPln       = C_PMNumVetoPln       ;
  PMNumVetoCh        = C_PMNumVetoCh        ;
  PMNumCh            = C_PMNumCh            ;
  PMNumCh_mid        = C_PMNumCh_mid        ;  
  PMNumCh_side       = C_PMNumCh_side       ;  
  PMNumCh1           = C_PMNumCh1           ; 
  PMPlnDist          = C_PMPlnDist          *mm; 
  PMPlnDist_First    = C_PMPlnDist_First    *mm; 
  PMPlnStart         = C_PMPlnStart         *mm; 
  PMChStart          = C_PMChStart          *mm; 
  PMScintiLength     = C_PMScintiLength     *mm;
  PMScintiWidth      = C_PMScintiWidth      *mm;
  PMScintiThick      = C_PMScintiThick      *mm;
  PMScintiHoleDia    = C_PMScintiHoleDia    *mm;
  PMScintiVertexX[0] = C_PMScintiVertexX[0] *mm;
  PMScintiVertexX[1] = C_PMScintiVertexX[1] *mm;
  PMScintiVertexX[2] = C_PMScintiVertexX[2] *mm;
  PMScintiVertexX[3] = C_PMScintiVertexX[3] *mm;
  PMScintiVertexX[4] = C_PMScintiVertexX[4] *mm;
  PMScintiVertexX[5] = C_PMScintiVertexX[5] *mm;
  PMScintiVertexX[6] = C_PMScintiVertexX[6] *mm;
  PMScintiVertexX[7] = C_PMScintiVertexX[7] *mm;
  PMScintiVertexY[0] = C_PMScintiVertexY[0] *mm;
  PMScintiVertexY[1] = C_PMScintiVertexY[1] *mm;
  PMScintiVertexY[2] = C_PMScintiVertexY[2] *mm;
  PMScintiVertexY[3] = C_PMScintiVertexY[3] *mm;
  PMScintiVertexY[4] = C_PMScintiVertexY[4] *mm;
  PMScintiVertexY[5] = C_PMScintiVertexY[5] *mm;
  PMScintiVertexY[6] = C_PMScintiVertexY[6] *mm;
  PMScintiVertexY[7] = C_PMScintiVertexY[7] *mm;
  PMLVetoLength      = C_PMLVetoLength      *mm; 
  PMLVetoWidth       = C_PMLVetoWidth       *mm; 
  PMLVetoThick       = C_PMLVetoThick       *mm; 
  PMVetoStartZ       = C_PMVetoStartZ       *mm;
  PMSVetoLength      = C_PMSVetoLength      *mm;
  PMSVetoThick       = C_PMSVetoThick       *mm;
  PMSVetoWidth       = C_PMSVetoWidth       *mm;
  PMVetoPos1X        = C_PMVetoPos1X        *mm;
  PMVetoPos1Y        = C_PMVetoPos1Y        *mm;
  PMVetoPos2X        = C_PMVetoPos2X        *mm;
  PMVetoPos2Y        = C_PMVetoPos2Y        *mm;
  PMVetoPos3X        = C_PMVetoPos3X        *mm;
  PMVetoPos3Y        = C_PMVetoPos3Y        *mm;
  PMVetoPos4X        = C_PMVetoPos4X        *mm;
  PMVetoPos4Y        = C_PMVetoPos4Y        *mm;



  // =============================================
  // =========      WAGASCI           ============
  // =============================================
  WMNumPln            = C_WMNumPln           ;
  WMNumCh             = C_WMNumCh            ;
  WMNumLayer          = C_WMNumLayer         ;
  WMNumXYLayerCh      = C_WMNumXYLayerCh     ;
  WMNumGridCh         = C_WMNumGridCh        ;    
  WMPlnDist           = C_WMPlnDist          *mm; 
  WMLayerDist         = C_WMLayerDist        *mm; 
  WMPlnStart          = C_WMPlnStart         *mm; 
  WMChStart           = C_WMChStart          *mm; 
  WMGridChStart       = C_WMGridChStart      *mm; 
  WMWaterTargetSizeX  = C_WMWaterTargetSizeX *mm;
  WMWaterTargetSizeY  = C_WMWaterTargetSizeY *mm;
  WMWaterTargetSizeZ  = C_WMWaterTargetSizeZ *mm;
  WMCHCellTargetSizeX = C_WMCHCellTargetSizeX*mm;
  WMCHCellTargetSizeY = C_WMCHCellTargetSizeY*mm;
  WMCHCellTargetSizeZ = C_WMCHCellTargetSizeZ*mm;
  WMScintiLength      = C_WMScintiLength     *mm;
  WMScintiWidth       = C_WMScintiWidth      *mm;
  WMScintiThick       = C_WMScintiThick      *mm;
  WMScintiHoleLength  = C_WMScintiHoleLength *mm; 
  WMScintiHoleWidth   = C_WMScintiHoleWidth  *mm;
  WMScintiHoleThick   = C_WMScintiHoleThick  *mm;
  WMScintiHoleShift1  = C_WMScintiHoleShift1 *mm;
  WMScintiHoleShift2  = C_WMScintiHoleShitt2 *mm;
  WMScintiSlitLength  = C_WMScintiSlitLength *mm; 
  WMScintiSlitLength2 = C_WMScintiSlitLength2*mm; 
  WMScintiSlitWidth   = C_WMScintiSlitWidth  *mm;
  WMScintiSlitThick   = C_WMScintiSlitThick  *mm;
  WMScintiSlitStep    = C_WMScintiSlitStep   *mm;
  WMScintiSlitShift   = C_WMScintiSlitShift  *mm;
  WMScintiVertexX[0]  = C_WMScintiVertexX[0] *mm;
  WMScintiVertexX[1]  = C_WMScintiVertexX[1] *mm;
  WMScintiVertexX[2]  = C_WMScintiVertexX[2] *mm;
  WMScintiVertexX[3]  = C_WMScintiVertexX[3] *mm;
  WMScintiVertexX[4]  = C_WMScintiVertexX[4] *mm;
  WMScintiVertexX[5]  = C_WMScintiVertexX[5] *mm;
  WMScintiVertexX[6]  = C_WMScintiVertexX[6] *mm;
  WMScintiVertexX[7]  = C_WMScintiVertexX[7] *mm;
  WMScintiVertexX2[0] = C_WMScintiVertexX2[0]*mm;
  WMScintiVertexX2[1] = C_WMScintiVertexX2[1]*mm;
  WMScintiVertexX2[2] = C_WMScintiVertexX2[2]*mm;
  WMScintiVertexX2[3] = C_WMScintiVertexX2[3]*mm;
  WMScintiVertexX2[4] = C_WMScintiVertexX2[4]*mm;
  WMScintiVertexX2[5] = C_WMScintiVertexX2[5]*mm;
  WMScintiVertexX2[6] = C_WMScintiVertexX2[6]*mm;
  WMScintiVertexX2[7] = C_WMScintiVertexX2[7]*mm;
  WMScintiVertexY[0]  = C_WMScintiVertexY[0] *mm;
  WMScintiVertexY[1]  = C_WMScintiVertexY[1] *mm;
  WMScintiVertexY[2]  = C_WMScintiVertexY[2] *mm;
  WMScintiVertexY[3]  = C_WMScintiVertexY[3] *mm;
  WMScintiVertexY[4]  = C_WMScintiVertexY[4] *mm;
  WMScintiVertexY[5]  = C_WMScintiVertexY[5] *mm;
  WMScintiVertexY[6]  = C_WMScintiVertexY[6] *mm;
  WMScintiVertexY[7]  = C_WMScintiVertexY[7] *mm;
  WMAlTankSizeX       = C_WMAlTankSizeX      *mm;
  WMAlTankSizeY       = C_WMAlTankSizeY      *mm;
  WMAlTankSizeZ       = C_WMAlTankSizeZ      *mm;

  // =============================================
  // =========           BG           ============
  // =============================================
  //-----Ceiling
  B2CeilingPosX = C_B2CeilingPosX *mm;
  B2CeilingPosY = C_B2CeilingPosY *mm;
  B2CeilingPosZ = C_B2CeilingPosZ *mm;
  B2TransPosX = C_B2TransPosX *mm;
  B2TransPosY = C_B2TransPosY *mm;
  B2TransPosZ = C_B2TransPosZ *mm;
  B2CeilingSizeX = C_B2CeilingSizeX *mm;
  B2CeilingSizeY = C_B2CeilingSizeY *mm;
  B2CeilingSizeZ = C_B2CeilingSizeZ *mm;
  B2VboxSizeX = C_B2VboxSizeX *mm;
  B2VboxSizeY = C_B2VboxSizeY *mm;
  B2VboxSizeZ = C_B2VboxSizeZ *mm;
  B2CeilingRadiusMin = C_B2CeilingRadiusMin *mm;
  B2CeilingRadiusMax = C_B2CeilingRadiusMax *mm;
  //-----Floor
  B2FloorPosX = C_B2FloorPosX *mm;
  B2FloorPosY = C_B2FloorPosY *mm;
  B2FloorPosZ = C_B2FloorPosZ *mm;
  B2FloorSizeX = C_B2FloorSizeX *mm;
  B2FloorSizeY = C_B2FloorSizeY *mm;
  B2FloorSizeZ = C_B2FloorSizeZ *mm;
  B2FloorRadiusMin = C_B2FloorRadiusMin *mm;
  B2FloorRadiusMax = C_B2FloorRadiusMax *mm;
  //-----Right pillar
  B2Right_PillarPosX = C_B2Right_PillarPosX *mm;
  B2Right_PillarPosY = C_B2Right_PillarPosY *mm;
  B2Right_PillarPosZ = C_B2Right_PillarPosZ *mm;
  B2Right_PillarSizeX = C_B2Right_PillarSizeX *mm;
  B2Right_PillarSizeY = C_B2Right_PillarSizeY *mm;
  B2Right_PillarSizeZ = C_B2Right_PillarSizeZ *mm;
  //-----Left pillar
  B2Left_PillarPosX = C_B2Left_PillarPosX *mm;
  B2Left_PillarPosY = C_B2Left_PillarPosY *mm;
  B2Left_PillarPosZ = C_B2Left_PillarPosZ *mm;
  B2Left_PillarSizeX = C_B2Left_PillarSizeX *mm;
  B2Left_PillarSizeY = C_B2Left_PillarSizeY *mm;
  B2Left_PillarSizeZ = C_B2Left_PillarSizeZ *mm;
  //-----Small left pillar
  B2Small_PillarPosX = C_B2Small_PillarPosX *mm;
  B2Small_PillarPosY = C_B2Small_PillarPosY *mm;
  B2Small_PillarPosZ = C_B2Small_PillarPosZ *mm;
  B2Small_PillarSizeX = C_B2Small_PillarSizeX *mm;
  B2Small_PillarSizeY = C_B2Small_PillarSizeY *mm;
  B2Small_PillarSizeZ = C_B2Small_PillarSizeZ *mm;
}
