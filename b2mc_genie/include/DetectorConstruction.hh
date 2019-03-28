#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

#include "DetectorDimension.hh"
#include "HLayerSD.hh"
#include "VetoSD.hh"
#include "VLayerSD.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4EllipticalTube.hh"
#include "G4ExtrudedSolid.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4RotationMatrix.hh"

#include "Const.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    void DefineMaterial();
    void DefineStructures();

    G4Material *Air,*Fe,*Scinti,*Water,*Concrete,*Al,*SUS304;

    char name[100];
    double x,y,z;
    G4double posX,posY,posZ;
    int mode;
    double sINGshiftX;
    double sINGshiftY;
    double sINGshiftZ;
    double sINGdeg   ;

    void PlaceINGRID(G4LogicalVolume *motherLV,int ingmod, int mod);
    void PlaceProtonModule(G4LogicalVolume *motherLV,int mod);
    void PlaceWaterModule(G4LogicalVolume *motherLV,G4LogicalVolume* targetLV, int mod);
    void PlaceCHModule(G4LogicalVolume *motherLV, G4LogicalVolume* targetLV, int mod);



    G4RotationMatrix *zrot,*zrot2,*xrot,*xrot2,*yrot,*yrot2,*rotgridY_v,*rotgridY_h,*rotgridX_v,*rotgridX_h;


    G4Box             *experimentalHall_box; 
    G4LogicalVolume   *worldLV;
    G4VPhysicalVolume *worldPV;
    G4Tubs            *HallDirtSLD;           
    G4LogicalVolume   *HallDirtLV;   
    G4VPhysicalVolume *HallDirtPV; 
    G4VisAttributes   *HallDirtVisAtt; 

    G4Box             *horizontalHall_box;        
    G4LogicalVolume   *horizontalLV;    
    G4VPhysicalVolume *horizontalPV;  

    G4Box             *verticalHall_box;       
    G4LogicalVolume   *verticalLV;   
    G4VPhysicalVolume *verticalPV; 

    G4Box             *protonHall_box;       
    G4LogicalVolume   *ProtonLV;   
    G4VPhysicalVolume *ProtonPV; 

    G4Box             *B2Hall_box;           
    G4LogicalVolume   *B2HallLV;   
    G4VPhysicalVolume *B2HallPV;

    G4Box *INGmodule_box;  
    G4Box *PMmodule_box; 
    G4Box *WMmodule_box; 
    G4LogicalVolume *moduleLV[NUMMAX_MODULE]; 


    //For B2BG
    //-----Floor
    G4Tubs *B2ceiling_tube;
    G4Box *B2Vbox;
    G4VSolid *B2Ceiling_tube;
    G4LogicalVolume *B2CeilingLV;
    G4VPhysicalVolume *B2CeilingPV;
    //-----Floor
    G4Tubs *B2Floor_tube;
    G4LogicalVolume *B2FloorLV;
    G4VPhysicalVolume *B2FloorPV;
    //-----Right pillar
    G4Box *B2Right_Pillar_box;
    G4LogicalVolume *B2Right_PillarLV;
    G4VPhysicalVolume *B2Right_PillarPV;
    //-----Left pillar
    G4Box *B2Left_Pillar_box;
    G4LogicalVolume *B2Left_PillarLV;
    G4VPhysicalVolume *B2Left_PillarPV;
    //-----Small left pillar
    G4Box *B2Small_Pillar_box;
    G4LogicalVolume *B2Small_PillarLV;
    G4VPhysicalVolume *B2Small_PillarPV;


    G4Box           *iron_block;
    G4VisAttributes *ironVisAtt;
    G4Box           *Water_box;
    G4LogicalVolume *WaterLV;
    G4VisAttributes *WaterVisAtt;
    G4LogicalVolume *B2WaterLV;
    G4VisAttributes *B2WaterVisAtt;
    G4Box           *Air_box; 
    G4LogicalVolume *B2AirLV;
    G4Box           *CH_box;
    G4VisAttributes *B2CHVisAtt; 
    G4Box           *Al_box;             
    G4LogicalVolume *AlLV;      
    G4VisAttributes *AlVisAtt;  


    G4ExtrudedSolid* vscint_tmp        ; 
    G4EllipticalTube* vsci_hole        ; 
    G4SubtractionSolid* vscint_int     ; 
    G4VisAttributes* vscint_intVisAtt  ; 
    G4ExtrudedSolid* hscint_tmp        ; 
    G4EllipticalTube* hsci_hole        ; 
    G4SubtractionSolid* hscint_int     ; 
    G4VisAttributes* hscint_intVisAtt  ; 
    G4ExtrudedSolid* vscint2_tmp       ;
    G4EllipticalTube* vsci_hole2       ;
    G4SubtractionSolid* vscint2_int    ;
    G4VisAttributes* vscint2_intVisAtt ;
    G4ExtrudedSolid* hscint2_tmp       ;
    G4EllipticalTube* hsci_hole2       ;
    G4SubtractionSolid* hscint2_int    ;
    G4VisAttributes* hscint2_intVisAtt ;
    G4ExtrudedSolid* vscint3_tmp       ;
    G4EllipticalTube* vsci_hole3       ;
    G4SubtractionSolid* vscint3_int    ;
    G4VisAttributes* vscint3_intVisAtt ;
    G4ExtrudedSolid* hscint3_tmp       ;
    G4EllipticalTube* hsci_hole3       ;
    G4SubtractionSolid* hscint3_int    ;
    G4VisAttributes* hscint3_intVisAtt ;
    G4ExtrudedSolid* vwscint_tmp           ; 
    G4Box* wmhole_v                        ; 
    G4SubtractionSolid* vwscint_int        ; 
    G4VisAttributes* vwscint_intVisAtt     ; 
    G4ExtrudedSolid* hwscint_tmp           ; 
    G4Box* wmhole_h                        ; 
    G4SubtractionSolid* hwscint_int        ; 
    G4VisAttributes* hwscint_intVisAtt     ; 
    G4VSolid* vwgridscint_tmp              ; 
    G4VSolid* vwgridscint_int              ; 
    G4VSolid* wmslit_v                     ; 
    G4VSolid* vwgridscint_tmp2             ; 
    G4VisAttributes* vwgridscint_intVisAtt ; 
    G4VSolid* hwgridscint_tmp              ; 
    G4VSolid* hwgridscint_int              ; 
    G4VSolid* wmslit_h                     ; 
    G4VSolid* hwgridscint_tmp2             ; 
    G4VisAttributes* hwgridscint_intVisAtt ; 
    G4Box* Lveto_box         ;
    G4Box* Sveto_box         ;
    G4Box* PLveto_box         ; 
    G4Box* PSveto_box         ; 

    //====================================================
    //====================================================
    //====================================================
    //====================================================
    G4LogicalVolume *ironLV;
    G4LogicalVolume *B2CHLV; 
    G4LogicalVolume* vscint_intLV      ; 
    G4LogicalVolume* hscint_intLV      ; 
    G4LogicalVolume* vscint2_intLV     ;
    G4LogicalVolume* hscint2_intLV     ;
    G4LogicalVolume* vscint3_intLV     ;
    G4LogicalVolume* hscint3_intLV     ;
    G4LogicalVolume* vwscint_intLV         ; 
    G4LogicalVolume* hwscint_intLV         ; 
    G4LogicalVolume* vwgridscint_intLV     ; 
    G4LogicalVolume* hwgridscint_intLV     ; 
    G4LogicalVolume *LvetoLV ;
    G4LogicalVolume *SvetoLV ;
    G4LogicalVolume *PLvetoLV ; 
    G4LogicalVolume *PSvetoLV ; 
    //====================================================
    //====================================================
    //====================================================


    VetoSD*   avetoSD;
    HLayerSD* ahlayerSD;  
    VLayerSD* avlayerSD; 


    //DetectorConstruction(int,double*);
    DetectorConstruction(int);
    DetectorConstruction(int,DetectorDimension*);
    ~DetectorConstruction();

    G4VPhysicalVolume* Construct();
    G4double WorldSizeX, WorldSizeY, WorldSizeZ;
    G4double INGHMotherSizeX, INGHMotherSizeY, INGHMotherSizeZ, INGVMotherSizeX, INGVMotherSizeY, INGVMotherSizeZ;
    G4double PMMotherSizeX, PMMotherSizeY, PMMotherSizeZ;
    G4double B2MotherSizeX, B2MotherSizeY, B2MotherSizeZ;
    G4double INGSizeX, INGSizeY, INGSizeZ;
    G4double PMSizeX, PMSizeY, PMSizeZ;
    G4double WMSizeX, WMSizeY, WMSizeZ;
    G4double WMWaterTargetSizeX, WMWaterTargetSizeY, WMWaterTargetSizeZ;
    G4double WMCHCellTargetSizeX, WMCHCellTargetSizeY, WMCHCellTargetSizeZ;
    G4double WMAlTankSizeX, WMAlTankSizeY, WMAlTankSizeZ;
    G4double HallDirtPosX, HallDirtPosY, HallDirtPosZ;
    G4double PillarPosX, PillarPosY, PillarPosZ;
    G4double B2floorPosY;
    G4double HallDirtRadiusMin, HallDirtRadiusMax, HallDirtHeight, HallDirtSPhi, HallDirtDPhi;
    G4double PillarX, PillarY, PillarZ;
    G4double INGHMotherPosX, INGHMotherPosY, INGHMotherPosZ, INGVMotherPosX, INGVMotherPosY, INGVMotherPosZ;
    G4double PMMotherPosX, PMMotherPosY, PMMotherPosZ;
    G4double B2MotherPosX, B2MotherPosY, B2MotherPosZ;
    G4double INGSpace, INGStart;
    G4double B2WMPosX, B2WMPosY, B2WMPosZ, B2CHPosX, B2CHPosY, B2CHPosZ, B2INGPosX, B2INGPosY, B2INGPosZ;
    G4double INGPlnDist, INGPlnStart, INGChStart, INGIronStart;
    G4double INGIronThick, INGIronXY, INGTotMassIron;
    G4double INGScintiLength, INGScintiWidth, INGScintiThick, INGScintiHoleDia_a, INGScintiHoleDia_b, INGTotMassScinti;
    G4double INGScintiVertexX[8], INGScintiVertexY[8];
    G4double INGLVetoLength, INGLVetoWidth, INGLVetoThick, INGVetoStartZ, INGSVetoLength, INGSVetoWidth, INGSVetoThick;
    G4double INGVetoPos1X, INGVetoPos1Y, INGVetoPos2X, INGVetoPos2Y, INGVetoPos3X, INGVetoPos3Y, INGVetoPos4X, INGVetoPos4Y;
    G4double PMPlnDist, PMPlnDist_First, PMPlnStart, PMChStart;
    G4double PMScintiLength, PMScintiWidth, PMScintiThick, PMScintiHoleDia, PMTotMassScinti, PMTotMassVetoSci;
    G4double PMScintiVertexX[8], PMScintiVertexY[8];
    G4double PMLVetoLength, PMLVetoWidth, PMLVetoThick, PMVetoStartZ, PMSVetoLength, PMSVetoWidth, PMSVetoThick;
    G4double PMVetoPos1X, PMVetoPos1Y, PMVetoPos2X, PMVetoPos2Y, PMVetoPos3X, PMVetoPos3Y, PMVetoPos4X, PMVetoPos4Y;
    G4double WMPlnDist, WMLayerDist, WMPlnStart, WMChStart, WMGridChStart;
    G4double WMScintiLength, WMScintiWidth, WMScintiThick, WMScintiHoleLength, WMScintiHoleWidth, WMScintiHoleThick;
    G4double WMScintiHoleShift1, WMScintiHoleShift2, WMScintiSlitLength, WMScintiSlitLength2, WMScintiSlitWidth, WMScintiSlitThick, WMScintiSlitStep, WMScintiSlitShift;
    G4double WMScintiVertexX[8], WMScintiVertexY[8];
    G4double WMScintiVertexX2[8];
    G4int INGNumPln, INGNumCh, INGNumVetoCh;
    G4int PMNumPln, PMNumVetoPln, PMNumVetoCh, PMNumCh, PMNumCh_mid, PMNumCh_side, PMNumCh1;
    G4int WMNumPln, WMNumCh, WMNumLayer, WMNumXYLayerCh, WMNumGridCh;

    //B2BG
    //-----Ceiling
    G4double B2CeilingPosX, B2CeilingPosY, B2CeilingPosZ;
    G4double B2TransPosX, B2TransPosY, B2TransPosZ;
    G4double B2CeilingSizeX, B2CeilingSizeY, B2CeilingSizeZ;
    G4double B2VboxSizeX, B2VboxSizeY, B2VboxSizeZ;
    G4double B2CeilingRadiusMax, B2CeilingRadiusMin;
    //-----Floor
    G4double B2FloorPosX, B2FloorPosY, B2FloorPosZ;
    G4double B2FloorSizeX, B2FloorSizeY, B2FloorSizeZ;
    G4double B2FloorRadiusMax, B2FloorRadiusMin;
    //-----Right pillar
    G4double B2Right_PillarPosX, B2Right_PillarPosY, B2Right_PillarPosZ;
    G4double B2Right_PillarSizeX, B2Right_PillarSizeY, B2Right_PillarSizeZ;
    //-----Left pillar
    G4double B2Left_PillarPosX, B2Left_PillarPosY, B2Left_PillarPosZ;
    G4double B2Left_PillarSizeX, B2Left_PillarSizeY, B2Left_PillarSizeZ;
    //-----Small left pillar
    G4double B2Small_PillarPosX, B2Small_PillarPosY, B2Small_PillarPosZ;
    G4double B2Small_PillarSizeX, B2Small_PillarSizeY, B2Small_PillarSizeZ;
  private:
    DetectorDimension *detdim;

};

#endif
