#ifndef ___CONSTH
#define ___CONSTH 1

//===========================================================================
//===========================================================================
//===========================================================================
//==========                  Scales of detectors                  ==========
//==========                   orrection Factors                   ==========
//==========      Length below are basically expressed in "mm"     ==========
//==========             Excemptions must be expressed.            ==========
//===========================================================================
//===========================================================================
//===========================================================================

//-------------------- Version information --------------------
// Ver.0: N. Chikuma (2016/06/26)
//         1) Set variables for WAGASCI
// Ver.1: K. Kin (2017/10/01 - 2017/10/19)
//         1) Add B2BG information
//         2) Adjust variables for Detector Response
//-------------------- Version information --------------------

#define MY_PATH "/home/t2k/nchikuma/b2_software"

//#define SUPPRESS_FOR_EXT_BG
//#define MASK_GRID_CHANNEL

//Run mode (option -m)
#define PROTON 2 //Proton Module @On-axis
#define HORIZONTAL 3
#define VERTICAL 4
#define INGWaterModule 5
#define INGWaterModuleBG 6
#define B2Water 7
#define B2CH 8 //Proton Module @B2
#define B2ING 9
#define B2Wall 10
#define B2Ceiling 11
#define B2Floor 12
#define B2Right_Pillar 13
#define B2Left_Pillar 14

//Module ID
#define MOD_INGRID_C 3
#define MOD_INGRID_H 6  //0-6
#define MOD_INGRID_V 13 //7-13
#define MOD_B2_INGRID 14
#define MOD_ONAXIS_WM 15
#define MOD_PM 16    //Proton Module Dectector Number, and Logical Volume Number for On-axis
#define MOD_B2_WM 21
#define MOD_B2_CH 22 //Logical Volume Number for B2 Proton Module/ CH WAGASCI

//ID for Top/Side view
#define SideView 0
#define TopView 1

//# of modules
#define NUMMAX_MODULE 30
//#define NUMINGMOD 14
#define NUMINGMOD 15 //Newly define of 14: B2-INGRID
#define NUMB2MOD 4

//Normalization of channle ID
#define NORMMOD 2000
#define NORMPLN 100

//Threshold of PE
#define PE_THRESHOLD 2.5
#define VETO_PE_THRESHOLD 2.5 //3.5

const double PI = 3.14159265358979323846;

//Set mod# put in B2MotherLV
// 0-6: H-INGRID
// 7-13: V-INGRID
// 14: B2-INGRID
// 16: Proton Module @B2
// 20: INGRIDWaterModule (On-axis)
// 21: WAGASCI water module

const int NotPlaceMod[4] = {24, 25, 26, 30};

//Detector Overlap checking
const bool check_MotherLV = true;
const bool check_moduleLV = true;
const bool check_ING      = false;
const bool check_PM       = false;
const bool check_WM       = false;
const bool check_CHM      = false;

const double C_WorldSizeX = 40000./2.; //mm //half of each side
const double C_WorldSizeY = 20000./2.; //mm
const double C_WorldSizeZ = 40000./2.; //mm

//=========================================
//==========       ND Hall       ==========
//=========================================
//------- Positions in World -------
const double C_HallDirtPosX      = -3222.; //mm
const double C_HallDirtPosY      =     0.; //mm
const double C_HallDirtPosZ      =  1700.; //mm
const double C_PillarPosX        = -8235.; //mm
const double C_PillarPosY        = -4443.; //mm
const double C_PillarPosZ        =  2631.; //mm
//----- Sizes ------
//Dimension of Hall dirt volume
const double C_HallDirtRadiusMin =  8950.; //mm
const double C_HallDirtRadiusMax = (8950. + 6000.); //mm
const double C_HallDirtHeight    = 10000.; //mm
const double C_HallDirtSPhi      = 0.;
const double C_HallDirtDPhi      = 2.*PI;
//Pillar
const double C_PillarX           = 1000./2.; //mm //half of each side
const double C_PillarY           = 3000./2.; //mm
const double C_PillarZ           = 4000./2.; //mm


//===========================================
//==========     Mother Volume     ==========
//===========================================
//--- B2 Floor ----
const double C_B2floorPosY =  -5943.; //mm
const double C_B2FloorPosX =  -3222.; //mm
const double C_B2FloorPosY = -7971.5; //mm
const double C_B2FloorPosZ =   1700.; //mm

//-------- Positions in World -------
//INGRID 7 horizontal modules
const double C_INGHMotherPosX =     0.; //mm
const double C_INGHMotherPosY =     0.; //mm
const double C_INGHMotherPosZ =     0.; //mm
//INGRID 7 vertical modules
const double C_INGVMotherPosX =     0.; //mm
const double C_INGVMotherPosY =     0.; //mm
const double C_INGVMotherPosZ = -4000.; //mm
//Proton Module & INGRIDWaterModule
const double C_PMMotherPosX   =     0.; //mm
const double C_PMMotherPosY   =     0.; //mm
const double C_PMMotherPosZ   = -1200.; //mm
//B2 modules 
const double C_B2MotherPosX   = -5735.; //mm
const double C_B2MotherPosY   = 1000.+C_B2floorPosY; //mm
const double C_B2MotherPosZ   =  2633.; //mm

//------- Sizes of box -------
//INGRID 7 horizontal modules 
const double C_INGHMotherSizeX = 10500./2.; //mm //half of each side
const double C_INGHMotherSizeY =  2000./2.; //mm
const double C_INGHMotherSizeZ =  1360./2.; //mm
//INGRID 7 vertical modules
const double C_INGVMotherSizeX =  2000./2.; //mm //half of each side
const double C_INGVMotherSizeY = 10500./2.; //mm
const double C_INGVMotherSizeZ =  1360./2.; //mm
//Proton Module & INGRIDWaterModule
const double C_PMMotherSizeX   =  2000./2.; //mm //half of each side
const double C_PMMotherSizeY   =  2000./2.; //mm
const double C_PMMotherSizeZ   =   960./2.; //mm
// B2 modules 
const double C_B2MotherSizeX   =  4000./2.; //mm //half of each side
const double C_B2MotherSizeY   =  2000./2.; //mm
const double C_B2MotherSizeZ   =  5000./2.; // 4000./2.; //mm


//===========================================
//==========     Module Volume     ==========
//===========================================
//-------- Positions -------
//Between each INGRID module in "INGH(V)Mother"
const double C_INGSpace  =  1500.; //mm | Between center of two adjacent modules
const double C_INGStart  = -4500.; //mm | Center of the first module
//Each sub-module in "B2Mother" 
const double C_B2WMPosX  =      0.; //mm
const double C_B2WMPosY  =    216.; //mm
const double C_B2WMPosZ  =   -998.; //-925.; //mm
const double C_B2CHPosX  =      0.; //mm
const double C_B2CHPosY  =    220.; //mm
const double C_B2CHPosZ  =  -1750.; //-1677.; //mm
const double C_B2INGPosX =      3.; //  0.; //mm
const double C_B2INGPosY =    117.; //113.; //mm
const double C_B2INGPosZ =    -50.; //mm
//B2 Backgrounds
//- Ceiling -
const double C_B2CeilingPosX = -3222.; //mm
const double C_B2CeilingPosY = -1573.; //mm
const double C_B2CeilingPosZ =  1700.; //mm
const double C_B2TransPosX   =  2500.; //mm
const double C_B2TransPosY   =  7000.; //mm
const double C_B2TransPosZ   =     0.; //mm
//- Right pillar -
const double C_B2Right_PillarPosX = -8235.; //mm
const double C_B2Right_PillarPosY = 3870./2.+C_B2floorPosY; //mm
const double C_B2Right_PillarPosZ = 2633.; //mm
//- Left pillar -
const double C_B2Left_PillarPosX = -3235.; //mm
const double C_B2Left_PillarPosY = 3870./2.+C_B2floorPosY; //mm
const double C_B2Left_PillarPosZ = 2633.; //mm
//- Small left pillar -
const double C_B2Small_PillarPosX = -3235.; //mm
const double C_B2Small_PillarPosY = 3870./2.+C_B2floorPosY; //mm
const double C_B2Small_PillarPosZ = -2000.; //mm

//------- Sizes of box -------
//INGRID Module
const double C_INGSizeX = 1500./2.; //mm
const double C_INGSizeY = 1500./2.; //mm
const double C_INGSizeZ = 1300./2.; //mm
//Proton Module (On-axis/B2)
const double C_PMSizeX = 1500./2.; //mm
const double C_PMSizeY = 1500./2.; //mm
const double C_PMSizeZ =  940./2.; //mm
//Water Module (On-axis/B2)
const double C_WMSizeX = 1500./2.; //mm
const double C_WMSizeY = 1500./2.; //mm
const double C_WMSizeZ =  550./2.; //600./2.; //mm
// B2 Backgrounds
//- Ceiling -
const double C_B2CeilingRadiusMin =    0.; //mm
const double C_B2CeilingRadiusMax = 8499.; //mm
const double C_B2CeilingSizeX     =  500.; //mm
const double C_B2CeilingSizeY     =    0.; //mm
const double C_B2CeilingSizeZ     = 2.*PI; //mm
const double C_B2VboxSizeX        = 2000.; //mm
const double C_B2VboxSizeY        = 2200.; //mm
const double C_B2VboxSizeZ        = 1250.; //mm
//-----Floor
const double C_B2FloorRadiusMin =     0.; //mm
const double C_B2FloorRadiusMax =  8499.; //mm
const double C_B2FloorSizeX     = 2028.5; //mm
const double C_B2FloorSizeY     =     0.; //mm
const double C_B2FloorSizeZ     =  2.*PI; //mm
//-----Right pillar
const double C_B2Right_PillarSizeX =     500.; //mm
const double C_B2Right_PillarSizeY = 3870./2.; //mm
const double C_B2Right_PillarSizeZ =    2000.; //mm
//-----Left pillar
const double C_B2Left_PillarSizeX =     500.; //mm
const double C_B2Left_PillarSizeY = 3870./2.; //mm
const double C_B2Left_PillarSizeZ =    2000.; //mm
//-----Small left pillar
const double C_B2Small_PillarSizeX =     500.; //mm
const double C_B2Small_PillarSizeY = 3870./2.; //mm
const double C_B2Small_PillarSizeZ =     500.; //mm


//====================================
//==========     INGRID     ==========
//====================================
const int    C_INGNumPln    =    11.;
const int    C_INGNumCh     =    24.;
const int    C_INGNumVetoCh =    22.;
const double C_INGPlnDist   =   107.; //mm | Between two adjacent planes
const double C_INGPlnStart  =  -540.; //mm | Center of 1st vertical scinti plane(pln#0)
const double C_INGChStart   =  -575.; //mm | Center of ch#0
const double C_INGIronStart = -481.5; //mm | Center of 1st iron
//------- Detector dimenstion -------
//Iron
const double C_INGIronThick   =   65.; //mm
const double C_INGIronXY      = 1240.; //mm
const double C_INGTotMassIron = 7.0605; //per plane ////99.54; //ton
//Scintillators
const double C_INGScintiLength    = 1200.; //mm
const double C_INGScintiWidth     =   50.; //mm
const double C_INGScintiThick     =   10.; //mm
const double C_INGScintiHoleDia_a =  2.57; //1.3 ; //mm
const double C_INGScintiHoleDia_b =  3.88; //1.95; //mm
const double C_INGTotMassScinti   =  0.3174; //per plane 4.01421; // 3.74; //ton
//8 Vertices of Scintillator Cross Section for Extruded Solid
const double C_INGScintiVertexX[8] = {-23.616, -24.389, -24.389, -23.616, 
				       23.616,  24.389,  24.389,  23.616}; //mm
const double C_INGScintiVertexY[8] = {-4.71, -0.5,  0.5,  4.71, 
				       4.71,  0.5, -0.5, -4.71 }; //mm
//Long veto planes
const double C_INGLVetoLength = 1300.; //mm
const double C_INGLVetoWidth  =   50.; //mm
const double C_INGLVetoThick  =   10.; //mm
const double C_INGVetoStartZ  = -525.; //mm | Center of 1st veto pln
//Short veto planes
const double C_INGSVetoLength = 1120.; //mm 
const double C_INGSVetoWidth  =   50.; //mm 
const double C_INGSVetoThick  =   10.; //mm 
//Veto Positions for INGRID
const double C_INGVetoPos1X =     59.; //mm | Top
const double C_INGVetoPos1Y =    684.; //mm
const double C_INGVetoPos2X =     -9.; //mm | Bottom
const double C_INGVetoPos2Y =   -659.; //mm
const double C_INGVetoPos3X =    709.; //mm | Left
const double C_INGVetoPos3Y =      3.; //mm
const double C_INGVetoPos4X = -705.75; //mm | Right
const double C_INGVetoPos4Y =    -37.; //mm


//===========================================
//==========     Proton Module     ==========
//===========================================
const int    C_PMNumPln        =    18.;
const int    C_PMNumVetoPln    =     4.;
const int    C_PMNumVetoCh     =    17.;
const int    C_PMNumCh         =    32.;
const int    C_PMNumCh_mid     =    16.; //Num of scibar type scinti
const int    C_PMNumCh_side    =     8.; //Num of INGRID type scinti for each side
const int    C_PMNumCh1        =    24.; //Num of ch for 1st pln
const double C_PMPlnDist       =    23.; //mm | Between two adjacent H pln & V pln
const double C_PMPlnDist_First =    27.; //mm | Between 1st V pln & 2nd H pln
const double C_PMPlnStart      = -404.5; //mm | Center of 1st horizontal plane(pln#0)
const double C_PMChStart       =  -575.; //mm | Center of ch#0
//------- Detector dimenstion -------
//Scintillators
const double C_PMScintiLength   =    1200.; //mm
const double C_PMScintiWidth    =      25.; //mm
const double C_PMScintiThick    =      13.; //mm
const double C_PMScintiHoleDia  =     1.8 ; // 0.9; //mm;
const double C_PMTotMassScinti  =  0.569205; //0.56904; //ton
const double C_PMTotMassVetoSci = 0.0288528; //0.028848; //ton
//8 Vertices of Scintillator Cross section for Extruded Solid
const double C_PMScintiVertexX[8] = {-11.672, -12.21, -12.21, -11.672, 
				      11.672,  12.21,  12.21,  11.672}; //mm
const double C_PMScintiVertexY[8] = {-6.17, -3.5,  3.5,  6.17, 
				      6.17,  3.5, -3.5, -6.17}; //mm
//Long veto planes
const double C_PMLVetoLength = 1250.; //mm 
const double C_PMLVetoWidth  =   50.; //mm 
const double C_PMLVetoThick  =   10.; //mm 
const double C_PMVetoStartZ  = -400.; //mm | Center of 1st veto pln
//Short veto planes
const double C_PMSVetoLength = 1200.; //mm 
const double C_PMSVetoWidth  =   50.; //mm 
const double C_PMSVetoThick  =   10.; //mm 
//Veto Positions for ProtomModule
const double C_PMVetoPos1X =   -5.; //mm | Top
const double C_PMVetoPos1Y =  655.; //mm
const double C_PMVetoPos2X =   15.; //mm | Bottom
const double C_PMVetoPos2Y = -655.; //mm
const double C_PMVetoPos3X =  655.; //mm | Left
const double C_PMVetoPos3Y =  -25.; //mm
const double C_PMVetoPos4X = -655.; //mm | Right
const double C_PMVetoPos4Y =    5.; //mm


//=====================================
//==========     WAGASCI     ==========
//=====================================
const int    C_WMNumPln       =      8.;
const int    C_WMNumCh        =     80.;
const int    C_WMNumLayer     =     16.;
const int    C_WMNumXYLayerCh =     40.;
const int    C_WMNumGridCh    =     20.; //#of ch of grid for each X or Y layer
const double C_WMPlnDist      =     57.; //mm | Between two adjacent X-layers (planes)
const double C_WMLayerDist    =    28.5; //mm | Between two adjacent X and Y-layers (layers)
const double C_WMPlnStart     =  -236.5; //-226.5; //mm | X center of 1st X-layer(pln#0 side view)
const double C_WMChStart      =  -512.5; //-537.5; //mm | Y center of ch#0
const double C_WMGridChStart  =  -474.9; //-524.9; //mm | Y center of ch#40 & ch#60
//------- Detector dimenstion -------
//Water target
const double C_WMWaterTargetSizeX = 1252.; //1256.; //mm
const double C_WMWaterTargetSizeY = 1180.; //1256.; //mm
const double C_WMWaterTargetSizeZ =  502.; // 482.; //466.; //mm
//CH target
const double C_WMCHCellTargetSizeX = 4.6; //mm
const double C_WMCHCellTargetSizeY = 4.6; //mm
const double C_WMCHCellTargetSizeZ = 2.5; //mm





//------- Detector components -------
//Scintillators
const double C_WMScintiLength = 1000.; //mm
const double C_WMScintiWidth  =   25.; //mm
const double C_WMScintiThick  =    3.; //mm
//Groove for fiber
const double C_WMScintiHoleLength = 1000.1; //mm 
const double C_WMScintiHoleWidth  =   1.21; //mm
const double C_WMScintiHoleThick  =   1.21; //mm
const double C_WMScintiHoleShift1 =    3.9; //mm
const double C_WMScintiHoleShitt2 =    0.9; //mm
//Slit for 3D-grid strudture
const double C_WMScintiSlitLength  =  6.5; //4.66; //mm, for WAGASCI
const double C_WMScintiSlitLength2 =  6.5; //3.50; //mm, for INGRID WM
const double C_WMScintiSlitWidth   = 13.1; //mm
const double C_WMScintiSlitThick   =  3.1; //mm
const double C_WMScintiSlitStep    =  50.; //mm
const double C_WMScintiSlitShift   =  5.5; //mm
//Margin from scintillator edge to bundle (MPPC)
const double C_WMFiberMargin = 200.; //mm
//8 Vertices of Scintillator Cross section for Extruded Solid
const double C_WMScintiVertexX[8] = {-10.4, -11.8, -11.8, -10.4, 
				      10.9,  11.8,  11.8,  10.9}; //mm
const double C_WMScintiVertexX2[8] = {-10.4, -11.8, -11.8, -10.4, 
				      10.9,  11.8,  11.8,  10.9}; //mm //For Grid Scintillator
const double C_WMScintiVertexY[8] = {-1.4, -0.75,  0.75,  1.4, 
				      1.4,  1.00, -1.00, -1.4};  //mm




//Water tank dimension (Aluminum)
const double C_WMAlTankSizeX = 1276.; //1456.; //mm //Outer scale
const double C_WMAlTankSizeY = 1204.; //1456.; //mm
const double C_WMAlTankSizeZ =  510.; // 606.; //mm

//==============================================
//==========     Detecor Response     ==========
//==============================================
//Energy deposit correction factors
const double C_Corr_Birks = 0.0208; //Quoted from SciBooNE MC
//Scintillator & Fiber related factors
const double C_ScintiAttLeng    =  10.46; //cm
const double C_WMScintiAttLeng  =    4.0; //cm
const double C_FiberAttLeng     =  241.7; //cm
const double C_TransTimeInFiber = 1./28.; //Calculation: 1cm / 2.8e10[cm/s] * 1e9 [ns]
const double C_INGRIDFactor     =     1.;//1.08; //P.E. factor for INGRID scintillator
const double C_ScibarFactor     =     1.;//1.77; //P.E. factor for SciBar scintillator
const double C_WMFactor         =     1.; //P.E. factor for WAGASCI scintillator
//MPPC
const double C_MeV2PE           =   45.9; //v3r4
const double C_WMMeV2PE         =   28.0; //31.5; //v3r4
const double C_MPPCPixel        =   667.; //v3
const double C_Eff_PDE          = -0.275; //@deltaV = 1.1 V
const double C_WMEff_PDE        =    -1.; //@deltaV = 1.1 V
const double C_ETC_CONST        =    5.0; //v3
const double C_rPDE             =    1.7; //@deltaV = 1.1 V
const double C_CrossAfterRate   =   0.18;  //0.09; //@deltaV = 1.1 V
const double C_WMCrossAfterRate =   0.027; //0.18; //@deltaV = 1.1 V, for IWM
//const double C_WMCrossAfterRate =   0.029; //@deltaV = 1.1 V
const double C_PixelGainVari    =  0.031; //Gain variation among pixels
const double C_WMPixelGainVari  =  0.040; //0.133; //For WAGASCI MPPC
//ADC
const double C_Pedestal          =    0.; //Pedeltal of ADC counts
const double C_WMPedestal        =    0.; //Pedeltal of ADC counts
const double C_Gain              =    10.;//Gain ADC counts of high gain channel
const double C_WMGain            =    10.;//Gain ADC counts of high gain channel
const double C_LowGain           =    1.; //Gain ADC counts of low gain channel
const double C_ADCtoCharge       = 135.5; //ADC to Charge
const double C_LowADCtoCharge    = 14.29; //ADC to Charge
const double C_LowGain_corr      = 14.29/13.55;
const double C_NonLinADCTh[4]    = {0.65, 3.2, 4.2, 14.};
const double C_NonLinADC1[4]     = {135.5, 217., 158.6, 5.1};
const double C_NonLinADC2[4]     = {0., -53., 133.9, 778.6};
const double C_ADCSaturation     = 850.;
const double C_WMADCSaturation   = 850.;
const double C_LowNonLinADCTh[4] = {7., 27., 35.5, 178.4};
const double C_LowNonLinADC1[4]  = {14.29, 26., 21.12, 0.7};
const double C_LowNonLinADC2[4]  = {0., -82., 48.24, 775.1};
const double C_LowADCSaturation  = 900.;
const double C_WMLowADCSaturation  = 900.;
const double C_ElecNoise         = 1.7; //Sigma of high gain electronics noise
const double C_LowElecNoise      = 1.2; //Sigma of low gain electronics noise

#endif
