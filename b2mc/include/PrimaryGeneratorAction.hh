#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleDefinition.hh"

#include "RunAction.hh"
#include "EventAction.hh"

#include "Neut.hh"
#include "Const.hh"

class G4ParticleGun;
class G4Event;

const double HallX             = C_HallDirtPosX/10.; //cm 
const double HallZ             = C_HallDirtPosZ/10.;    //cm
const double HallRadiusMax     = C_HallDirtRadiusMax/10.;//cm
const double HallRadiusMin     = C_HallDirtRadiusMin/10.;//cm

const double INGRID_VMODOFFSET = C_INGVMotherPosZ/10.; //cm
const double PM_OFFSET         = C_PMMotherPosZ/10.; //cm
const double WM_OFFSET         = C_PMMotherPosZ/10.; //cm //same as ProtonModule
const double WMBG_OFFSET       = C_PMMotherPosZ/10.; //cm //same as ProtonModule
const double PM_OFFSETX        = C_PMMotherPosX/10.; //cm
const double PM_OFFSETY        = C_PMMotherPosY/10.; //cm
const double PM_OFFSETZ        = C_PMMotherPosZ/10.; //cm
const double B2WM_OFFSETX      = (C_B2MotherPosX+C_B2WMPosX)/10. ; //cm
const double B2WM_OFFSETY      = (C_B2MotherPosY+C_B2WMPosY)/10. ; //cm
const double B2WM_OFFSETZ      = (C_B2MotherPosZ+C_B2WMPosZ)/10. ; //cm
const double B2CH_OFFSETX      = (C_B2MotherPosX+C_B2CHPosX)/10. ; //cm
const double B2CH_OFFSETY      = (C_B2MotherPosY+C_B2CHPosY)/10. ; //cm
const double B2CH_OFFSETZ      = (C_B2MotherPosZ+C_B2CHPosZ)/10. ; //cm
const double B2ING_OFFSETX     = (C_B2MotherPosX+C_B2INGPosX)/10. ; //cm
const double B2ING_OFFSETY     = (C_B2MotherPosY+C_B2INGPosY)/10. ; //cm
const double B2ING_OFFSETZ     = (C_B2MotherPosZ+C_B2INGPosZ)/10. ; //cm

//BG
const double B2Wall_OFFSETY    = 17.38; //cm
const double B2Wall_OFFSETZ    = 680.8; //cm

//INGRID
const double INGRID_MODSpace     = C_INGSpace/10.; //cm
const double total_mass_fe       = C_INGTotMassIron; //ton
const double total_mass_sci      = C_INGTotMassScinti; //ton
const double scinti_start        = C_INGPlnStart/10.; //cm 
                                    // the center of the first scintillator
const double iron_start          = (C_INGIronStart)/10.; //cm 
                                   // the center fo the first iron plate
const double width_fe            = C_INGIronThick/10.; // cm
const double width_sci           = C_INGScintiThick/10.; // cm
const double GAP                 = C_INGPlnDist/10.; // 6.5 + 2 + 2.2
const double scintileng_ing      = C_INGScintiLength/10.; //cm
const double ironxyleng_ing      = C_INGIronXY/10.; //cm
const double iron_density        = 7.874; //g/cm3
const double scinti_density      = 1.03; //g/cm3
//for Mod03
const double scinti_start03      = scinti_start + 1.0;
const double iron_start03        = iron_start + 0.9;
const double GAP03               = GAP - 0.2;


//for Proton Module
const double total_mass_sci_pm   = C_PMTotMassScinti;//ton (total mass)
const double total_mass_front_pm = C_PMTotMassVetoSci;//ton (total mass of front veto)
const double ingrid_width        = C_INGScintiThick/10.;//(cm) (total mass of ingrid type)
const double scibar_width        = C_PMScintiThick/10.; //(cm) (total mass of scibar type)
const double scibar_region       = C_PMScintiThick/10.*2.*(C_PMNumPln-1);//44.2cm
const double sciing_region       = (C_INGScintiThick+C_PMScintiThick)/10.*(C_PMNumPln-1);//39.1cm
const double ingrid_region       = C_INGScintiThick/10.*2.*(C_PMNumPln-1); //34cm
const double Pscinti_start       = C_PMPlnStart/10.;//(cm)
const double distance_first      = C_PMPlnDist_First/10.;//(cm)
const double distance_pln        = C_PMPlnDist/10.;//(cm)
const double scintileng_pm       = C_PMScintiLength/10.; //cm

//for Water Module
const double watertank_x         = C_WMWaterTargetSizeX/10.; //cm
const double watertank_y         = C_WMWaterTargetSizeY/10.; //cm
const double watertank_z         = C_WMWaterTargetSizeZ/10.; //cm
const double altank_x            = C_WMAlTankSizeX/10.; //cm
const double altank_y            = C_WMAlTankSizeY/10.; //cm
const double altank_z            = C_WMAlTankSizeZ/10.; //cm




class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction(Neut*, RunAction*, EventAction*, int, int);
  ~PrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event* anEvent);

private:
  G4ParticleGun* particleGun;
  G4ParticleTable* particleTable;

  int module_mode;
  int neutrino_flavor;

  RunAction* runaction;
  EventAction* eventaction;
  
  Neut *neut_file;
};

#endif
