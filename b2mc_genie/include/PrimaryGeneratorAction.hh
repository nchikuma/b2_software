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



class GENIE_parameters
{
  public:
    int    Iev         ; // Event number
    int    Neutrino    ; // Neutrino pdg code
    int    FSPrimLept  ; // Final state primary lepton pdg code
    int    Target      ; // Nuclear target pdg code (10LZZZAAAI)
    int    TargetZ     ; // Nuclear target Z (extracted from pdg code above)
    int    TargetA     ; // Nuclear target A (extracted from pdg code above)
    int    HitNuc      ; // Hit nucleon pdg code      (not set for COH,IMD and NuEL events)
    int    HitQrk      ; // Hit quark pdg code        (set for DIS events only)
    bool   FromSea     ; // Hit quark is from sea     (set for DIS events only)
    int    ResId       ; // Produced baryon resonance (set for resonance events only)
    bool   IsQel       ; // Is QEL?
    bool   IsRes       ; // Is RES?
    bool   IsDis       ; // Is DIS?
    bool   IsCoh       ; // Is Coherent?
    bool   IsMec       ; // Is MEC?
    bool   IsDfr       ; // Is Diffractive?
    bool   IsImd       ; // Is IMD?
    bool   IsImdAnh    ; // Is IMD annihilation?
    bool   IsNuEL      ; // Is ve elastic?
    bool   IsEM        ; // Is EM process?
    bool   IsCC        ; // Is Weak CC process?
    bool   IsNC        ; // Is Weak NC process?
    bool   IsCharmPro  ; // Produces charm?
    int    CodeNeut    ; // The equivalent NEUT reaction code (if any)
    int    CodeNuance  ; // The equivalent NUANCE reaction code (if any)
    double Weight      ; // Event weight
    double KineXs      ; // Bjorken x as was generated during kinematical selection                                          ; takes fermi momentum / off-shellness into account
    double KineYs      ; // Inelasticity y as was generated during kinematical selection                                     ; takes fermi momentum / off-shellness into account
    double KineTs      ; // Energy transfer to nucleus at COH events as was generated during kinematical selection
    double KineQ2s     ; // Momentum transfer Q^2 as was generated during kinematical selection                              ; takes fermi momentum / off-shellness into account
    double KineWs      ; // Hadronic invariant mass W as was generated during kinematical selection                          ; takes fermi momentum / off-shellness into account
    double KineX       ; // Experimental-like Bjorken x                                                                      ; neglects fermi momentum / off-shellness
    double KineY       ; // Experimental-like inelasticity y                                                                 ; neglects fermi momentum / off-shellness
    double KineT       ; // Experimental-like energy transfer to nucleus at COH events
    double KineQ2      ; // Experimental-like momentum transfer Q^2                                                          ; neglects fermi momentum / off-shellness
    double KineW       ; // Experimental-like hadronic invariant mass W                                                      ; neglects fermi momentum / off-shellness
    double EvRF        ; // Neutrino energy @ the rest-frame of the hit-object (eg nucleon for CCQE, e- for ve- elastic,...)
    double Ev          ; // Neutrino energy @ LAB
    double Pxv         ; // Neutrino px @ LAB
    double Pyv         ; // Neutrino py @ LAB
    double Pzv         ; // Neutrino pz @ LAB
    double En          ; // Initial state hit nucleon energy @ LAB
    double Pxn         ; // Initial state hit nucleon px @ LAB
    double Pyn         ; // Initial state hit nucleon py @ LAB
    double Pzn         ; // Initial state hit nucleon pz @ LAB
    double El          ; // Final state primary lepton energy @ LAB
    double Pxl         ; // Final state primary lepton px @ LAB
    double Pyl         ; // Final state primary lepton py @ LAB
    double Pzl         ; // Final state primary lepton pz @ LAB
    double Pl          ; // Final state primary lepton p  @ LAB
    double Costhl      ; // Final state primary lepton cos(theta) wrt to neutrino direction
    int    NfP         ; // Nu. of final state p's + \bar{p}'s (after intranuclear rescattering)
    int    NfN         ; // Nu. of final state n's + \bar{n}'
    int    NfPip       ; // Nu. of final state pi+'s
    int    NfPim       ; // Nu. of final state pi-'s
    int    NfPi0       ; // Nu. of final state pi0's (
    int    NfKp        ; // Nu. of final state K+'s
    int    NfKm        ; // Nu. of final state K-'s
    int    NfK0        ; // Nu. of final state K0's + \bar{K0}'s
    int    NfEM        ; // Nu. of final state gammas and e-/e+
    int    NfOther     ; // Nu. of heavier final state hadrons (D+/-,D0,Ds+/-,Lamda,Sigma,Lamda_c,Sigma_c,...)
    int    NiP         ; // Nu. of `primary' (: before intranuclear rescattering) p's + \bar{p}'s
    int    NiN         ; // Nu. of `primary' n's + \bar{n}'s
    int    NiPip       ; // Nu. of `primary' pi+'s
    int    NiPim       ; // Nu. of `primary' pi-'s
    int    NiPi0       ; // Nu. of `primary' pi0's
    int    NiKp        ; // Nu. of `primary' K+'s
    int    NiKm        ; // Nu. of `primary' K-'s
    int    NiK0        ; // Nu. of `primary' K0's + \bar{K0}'s
    int    NiEM        ; // Nu. of `primary' gammas and e-/e+
    int    NiOther     ; // Nu. of other `primary' hadron shower particles
    int    Nf          ; // Nu. of final state particles in hadronic system
    int    Pdgf  [250] ; // Pdg code of k^th final state particle in hadronic system
    double Ef    [250] ; // Energy     of k^th final state particle in hadronic system @ LAB
    double Pxf   [250] ; // Px         of k^th final state particle in hadronic system @ LAB
    double Pyf   [250] ; // Py         of k^th final state particle in hadronic system @ LAB
    double Pzf   [250] ; // Pz         of k^th final state particle in hadronic system @ LAB
    double Pf    [250] ; // P          of k^th final state particle in hadronic system @ LAB
    double Costhf[250] ; // cos(theta) of k^th final state particle in hadronic system @ LAB wrt to neutrino direction
    int    Ni          ; // Nu. of particles in 'primary' hadronic system (before intranuclear rescattering)
    int    Pdgi[250]   ; // Pdg code of k^th particle in 'primary' hadronic system
    int    Resc[250]   ; // FSI code of k^th particle in 'primary' hadronic system
    double Ei  [250]   ; // Energy   of k^th particle in 'primary' hadronic system @ LAB
    double Pxi [250]   ; // Px       of k^th particle in 'primary' hadronic system @ LAB
    double Pyi [250]   ; // Py       of k^th particle in 'primary' hadronic system @ LAB
    double Pzi [250]   ; // Pz       of k^th particle in 'primary' hadronic system @ LAB
    double VtxX        ; // Vertex x in detector coord system (SI)
    double VtxY        ; // Vertex y in detector coord system (SI)
    double VtxZ        ; // Vertex z in detector coord system (SI)
    double VtxT        ; // Vertex t in detector coord system (SI)
    double SumKEf      ; // Sum of kinetic energies of all final state particles
    double CalResp0    ; // Approximate calorimetric response to the hadronic system computed as sum of

};

void SetBranchGENIE(TTree* tree,GENIE_parameters& genie)
{
  tree->SetBranchAddress("iev"         ,          & (genie.Iev        )); 
  tree->SetBranchAddress("neu"         ,          & (genie.Neutrino   )); 
  tree->SetBranchAddress("fspl"        ,          & (genie.FSPrimLept )); 
  tree->SetBranchAddress("tgt"         ,          & (genie.Target     )); 
  tree->SetBranchAddress("Z"           ,          & (genie.TargetZ    )); 
  tree->SetBranchAddress("A"           ,          & (genie.TargetA    )); 
  tree->SetBranchAddress("hitnuc"      ,          & (genie.HitNuc     )); 
  tree->SetBranchAddress("hitqrk"      ,          & (genie.HitQrk     )); 
  tree->SetBranchAddress("resid"       ,          & (genie.ResId      )); 
  tree->SetBranchAddress("sea"         ,          & (genie.FromSea    )); 
  tree->SetBranchAddress("qel"         ,          & (genie.IsQel      )); 
  tree->SetBranchAddress("mec"         ,          & (genie.IsMec      )); 
  tree->SetBranchAddress("res"         ,          & (genie.IsRes      )); 
  tree->SetBranchAddress("dis"         ,          & (genie.IsDis      )); 
  tree->SetBranchAddress("coh"         ,          & (genie.IsCoh      )); 
  tree->SetBranchAddress("dfr"         ,          & (genie.IsDfr      )); 
  tree->SetBranchAddress("imd"         ,          & (genie.IsImd      )); 
  tree->SetBranchAddress("imdanh"      ,          & (genie.IsImdAnh   )); 
  tree->SetBranchAddress("nuel"        ,          & (genie.IsNuEL     )); 
  tree->SetBranchAddress("em"          ,          & (genie.IsEM       )); 
  tree->SetBranchAddress("cc"          ,          & (genie.IsCC       )); 
  tree->SetBranchAddress("nc"          ,          & (genie.IsNC       )); 
  tree->SetBranchAddress("charm"       ,          & (genie.IsCharmPro )); 
  tree->SetBranchAddress("neut_code"   ,          & (genie.CodeNeut   )); 
  tree->SetBranchAddress("nuance_code" ,          & (genie.CodeNuance )); 
  tree->SetBranchAddress("wght"        ,          & (genie.Weight     )); 
  tree->SetBranchAddress("xs"          ,          & (genie.KineXs     )); 
  tree->SetBranchAddress("ys"          ,          & (genie.KineYs     )); 
  tree->SetBranchAddress("ts"          ,          & (genie.KineTs     )); 
  tree->SetBranchAddress("Q2s"         ,          & (genie.KineQ2s    )); 
  tree->SetBranchAddress("Ws"          ,          & (genie.KineWs     )); 
  tree->SetBranchAddress("x"           ,          & (genie.KineX      )); 
  tree->SetBranchAddress("y"           ,          & (genie.KineY      )); 
  tree->SetBranchAddress("t"           ,          & (genie.KineT      )); 
  tree->SetBranchAddress("Q2"          ,          & (genie.KineQ2     )); 
  tree->SetBranchAddress("W"           ,          & (genie.KineW      )); 
  tree->SetBranchAddress("EvRF"        ,          & (genie.EvRF       )); 
  tree->SetBranchAddress("Ev"          ,          & (genie.Ev         )); 
  tree->SetBranchAddress("pxv"         ,          & (genie.Pxv        )); 
  tree->SetBranchAddress("pyv"         ,          & (genie.Pyv        )); 
  tree->SetBranchAddress("pzv"         ,          & (genie.Pzv        )); 
  tree->SetBranchAddress("En"          ,          & (genie.En         )); 
  tree->SetBranchAddress("pxn"         ,          & (genie.Pxn        )); 
  tree->SetBranchAddress("pyn"         ,          & (genie.Pyn        )); 
  tree->SetBranchAddress("pzn"         ,          & (genie.Pzn        )); 
  tree->SetBranchAddress("El"          ,          & (genie.El         )); 
  tree->SetBranchAddress("pxl"         ,          & (genie.Pxl        )); 
  tree->SetBranchAddress("pyl"         ,          & (genie.Pyl        )); 
  tree->SetBranchAddress("pzl"         ,          & (genie.Pzl        )); 
  tree->SetBranchAddress("pl"          ,          & (genie.Pl         )); 
  tree->SetBranchAddress("cthl"        ,          & (genie.Costhl     )); 
  tree->SetBranchAddress("nfp"         ,          & (genie.NfP        )); 
  tree->SetBranchAddress("nfn"         ,          & (genie.NfN        )); 
  tree->SetBranchAddress("nfpip"       ,          & (genie.NfPip      )); 
  tree->SetBranchAddress("nfpim"       ,          & (genie.NfPim      )); 
  tree->SetBranchAddress("nfpi0"       ,          & (genie.NfPi0      )); 
  tree->SetBranchAddress("nfkp"        ,          & (genie.NfKp       )); 
  tree->SetBranchAddress("nfkm"        ,          & (genie.NfKm       )); 
  tree->SetBranchAddress("nfk0"        ,          & (genie.NfK0       )); 
  tree->SetBranchAddress("nfem"        ,          & (genie.NfEM       )); 
  tree->SetBranchAddress("nfother"     ,          & (genie.NfOther    )); 
  tree->SetBranchAddress("nip"         ,          & (genie.NiP        )); 
  tree->SetBranchAddress("nin"         ,          & (genie.NiN        )); 
  tree->SetBranchAddress("nipip"       ,          & (genie.NiPip      )); 
  tree->SetBranchAddress("nipim"       ,          & (genie.NiPim      )); 
  tree->SetBranchAddress("nipi0"       ,          & (genie.NiPi0      )); 
  tree->SetBranchAddress("nikp"        ,          & (genie.NiKp       )); 
  tree->SetBranchAddress("nikm"        ,          & (genie.NiKm       )); 
  tree->SetBranchAddress("nik0"        ,          & (genie.NiK0       )); 
  tree->SetBranchAddress("niem"        ,          & (genie.NiEM       )); 
  tree->SetBranchAddress("niother"     ,          & (genie.NiOther    )); 
  tree->SetBranchAddress("ni"          ,          & (genie.Ni         )); 
  tree->SetBranchAddress("pdgi"        ,            (genie.Pdgi       )); 
  tree->SetBranchAddress("resc"        ,            (genie.Resc       )); 
  tree->SetBranchAddress("Ei"          ,            (genie.Ei         )); 
  tree->SetBranchAddress("pxi"         ,            (genie.Pxi        )); 
  tree->SetBranchAddress("pyi"         ,            (genie.Pyi        )); 
  tree->SetBranchAddress("pzi"         ,            (genie.Pzi        )); 
  tree->SetBranchAddress("nf"          ,          & (genie.Nf         )); 
  tree->SetBranchAddress("pdgf"        ,            (genie.Pdgf       )); 
  tree->SetBranchAddress("Ef"          ,            (genie.Ef         )); 
  tree->SetBranchAddress("pxf"         ,            (genie.Pxf        )); 
  tree->SetBranchAddress("pyf"         ,            (genie.Pyf        )); 
  tree->SetBranchAddress("pzf"         ,            (genie.Pzf        )); 
  tree->SetBranchAddress("pf"          ,            (genie.Pf         )); 
  tree->SetBranchAddress("cthf"        ,            (genie.Costhf     )); 
  tree->SetBranchAddress("vtxx"        ,          & (genie.VtxX       )); 
  tree->SetBranchAddress("vtxy"        ,          & (genie.VtxY       )); 
  tree->SetBranchAddress("vtxz"        ,          & (genie.VtxZ       )); 
  tree->SetBranchAddress("vtxt"        ,          & (genie.VtxT       )); 
  tree->SetBranchAddress("sumKEf"      ,          & (genie.SumKEf     )); 
  tree->SetBranchAddress("calresp0"    ,          & (genie.CalResp0   )); 

}


class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  //PrimaryGeneratorAction(Neut*, RunAction*, EventAction*, int, int);
  PrimaryGeneratorAction(RunAction*, EventAction*, int, int);
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
  
  //Neut *neut_file;
  TTree *genie_tree;
  GENIE_parameters genie_para;

  int iEvt;
  int NumEvt;

};

#endif




