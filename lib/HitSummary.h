#ifndef HITSUMMARY_H
#define HITSUMMARY_H
#include <iostream>
#include "TObject.h"

#include "TRefArray.h"
#include "TRef.h"

#include "SimHitSummary.h"

#define HIT_MAXSIMHITS 1

//......................................................................

class HitSummary : public TObject
{
public:
    HitSummary();
    HitSummary(const HitSummary& hit);
    virtual ~HitSummary();
    
    void Clear   (Option_t* option="");
    void Print();

public:

    int mod;                  // Module ID
                              // 0~6:Horizontal, 7~13:Vertical
                              // 14, 15: Off-axis  module
                              // 16: proton module??
    int cyc;                  // Integration number 0~22
                              // usually 1st beam bunch <-> 4 cycle

    int view;                 // xz view (0), yz view (1)
    int pln;                  // plane number
    int ch;                   // strip number(0~23)
    int adc;                  // high gain ADC value
    int loadc;                // low  gain ADC value
    float pe;                 // number of photoelectrons, without correction
                              // for light attenuation in fiber,
                              //     MPPC crosstalk&afterpulse,
                              //     MPPC liniarity
                              //     Elec. linearity
    float lope;               // number of photoelectrons, without correction
                              // for light attenuation in fiber,
                              //     MPPC crosstalk&afterpulse,
                              //     MPPC liniarity
                              //     Elec. linearity
    float Pe(){ if(pe<50)return pe; else return lope;}
    float pecorr;             // number of photoelectrons, with correction
                              // after correction
                              // (needs 3D reconstruction)
    float vise;               // visible energy (MeV), without correction
                              // before correction
    float visecorr;           // visible energy (MeV), with correction
                              // after  correction
    long   tdc;                // raw TDC value
    float time;               // hit time (ns).
                              // before correction 
                              // for light propagation in fiber
    float tnearhit;           // minumum value of hit time difference (ns).
                              // comparing other hit

    float timecorr;           // hit time (ns). Only the first hit inside 
                              // after  correction
    float xy;                 // transverse position (cm), x or y
                              // depending on view
    float z;                  // position (cm) along beam direction
    int   addbasicrecon;      // This hit is member of basic recon or not
    int   dummy;              // this is dummy(study for MPPC noise)

    bool  gocosmic;           // For efficiency study with cosmic, it is denominator for efficiency
    bool  hitcosmic;          // For efficiency study with cosmic, it is numerator for efficiency
    bool  isohit;

    int gridcell_id_x1; //for cross talk study 2016/1/13
    int gridcell_id_x2; //for cross talk study 2016/1/13
    int gridcell_id_y1; //for cross talk study 2016/1/13
    int gridcell_id_y2; //for cross talk study 2016/1/13
    float pe_cross; //for cross talk study 2016/1/13
    
    float pathlength; // for PE correction, 2018/04/30
    float pe_permm;   // for PE correction, 2018/04/30

    int NSimHits() const {return nsimhits;}
    // number of sim hits associated to this reco hit

    void AddSimHit(SimHitSummary* hit);
    SimHitSummary* GetSimHit(int i) const;

private:

    int nsimhits;
    TRef fSimHit[HIT_MAXSIMHITS];

    ClassDef(HitSummary, 7) //  Hit Summary
};

#endif // HITSUMMARY_H
////////////////////////////////////////////////////////////////////////
