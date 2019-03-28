#ifndef FIRSTREDUCSUMMARY_H
#define FIRSTREDUCSUMMARY_H
#include <iostream>
#include "TObject.h"
#include "TClonesArray.h"

#include "TRefArray.h"
#include "TRef.h"
#include "vector"
using namespace std;

#include "HitSummary.h"
#define HIT_MAX_1STREDUC 1000 //## Maximum # of HitSummaries 

//......................................................................
class FirstReducSummary : public TObject{
public:
    int            hitmod;  //### hit module
    int            hitcyc;  //### hit cycle
    unsigned int  xhitbit;  //### hit X layer 11(TPL) bits
    unsigned int  yhitbit;  //### hit Y layer 11(TPL) bits
    int          nhitxlyr;  //### number of hit X layers
    int          nhitylyr;  //### number of hit Y layers
    float          xtotpe;  //### total p.e. in X layers
    float          ytotpe;  //### total p.e. in Y layers
    bool       xtracklike;  //### definition of track-like will be changed or removed
    bool       ytracklike;  //### now(2010/2/12), 3 layers hit in a row 

    vector<int>      hity;  //### temporary for stdudying external(?) redioactivity
    vector<int>      hitx;
    vector<int>      hitxz;
    vector<int>      hityz;


    FirstReducSummary();
    FirstReducSummary(const FirstReducSummary& basicsum);
    virtual ~FirstReducSummary();
    void Clear   (Option_t* option="");
    void Print();
    void AddHit(HitSummary* hit);
    HitSummary* GetHit(int i) const;
    int  nhits;   
private:
    TRef fHit[HIT_MAX_1STREDUC];
    ClassDef(FirstReducSummary, 4)
};

#endif // HITSUMMARY_H
////////////////////////////////////////////////////////////////////////
