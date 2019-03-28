#include "FirstReducSummary.h"


//......................................................................

FirstReducSummary::FirstReducSummary():
  nhits(0)
{ 
  hitmod   =    -1;
  hitcyc   =    -1;
  xhitbit  = 0x000;
  xhitbit  = 0x000;
  nhitxlyr =    -1;
  nhitylyr =    -1;
  xtotpe   = -1e-5;
  ytotpe   = -1e-5;
  xtracklike = false;
  ytracklike = false;

  hitx.clear();
  hity.clear();
  hitxz.clear();
  hityz.clear();
}
//......................................................................
FirstReducSummary::FirstReducSummary(const FirstReducSummary& reduc) :
  nhits(0)
{ 
  hitmod    =  reduc.hitmod;
  hitcyc    =  reduc.hitcyc;
  xhitbit   =  reduc.xhitbit;
  yhitbit   =  reduc.yhitbit;
  nhitxlyr  =  reduc.nhitxlyr;
  nhitylyr  =  reduc.nhitylyr;
  xtotpe    =  reduc.xtotpe;
  ytotpe    =  reduc.ytotpe;
  xtracklike=  reduc.xtracklike;
  ytracklike=  reduc.ytracklike;
  nhits     = 0;
  for (int i=0; i < HIT_MAX_1STREDUC; ++i) {
    fHit[i] = TRef(NULL);
  }
  for (int i=0; i < reduc.nhits; ++i) 
    AddHit(reduc.GetHit(i));

  hitx  = reduc.hitx;
  hity  = reduc.hity;
  hitxz = reduc.hitxz;
  hityz = reduc.hityz;
  //for (int i=0; i < reduc.hitx.size(); ++i)
  //hitx.push_back( reduc.hitx[i] );



}

//......................................................................



FirstReducSummary::~FirstReducSummary() 
{

}

//......................................................................


void FirstReducSummary::Clear(Option_t* option)
{
  for (int i=0; i<HIT_MAX_1STREDUC; ++i)
    fHit[i] = TRef(NULL);
  nhits = 0;

}

//......................................................................

void FirstReducSummary::Print()
{

}
//......................................................................
void FirstReducSummary::AddHit(HitSummary* hit) 
{

    if (nhits < HIT_MAX_1STREDUC) {
        fHit[nhits] = TRef((HitSummary*) hit);
        ++nhits;
    }

}
//......................................................................


HitSummary* FirstReducSummary::GetHit(int i) const
{ 

  return (HitSummary*)fHit[i].GetObject();
}


//......................................................................


ClassImp(FirstReducSummary)

////////////////////////////////////////////////////////////////////////
