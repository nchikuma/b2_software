#include "HitSummary.h"

//......................................................................

HitSummary::HitSummary():
  nsimhits(0)
{ 
  view     = -1;
  mod      = -1;
  pln      = -1;
  ch       = -1;
  adc      = -1;
  loadc    = -1;
  cyc      = -1;
  pe       = -1.e5;
  lope     = -1.e5;
  pecorr   = -1.e5;
  vise     = -1.e5;
  visecorr = -1.e5;
  tdc      = -1.e5;
  time     = -1.e5;
  tnearhit = -1.e5;
  timecorr = -1.e5;
  xy       = -1.e5;
  z        = -1.e5;
  addbasicrecon  = 0;
  dummy          = 0;
  gocosmic = false;
  hitcosmic= false;
  isohit   = false;
  gridcell_id_x1=0; //for cross talk study 2016/1/13
  gridcell_id_x2=0; //for cross talk study 2016/1/13
  gridcell_id_y1=0; //for cross talk study 2016/1/13
  gridcell_id_y2=0; //for cross talk study 2016/1/13
  pe_cross=0; //for cross talk study 2016/1/13
  
  pathlength = 0; // for PE correction, 2018/04/30
  pe_permm = 0;   // for PE correction, 2018/04/30
}

//......................................................................

HitSummary::HitSummary(const HitSummary& hit) :
  nsimhits(0)
{ 

  mod      = hit.mod;
  view     = hit.view;
  pln      = hit.pln;
  ch       = hit.ch;
  cyc      = hit.cyc;
  adc      = hit.adc;
  loadc    = hit.loadc;
  pe       = hit.pe;
  lope     = hit.lope;
  pecorr   = hit.pecorr;
  vise     = hit.vise;
  visecorr = hit.visecorr;
  tdc      = hit.tdc;
  time     = hit.time;
  tnearhit = hit.tnearhit;
  timecorr = hit.timecorr;
  xy       = hit.xy;
  z        = hit.z;
  addbasicrecon  = hit.addbasicrecon;
  dummy          = hit.dummy;
  gocosmic       = hit.gocosmic;
  hitcosmic      = hit.hitcosmic;
  isohit         = hit.isohit;
  gridcell_id_x1 = hit.gridcell_id_x1; 
  gridcell_id_x2 = hit.gridcell_id_x2; 
  gridcell_id_y1 = hit.gridcell_id_y1; 
  gridcell_id_y2 = hit.gridcell_id_y2; 
  pe_cross       = hit.pe_cross;
  pathlength     = hit.pathlength;
  pe_permm       = hit.pe_permm;

  for (int i=0; i<HIT_MAXSIMHITS; ++i) {
    fSimHit[i] = TRef(NULL);
  }

  for (int i=0; i < hit.nsimhits; ++i) 
    AddSimHit(hit.GetSimHit(i));
}

//......................................................................

HitSummary::~HitSummary() 
{
}

//......................................................................


void HitSummary::Clear(Option_t* option)
{
  for (int i=0; i<HIT_MAXSIMHITS; ++i)
    fSimHit[i] = TRef(NULL);
  nsimhits = 0;
}

//......................................................................

void HitSummary::Print()
{
  std::cout << "-------------------" <<std::endl
    << "Module# = " << mod   <<std::endl
    << "View    = " << view  <<std::endl
    << "Plane # = " << pln   <<std::endl
    << "Ch. #   = " << ch    <<std::endl
    << "ADC     = " << adc   <<std::endl
    << "TDC     = " << tdc   <<std::endl
    << "xy      = " << xy    <<std::endl
    << "z       = " << z     <<std::endl
    << std::endl;
}


//......................................................................

void HitSummary::AddSimHit(SimHitSummary* sbhitsum) 
{
  if (nsimhits < HIT_MAXSIMHITS) {
    fSimHit[nsimhits] = TRef((SimHitSummary*) sbhitsum);
    ++nsimhits;
  }
}
//......................................................................


SimHitSummary* HitSummary::GetSimHit(int i) const
{ 
  return (SimHitSummary*)fSimHit[i].GetObject();
}

//......................................................................


ClassImp(HitSummary)

  ////////////////////////////////////////////////////////////////////////
