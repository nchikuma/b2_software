#include <iostream>
#include "wgGetTree.h"
#include "wgErrorCode.h"
#include "Const.h"
#include "TH1F.h"

using namespace std;

TFile* wgGetTree::finput;    
string wgGetTree::finputname;
TTree* wgGetTree::tree;
TTree* wgGetTree::info;

TFile* wgGetTree::foutput;    
string wgGetTree::foutputname;
TTree* wgGetTree::tree_out;

TFile* wgGetTree::ingfinput;    
string wgGetTree::ingfinputname;
TTree* wgGetTree::ingtree;

//************************************************************************
wgGetTree::wgGetTree(){
}

//************************************************************************
wgGetTree::wgGetTree(string& str){
  this->Open(str);
}

//************************************************************************
wgGetTree::wgGetTree(string& str,IngRecon_t& ingrecon){
  cout << "INGRID tree date to be opened." << endl;
  CheckExist *Check = new CheckExist;
  if(!Check->RootFile(str)){
    cout << "ERROR!! FAIL TO SET TREEFILE (IngRecon_t)" <<endl;
    cout << "Reading file : " << str << endl;
  }
  else{
    wgGetTree::ingfinputname = str;
    if(wgGetTree::ingfinput){
      cout << "WARNING!! TFile is overwrited!" <<endl;
      this->Close();
    }
    wgGetTree::ingfinput = new TFile(str.c_str(),"open");
    wgGetTree::ingtree = (TTree*)wgGetTree::ingfinput->Get("tree");
    this->SetTreeFile(ingrecon);
  }
  delete Check;
}

//************************************************************************
bool wgGetTree::Open(string& str){
  CheckExist *Check = new CheckExist;
  wgGetTree::finputname = str;
  if(!Check->RootFile(str)){
    cout << "ERROR!! FAIL TO OPEN ROOT FILE" <<endl;
    cout << "Reading file : " << str << endl;
    return false;
  }
  delete Check;
  if(wgGetTree::finput){
    cout << "WARNING!! TFile is overwrited!" <<endl;
    this->Close();
  }
  wgGetTree::finput = new TFile(str.c_str(),"open");
  wgGetTree::tree = (TTree*)wgGetTree::finput->Get("tree");
  wgGetTree::info = (TTree*)wgGetTree::finput->Get("info");
  return true;
}

//************************************************************************
void wgGetTree::Close(){
  if(this->finput){
    //finput->cd();
    //wgGetTree::finput->Close();
    delete wgGetTree::finput;
  }
  /*
     if(wgGetTree::bsd.infile){
     wgGetTree::bsd.infile->cd();
     wgGetTree::bsd.infile->Close();
     delete wgGetTree::bsd.infile;
     }
     */
}

//************************************************************************
wgGetTree::~wgGetTree(){
  //this->Close();
}

//************************************************************************
bool wgGetTree::SetTreeFile(Raw_t& rdin){
  if(!wgGetTree::finput){
    cout << "WARNING!! TFile is not open!" <<endl;
    return false;
  }
  tree->SetBranchAddress("spill",&rdin.spill);
  tree->SetBranchAddress("spill_flag",&rdin.spill_flag);
  tree->SetBranchAddress("spill_mode",&rdin.spill_mode);
  tree->SetBranchAddress("spill_count",&rdin.spill_count);
  tree->SetBranchAddress("bcid",rdin.bcid);
  tree->SetBranchAddress("charge",rdin.charge);
  tree->SetBranchAddress("time",rdin.time);
  tree->SetBranchAddress("gs",rdin.gs);
  tree->SetBranchAddress("hit",rdin.hit);
  tree->SetBranchAddress("chipid",rdin.chipid);
  tree->SetBranchAddress("col",rdin.col);
  tree->SetBranchAddress("chipch",rdin.chipch);
  tree->SetBranchAddress("chip",rdin.chip);
  tree->SetBranchAddress("debug",rdin.debug);
  tree->SetBranchAddress("view",&rdin.view);
  tree->SetBranchAddress("pln",rdin.pln);
  tree->SetBranchAddress("ch",rdin.ch);
  tree->SetBranchAddress("grid",rdin.grid);
  tree->SetBranchAddress("x", rdin.x);
  tree->SetBranchAddress("y", rdin.y);
  tree->SetBranchAddress("z", rdin.z);
  tree->SetBranchAddress("time_ns",rdin.time_ns);
  tree->SetBranchAddress("pe",rdin.pe);
  tree->SetBranchAddress("gain",rdin.gain);
  tree->SetBranchAddress("pedestal",rdin.pedestal);
  tree->SetBranchAddress("ped_nohit",rdin.ped_nohit);
  return true;
}

//************************************************************************
bool wgGetTree::SetTreeFile(Hit_t& hit){
  if(!wgGetTree::finput){
    cout << "WARNING!! TFile is not open!" <<endl;
    return false;
  }
  tree->SetBranchAddress("num_hits"         ,&hit.num_hits          );
  tree->SetBranchAddress("num_bcid_cluster" ,&hit.num_bcid_cluster  );
  tree->SetBranchAddress("num_bcid_hits"    , hit.num_bcid_hits     );
  tree->SetBranchAddress("clustered_bcid"   , hit.clustered_bcid    );
  tree->SetBranchAddress("clustered_view"   , hit.clustered_view    );
  tree->SetBranchAddress("clustered_hitid"  , hit.clustered_hitid   );
  tree->SetBranchAddress("hit_bcid"         , hit.hit_bcid          );
  tree->SetBranchAddress("hit_view"         , hit.hit_view          );
  tree->SetBranchAddress("hit_pln"          , hit.hit_pln           );
  tree->SetBranchAddress("hit_ch"           , hit.hit_ch            );
  tree->SetBranchAddress("hit_dif"          , hit.hit_dif           );
  tree->SetBranchAddress("hit_chip"         , hit.hit_chip          );
  tree->SetBranchAddress("hit_chipch"       , hit.hit_chipch        );
  tree->SetBranchAddress("hit_sca"          , hit.hit_sca           );
  tree->SetBranchAddress("hit_adc"          , hit.hit_adc           );
  tree->SetBranchAddress("hit_gs"           , hit.hit_gs            );
  tree->SetBranchAddress("hit_tdc"          , hit.hit_tdc           );
  tree->SetBranchAddress("hit_tdc_mod"      , hit.hit_tdc_mod       );
  tree->SetBranchAddress("hit_pe"           , hit.hit_pe            );
  tree->SetBranchAddress("hit_time"         , hit.hit_time          );
  tree->SetBranchAddress("hit_pe_permm"     , hit.hit_pe_permm      );
  tree->SetBranchAddress("hit_pathlength"   , hit.hit_pathlength    );
  tree->SetBranchAddress("hit_ontrack"      , hit.hit_ontrack       );
  tree->SetBranchAddress("hit_grid"         , hit.hit_grid          );
  return true;
}

//************************************************************************
bool wgGetTree::MakeTreeFile(string& str, Hit_t& hit){
  if(!wgGetTree::foutput){
    wgGetTree::foutputname = str;
    wgGetTree::foutput  = new TFile(str.c_str(),"recreate");
    wgGetTree::tree_out = new TTree("tree","tree");
  }
  tree_out->Branch("num_hits"        ,&hit.num_hits  ,"num_hits/I");
  tree_out->Branch("num_bcid_cluster",&hit.num_bcid_cluster,"num_bcid_cluster/I");
  tree_out->Branch("num_time_cluster",&hit.num_time_cluster,"num_time_cluster/I");
  tree_out->Branch("num_bcid_hits"   , hit.num_bcid_hits  ,Form("num_bcid_hits[%d]/I"  ,MAX_NUM_BCID_CLUSTER));
  tree_out->Branch("clustered_bcid"  , hit.clustered_bcid ,Form("clustered_bcid[%d]/I" ,MAX_NUM_BCID_CLUSTER));
  tree_out->Branch("clustered_view"  , hit.clustered_view ,Form("clustered_view[%d]/I" ,MAX_NUM_BCID_CLUSTER));
  tree_out->Branch("clustered_hitid" , hit.clustered_hitid,Form("clustered_hitid[%d][%d]/I",MAX_NUM_BCID_CLUSTER,MAX_NUM_HIT));
  tree_out->Branch("hit_bcid"        , hit.hit_bcid       ,Form("hit_bcid[%d]/I"       ,MAX_NUM_HIT));
  tree_out->Branch("hit_view"        , hit.hit_view       ,Form("hit_view[%d]/I"       ,MAX_NUM_HIT));
  tree_out->Branch("hit_pln"         , hit.hit_pln        ,Form("hit_pln[%d]/I"        ,MAX_NUM_HIT));
  tree_out->Branch("hit_ch"          , hit.hit_ch         ,Form("hit_ch[%d]/I"         ,MAX_NUM_HIT));
  tree_out->Branch("hit_dif"         , hit.hit_dif        ,Form("hit_dif[%d]/I"        ,MAX_NUM_HIT));
  tree_out->Branch("hit_chip"        , hit.hit_chip       ,Form("hit_chip[%d]/I"       ,MAX_NUM_HIT));
  tree_out->Branch("hit_chipch"      , hit.hit_chipch     ,Form("hit_chipch[%d]/I"     ,MAX_NUM_HIT));
  tree_out->Branch("hit_sca"         , hit.hit_sca        ,Form("hit_sca[%d]/I"        ,MAX_NUM_HIT));
  tree_out->Branch("hit_adc"         , hit.hit_adc        ,Form("hit_adc[%d]/I"        ,MAX_NUM_HIT));
  tree_out->Branch("hit_gs"          , hit.hit_gs         ,Form("hit_gs[%d]/I"         ,MAX_NUM_HIT));
  tree_out->Branch("hit_tdc"         , hit.hit_tdc        ,Form("hit_tdc[%d]/I"        ,MAX_NUM_HIT));
  tree_out->Branch("hit_tdc_mod"     , hit.hit_tdc_mod    ,Form("hit_tdc_mod[%d]/I"    ,MAX_NUM_HIT));
  tree_out->Branch("hit_pe"          , hit.hit_pe         ,Form("hit_pe[%d]/D"         ,MAX_NUM_HIT));
  tree_out->Branch("hit_time"        , hit.hit_time       ,Form("hit_time[%d]/D"       ,MAX_NUM_HIT));
  tree_out->Branch("hit_pe_permm"    , hit.hit_pe_permm   ,Form("hit_pe_permm[%d]/D"   ,MAX_NUM_HIT));
  tree_out->Branch("hit_cos"         , hit.hit_cos        ,Form("hit_cos[%d]/D"        ,MAX_NUM_HIT));
  tree_out->Branch("hit_pathlength"  , hit.hit_pathlength ,Form("hit_pathlength[%d]/D" ,MAX_NUM_HIT));
  tree_out->Branch("hit_ontrack"     , hit.hit_ontrack    ,Form("hit_ontrack[%d]/O"    ,MAX_NUM_HIT));
  tree_out->Branch("hit_grid"        , hit.hit_grid       ,Form("hit_grid[%d]/O"       ,MAX_NUM_HIT));
  tree_out->Branch("hit_numtrack"    , hit.hit_numtrack   ,Form("hit_numtrack[%d]/I"   ,MAX_NUM_HIT));
  tree_out->Branch("hit_cluster_pe"  , hit.hit_cluster_pe ,Form("hit_cluster_pe[%d]/D" ,MAX_NUM_HIT));


  return true;
}

//************************************************************************
bool wgGetTree::SetTreeFile(Recon_t& recon){
  if(!wgGetTree::finput){
    cout << "WARNING!! TFile is not open!" <<endl;
    return false;
  }
  tree->SetBranchAddress("num_recon"       ,&recon.num_recon       );
  tree->SetBranchAddress("num_recon_hits"  , recon.num_recon_hits  );
  tree->SetBranchAddress("recon_hits_hitid", recon.recon_hits_hitid);
  tree->SetBranchAddress("recon_start_z"   , recon.recon_start_z   );
  tree->SetBranchAddress("recon_stop_z"    , recon.recon_stop_z    );
  tree->SetBranchAddress("recon_start_pln" , recon.recon_start_pln );
  tree->SetBranchAddress("recon_stop_pln"  , recon.recon_stop_pln  );
  tree->SetBranchAddress("recon_start_ch"  , recon.recon_start_ch  );
  tree->SetBranchAddress("recon_stop_ch"   , recon.recon_stop_ch   );
  tree->SetBranchAddress("recon_start_xy"  , recon.recon_start_xy  );
  tree->SetBranchAddress("recon_stop_xy"   , recon.recon_stop_xy   );
  tree->SetBranchAddress("recon_slope"     , recon.recon_slope     );
  tree->SetBranchAddress("recon_intercept" , recon.recon_intercept );
  tree->SetBranchAddress("recon_pathlength", recon.recon_pathlength);
  tree->SetBranchAddress("recon_total_pe"  , recon.recon_total_pe  );
  tree->SetBranchAddress("recon_mean_dedx" , recon.recon_mean_dedx );
  tree->SetBranchAddress("recon_mean_time" , recon.recon_mean_time );
  tree->SetBranchAddress("recon_view"      , recon.recon_view      );
  tree->SetBranchAddress("recon_bcid"      , recon.recon_bcid      );
  tree->SetBranchAddress("recon_bcid_id"   , recon.recon_bcid_id    );
  tree->SetBranchAddress("recon_veto"      , recon.recon_veto      );
  tree->SetBranchAddress("recon_sideescape", recon.recon_sideescape);
  tree->SetBranchAddress("recon_cos"       , recon.recon_cos       );
  tree->SetBranchAddress("recon_len"       , recon.recon_len       );
  tree->SetBranchAddress("recon_chi2"      , recon.recon_chi2      );
  tree->SetBranchAddress("recon_connected" , recon.recon_connected );
  tree->SetBranchAddress("recon_miss_paired" , recon.recon_miss_paired );
  return true;
}


//************************************************************************
bool wgGetTree::SetTreeFile(BSD_t& bsd){
  if(!wgGetTree::finput){
    cout << "WARNING!! TFile is not open!" <<endl;
    return false;
  }
  tree->SetBranchAddress("spill"          ,&bsd.spill          );
  tree->SetBranchAddress("spill_mode"     ,&bsd.spill_mode     );
  tree->SetBranchAddress("spill_count"    ,&bsd.spill_count    );
  tree->SetBranchAddress("spill_match"    ,&bsd.spill_match    );
  tree->SetBranchAddress("spillbit_match" ,&bsd.spillbit_match );
  tree->SetBranchAddress("missing_spill"  ,&bsd.missing_spill  );
  tree->SetBranchAddress("nurun"          ,&(bsd.nurun));
  tree->SetBranchAddress("midas_event"    ,&(bsd.midas_event));
  tree->SetBranchAddress("mrrun"          ,&(bsd.mrrun));
  tree->SetBranchAddress("mrshot"         ,&(bsd.mrshot));
  tree->SetBranchAddress("spillnum"       ,&(bsd.spillnum));
  tree->SetBranchAddress("trg_sec"        , (bsd.trg_sec));
  tree->SetBranchAddress("trg_nano"       , (bsd.trg_nano));
  tree->SetBranchAddress("gpsstat"        , (bsd.gpsstat));
  tree->SetBranchAddress("ct_np"          , (bsd.ct_np));
  tree->SetBranchAddress("beam_time"      , (bsd.beam_time));
  tree->SetBranchAddress("beam_flag"      , (bsd.beam_flag));
  tree->SetBranchAddress("hct"            , (bsd.hct));
  tree->SetBranchAddress("tpos"           , (bsd.tpos));
  tree->SetBranchAddress("tdir"           , (bsd.tdir));
  tree->SetBranchAddress("tsize"          , (bsd.tsize));
  tree->SetBranchAddress("mumon"          , (bsd.mumon));
  tree->SetBranchAddress("otr"            , (bsd.otr));
  tree->SetBranchAddress("good_gps_flag"  ,&(bsd.good_gps_flag));
  tree->SetBranchAddress("trigger_flag"   ,&(bsd.trigger_flag));
  tree->SetBranchAddress("spill_flag"     ,&(bsd.spill_flag));
  tree->SetBranchAddress("good_spill_flag",&(bsd.good_spill_flag));
  tree->SetBranchAddress("target_eff"     , (bsd.target_eff));
  tree->SetBranchAddress("run_type"       ,&(bsd.run_type));
  tree->SetBranchAddress("magset_id"      ,&(bsd.magset_id));
  tree->SetBranchAddress("hctx"           , (bsd.hctx));
  tree->SetBranchAddress("htrans"         , (bsd.htrans));
  tree->SetBranchAddress("hps"            , (bsd.hps));
  return true;
}


//************************************************************************
bool wgGetTree::MakeTreeFile(string& str, Recon_t& recon){
  if(!wgGetTree::foutput){
    wgGetTree::foutputname = str;
    wgGetTree::foutput  = new TFile(str.c_str(),"recreate");
    wgGetTree::tree_out = new TTree("tree","tree");
  }
  tree_out->Branch("num_recon"        ,&recon.num_recon        ,"num_recon/I");
  tree_out->Branch("num_recon_hits"   , recon.num_recon_hits   ,Form("num_recon_hits[%d]/I"   ,MAX_NUM_TRACK));
  tree_out->Branch("recon_hits_hitid" , recon.recon_hits_hitid ,Form("recon_hits_hitid[%d][%d]/I",MAX_NUM_TRACK,MAX_NUM_TRACKHIT)); 
  tree_out->Branch("recon_start_z"    , recon.recon_start_z    ,Form("recon_start_z[%d]/D"    ,MAX_NUM_TRACK)); 
  tree_out->Branch("recon_stop_z"     , recon.recon_stop_z     ,Form("recon_stop_z[%d]/D"     ,MAX_NUM_TRACK)); 
  tree_out->Branch("recon_start_pln"  , recon.recon_start_pln  ,Form("recon_start_pln[%d]/I"  ,MAX_NUM_TRACK)); 
  tree_out->Branch("recon_stop_pln"   , recon.recon_stop_pln   ,Form("recon_stop_pln[%d]/I"   ,MAX_NUM_TRACK)); 
  tree_out->Branch("recon_start_ch"   , recon.recon_start_ch   ,Form("recon_start_ch[%d]/I"   ,MAX_NUM_TRACK)); 
  tree_out->Branch("recon_stop_ch"    , recon.recon_stop_ch    ,Form("recon_stop_ch[%d]/I"    ,MAX_NUM_TRACK)); 
  tree_out->Branch("recon_start_xy"   , recon.recon_start_xy   ,Form("recon_start_xy[%d]/D"   ,MAX_NUM_TRACK)); 
  tree_out->Branch("recon_stop_xy"    , recon.recon_stop_xy    ,Form("recon_stop_xy[%d]/D"    ,MAX_NUM_TRACK)); 
  tree_out->Branch("recon_slope"      , recon.recon_slope      ,Form("recon_slope[%d]/D"      ,MAX_NUM_TRACK));   
  tree_out->Branch("recon_intercept"  , recon.recon_intercept  ,Form("recon_intercept [%d]/D" ,MAX_NUM_TRACK)); 
  tree_out->Branch("recon_pathlength" , recon.recon_pathlength ,Form("recon_pathlength[%d]/D" ,MAX_NUM_TRACK)); 
  tree_out->Branch("recon_total_pe"   , recon.recon_total_pe   ,Form("recon_total_pe[%d]/D"   ,MAX_NUM_TRACK));  
  tree_out->Branch("recon_mean_dedx"  , recon.recon_mean_dedx  ,Form("recon_mean_dedx[%d]/D"  ,MAX_NUM_TRACK));  
  tree_out->Branch("recon_mean_time"  , recon.recon_mean_time  ,Form("recon_mean_time[%d]/D"  ,MAX_NUM_TRACK));  
  tree_out->Branch("recon_view"       , recon.recon_view       ,Form("recon_view[%d]/I"       ,MAX_NUM_TRACK));  
  tree_out->Branch("recon_bcid"       , recon.recon_bcid       ,Form("recon_bcid[%d]/I"       ,MAX_NUM_TRACK)); 
  tree_out->Branch("recon_bcid_id"    , recon.recon_bcid_id    ,Form("recon_bcid_id[%d]/I"    ,MAX_NUM_TRACK)); 
  tree_out->Branch("recon_veto"       , recon.recon_veto       ,Form("recon_veto[%d]/I"       ,MAX_NUM_TRACK)); 
  tree_out->Branch("recon_sideescape" , recon.recon_sideescape ,Form("recon_sideescape[%d]/I" ,MAX_NUM_TRACK)); 
  tree_out->Branch("recon_cos"        , recon.recon_cos        ,Form("recon_cos[%d]/D"        ,MAX_NUM_TRACK));   
  tree_out->Branch("recon_len"        , recon.recon_len        ,Form("recon_len[%d]/D"        ,MAX_NUM_TRACK));   
  tree_out->Branch("recon_chi2"       , recon.recon_chi2       ,Form("recon_chi2[%d]/D"       ,MAX_NUM_TRACK));   
  tree_out->Branch("recon_connected"  , recon.recon_connected  ,Form("recon_connected[%d]/O"  ,MAX_NUM_TRACK));   
  tree_out->Branch("recon_miss_paired", recon.recon_miss_paired,Form("recon_miss_paired[%d]/O",MAX_NUM_TRACK));   

  return true;
}

//************************************************************************
bool wgGetTree::SetTreeFile(Track_t& track){
  if(!wgGetTree::finput){
    cout << "WARNING!! TFile is not open!" <<endl;
    return false;
  }
  tree->SetBranchAddress("num_track"       ,&track.num_track       );
  tree->SetBranchAddress("num_trackid"     ,&track.num_trackid     );
  tree->SetBranchAddress("num_track_hits"  , track.num_track_hits  );
  tree->SetBranchAddress("track_hits_hitid", track.track_hits_hitid);
  tree->SetBranchAddress("track_recon_id"  , track.track_recon_id  );
  tree->SetBranchAddress("track_start_z"   , track.track_start_z   );
  tree->SetBranchAddress("track_stop_z"    , track.track_stop_z    );
  tree->SetBranchAddress("track_start_pln" , track.track_start_pln );
  tree->SetBranchAddress("track_stop_pln"  , track.track_stop_pln  );
  tree->SetBranchAddress("track_start_ch"  , track.track_start_ch  );
  tree->SetBranchAddress("track_stop_ch"   , track.track_stop_ch   );
  tree->SetBranchAddress("track_start_xy"  , track.track_start_xy  );
  tree->SetBranchAddress("track_stop_xy"   , track.track_stop_xy   );
  tree->SetBranchAddress("track_slope"     , track.track_slope     );
  tree->SetBranchAddress("track_intercept" , track.track_intercept );
  tree->SetBranchAddress("track_pathlength", track.track_pathlength);
  tree->SetBranchAddress("track_total_pe"  , track.track_total_pe  );
  tree->SetBranchAddress("track_mean_dedx" , track.track_mean_dedx );
  tree->SetBranchAddress("track_mean_time" , track.track_mean_time );
  tree->SetBranchAddress("track_view"      , track.track_view      );
  tree->SetBranchAddress("track_bcid"      , track.track_bcid      );
  tree->SetBranchAddress("track_bcid_id"   , track.track_bcid_id   );
  tree->SetBranchAddress("track_veto"      , track.track_veto      );
  tree->SetBranchAddress("track_sideescape", track.track_sideescape);
  tree->SetBranchAddress("track_cos_zen"   , track.track_cos_zen   );
  tree->SetBranchAddress("track_cos_azi"   , track.track_cos_azi   );
  tree->SetBranchAddress("track_len"       , track.track_len       );
  tree->SetBranchAddress("track_id"        , track.track_id        );
  tree->SetBranchAddress("track_likelihood", track.track_likelihood);
  tree->SetBranchAddress("track_ing_match" , track.track_ing_match );
  tree->SetBranchAddress("track_ingtrk_id" , track.track_ingtrk_id );
  tree->SetBranchAddress("track_ing_diff"  , track.track_ing_diff  );
  tree->SetBranchAddress("track_pm_match"  , track.track_pm_match  );
  tree->SetBranchAddress("track_pmtrk_id"  , track.track_pmtrk_id  );
  tree->SetBranchAddress("track_pm_diff"   , track.track_pm_diff   );
  tree->SetBranchAddress("track_ingcyc"    , track.track_ingcyc    );
  tree->SetBranchAddress("track_ing_stopz" , track.track_ing_stopz );
  tree->SetBranchAddress("track_ing_stopxy", track.track_ing_stopxy);

  return true;
}

//************************************************************************
bool wgGetTree::MakeTreeFile(string& str, Track_t& track){
  if(!wgGetTree::foutput){
    wgGetTree::foutputname = str;
    wgGetTree::foutput  = new TFile(str.c_str(),"recreate");
    wgGetTree::tree_out = new TTree("tree","tree");
  }
  tree_out->Branch("num_track"       ,&track.num_track       ,"num_track/I"); 
  tree_out->Branch("num_trackid"     ,&track.num_trackid     ,"num_trackid/I");
  tree_out->Branch("num_track_hits"  , track.num_track_hits  ,Form("num_track_hits[%d]/I"  ,MAX_NUM_TRACK));
  tree_out->Branch("track_hits_hitid", track.track_hits_hitid,Form("track_hits_hitid[%d][%d]/I",MAX_NUM_TRACK,MAX_NUM_TRACKHIT));
  tree_out->Branch("track_recon_id"  , track.track_recon_id  ,Form("track_recon_id[%d][%d]/I"  ,MAX_NUM_TRACK,MAX_NUM_TRACK)); 
  tree_out->Branch("track_start_z"   , track.track_start_z   ,Form("track_start_z[%d]/D"   ,MAX_NUM_TRACK)); 
  tree_out->Branch("track_stop_z"    , track.track_stop_z    ,Form("track_stop_z[%d]/D"    ,MAX_NUM_TRACK)); 
  tree_out->Branch("track_start_xy"  , track.track_start_xy  ,Form("track_start_xy[%d]/D"  ,MAX_NUM_TRACK)); 
  tree_out->Branch("track_stop_xy"   , track.track_stop_xy   ,Form("track_stop_xy[%d]/D"   ,MAX_NUM_TRACK)); 
  tree_out->Branch("track_start_pln" , track.track_start_pln ,Form("track_start_pln[%d]/I" ,MAX_NUM_TRACK)); 
  tree_out->Branch("track_stop_pln"  , track.track_stop_pln  ,Form("track_stop_pln[%d]/I"  ,MAX_NUM_TRACK)); 
  tree_out->Branch("track_start_ch"  , track.track_start_ch  ,Form("track_start_ch[%d]/I"  ,MAX_NUM_TRACK)); 
  tree_out->Branch("track_stop_ch"   , track.track_stop_ch   ,Form("track_stop_ch[%d]/I"   ,MAX_NUM_TRACK)); 
  tree_out->Branch("track_slope"     , track.track_slope     ,Form("track_slope[%d]/D"     ,MAX_NUM_TRACK));   
  tree_out->Branch("track_intercept" , track.track_intercept ,Form("track_intercept [%d]/D",MAX_NUM_TRACK)); 
  tree_out->Branch("track_pathlength", track.track_pathlength,Form("track_pathlength[%d]/D",MAX_NUM_TRACK)); 
  tree_out->Branch("track_total_pe"  , track.track_total_pe  ,Form("track_total_pe[%d]/D"  ,MAX_NUM_TRACK));  
  tree_out->Branch("track_mean_dedx" , track.track_mean_dedx ,Form("track_mean_dedx[%d]/D" ,MAX_NUM_TRACK));  
  tree_out->Branch("track_mean_time" , track.track_mean_time ,Form("track_mean_time[%d]/D" ,MAX_NUM_TRACK));  
  tree_out->Branch("track_view"      , track.track_view      ,Form("track_view[%d]/I"      ,MAX_NUM_TRACK));  
  tree_out->Branch("track_bcid"      , track.track_bcid      ,Form("track_bcid[%d]/I"      ,MAX_NUM_TRACK)); 
  tree_out->Branch("track_bcid_id"   , track.track_bcid_id   ,Form("track_bcid_id[%d]/I"   ,MAX_NUM_TRACK)); 
  tree_out->Branch("track_veto"      , track.track_veto      ,Form("track_veto[%d]/I"      ,MAX_NUM_TRACK)); 
  tree_out->Branch("track_sideescape", track.track_sideescape,Form("track_sideescape[%d]/I",MAX_NUM_TRACK)); 
  tree_out->Branch("track_cos_zen"   , track.track_cos_zen   ,Form("track_cos_zen[%d]/D"   ,MAX_NUM_TRACK));   
  tree_out->Branch("track_cos_azi"   , track.track_cos_azi   ,Form("track_cos_azi[%d]/D"   ,MAX_NUM_TRACK));   
  tree_out->Branch("track_len"       , track.track_len       ,Form("track_len[%d]/D"       ,MAX_NUM_TRACK));   
  tree_out->Branch("track_id"        , track.track_id        ,Form("track_id[%d]/I"        ,MAX_NUM_TRACK));   
  tree_out->Branch("track_likelihood", track.track_likelihood,Form("track_likelihood[%d]/I",MAX_NUM_TRACK));
  tree_out->Branch("track_ing_match" , track.track_ing_match ,Form("track_ing_match[%d]/O" ,MAX_NUM_TRACK));
  tree_out->Branch("track_ingtrk_id" , track.track_ingtrk_id ,Form("track_ingtrk_id[%d]/I" ,MAX_NUM_TRACK));
  tree_out->Branch("track_ing_diff"  , track.track_ing_diff  ,Form("track_ing_diff[%d]/D"  ,MAX_NUM_TRACK));
  tree_out->Branch("track_pm_match"  , track.track_pm_match  ,Form("track_pm_match[%d]/O"  ,MAX_NUM_TRACK));
  tree_out->Branch("track_pmtrk_id"  , track.track_pmtrk_id  ,Form("track_pmtrk_id[%d]/I"  ,MAX_NUM_TRACK));
  tree_out->Branch("track_pm_diff"   , track.track_pm_diff   ,Form("track_pm_diff[%d]/D"   ,MAX_NUM_TRACK));
  tree_out->Branch("track_ingcyc"    , track.track_ingcyc    ,Form("track_ingcyc[%d]/I"    ,MAX_NUM_TRACK));
  tree_out->Branch("track_ing_stopz" , track.track_ing_stopz ,Form("track_ing_stopz[%d]/I" ,MAX_NUM_TRACK));
  tree_out->Branch("track_ing_stopxy", track.track_ing_stopxy,Form("track_ing_stopxy[%d]/D",MAX_NUM_TRACK));

  return true;
}
//************************************************************************
//
bool wgGetTree::MakeTreeFile(string& str, IngRecon_t& IngRecon){
  if(!wgGetTree::foutput){
    wgGetTree::foutputname = str;
    wgGetTree::foutput  = new TFile(str.c_str(),"recreate");
    wgGetTree::tree_out = new TTree("tree","tree");
  }

  tree_out-> Branch("unixtime"    ,&IngRecon.unixtime   ,"unixtime/I"    );
  tree_out-> Branch("nd280nspill" ,&IngRecon.nd280nspill,"nd280nspill/I" );

  tree_out-> Branch("num_inghit"  ,IngRecon.num_inghit  ,Form("num_inghit[%d]/I"  ,NUM_CYC));
  tree_out-> Branch("num_pmhit"   ,IngRecon.num_pmhit   ,Form("num_pmhit[%d]/I"   ,NUM_CYC));
  tree_out-> Branch("num_ingrecon",IngRecon.num_ingrecon,Form("num_ingrecon[%d]/I",NUM_CYC));
  tree_out-> Branch("num_pmrecon" ,IngRecon.num_pmrecon ,Form("num_pmrecon[%d]/I" ,NUM_CYC));

  tree_out-> Branch("inghit_cyc"  , IngRecon.inghit_cyc ,Form("inghit_cyc[%d]/I" ,MAX_NUM_INGHIT));
  tree_out-> Branch("inghit_view" , IngRecon.inghit_view,Form("inghit_view[%d]/I",MAX_NUM_INGHIT));
  tree_out-> Branch("inghit_pln"  , IngRecon.inghit_pln ,Form("inghit_pln[%d]/I" ,MAX_NUM_INGHIT));
  tree_out-> Branch("inghit_ch"   , IngRecon.inghit_ch  ,Form("inghit_ch[%d]/I"  ,MAX_NUM_INGHIT));
  tree_out-> Branch("inghit_pe"   , IngRecon.inghit_pe  ,Form("inghit_pe[%d]/D"  ,MAX_NUM_INGHIT));
  tree_out-> Branch("inghit_time" , IngRecon.inghit_time,Form("inghit_time[%d]/D",MAX_NUM_INGHIT));
  tree_out-> Branch("pmhit_cyc"   , IngRecon.pmhit_cyc  ,Form("pmhit_cyc[%d]/I"  ,MAX_NUM_INGHIT));
  tree_out-> Branch("pmhit_view"  , IngRecon.pmhit_view ,Form("pmhit_view[%d]/I" ,MAX_NUM_INGHIT));
  tree_out-> Branch("pmhit_pln"   , IngRecon.pmhit_pln  ,Form("pmhit_pln[%d]/I"  ,MAX_NUM_INGHIT));
  tree_out-> Branch("pmhit_ch"    , IngRecon.pmhit_ch   ,Form("pmhit_ch[%d]/I"   ,MAX_NUM_INGHIT));
  tree_out-> Branch("pmhit_pe"    , IngRecon.pmhit_pe   ,Form("pmhit_pe[%d]/D"   ,MAX_NUM_INGHIT));
  tree_out-> Branch("pmhit_time"  , IngRecon.pmhit_time ,Form("pmhit_time[%d]/D" ,MAX_NUM_INGHIT));

  tree_out-> Branch("ingrecon_clstime"  , IngRecon.ingrecon_clstime ,Form("ingrecon_clstime[%d]/D" ,MAX_NUM_INGRECON));
  tree_out-> Branch("ingrecon_view"     , IngRecon.ingrecon_view    ,Form("ingrecon_view[%d]/I"    ,MAX_NUM_INGRECON));
  tree_out-> Branch("ingrecon_slope"    , IngRecon.ingrecon_slope   ,Form("ingrecon_slope[%d]/D"   ,MAX_NUM_INGRECON));
  tree_out-> Branch("ingrecon_angle"    , IngRecon.ingrecon_angle   ,Form("ingrecon_angle[%d]/D"   ,MAX_NUM_INGRECON));
  tree_out-> Branch("ingrecon_startz"   , IngRecon.ingrecon_startz  ,Form("ingrecon_startz[%d]/D"  ,MAX_NUM_INGRECON));
  tree_out-> Branch("ingrecon_startpln" , IngRecon.ingrecon_startpln,Form("ingrecon_startpln[%d]/I",MAX_NUM_INGRECON));
  tree_out-> Branch("ingrecon_startxy"  , IngRecon.ingrecon_startxy ,Form("ingrecon_startxy[%d]/D" ,MAX_NUM_INGRECON));
  tree_out-> Branch("ingrecon_endz"     , IngRecon.ingrecon_endz    ,Form("ingrecon_endz[%d]/D"    ,MAX_NUM_INGRECON));
  tree_out-> Branch("ingrecon_endpln"   , IngRecon.ingrecon_endpln  ,Form("ingrecon_endpln[%d]/I"  ,MAX_NUM_INGRECON));
  tree_out-> Branch("ingrecon_endxy"    , IngRecon.ingrecon_endxy   ,Form("ingrecon_endxy[%d]/D"   ,MAX_NUM_INGRECON));
  tree_out-> Branch("ingrecon_wg_match" , IngRecon.ingrecon_wg_match,Form("ingrecon_wg_match[%d]/O",MAX_NUM_INGRECON));
  tree_out-> Branch("ingrecon_wgtrk_id" , IngRecon.ingrecon_wgtrk_id,Form("ingrecon_wgtrk_id[%d]/I",MAX_NUM_INGRECON));

  tree_out-> Branch("pmrecon_clstime"   , IngRecon.pmrecon_clstime ,Form("pmrecon_clstime[%d]/D" ,MAX_NUM_INGRECON));
  tree_out-> Branch("pmrecon_view"      , IngRecon.pmrecon_view    ,Form("pmrecon_view[%d]/I"    ,MAX_NUM_INGRECON));
  tree_out-> Branch("pmrecon_slope"     , IngRecon.pmrecon_slope   ,Form("pmrecon_slope[%d]/D"   ,MAX_NUM_INGRECON));
  tree_out-> Branch("pmrecon_angle"     , IngRecon.pmrecon_angle   ,Form("pmrecon_angle[%d]/D"   ,MAX_NUM_INGRECON));
  tree_out-> Branch("pmrecon_startz"    , IngRecon.pmrecon_startz  ,Form("pmrecon_startz[%d]/D"  ,MAX_NUM_INGRECON));
  tree_out-> Branch("pmrecon_startpln"  , IngRecon.pmrecon_startpln,Form("pmrecon_startpln[%d]/I",MAX_NUM_INGRECON));
  tree_out-> Branch("pmrecon_startxy"   , IngRecon.pmrecon_startxy ,Form("pmrecon_startxy[%d]/D" ,MAX_NUM_INGRECON));
  tree_out-> Branch("pmrecon_endz"      , IngRecon.pmrecon_endz    ,Form("pmrecon_endz[%d]/D"    ,MAX_NUM_INGRECON));
  tree_out-> Branch("pmrecon_endpln"    , IngRecon.pmrecon_endpln  ,Form("pmrecon_endpln[%d]/I"  ,MAX_NUM_INGRECON));
  tree_out-> Branch("pmrecon_endxy"     , IngRecon.pmrecon_endxy   ,Form("pmrecon_endxy[%d]/D"   ,MAX_NUM_INGRECON));
  tree_out-> Branch("pmrecon_wg_match"  , IngRecon.pmrecon_wg_match,Form("pmrecon_wg_match[%d]/O",MAX_NUM_INGRECON));
  tree_out-> Branch("pmrecon_wgtrk_id"  , IngRecon.pmrecon_wgtrk_id,Form("pmrecon_wgtrk_id[%d]/I",MAX_NUM_INGRECON));
  return true;
}

//************************************************************************
//
bool wgGetTree::MakeTreeFile(string& str,BSD_t& Bsd){
  if(!wgGetTree::foutput){
    wgGetTree::foutputname = str;
    wgGetTree::foutput  = new TFile(str.c_str(),"recreate");
    wgGetTree::tree_out = new TTree("tree","tree");
  }

  tree_out->Branch("nurun"          , &Bsd.nurun           , "nurun/I"          );
  tree_out->Branch("midas_event"    , &Bsd.midas_event     , "midas_event/I"    );
  tree_out->Branch("mrrun"          , &Bsd.mrrun           , "mrrun/I"          );
  tree_out->Branch("mrshot"         , &Bsd.mrshot          , "mrshot/I"         );
  tree_out->Branch("spillnum"       , &Bsd.spillnum        , "spillnum/I"       );
  tree_out->Branch("trg_sec"        ,  Bsd.trg_sec         , "trg_sec[3]/I"     );
  tree_out->Branch("trg_nano"       ,  Bsd.trg_nano        , "trg_nano[3]/I"    );
  tree_out->Branch("gpsstat"        ,  Bsd.gpsstat         , "gpsstat[2]/I"     );
  tree_out->Branch("ct_np"          ,  Bsd.ct_np           , "ct_np[5][9]/D"    );
  tree_out->Branch("beam_time"      ,  Bsd.beam_time       , "beam_time[5][9]/D");
  tree_out->Branch("beam_flag"      ,  Bsd.beam_flag       , "beam_flag[5][9]/I");
  tree_out->Branch("hct"            ,  Bsd.hct             , "hct[3][5]/D"      );
  tree_out->Branch("tpos"           ,  Bsd.tpos            , "tpos[2]/D"        );
  tree_out->Branch("tdir"           ,  Bsd.tdir            , "tdir[2]/D"        );
  tree_out->Branch("tsize"          ,  Bsd.tsize           , "tsize[2]/D"       );
  tree_out->Branch("mumon"          ,  Bsd.mumon           , "mumon[12]/D"      );
  tree_out->Branch("otr"            ,  Bsd.otr             , "otr[13]/D"        );
  tree_out->Branch("good_gps_flag"  , &Bsd.good_gps_flag   , "good_gps_flag/I"  );
  tree_out->Branch("trigger_flag"   , &Bsd.trigger_flag    , "trigger_flag/I"   );
  tree_out->Branch("spill_flag"     , &Bsd.spill_flag      , "spill_flag/I"     );
  tree_out->Branch("good_spill_flag", &Bsd.good_spill_flag , "good_spill_flag/I");
  tree_out->Branch("target_eff"     ,  Bsd.target_eff      , "target_eff[3]/D"  );
  tree_out->Branch("run_type"       , &Bsd.run_type        , "run_type/I"       );
  tree_out->Branch("magset_id"      , &Bsd.magset_id       , "magset_id/I"      );
  tree_out->Branch("hctx"           ,  Bsd.hctx            , "hctx[2][5]/D"     );
  tree_out->Branch("htrans"         ,  Bsd.htrans          , "htrans[2]/D"      );
  tree_out->Branch("hps"            ,  Bsd.hps             , "hps[2]/D"         );

  return true;
}


//************************************************************************
bool wgGetTree::SetTreeFile(IngRecon_t& IngRecon,bool ingonly){
  if(ingonly){
    if(!wgGetTree::ingfinput){
      cout << "WARNING!! INGRID TFile is not open!" <<endl;
      return false;
    }
    ingtree-> SetBranchAddress("unixtime"    ,&IngRecon.unixtime   );
    ingtree-> SetBranchAddress("nd280nspill" ,&IngRecon.nd280nspill);

    ingtree-> SetBranchAddress("num_inghit"  ,IngRecon.num_inghit  );
    ingtree-> SetBranchAddress("num_pmhit"   ,IngRecon.num_pmhit   );
    ingtree-> SetBranchAddress("num_ingrecon",IngRecon.num_ingrecon);
    ingtree-> SetBranchAddress("num_pmrecon" ,IngRecon.num_pmrecon );

    ingtree-> SetBranchAddress("inghit_cyc"  , IngRecon.inghit_cyc);
    ingtree-> SetBranchAddress("inghit_view" , IngRecon.inghit_view);
    ingtree-> SetBranchAddress("inghit_pln"  , IngRecon.inghit_pln);
    ingtree-> SetBranchAddress("inghit_ch"   , IngRecon.inghit_ch);
    ingtree-> SetBranchAddress("inghit_pe"   , IngRecon.inghit_pe);
    ingtree-> SetBranchAddress("inghit_time" , IngRecon.inghit_time);
    ingtree-> SetBranchAddress("pmhit_cyc"   , IngRecon.pmhit_cyc);
    ingtree-> SetBranchAddress("pmhit_view"  , IngRecon.pmhit_view);
    ingtree-> SetBranchAddress("pmhit_pln"   , IngRecon.pmhit_pln);
    ingtree-> SetBranchAddress("pmhit_ch"    , IngRecon.pmhit_ch);
    ingtree-> SetBranchAddress("pmhit_pe"    , IngRecon.pmhit_pe);
    ingtree-> SetBranchAddress("pmhit_time"  , IngRecon.pmhit_time);

    ingtree-> SetBranchAddress("ingrecon_clstime"  , IngRecon.ingrecon_clstime);
    ingtree-> SetBranchAddress("ingrecon_view"     , IngRecon.ingrecon_view);
    ingtree-> SetBranchAddress("ingrecon_slope"    , IngRecon.ingrecon_slope);
    ingtree-> SetBranchAddress("ingrecon_angle"    , IngRecon.ingrecon_angle);
    ingtree-> SetBranchAddress("ingrecon_startz"   , IngRecon.ingrecon_startz);
    ingtree-> SetBranchAddress("ingrecon_startpln" , IngRecon.ingrecon_startpln);
    ingtree-> SetBranchAddress("ingrecon_startxy"  , IngRecon.ingrecon_startxy);
    ingtree-> SetBranchAddress("ingrecon_endz"     , IngRecon.ingrecon_endz);
    ingtree-> SetBranchAddress("ingrecon_endpln"   , IngRecon.ingrecon_endpln);
    ingtree-> SetBranchAddress("ingrecon_endxy"    , IngRecon.ingrecon_endxy);
    ingtree-> SetBranchAddress("ingrecon_wg_match" , IngRecon.ingrecon_wg_match);
    ingtree-> SetBranchAddress("ingrecon_wgtrk_id" , IngRecon.ingrecon_wgtrk_id);

    ingtree-> SetBranchAddress("pmrecon_clstime"   , IngRecon.pmrecon_clstime);
    ingtree-> SetBranchAddress("pmrecon_view"      , IngRecon.pmrecon_view);
    ingtree-> SetBranchAddress("pmrecon_slope"     , IngRecon.pmrecon_slope);
    ingtree-> SetBranchAddress("pmrecon_angle"     , IngRecon.pmrecon_angle);
    ingtree-> SetBranchAddress("pmrecon_startz"    , IngRecon.pmrecon_startz);
    ingtree-> SetBranchAddress("pmrecon_startpln"  , IngRecon.pmrecon_startpln);
    ingtree-> SetBranchAddress("pmrecon_startxy"   , IngRecon.pmrecon_startxy);
    ingtree-> SetBranchAddress("pmrecon_endz"      , IngRecon.pmrecon_endz);
    ingtree-> SetBranchAddress("pmrecon_endpln"    , IngRecon.pmrecon_endpln);
    ingtree-> SetBranchAddress("pmrecon_endxy"     , IngRecon.pmrecon_endxy);
    ingtree-> SetBranchAddress("pmrecon_wg_match"  , IngRecon.pmrecon_wg_match);
    ingtree-> SetBranchAddress("pmrecon_wgtrk_id"  , IngRecon.pmrecon_wgtrk_id);
  }
  else{
    if(!wgGetTree::finput){
      cout << "WARNING!! TFile is not open!" <<endl;
      return false;
    }
    tree-> SetBranchAddress("unixtime"    ,&IngRecon.unixtime   );
    tree-> SetBranchAddress("nd280nspill" ,&IngRecon.nd280nspill);

    tree-> SetBranchAddress("num_inghit"  ,IngRecon.num_inghit  );
    tree-> SetBranchAddress("num_pmhit"   ,IngRecon.num_pmhit   );
    tree-> SetBranchAddress("num_ingrecon",IngRecon.num_ingrecon);
    tree-> SetBranchAddress("num_pmrecon" ,IngRecon.num_pmrecon );

    tree-> SetBranchAddress("inghit_cyc"  , IngRecon.inghit_cyc);
    tree-> SetBranchAddress("inghit_view" , IngRecon.inghit_view);
    tree-> SetBranchAddress("inghit_pln"  , IngRecon.inghit_pln);
    tree-> SetBranchAddress("inghit_ch"   , IngRecon.inghit_ch);
    tree-> SetBranchAddress("inghit_pe"   , IngRecon.inghit_pe);
    tree-> SetBranchAddress("inghit_time" , IngRecon.inghit_time);
    tree-> SetBranchAddress("pmhit_cyc"   , IngRecon.pmhit_cyc);
    tree-> SetBranchAddress("pmhit_view"  , IngRecon.pmhit_view);
    tree-> SetBranchAddress("pmhit_pln"   , IngRecon.pmhit_pln);
    tree-> SetBranchAddress("pmhit_ch"    , IngRecon.pmhit_ch);
    tree-> SetBranchAddress("pmhit_pe"    , IngRecon.pmhit_pe);
    tree-> SetBranchAddress("pmhit_time"  , IngRecon.pmhit_time);

    tree-> SetBranchAddress("ingrecon_clstime"  , IngRecon.ingrecon_clstime);
    tree-> SetBranchAddress("ingrecon_view"     , IngRecon.ingrecon_view);
    tree-> SetBranchAddress("ingrecon_slope"    , IngRecon.ingrecon_slope);
    tree-> SetBranchAddress("ingrecon_angle"    , IngRecon.ingrecon_angle);
    tree-> SetBranchAddress("ingrecon_startz"   , IngRecon.ingrecon_startz);
    tree-> SetBranchAddress("ingrecon_startpln" , IngRecon.ingrecon_startpln);
    tree-> SetBranchAddress("ingrecon_startxy"  , IngRecon.ingrecon_startxy);
    tree-> SetBranchAddress("ingrecon_endz"     , IngRecon.ingrecon_endz);
    tree-> SetBranchAddress("ingrecon_endpln"   , IngRecon.ingrecon_endpln);
    tree-> SetBranchAddress("ingrecon_endxy"    , IngRecon.ingrecon_endxy);
    tree-> SetBranchAddress("ingrecon_wg_match" , IngRecon.ingrecon_wg_match);
    tree-> SetBranchAddress("ingrecon_wgtrk_id" , IngRecon.ingrecon_wgtrk_id);

    tree-> SetBranchAddress("pmrecon_clstime"   , IngRecon.pmrecon_clstime);
    tree-> SetBranchAddress("pmrecon_view"      , IngRecon.pmrecon_view);
    tree-> SetBranchAddress("pmrecon_slope"     , IngRecon.pmrecon_slope);
    tree-> SetBranchAddress("pmrecon_angle"     , IngRecon.pmrecon_angle);
    tree-> SetBranchAddress("pmrecon_startz"    , IngRecon.pmrecon_startz);
    tree-> SetBranchAddress("pmrecon_startpln"  , IngRecon.pmrecon_startpln);
    tree-> SetBranchAddress("pmrecon_startxy"   , IngRecon.pmrecon_startxy);
    tree-> SetBranchAddress("pmrecon_endz"      , IngRecon.pmrecon_endz);
    tree-> SetBranchAddress("pmrecon_endpln"    , IngRecon.pmrecon_endpln);
    tree-> SetBranchAddress("pmrecon_endxy"     , IngRecon.pmrecon_endxy);
    tree-> SetBranchAddress("pmrecon_wg_match"  , IngRecon.pmrecon_wg_match);
    tree-> SetBranchAddress("pmrecon_wgtrk_id"  , IngRecon.pmrecon_wgtrk_id);
  }

  return true;
}

//************************************************************************
bool wgGetTree::SetInfo(Info_t& inform){
  if(!wgGetTree::finput){
    cout << "WARNING!! TFile is not open!" <<endl;
    return false;
  }
  if(!wgGetTree::info){
    cout << "WARNING!! Info tree is not open!" <<endl;
    return false;
  }
  info->SetBranchAddress("start_time"    ,&inform.start_time   );
  info->SetBranchAddress("stop_time"     ,&inform.stop_time    );
  info->SetBranchAddress("nb_data_pkts"  ,&inform.nb_data_pkts );
  info->SetBranchAddress("nb_lost_pkts"  ,&inform.nb_lost_pkts );

  return true;
}

//************************************************************************
void wgGetTree::GetEntry(int i){
  wgGetTree::tree->GetEntry(i);
}


//***************************************
bool wgGetTree::OpenBsdFile(string filename)
{
  wgGetTree::bsd.infile = new TFile(filename.c_str(),"open");
  if(wgGetTree::bsd.infile->IsZombie()){
    cout << "No such a ROOT file: " << filename << endl;
    return false;
  }

  bsd.version_name = (TNamed*) wgGetTree::bsd.infile->Get("version");
  if(bsd.version_name==NULL){
    cerr << "[BSD] Error! Empty version name! " << endl;
    return false;
  }

  bsd.version = (string) bsd.version_name->GetTitle();
  if     (wgGetTree::bsd.version=="p06") wgGetTree::bsd.is_bsd = 0;
  else if(wgGetTree::bsd.version=="v01") wgGetTree::bsd.is_bsd = 1;
  else                           wgGetTree::bsd.is_bsd = -1;

  wgGetTree::bsd.bsd = (TTree*)wgGetTree::bsd.infile->Get("bsd");
  if(wgGetTree::bsd.bsd==NULL){
    cerr << "[BSD] Error! Empty bsd file! " << endl;
    return false;
  }
  wgGetTree::bsd.bsd->SetBranchAddress("nurun"          ,&(wgGetTree::bsd.nurun));
  wgGetTree::bsd.bsd->SetBranchAddress("midas_event"    ,&(wgGetTree::bsd.midas_event));
  wgGetTree::bsd.bsd->SetBranchAddress("mrrun"          ,&(wgGetTree::bsd.mrrun));
  wgGetTree::bsd.bsd->SetBranchAddress("mrshot"         ,&(wgGetTree::bsd.mrshot));
  wgGetTree::bsd.bsd->SetBranchAddress("spillnum"       ,&(wgGetTree::bsd.spillnum));
  wgGetTree::bsd.bsd->SetBranchAddress("trg_sec"        , (wgGetTree::bsd.trg_sec));
  wgGetTree::bsd.bsd->SetBranchAddress("trg_nano"       , (wgGetTree::bsd.trg_nano));
  wgGetTree::bsd.bsd->SetBranchAddress("gpsstat"        , (wgGetTree::bsd.gpsstat));
  wgGetTree::bsd.bsd->SetBranchAddress("ct_np"          , (wgGetTree::bsd.ct_np));
  wgGetTree::bsd.bsd->SetBranchAddress("beam_time"      , (wgGetTree::bsd.beam_time));
  wgGetTree::bsd.bsd->SetBranchAddress("beam_flag"      , (wgGetTree::bsd.beam_flag));
  wgGetTree::bsd.bsd->SetBranchAddress("hct"            , (wgGetTree::bsd.hct));
  wgGetTree::bsd.bsd->SetBranchAddress("tpos"           , (wgGetTree::bsd.tpos));
  wgGetTree::bsd.bsd->SetBranchAddress("tdir"           , (wgGetTree::bsd.tdir));
  wgGetTree::bsd.bsd->SetBranchAddress("tsize"          , (wgGetTree::bsd.tsize));
  wgGetTree::bsd.bsd->SetBranchAddress("mumon"          , (wgGetTree::bsd.mumon));
  wgGetTree::bsd.bsd->SetBranchAddress("otr"            , (wgGetTree::bsd.otr));
  wgGetTree::bsd.bsd->SetBranchAddress("good_gps_flag"  ,&(wgGetTree::bsd.good_gps_flag));
  wgGetTree::bsd.bsd->SetBranchAddress("trigger_flag"   ,&(wgGetTree::bsd.trigger_flag));
  wgGetTree::bsd.bsd->SetBranchAddress("spill_flag"     ,&(wgGetTree::bsd.spill_flag));
  wgGetTree::bsd.bsd->SetBranchAddress("good_spill_flag",&(wgGetTree::bsd.good_spill_flag));
  wgGetTree::bsd.bsd->SetBranchAddress("target_eff"     , (wgGetTree::bsd.target_eff));
  wgGetTree::bsd.bsd->SetBranchAddress("run_type"       ,&(wgGetTree::bsd.run_type));
  wgGetTree::bsd.bsd->SetBranchAddress("magset_id"      ,&(wgGetTree::bsd.magset_id));
  wgGetTree::bsd.bsd->SetBranchAddress("hctx"           , (wgGetTree::bsd.hctx));
  wgGetTree::bsd.bsd->SetBranchAddress("htrans"         , (wgGetTree::bsd.htrans));
  wgGetTree::bsd.bsd->SetBranchAddress("hps"            , (wgGetTree::bsd.hps));

  return true;
};

//***************************************
void wgGetTree::GetBSDEntry(int i)
{
  wgGetTree::bsd.bsd->GetEntry(i);
};


//************************************************************************
double wgGetTree::GetStartTime(){
  TH1F* h;
  h=(TH1F*)wgGetTree::finput->Get("start_time");
  double ret = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
  delete h;
  return ret;
}

//************************************************************************
double wgGetTree::GetStopTime(){
  TH1F* h;
  h=(TH1F*)wgGetTree::finput->Get("stop_time");
  double ret = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
  delete h;
  return ret;
}

//************************************************************************
double wgGetTree::GetDataPacket(){
  TH1F* h;
  h=(TH1F*)wgGetTree::finput->Get("nb_data_pkts");
  double ret = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
  delete h;
  return ret;
}

//************************************************************************
double wgGetTree::GetLostPacket(){
  TH1F* h;
  h=(TH1F*)wgGetTree::finput->Get("nb_lost_pkts");
  double ret = h->GetXaxis()->GetBinCenter(h->GetMaximumBin());
  delete h;
  return ret;
}

//************************************************************************
TH1F* wgGetTree::GetHist_StartTime(){
  TH1F* h;
  h=(TH1F*)wgGetTree::finput->Get("start_time");
  return h;
}

//************************************************************************
TH1F* wgGetTree::GetHist_StopTime(){
  TH1F* h;
  h=(TH1F*)wgGetTree::finput->Get("stop_time");
  return h;
}

//************************************************************************
TH1F* wgGetTree::GetHist_DataPacket(){
  TH1F* h;
  h=(TH1F*)wgGetTree::finput->Get("nb_data_pkts");
  return h;
}

//************************************************************************
TH1F* wgGetTree::GetHist_LostPacket(){
  TH1F* h;
  h=(TH1F*)wgGetTree::finput->Get("nb_lost_pkts");
  return h;
}
