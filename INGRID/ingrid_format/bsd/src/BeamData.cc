// -*- C++ -*-

#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include "BeamData/BeamData.h"


BeamData* BeamData::s_instance = 0;

static BEAM_DATA s_dummy_data;  // dummy


BeamData& BeamData::instance()
{
  if( s_instance == 0 ) {
    s_instance = new BeamData;
    if( s_instance == 0 ) { // failed
      std::cerr << "BeamData::instance() " 
        << "cannot create BeamData object"
        << std::endl;
      perror( "new" );
      abort();
    }
  }
  return( *s_instance );
}


// constructor
  BeamData::BeamData()
: m_mrrun_begin(27),
  m_mrrun_end(100),
  m_db_nelem(0),
  m_current_spill_nd280(0),
  m_current_unixtime(0),
  m_init(false),
  m_top_index(0),
  m_start_index(0),
  m_last_index(0),
  m_nentries(0),
  m_top_unixtime(0),
  m_last_unixtime(0)
{
  m_chain = new TChain("bsd");
  m_data = s_dummy_data; // clear
  m_db.clear();

  m_dbdir  .Form("/home/t2k/nchikuma/b2_data/bsd/spilldb");
  m_version.Form("v01");
  m_basedir.Form("/home/t2k/nchikuma/b2_data/bsd");
}


// destractor
BeamData::~BeamData()
{
}


void BeamData::init_db( const char* dbdir, const char* version, const char* datadir, int mrrun_begin, int mrrun_end )
{
  init_db( TString(dbdir), TString(version), TString(datadir), mrrun_begin, mrrun_end );
}


void BeamData::init_db( const TString& dbdir, const TString& version, const TString& datadir, int mrrun_begin, int mrrun_end )
{
  m_dbdir       = dbdir;
  m_version     = version;
  m_basedir     = datadir;
  m_mrrun_begin = mrrun_begin;
  m_mrrun_end   = mrrun_end;

  init_db();
}

void BeamData::init_db()
{
  std::cout << "BeamData::init_db() : "  << m_dbdir << " " << m_basedir << std::endl;

  for(int r=m_mrrun_begin; r<=m_mrrun_end; r++) {
    TString dbfname;
    dbfname.Form("/run%07d.db", r);
    dbfname = m_dbdir + dbfname;

    if( access(dbfname.Data(),R_OK) == 0 ){
      std::cout << "reading " << dbfname.Data() << std::endl;

      std::ifstream dbfile( dbfname.Data() );
      if( dbfile.good() ) {
        while(1) {
          std::string fname;
          unsigned int  runnum;
          unsigned int  start_unixtime;
          unsigned int  stop_unixtime;
          unsigned int  start_spillnum;
          unsigned int  stop_spillnum;

          dbfile >> fname >> runnum >> start_unixtime >> stop_unixtime >> start_spillnum >> stop_spillnum;
          if( !dbfile.good() )
            break;
          TString fpathname;
          fpathname.Form("%s/%s", m_basedir.Data(), fname.c_str());

          add_db( fpathname, runnum, start_unixtime, stop_unixtime, start_spillnum, stop_spillnum);
        }
      }
    }
  }

  if( m_db.empty() ) {
    std::cerr << "No database for spillnum !! Please check database files in " << m_dbdir << ", version=" << m_version << std::endl;
  }

  for(unsigned int i=0; i<m_db.size(); i++) {
  }

}

void BeamData::add_db( const TString& fname, int runnum, int start_time, int stop_time, int start_spillnum, int stop_spillnum )
{
  SPILL_DB db;

  db.fname          = fname;
  db.runnum         = runnum;
  db.start_time     = start_time;
  db.stop_time      = stop_time;
  db.start_spillnum = start_spillnum;
  db.stop_spillnum  = stop_spillnum;

  db.is_load        = false;

  db.nevent = stop_spillnum - start_spillnum;
  if( db.nevent < 0 ) {
    db.nevent  = 0xFFFFFFFF - start_spillnum + stop_spillnum;
  }

  db.top_index = 0;
  m_db.push_back( db );


  m_db_nelem++;
}




bool BeamData::get_spill( int spill_nd280, int unixtime, int dt_cut )
{
  if( m_db.empty() ) return false;

  // check same event or not
  if( (spill_nd280 == m_current_spill_nd280) && (unixtime == m_current_unixtime) ) {
    return true;
  }

  // check same run_file or not
  if( m_top_unixtime+10 <= unixtime && unixtime <= m_last_unixtime-10 ) { 
    if( (spill_nd280 < m_current_spill_nd280) && (unixtime < m_current_unixtime) ) {
      // this may be previous event
      m_start_index = m_top_index;
    }
  }
  else {
    m_top_index     = 0;
    m_top_unixtime  = 0;
    m_last_index    = 0;
    m_last_unixtime = 0;

    int nfile = 0;
    for(int i=0; i<m_db_nelem; i++) {
      if( m_db[i].start_time-10 <= unixtime && unixtime <= m_db[i].stop_time+10 ) { // select by wide range 
        if( ! m_db[i].is_load ) {
          m_chain->Add( m_db[i].fname.Data() );

          m_db[i].is_load   = true;
          m_db[i].top_index = m_nentries;

          m_nentries = m_chain->GetEntries();
        }

        if( nfile == 0 ) {
          m_top_index        = m_db[i].top_index;
          m_top_unixtime     = m_db[i].start_time;
        }

        m_last_index       = m_db[i].top_index + m_db[i].nevent;
        m_last_unixtime    = m_db[i].stop_time;

        nfile++;
      }
    }
    m_start_index = m_top_index;
  }
  set_branch();

  bool is_good = false;
  for(int i=m_start_index; i<=m_last_index; i++) {
    int ret = m_chain->GetEntry(i);
    if( ret < 0 ) {
      std::cerr << "Error in reading file" << std::endl;
    }
    int spill_nd280_this = 1 + (m_data.spillnum)&0xFFFF;
    int dt               = TMath::Abs( m_data.trg_sec[0] - unixtime );

    if( (spill_nd280_this == spill_nd280) && (dt < dt_cut) ) {
      is_good = true;
      m_current_spill_nd280 = spill_nd280;
      m_current_unixtime    = unixtime;

      m_start_index = i;

      break;
    }
  }
  if( ! is_good ) {
    m_data = s_dummy_data;
  }
  return is_good;
}

void BeamData::set_branch()
{
  if( ! m_init ) {
    m_chain->SetBranchAddress("nurun",            &(m_data.nurun),            &b_nurun);
    m_chain->SetBranchAddress("midas_event",      &(m_data.midas_event),      &b_midas_event);
    m_chain->SetBranchAddress("mrrun",            &(m_data.mrrun),            &b_mrrun);
    m_chain->SetBranchAddress("mrshot",           &(m_data.mrshot),           &b_mrshot);
    m_chain->SetBranchAddress("spillnum",         &(m_data.spillnum),         &b_spillnum);
    m_chain->SetBranchAddress("trg_sec",          &(m_data.trg_sec),          &b_trig_sec);
    m_chain->SetBranchAddress("trg_nano",         &(m_data.trg_nano),         &b_trig_nano);
    m_chain->SetBranchAddress("gpsstat",          &(m_data.gpsstat),          &b_gpsstat);
    m_chain->SetBranchAddress("ct_np",            &(m_data.ct_np),            &b_ct_np);
    m_chain->SetBranchAddress("beam_time",        &(m_data.beam_time),        &b_beam_time);
    m_chain->SetBranchAddress("beam_flag",        &(m_data.beam_flag),        &b_beam_flag);
    m_chain->SetBranchAddress("hct",              &(m_data.hct),              &b_hct);
    m_chain->SetBranchAddress("hctx",             &(m_data.hctx),             &b_hctx);
    m_chain->SetBranchAddress("htrans",           &(m_data.htrans),           &b_htrans);
    m_chain->SetBranchAddress("hps",              &(m_data.hps),              &b_hps);
    m_chain->SetBranchAddress("tpos",             &(m_data.tpos),             &b_tpos);
    m_chain->SetBranchAddress("tdir",             &(m_data.tdir),             &b_tdir);
    m_chain->SetBranchAddress("tsize",            &(m_data.tsize),            &b_tsize);
    m_chain->SetBranchAddress("mumon",            &(m_data.mumon),            &b_mumon);
    m_chain->SetBranchAddress("otr",              &(m_data.otr),              &b_otr);
    m_chain->SetBranchAddress("good_gps_flag",    &(m_data.good_gps_flag),    &b_good_gps_flag);
    m_chain->SetBranchAddress("trigger_flag",     &(m_data.trigger_flag),     &b_trigger_flag);
    m_chain->SetBranchAddress("spill_flag",       &(m_data.spill_flag),       &b_spill_flag);
    m_chain->SetBranchAddress("good_spill_flag",  &(m_data.good_spill_flag),  &b_good_spill_flag);
    m_chain->SetBranchAddress("target_eff",       &(m_data.target_eff),       &b_target_eff);
    m_chain->SetBranchAddress("run_type",         &(m_data.run_type),         &b_run_type);
    m_chain->SetBranchAddress("magset_id",        &(m_data.magset_id),        &b_magset_id);


    m_init = true;
  }

}
