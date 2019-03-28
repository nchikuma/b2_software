// -*- C++ -*-
#ifndef BEAM_DATA_INCLUDED
#define BEAM_DATA_INCLUDED

#include <iostream>
#include <vector>

#include <TROOT.h>
#include <TChain.h>

//#define DEBUG_BEAMDATA

//
typedef struct 
{
    TString fname;
    unsigned int runnum;
    unsigned int start_time;
    unsigned int stop_time;
    unsigned int start_spillnum;
    unsigned int stop_spillnum;

    bool is_load;
    unsigned long  nevent;
    unsigned long  top_index;


} SPILL_DB;


typedef struct 
{
    //
    Int_t nurun;
    Int_t midas_event;

    Int_t mrrun;
    Int_t mrshot;

    Int_t spillnum;
    Int_t trg_sec[3];
    Int_t trg_nano[3];

    Int_t gpsstat[2];

    Double_t ct_np[5][9];  // [][0] is spill info, [][1..8] is bunch info
    Double_t beam_time[5][9];
    Int_t    beam_flag[5][9];

    Double_t hct[3][5]; // [][0] is total current, [][1..4] is individual ch
    Double_t hctx[2][5]; 
    Double_t htrans[2]; 
    Double_t hps[2]; 

    //
    Double_t tpos[2]; // position on target (x,y)
    Double_t tdir[2]; // direction on target (x',y')
    Double_t tsize[2];  // beam size on target (x,y)

    Double_t mumon[12]; // mumon fit results
    Double_t otr[13];   // OTR fit results

    Int_t  good_gps_flag;
    Int_t  trigger_flag;
    Int_t  spill_flag;
    Int_t  good_spill_flag;

    Double_t target_eff[3];  // new, since 2010/Feb/23
    
    Int_t  run_type;    // new, since 2010/Feb/23
    Int_t  magset_id;   // new, since 2010/Feb/23


} BEAM_DATA;


class BeamData
{
public:
    // instance
    static BeamData& instance();

    
    // method
    void init_db();
    void init_db( const TString& dbdir, const TString& version, const TString& datadir, int mrrun_begin=29, int mrrun_end=100 );
    void init_db( const char* dbdir, const char* version, const char* datadir, int mrrun_begin=29, int mrrun_end=100 );
    bool get_spill(int spillnum, int unixtime, int dt_cut=5);
    

    // extractor
    const BEAM_DATA&   data() const { return m_data; }

    //
    Int_t beam_runnumber() const { return m_data.nurun; }
    Int_t MR_runnumber()   const { return m_data.mrrun; }
    Int_t spill_number()   const { return m_data.spillnum; }

    
    Int_t trigger_time_sec( int igps )  const { return (0<=igps && igps<3) ? m_data.trg_sec[igps] : 0; }
    Int_t trigger_time_nano( int igps ) const { return (0<=igps && igps<3) ? m_data.trg_nano[igps]: 0; }
 
    Int_t gps_status( int igps ) const { return (0<=igps && igps<2) ? m_data.gpsstat[igps] : -1; }

    Double_t ct_np( int ch )               const { return (0<=ch && ch<5) ? m_data.ct_np[ch][0] : 0; }
    Double_t ct_np_bunch( int ch, int b )  const { return (ch<=ch && ch<5 && 1<=b && b<=8) ? m_data.ct_np[ch][b] : 0; }

    Double_t beam_timing( int ch )         const { return (0<=ch && ch<5) ? m_data.beam_time[ch][0] : 0; }
    Double_t bunch_timing( int ch, int b ) const { return (ch<=ch && ch<5 && 1<=b && b<=8) ? m_data.beam_time[ch][b] : 0; }

    Double_t beam_flag( int ch )           const { return (0<=ch && ch<5) ? m_data.beam_flag[ch][0] : 0; }
    Double_t bunch_flag( int ch, int b )   const { return (ch<=ch && ch<5 && 1<=b && b<=8) ? m_data.beam_flag[ch][b] : 0; }


    Double_t horn_current( int ihorn )                const { return (0<=ihorn && ihorn<3) ? m_data.hct[ihorn][0] : 0; }
    Double_t horn_busbar_current( int ihorn, int ch ) const { return (0<=ihorn && ihorn<3 && 1<=ch && ch<4) ? m_data.hct[ihorn][ch] : 0; }

    Double_t beam_pos_target( int ixy )   const { return (0<ixy && ixy<2) ? m_data.tpos[ixy] : -9999; }
    Double_t beam_dir_target( int ixy )   const { return (0<ixy && ixy<2) ? m_data.tdir[ixy] : -9999; }
    Double_t beam_size_target( int ixy )  const { return (0<ixy && ixy<2) ? m_data.tsize[ixy] : -9999; }

    Double_t mumon_si_totq() const { return m_data.mumon[0]; }
    Double_t mumon_si_peak() const { return m_data.mumon[1]; }
    Double_t mumon_si_x()    const { return m_data.mumon[2]; }
    Double_t mumon_si_wx()   const { return m_data.mumon[3]; }
    Double_t mumon_si_y()    const { return m_data.mumon[4]; }
    Double_t mumon_si_wy()   const { return m_data.mumon[5]; }

    Double_t mumon_ic_totq() const { return m_data.mumon[6]; }
    Double_t mumon_ic_peak() const { return m_data.mumon[7]; }
    Double_t mumon_ic_x()    const { return m_data.mumon[8]; }
    Double_t mumon_ic_wx()   const { return m_data.mumon[9]; }
    Double_t mumon_ic_y()    const { return m_data.mumon[10]; }
    Double_t mumon_ic_wy()   const { return m_data.mumon[11]; }


    Double_t otr_x()      const { return m_data.otr[0]; }
    Double_t otr_y()      const { return m_data.otr[1]; }
    Double_t otr_wx()     const { return m_data.otr[2]; }
    Double_t otr_wy()     const { return m_data.otr[3]; }
    Double_t otr_x_err()  const { return m_data.otr[4]; }
    Double_t otr_y_err()  const { return m_data.otr[5]; }
    Double_t otr_wx_err() const { return m_data.otr[8]; }
    Double_t otr_wy_err() const { return m_data.otr[9]; }
    Double_t otr_lyield() const { return m_data.otr[12]; }
    

    Int_t good_gps_flag() const { return m_data.good_gps_flag; }
    Int_t trigger_flag()  const { return m_data.trigger_flag; }
    Int_t spill_flag()    const { return m_data.spill_flag; }
    Int_t good_spill()    const { return m_data.good_spill_flag; }

    Double_t target_eff(int index ) const { return m_data.target_eff[index]; }

    Int_t run_type()  const { return m_data.run_type; }
    Int_t magset_id() const { return m_data.magset_id; }

    
protected:
    // constructor
    BeamData();
    
    // destructor
    ~BeamData();

private:
    // 
    static BeamData* s_instance;


    //
    TString  m_dbdir;
    TString  m_basedir;
    TString  m_version;

    //
    TChain*  m_chain;

    int      m_mrrun_begin;
    int      m_mrrun_end;

    
    std::vector<SPILL_DB> m_db;
    int      m_db_nelem;

    int      m_current_spill_nd280;
    int      m_current_unixtime;

    bool     m_init;

    long     m_top_index;
    long     m_start_index;
    long     m_last_index;
    long     m_nentries;

    int      m_top_unixtime;
    int      m_last_unixtime;


    //
    // data
    // Declaration of leaf types
    BEAM_DATA      m_data;
    
    // List of branches
    TBranch        *b_nurun;   //!
    TBranch        *b_midas_event;   //!
    TBranch        *b_mrrun;   //!
    TBranch        *b_mrshot;   //!
    TBranch        *b_spillnum;   //!
    TBranch        *b_trig_sec;   //!
    TBranch        *b_trig_nano;   //!
    TBranch        *b_gpsstat;   //!
    TBranch        *b_ct_np;   //!
    TBranch        *b_beam_time;   //!
    TBranch        *b_beam_flag;   //!
    TBranch        *b_hct;   //!
    TBranch        *b_hctx;   //!
    TBranch        *b_htrans;   //!
    TBranch        *b_hps;   //!
    TBranch        *b_tpos;   //!
    TBranch        *b_tdir;   //!
    TBranch        *b_tsize;   //!
    TBranch        *b_mumon;   //!
    TBranch        *b_otr;   //!
    TBranch        *b_good_gps_flag;   //!
    TBranch        *b_trigger_flag;   //!
    TBranch        *b_spill_flag;   //!
    TBranch        *b_good_spill_flag;   //!
    TBranch        *b_target_eff;   //!
    TBranch        *b_run_type;   //!
    TBranch        *b_magset_id;   //!


    
    // method
    void add_db(  const TString& fname, int runnum, int start_time, int stop_time, int start_spillnum, int stop_spillnum );
    void set_branch();

};


#endif // BEAM_DATA_INCLUDED
