#!/usr/bin/python
# -*- coding: utf-8 -*-

############################################################################################
#  Here are all of the fixed parameters used in the slow monitor system written in python. #
#  These must be identical to the other setting file "Setup.h", that is used to compile   #
#  C++ programs and to run SH scripts.                                                     #
#                                                                                          #
#  2017/09/23                                                                              #
#  Naruhiro Chikuma                                                                        #
#  The University of Tokyo                                                                 #
#                                                                                          #
############################################################################################

import os, subprocess,datetime

class Setup:
  def __init__(self):


    ##################
    # Directories
    ##################
    self.MAIN_DATA_DIR   = "/home/t2k/nchikuma/b2_data"
    self.MAIN_DATA_DIR2  = "/home/t2k/nchikuma/b2_data2"
    self.MAIN_SOFT_DIR   = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/wagasci_software"
    self.ANALYSIS_DIR    = "{0}/Analysis"        .format(self.MAIN_SOFT_DIR)

    ########################
    # Data copy from DAQ PC
    ########################
    self.RUNNAME       = "run"
    self.BACKUPDATA_DIR = "{0}/daqdata".format(self.MAIN_DATA_DIR)
    self.ID_DIR         = "{0}/runid"  .format(self.MAIN_DATA_DIR2)

    ########################
    # Auto process
    ########################
    self.DECODE_DATA_DIR   = "{0}/decode"     .format(self.MAIN_DATA_DIR)
    self.HIST_DATA_DIR     = "{0}/hist"       .format(self.MAIN_DATA_DIR2)
    self.RECON_DATA_DIR    = "{0}/recon"      .format(self.MAIN_DATA_DIR2)
    self.CALIB_DIR         = "{0}/calibration".format(self.MAIN_DATA_DIR)
    self.DQ_MERGE_DIR      = "{0}/dq_merge"   .format(self.MAIN_DATA_DIR2)
    self.MERGE_DIR         = "{0}/merge"      .format(self.MAIN_DATA_DIR2)
    self.BSD_DIR           = "{0}/bsd"        .format(self.MAIN_DATA_DIR)
    self.XML_DATA_DIR      = "{0}/xmlfile"    .format(self.MAIN_DATA_DIR2)
    self.DQ_HISTORY_DIR    = "{0}/dq_history" .format(self.MAIN_DATA_DIR2)
    self.CALIB_ID_FILE     = "{0}/calib_id.txt"               .format(self.ID_DIR)
    self.AUTO_RUNID_LIST   = "{0}/auto_process_runid_list.txt".format(self.ID_DIR)
    self.PROCESS_DECODE    = "{0}/bin/Decoder"         .format(self.ANALYSIS_DIR)
    self.PROCESS_HIST      = "{0}/bin/wgMakeHist"      .format(self.ANALYSIS_DIR)
    self.PROCESS_RECON     = "{0}/bin/wgRecon"         .format(self.ANALYSIS_DIR)
    self.PROCESS_MERGEING  = "{0}/bin/wgMerge_inglib"  .format(self.ANALYSIS_DIR)
    self.PROCESS_ANAHIST   = "{0}/bin/wgAnaHist"       .format(self.ANALYSIS_DIR)
    self.PROCESS_ANAHISTSUM= "{0}/bin/wgAnaHistSummary".format(self.ANALYSIS_DIR)
    self.PROCESS_DQCHECK   = "{0}/bin/wgDQCheck"       .format(self.ANALYSIS_DIR)
    self.PROCESS_DQHISTORY = "{0}/bin/wgDQHistory"     .format(self.ANALYSIS_DIR)
    self.BSD_SPILLCHECK    = "{0}/bin/wgBsdSpillCheck" .format(self.ANALYSIS_DIR)
    self.SPILL_CHECK       = "{0}/bin/wgSpillCheck"    .format(self.ANALYSIS_DIR)
    self.SPILL_EFF         = "{0}/bin/wgSpillEff"      .format(self.ANALYSIS_DIR)

  def insertstr(self,s="",pos=-1,x=""):
    return x.join([s[:pos],s[pos:]])

  def str2ascii(self,string=""):
    tmp = binascii.hexlify(string)
    for i in xrange(len(tmp)-2, 0, -2):
      tmp = self.insertstr(tmp,i,'\\x')
    return str('\\x'+tmp)

  def set_process_state(self,runid=-1,acqid=-1,ini=-1,fin=-1):
    cmd = "sed -i 's/%05d %03d %d/%05d %03d %d/' %s"%(
        runid,acqid,ini,
        runid,acqid,fin,
        self.AUTO_RUNID_LIST)
    subprocess.call(cmd,shell=True)

  def get_process_state(self,runid=-1,acqid=-1):
    cmd = "cat %s | grep '%05d %03d'"%(self.AUTO_RUNID_LIST,runid,acqid)
    res = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    result = str(res.communicate()[0].replace("\n","")).split()
    if len(result)!=4 or (not result[2].isdigit()):
      return -1
    else:
      return int(result[2])

  def get_calibname(self,runid=-1,acqid=-1):
    cmd = "cat %s | grep '%05d %03d'"%(self.AUTO_RUNID_LIST,runid,acqid)
    res = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    result = str(res.communicate()[0].replace("\n","")).split()
    if len(result)!=4:
      return None
    else:
      return result[3]


  def __enter__(self):
    return self

  def __exit__(self,exc_type,exc_value,traceback):
    return True
