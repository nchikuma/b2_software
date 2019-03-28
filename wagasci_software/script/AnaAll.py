#!/usr/bin/python
# -*- coding: utf-8 -*-


import sys, os, time, subprocess, datetime, math, glob
import threading
import Setup

############################################################################################

def Decode(rawfile1,rawfile2,calibfile,mode):
  setup = Setup.Setup()
  if os.path.exists(rawfile1) and os.path.exists(rawfile2):

    cmd = "%s -rf %s -i %s"%(setup.PROCESS_DECODE,rawfile1,calibfile)
    print cmd
    sys.stdout.flush()
    subprocess.call(cmd,shell=True)

    cmd = "%s -rf %s -i %s"%(setup.PROCESS_DECODE,rawfile2,calibfile)
    print cmd
    sys.stdout.flush()
    subprocess.call(cmd,shell=True)

    return True
  else:
    msg ="No such raw files. %s %s"%(rawfile1, rawfile2)
    print msg
    sys.stdout.flush()

    return False

############################################################################################

def MergeINGlib(decodefile1,decodefile2,rundi,acqid,mode="default"):
  setup = Setup.Setup()
  if os.path.exists(decodefile1) and os.path.exists(decodefile2):

    cmd = "%s -r %d -s %d"%(setup.PROCESS_MERGEING,runid,acqid)
    print cmd
    sys.stdout.flush()

    subprocess.call(cmd,shell=True)

    cmd = "%s -r %d -s %d -c"%(setup.PROCESS_MERGEING,runid,acqid)
    print cmd
    sys.stdout.flush()
    subprocess.call(cmd,shell=True)

    return True
  else:
    msg =  "No decoded files. %s %s"%(decodefile1, decodefile2)
    print msg
    sys.stdout.flush()

    return False


############################################################################################

def analysis_all(runid=-1,acqid=-1,mode="default"):
  setup = Setup.Setup()

  calibname = setup.get_calibname(runid,acqid)
  if calibname==None:
    msg = "No such runid/acqid, %08d, %03d"%(runid, acqid)
    print msg
    sys.stdout.flush()
    return

  runidname = "%05d"%(runid)
  acqidname = "%03d"%(acqid)
  
  rawfile1 = glob.glob("%s/%s_%s*/%s_%s_%s*_dif_1_1_1.raw"%(
        setup.BACKUPDATA_DIR,
        setup.RUNNAME,runidname,
        setup.RUNNAME,runidname,acqidname))
  rawfile2 = glob.glob("%s/%s_%s*/%s_%s_%s*_dif_1_1_2.raw"%(
        setup.BACKUPDATA_DIR,
        setup.RUNNAME,runidname,
        setup.RUNNAME,runidname,acqidname))

  if (not len(rawfile1)==1) or (not len(rawfile2)==1): 
    msg =  "There is no such files, otherwise too many files : %s ; %s"%(rawfile1 , rawfile2)
    print msg
    sys.stdout.flush()
    return

  rawfile1 = rawfile1[0]
  rawfile2 = rawfile2[0]
  suffix = rawfile1.replace("%s/%s_%s"%(
        setup.BACKUPDATA_DIR,
        setup.RUNNAME,runidname),"").split("/")[0]
  
  rawfile1 = "%s/%s_%s%s/%s_%s_%s%s_dif_1_1_1.raw"%(
        setup.BACKUPDATA_DIR,
        setup.RUNNAME,runidname,suffix,
        setup.RUNNAME,runidname,acqidname,suffix)
  rawfile2 = "%s/%s_%s%s/%s_%s_%s%s_dif_1_1_2.raw"%(
        setup.BACKUPDATA_DIR,
        setup.RUNNAME,runidname,suffix,
        setup.RUNNAME,runidname,acqidname,suffix)
  configfile = "%s/%s_%s%s/%s_%s_%s%s.xml"%(
        setup.BACKUPDATA_DIR,
        setup.RUNNAME,runidname,suffix,
        setup.RUNNAME,runidname,acqidname,suffix)
  calibfile = "%s/%s/calib_result.xml"%(setup.CALIB_DIR,calibname)
  decodefile1 = "%s/%s_%s_%s%s_dif_1_1_1_tree.root"%(
        setup.DECODE_DATA_DIR,setup.RUNNAME,runidname,acqidname,suffix)
  decodefile2 = "%s/%s_%s_%s%s_dif_1_1_2_tree.root"%(
        setup.DECODE_DATA_DIR,setup.RUNNAME,runidname,acqidname,suffix)
  histfile1 = "%s/%s_%s_%s%s_dif_1_1_1_hist.root"%(
        setup.HIST_DATA_DIR,setup.RUNNAME,runidname,acqidname,suffix)
  histfile2 = "%s/%s_%s_%s%s_dif_1_1_2_hist.root"%(
        setup.HIST_DATA_DIR,setup.RUNNAME,runidname,acqidname,suffix)
  xmldir1 = "%s/%s_%s_%s%s_dif_1_1_1"%(
        setup.XML_DATA_DIR,setup.RUNNAME,runidname,acqidname,suffix)
  xmldir2 = "%s/%s_%s_%s%s_dif_1_1_2"%(
        setup.XML_DATA_DIR,setup.RUNNAME,runidname,acqidname,suffix)

  decodeacqname = "%s/%s_%s_%s%s"%(
        setup.DECODE_DATA_DIR,setup.RUNNAME,runidname,acqidname,suffix)
  acqname = "%s_%s_%s%s"%(
        setup.RUNNAME,runidname,acqidname,suffix)
  reconfile = "%s/%s_%s_%s%s_recon.root"%(
        setup.RECON_DATA_DIR,setup.RUNNAME,runidname,acqidname,suffix)


  ###################################################
  # Decode

  if Decode(rawfile1,rawfile2,calibfile,mode):
    msg = "OK. Decode process is done."
    print msg
    sys.stdout.flush()
  else:
    msg = "ERROR: Decode process is failed."
    print msg
    sys.stdout.flush()
        
  ###################################################
  # Merge INGRID library

  if MergeINGlib(decodefile1,decodefile2,runid,acqid,mode):
    msg = "OK: MergeINGlib process is done."
    print msg
    sys.stdout.flush()
  else:
    msg = "ERROR: MergeINGlib process is failed."
    print msg
    sys.stdout.flush()


def PrintUsage():
  print "Usage1: {0} <runid> <acqid>".format(sys.argv[0])

############################################################################################

if __name__ == '__main__':

  setup = Setup.Setup()

  if len(sys.argv)<3:
    PrintUsage()
    exit(0)

  if (not sys.argv[1].isdigit()) or (not sys.argv[2].isdigit()):
    PrintUsage()
    exit(0)
  runid   = int(sys.argv[1])
  acqid   = int(sys.argv[2])
  mode    = "default"
  analysis_all(runid,acqid,mode)
