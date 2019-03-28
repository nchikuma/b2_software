#!/usr/bin/python
# coding: utf-8


import sys,subprocess,os,time

def get_ingrid_list(filename):
  if not os.path.isfile(filename):
    print "No such a file:", filename
    return None
  
  ingridlist = []
  f = open(filename)
  line = f.readline()
  while line:
    ingridlist.append(line.strip().split())
    line = f.readline()
  return ingridlist
  f.close()


if __name__ == '__main__':

  
  #runlist1 = get_ingrid_list("ingrid_list.txt")
  runlist1 = get_ingrid_list("ingrid_list_t2krun9.txt")
  run1 = []
  for run in runlist1:
    tmp = [int(run[0]),int(run[1])]
    run1.append(tmp)
  print len(runlist1)
  print len(run1)

  runlist2 = get_ingrid_list("wagasci_runid2.txt")
  run2 = []
  for run in runlist2:
    tmp = [int(run[0]),int(run[1])]
    run2.append(tmp)
  print len(runlist2)
  print len(run2)
  
  data_dir  = "/home/t2k/nchikuma/b2_data/data_dst"
  data_dir2 = "/home/t2k/nchikuma/b2_data2/data_dst"
  #data_dir = "/home/t2k/nchikuma/b2_data/data_cosmic"
  numNoFile=0
  #for run in run2:
  for run in run1:
    runid  = int(run[0])
    srunid = int(run[1])
    tmp = [runid,srunid]
    decodefile = "%s/ingrid_%08d_%04d.root"      %(data_dir ,runid,srunid)
    calibfile  = "%s/ingrid_%08d_%04d_calib.root"%(data_dir ,runid,srunid)
    reconfile  = "%s/ingrid_%08d_%04d_recon.root"%(data_dir2,runid,srunid)
    anafile    = "%s/ingrid_%08d_%04d_anas1.root"%(data_dir2,runid,srunid)
    bsdfile    = "%s/ingrid_%08d_%04d_bsd1.root" %(data_dir2,runid,srunid)
    plotfile   = "%s/ingrid_%08d_%04d_plot.root" %(data_dir2,runid,srunid)
    if (not run in run2) and os.path.exists(bsdfile):
      cmd = "rm -f %s"%(bsdfile)
      subprocess.call(cmd,shell=True)
      numNoFile += 1
    #if not os.path.exists(plotfile):
    #  #cmd = "rm -f %s"%(plotfile)
    #  #subprocess.call(cmd,shell=True)
    #  print plotfile
    #  numNoFile += 1
    #  #print tmp
    #if not os.path.exists(decodefile):
    #  print "No Decode File: ", decodefile
    #  numNoFile += 1
    #  next
    #if not os.path.exists(calibfile):
    #  print "No Calib File: ", calibfile
    #  numNoFile += 1
    #  next
    #if not os.path.exists(reconfile):
    #  print "No Recon File: ", reconfile
    #  numNoFile += 1
    #  next
  if numNoFile==0:
    print "All exist"
  else:
    print "Non existing file :",numNoFile
