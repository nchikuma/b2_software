#!/usr/bin/python
# coding: utf-8


import sys,subprocess,os,time

MAXPRINT=10
NUM_JOBS_SUMBIT=10
NUM_MC_FILE=1000
QUEUE = "s"


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


jobname = [
  "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/shell/decode_INGRID.sh",
  "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/shell/decode_WAGASCI.sh",
  "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/shell/analysis_B2beam.sh",
  "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/shell/analysis_B2cosmic_ing.sh",
  "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/shell/analysis_B2cosmic_wg.sh",
  "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/shell/hiteff.sh",
  "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/shell/hiteff_wgcosmic.sh",
  "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/shell/trackeff.sh",
  "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/shell/hiteff_mc.sh",
  "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/shell/trackeff_mc.sh",
  "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/shell/make_Calib.sh",
  "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/shell/mc_neut.sh",
  "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/shell/trackeff_neut.sh",
  "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/shell/mc_cosmic.sh",
  "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/shell/mc_sandmu.sh",
  "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/shell/pe_check.sh",
  "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/shell/pe_check_wgcosmic.sh",
  "/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software/shell/pe_check_mc.sh"
  ]
listfilename = [
  "ingrid_list.txt",
  "wagasci_list.txt",
  "ingrid_list.txt",
  "ingrid_list.txt",
  "wagasci_list.txt",
  "ingrid_list.txt",
  "wagasci_list.txt",
  "ingrid_list.txt",
  "",
  "",
  "ingrid_list.txt",
  "neut_list.txt",
  "neut_list.txt",
  "",
  "",
  "ingrid_list.txt",
  "wagasci_list.txt",
  ""
  ]
submitname = [
  "submit_ing",
  "submit_wg",
  "submit_ana",
  "submit_cosmic_ing",
  "submit_cosmic_wg",
  "submit_hiteff",
  "submit_hiteff_wg",
  "submit_trkeff",
  "submit_hiteff_mc",
  "submit_trkeff_mc",
  "submit_makecalib",
  "submit_neut",
  "submit_trkeff_neut",
  "submit_mc_cosmic",
  "submit_mc_sandmu",
  "submit_pe",
  "submit_pe_wg",
  "submit_pe_mc"
  ]

def PrintJoblist():

  numjob  = len(jobname)
  numlist = len(listfilename)
  numname = len(submitname)

  nummax = numjob
  if nummax<numlist: nummax = numlist
  if nummax<numname: nummax = numname
  print ""
  for i in range(nummax):
    print "=== Job ID :",i,"==="
    if i<numjob : print "  jobname     :",jobname     [i]
    if i<numlist: print "  listfilename:",listfilename[i]
    if i<numname: print "  submitname  :",submitname  [i]
  print ""

  print "./job_submit.py <jobid>"
  print "./job_submit.py <jobid> f  : To force submission."



if __name__ == '__main__':

  if len(jobname)!=len(listfilename) or len(jobname)!=len(submitname):
    print "Correct the job list."
    PrintJoblist()
    exit(0)

  if (len(sys.argv)!=2 and len(sys.argv)!=3) or not sys.argv[1].isdigit():
    print "Put an argument for selecting a submitted job."
    PrintJoblist()
    exit(0)
  
  force_mode = False
  if len(sys.argv)==3:
    if sys.argv[2]=="f":
      force_mode = True
    else:
      PrintJoblist()
      print "ERROR : This mode is not available." 
      exit(0)
  

  jobid = int(sys.argv[1])
  if jobid>=len(jobname):
    PrintJoblist()
    print "ERROR : This job ID is not listed in the script."
    exit(0)

  jobname      = jobname     [jobid]
  listfilename = listfilename[jobid]
  submitname   = submitname  [jobid]


  runlist = []
  if listfilename=="":
    for i in range(NUM_MC_FILE):
      tmp = [0,i+1]
      runlist.append(tmp)
  else:
    runlist = get_ingrid_list(listfilename)

  print "================"
  print "This script might work for submitting jobs with a type of;"
  print "$<jobname> <runid> <acqid>"
  print "JOB:",jobname
  print "runlist:",listfilename
  print "queue:",QUEUE
  print "=== RUN LIST ==="
  count = 0
  num_runs=len(runlist)
  for run in runlist:
    if count<MAXPRINT:
      print run
    elif count==MAXPRINT:
      print "More..."
      break
    count+=1
  print "In total,",num_runs," runs are listed."
  print "NUM_JOBS_SUMBIT=",NUM_JOBS_SUMBIT

  ans = ""
  if not force_mode:
    print "Are you sure to submit jobs? Yes(yes/y/Y) to submit."
    ans = sys.stdin.readline().strip()
  if ans in ["Yes","yes","y","Y"] or force_mode:
    num_submit = 0
    count = 0
    for run in runlist:
      runid  = int(run[0])
      srunid = int(run[1])
      #if jobid==11 or jobid==13 or jobid==14:
      #  mc_outputdir = "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/mc_neut/"
      #  mc_outputname= "ingrid_%08d_%04d_calib.root"%(runid,srunid)
      #  mc_name = mc_outputdir + mc_outputname
      #  if os.path.isfile(mc_name):
      #    print "The MC file already exists : ", mc_name
      #    continue
      if count==0:
        jobfile = "./submit/%s_%04d.sh"%(submitname,num_submit)
        f = open(jobfile,'w')
        f.write("#!/bin/sh\n")
        f.write("\n")
      cmd = "%s %d %d\n"%(jobname,runid,srunid)
      f.write(cmd)
      count+=1
      if count>=NUM_JOBS_SUMBIT or runlist.index(run)==num_runs-1:
        count = 0
        num_submit+=1
        f.close()
        #time.sleep(1)
        os.chmod(jobfile,0755)
        cmd = "bsub -q %s %s"%(QUEUE,jobfile)
        print cmd
        subprocess.call(cmd,shell=True)
  else:
    print "Nothing has been done."
