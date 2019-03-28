#!/bin/sh

MCSIM=0
TWODRECON=0
THRDRECON=1
PLOT_ALL=0
PLOT_MODE=0
PLOT_REWEI=0
DETRES=0
DETRES2=0
TRKEFF=0
PECHECK=0
HITEFF=0

XSEC_RATIO=0


OPT_RATIO=""
if test $XSEC_RATIO -eq 1
then
  OPT_RATIO="-d"
fi


if [ $# -ne 2 ]; then
  echo "Two arguments are required.            "
  echo "==============================================="
  echo "Usage:                                         "
  echo "  ./mc_neut.sh <run_id> <srun_id>"
  echo "==============================================="
  exit
fi

RUNID=$1;
SRUNID=$2;

expr "${RUNID}"  + 1 > /dev/null 2>&1;status1=$?
expr "${SRUNID}" + 1 > /dev/null 2>&1;status2=$?
if [ $status1 -ne 0 -o $status2 -ne 0 ] ; then
  echo "run_id must be an integer."
  exit
fi


WORK_DIR="/gpfs/fs03/t2k/beam/work/nchikuma/B2/b2_software"
WAGASCI_SOFT_DIR="$WORK_DIR/wagasci_software/Analysis/bin"
INGRID_SOFT_DIR="$WORK_DIR/INGRID/INGRID/v1r1"  
INGRID_FORM_DIR="$WORK_DIR/INGRID/ingrid_format"
ANALYSIS_DIR="$WORK_DIR/analysis/app"
MC_DIR="$WORK_DIR/b2mc"
NEUT_DIR="/hsm/nu/wagasci/neut/rhc/5.3.2.1"
OUTPUT_DIR="/home/t2k/nchikuma/b2_data2/mc_neut"
OUTPUT_DIR2_0="/home/t2k/nchikuma/b2_data2/reweight/rfg"
OUTPUT_DIR2_1="/home/t2k/nchikuma/b2_data2/reweight/rpa"
REWEIGHT_DIR="/hsm/t2k/t2k_JB/t2k_beam/kenichi/reweight/rfg"
REWEIGHT_DIR2="/hsm/t2k/t2k_JB/t2k_beam/kenichi/reweight/rfg_rpa"
DETRES_DIR="/home/t2k/nchikuma/b2_data3/sysE_detRes/data"
DETRES2_DIR="/home/t2k/nchikuma/b2_data3/sysE_detRes2/data"


# File
NEUTFILE1=`printf "%s/nd7/b2_numubar_nd7_h2o_%05d.nt"   $NEUT_DIR $SRUNID`
NEUTFILE2=`printf "%s/nd7/b2_numu_nd7_h2o_%05d.nt"      $NEUT_DIR $SRUNID`

ALLFILE1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_*.root" $OUTPUT_DIR $RUNID $SRUNID`
ALLFILE2=`printf "%s/h2o/numu/ingrid_%08d_%04d_*.root"    $OUTPUT_DIR $RUNID $SRUNID`

OUTPUTFILE1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_calib.root" $OUTPUT_DIR $RUNID $SRUNID`
OUTPUTFILE2=`printf "%s/h2o/numu/ingrid_%08d_%04d_calib.root"    $OUTPUT_DIR $RUNID $SRUNID`

RECONFILE1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_recon.root" $OUTPUT_DIR $RUNID $SRUNID`
RECONFILE2=`printf "%s/h2o/numu/ingrid_%08d_%04d_recon.root"    $OUTPUT_DIR $RUNID $SRUNID`

ANAFILE1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_anas1.root" $OUTPUT_DIR $RUNID $SRUNID`
ANAFILE2=`printf "%s/h2o/numu/ingrid_%08d_%04d_anas1.root"    $OUTPUT_DIR $RUNID $SRUNID`

PLOTFILE1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_plot.root" $OUTPUT_DIR $RUNID $SRUNID`
PLOTFILE2=`printf "%s/h2o/numu/ingrid_%08d_%04d_plot.root"    $OUTPUT_DIR $RUNID $SRUNID`

REWEIGHTFILE1_0=`printf "%s/h2o/numubar/ingrid_%08d_%04d_reweight.root" $REWEIGHT_DIR $RUNID $SRUNID`
REWEIGHTFILE2_0=`printf "%s/h2o/numu/ingrid_%08d_%04d_reweight.root"    $REWEIGHT_DIR $RUNID $SRUNID`

TRKEFFFILE1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_trkeff.root" $OUTPUT_DIR $RUNID $SRUNID`
TRKEFFFILE2=`printf "%s/h2o/numu/ingrid_%08d_%04d_trkeff.root"    $OUTPUT_DIR $RUNID $SRUNID`

PEFILE1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_pe.root" $OUTPUT_DIR $RUNID $SRUNID`
PEFILE2=`printf "%s/h2o/numu/ingrid_%08d_%04d_pe.root"    $OUTPUT_DIR $RUNID $SRUNID`

HITEFFFILE1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_hiteff.root" $OUTPUT_DIR $RUNID $SRUNID`
HITEFFFILE2=`printf "%s/h2o/numu/ingrid_%08d_%04d_hiteff.root"    $OUTPUT_DIR $RUNID $SRUNID`


TUNEFILE="~/b2_data/jnubeam/tunefile/tune_nd7_8_9.root"


#Set up
cmd="source ${WORK_DIR}/Run_At_Start_B2.sh"
echo -e "---------\n ${cmd}\n"; eval $cmd
cmd="source ${INGRID_SOFT_DIR}/cmt/setup.sh";
echo -e "---------\n ${cmd}\n"; eval $cmd
cmd="source ${MC_DIR}/setup.sh";
echo -e "---------\n ${cmd}\n"; eval $cmd


if test $MCSIM -eq 1
then
  echo "MC simulation"
  #Remove All Files
  cmd="rm -f $ALLFILE1"
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd="rm -f $ALLFILE2"
  echo -e "---------\n ${cmd}\n"; eval $cmd

  #MC Job
  cmd=`printf "%s/bin/Linux-g++/b2mc -m7 -f2 -i %s -o %s" $MC_DIR $NEUTFILE1 $OUTPUTFILE1`
  echo -e "---------\n ${cmd}\n"; eval $cmd

fi

if test $TWODRECON -eq 1
then
  echo "2D reconstruction"

  #Recon
  cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE1} ${RECONFILE1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd

fi

if test $THRDRECON -eq 1
then
  echo "3D reconstruction"

  #Three Dim Recon
  cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x " ${ANALYSIS_DIR} ${RECONFILE1} ${ANAFILE1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd

fi


if test $PLOT_ALL -eq 1
then
  echo "plot all"

  ##Plot 
  cmd=`printf "%s/B2Plot -i %s -o %s -m 7 -t %s %s" ${ANALYSIS_DIR} ${ANAFILE1} ${PLOTFILE1} ${TUNEFILE} ${OPT_RATIO}`
  echo -e "---------\n ${cmd}\n"; eval $cmd

fi
