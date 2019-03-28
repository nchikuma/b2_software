#!/bin/sh

MCSIM=0
TWODRECON=0
THRDRECON=0
PLOT_ALL=0
PLOT_REWEI=1


XSEC_RATIO=1
USENO3DTRK=1

ONLYNSEL=1
ONLY_NOMINAL=1


OPT_RATIO=""
if test $XSEC_RATIO -eq 1
then
  OPT_RATIO="-d"
fi

OPT_NOT3DRRK=""
if test $USENO3DTRK -eq 1
then
  OPT_NOT3DRRK=" -t"
fi

OPT_ONLYNSEL=""
if test $ONLYNSEL -eq 1
then
  OPT_ONLYNSEL="-p"
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

#NEUT_DIR="/hsm/nu/wagasci/neut/rhc/5.3.2.1"
#NEUT_DIR="/hsm/t2k/t2k_JB/t2k_beam/kenichi/neut"
#NEUT_DIR="/home/t2k/nchikuma/neut_nfsi"
#EUT_DIR="/home/t2k/kenichi/5401"
NEUT_DIR="/hsm/nu/wagasci/neut/rhc/5.4.0.1"

#OUTPUT_DIR="/home/t2k/nchikuma/b2_data2/mc_neut"
#OUTPUT_DIR="/home/t2k/nchikuma/b2_data2/mc_neut2"
OUTPUT_DIR="/home/t2k/nchikuma/b2_data2/mc_neut3"

REWEIGHT_DIR="/hsm/t2k/t2k_JB/t2k_beam/kenichi/reweight3/rfg"
REWEIGHT_DIR2="/hsm/t2k/t2k_JB/t2k_beam/kenichi/reweight3/rfg_rpa"


# File
NEUTFILE1=`printf "%s/nd7/b2_numubar_nd7_h2o_%05d.nt"   $NEUT_DIR $SRUNID`
NEUTFILE2=`printf "%s/nd7/b2_numu_nd7_h2o_%05d.nt"      $NEUT_DIR $SRUNID`
NEUTFILE3=`printf "%s/nd8/b2_numubar_nd8_ch_%05d.nt"    $NEUT_DIR $SRUNID`
NEUTFILE4=`printf "%s/nd8/b2_numu_nd8_ch_%05d.nt"       $NEUT_DIR $SRUNID`
NEUTFILE5=`printf "%s/nd9/b2_numubar_nd9_fe_%05d.nt"    $NEUT_DIR $SRUNID`
NEUTFILE6=`printf "%s/nd9/b2_numu_nd9_fe_%05d.nt"       $NEUT_DIR $SRUNID`

ALLFILE1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_*.root" $OUTPUT_DIR $RUNID $SRUNID`
ALLFILE2=`printf "%s/h2o/numu/ingrid_%08d_%04d_*.root"    $OUTPUT_DIR $RUNID $SRUNID`
ALLFILE3=`printf "%s/ch/numubar/ingrid_%08d_%04d_*.root"  $OUTPUT_DIR $RUNID $SRUNID`
ALLFILE4=`printf "%s/ch/numu/ingrid_%08d_%04d_*.root"     $OUTPUT_DIR $RUNID $SRUNID`
ALLFILE5=`printf "%s/fe/numubar/ingrid_%08d_%04d_*.root"  $OUTPUT_DIR $RUNID $SRUNID`
ALLFILE6=`printf "%s/fe/numu/ingrid_%08d_%04d_*.root"     $OUTPUT_DIR $RUNID $SRUNID`

OUTPUTFILE1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_calib.root" $OUTPUT_DIR $RUNID $SRUNID`
OUTPUTFILE2=`printf "%s/h2o/numu/ingrid_%08d_%04d_calib.root"    $OUTPUT_DIR $RUNID $SRUNID`
OUTPUTFILE3=`printf "%s/ch/numubar/ingrid_%08d_%04d_calib.root"  $OUTPUT_DIR $RUNID $SRUNID`
OUTPUTFILE4=`printf "%s/ch/numu/ingrid_%08d_%04d_calib.root"     $OUTPUT_DIR $RUNID $SRUNID`
OUTPUTFILE5=`printf "%s/fe/numubar/ingrid_%08d_%04d_calib.root"  $OUTPUT_DIR $RUNID $SRUNID`
OUTPUTFILE6=`printf "%s/fe/numu/ingrid_%08d_%04d_calib.root"     $OUTPUT_DIR $RUNID $SRUNID`

RECONFILE1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_recon.root" $OUTPUT_DIR $RUNID $SRUNID`
RECONFILE2=`printf "%s/h2o/numu/ingrid_%08d_%04d_recon.root"    $OUTPUT_DIR $RUNID $SRUNID`
RECONFILE3=`printf "%s/ch/numubar/ingrid_%08d_%04d_recon.root"  $OUTPUT_DIR $RUNID $SRUNID`
RECONFILE4=`printf "%s/ch/numu/ingrid_%08d_%04d_recon.root"     $OUTPUT_DIR $RUNID $SRUNID`
RECONFILE5=`printf "%s/fe/numubar/ingrid_%08d_%04d_recon.root"  $OUTPUT_DIR $RUNID $SRUNID`
RECONFILE6=`printf "%s/fe/numu/ingrid_%08d_%04d_recon.root"     $OUTPUT_DIR $RUNID $SRUNID`

ANAFILE1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_anas1.root" $OUTPUT_DIR $RUNID $SRUNID`
ANAFILE2=`printf "%s/h2o/numu/ingrid_%08d_%04d_anas1.root"    $OUTPUT_DIR $RUNID $SRUNID`
ANAFILE3=`printf "%s/ch/numubar/ingrid_%08d_%04d_anas1.root"  $OUTPUT_DIR $RUNID $SRUNID`
ANAFILE4=`printf "%s/ch/numu/ingrid_%08d_%04d_anas1.root"     $OUTPUT_DIR $RUNID $SRUNID`
ANAFILE5=`printf "%s/fe/numubar/ingrid_%08d_%04d_anas1.root"  $OUTPUT_DIR $RUNID $SRUNID`
ANAFILE6=`printf "%s/fe/numu/ingrid_%08d_%04d_anas1.root"     $OUTPUT_DIR $RUNID $SRUNID`

PLOTFILE1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_plot.root" $OUTPUT_DIR $RUNID $SRUNID`
PLOTFILE2=`printf "%s/h2o/numu/ingrid_%08d_%04d_plot.root"    $OUTPUT_DIR $RUNID $SRUNID`
PLOTFILE3=`printf "%s/ch/numubar/ingrid_%08d_%04d_plot.root"  $OUTPUT_DIR $RUNID $SRUNID`
PLOTFILE4=`printf "%s/ch/numu/ingrid_%08d_%04d_plot.root"     $OUTPUT_DIR $RUNID $SRUNID`
PLOTFILE5=`printf "%s/fe/numubar/ingrid_%08d_%04d_plot.root"  $OUTPUT_DIR $RUNID $SRUNID`
PLOTFILE6=`printf "%s/fe/numu/ingrid_%08d_%04d_plot.root"     $OUTPUT_DIR $RUNID $SRUNID`

REWEIGHTFILE1_0=`printf "%s/h2o/numubar/ingrid_%08d_%04d_reweight.root" $REWEIGHT_DIR $RUNID $SRUNID`
REWEIGHTFILE2_0=`printf "%s/h2o/numu/ingrid_%08d_%04d_reweight.root"    $REWEIGHT_DIR $RUNID $SRUNID`
REWEIGHTFILE3_0=`printf "%s/ch/numubar/ingrid_%08d_%04d_reweight.root"  $REWEIGHT_DIR $RUNID $SRUNID`
REWEIGHTFILE4_0=`printf "%s/ch/numu/ingrid_%08d_%04d_reweight.root"     $REWEIGHT_DIR $RUNID $SRUNID`
REWEIGHTFILE1_1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_reweight.root" $REWEIGHT_DIR2 $RUNID $SRUNID`
REWEIGHTFILE2_1=`printf "%s/h2o/numu/ingrid_%08d_%04d_reweight.root"    $REWEIGHT_DIR2 $RUNID $SRUNID`
REWEIGHTFILE3_1=`printf "%s/ch/numubar/ingrid_%08d_%04d_reweight.root"  $REWEIGHT_DIR2 $RUNID $SRUNID`
REWEIGHTFILE4_1=`printf "%s/ch/numu/ingrid_%08d_%04d_reweight.root"     $REWEIGHT_DIR2 $RUNID $SRUNID`


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


#  #Remove All Files
#  cmd="rm -f $ALLFILE1"
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd="rm -f $ALLFILE2"
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd="rm -f $ALLFILE3"
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd="rm -f $ALLFILE4"
#  echo -e "---------\n ${cmd}\n"; eval $cmd

  #MC Job
  cmd=`printf "%s/bin/Linux-g++/b2mc -m7 -f2 -i %s -o %s" $MC_DIR $NEUTFILE1 $OUTPUTFILE1`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/bin/Linux-g++/b2mc -m7 -f1 -i %s -o %s" $MC_DIR $NEUTFILE2 $OUTPUTFILE2`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/bin/Linux-g++/b2mc -m8 -f2 -i %s -o %s" $MC_DIR $NEUTFILE3 $OUTPUTFILE3`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/bin/Linux-g++/b2mc -m8 -f1 -i %s -o %s" $MC_DIR $NEUTFILE4 $OUTPUTFILE4`
  echo -e "---------\n ${cmd}\n"; eval $cmd

fi

if test $TWODRECON -eq 1
then
  echo "2D reconstruction"

  #Recon
  cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE1} ${RECONFILE1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE2} ${RECONFILE2}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE3} ${RECONFILE3}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE4} ${RECONFILE4}`
  echo -e "---------\n ${cmd}\n"; eval $cmd

fi

if test $THRDRECON -eq 1
then
  echo "3D reconstruction"

  #Three Dim Recon
  cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${RECONFILE1} ${ANAFILE1} ${OPT_NOT3DRRK}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${RECONFILE2} ${ANAFILE2} ${OPT_NOT3DRRK}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${RECONFILE3} ${ANAFILE3} ${OPT_NOT3DRRK}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${RECONFILE4} ${ANAFILE4} ${OPT_NOT3DRRK}`
  echo -e "---------\n ${cmd}\n"; eval $cmd

fi


if test $PLOT_ALL -eq 1
then
  echo "plot all"

  ##Plot 
  cmd=`printf "%s/B2Plot -i %s -o %s -m 7 -t %s %s %s" ${ANALYSIS_DIR} ${ANAFILE1} ${PLOTFILE1} ${TUNEFILE} ${OPT_RATIO} ${OPT_ONLYNSEL}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2Plot -i %s -o %s -m 7 -t %s %s %s" ${ANALYSIS_DIR} ${ANAFILE2} ${PLOTFILE2} ${TUNEFILE} ${OPT_RATIO} ${OPT_ONLYNSEL}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2Plot -i %s -o %s -m 8 -t %s %s %s" ${ANALYSIS_DIR} ${ANAFILE3} ${PLOTFILE3} ${TUNEFILE} ${OPT_RATIO} ${OPT_ONLYNSEL}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2Plot -i %s -o %s -m 8 -t %s %s %s" ${ANALYSIS_DIR} ${ANAFILE4} ${PLOTFILE4} ${TUNEFILE} ${OPT_RATIO} ${OPT_ONLYNSEL}`
  echo -e "---------\n ${cmd}\n"; eval $cmd

fi

if test $PLOT_REWEI -eq 1
then
  echo "Plot with t2kReweight"


  for l in `seq 0 1`
  do
    if test $l -eq 0
    then
      REWEIGHTFILE1=${REWEIGHTFILE1_0}
      REWEIGHTFILE2=${REWEIGHTFILE2_0}
      REWEIGHTFILE3=${REWEIGHTFILE3_0}
      REWEIGHTFILE4=${REWEIGHTFILE4_0}
    else 
      REWEIGHTFILE1=${REWEIGHTFILE1_1}
      REWEIGHTFILE2=${REWEIGHTFILE2_1}
      REWEIGHTFILE3=${REWEIGHTFILE3_1}
      REWEIGHTFILE4=${REWEIGHTFILE4_1}
    fi

    dial=10

    echo $dial

    OPTPLOT=""
    if test $dial -eq 10
    then
      OPTPLOT=$OPT_ONLYNSEL
    else
      OPTPLOT="-p"
    fi
    
    PLOTFILE1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_plot_%d.root" $OUTPUT_DIR $RUNID $SRUNID $l`
    PLOTFILE2=`printf "%s/h2o/numu/ingrid_%08d_%04d_plot_%d.root"    $OUTPUT_DIR $RUNID $SRUNID $l`
    PLOTFILE3=`printf "%s/ch/numubar/ingrid_%08d_%04d_plot_%d.root"  $OUTPUT_DIR $RUNID $SRUNID $l`
    PLOTFILE4=`printf "%s/ch/numu/ingrid_%08d_%04d_plot_%d.root"     $OUTPUT_DIR $RUNID $SRUNID $l`

    cmd=`printf "rm -f %s" ${PLOTFILE1}` echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "rm -f %s" ${PLOTFILE2}` echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "rm -f %s" ${PLOTFILE3}` echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "rm -f %s" ${PLOTFILE4}` echo -e "---------\n ${cmd}\n"; eval $cmd


    cmd=`printf "%s/B2Plot -i %s -o %s -m 7 -t %s -w %d -r %s %s %s" ${ANALYSIS_DIR} ${ANAFILE1} ${PLOTFILE1} ${TUNEFILE} ${dial} ${REWEIGHTFILE1} ${OPT_RATIO} ${OPTPLOT}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/B2Plot -i %s -o %s -m 7 -t %s -w %d -r %s %s %s" ${ANALYSIS_DIR} ${ANAFILE2} ${PLOTFILE2} ${TUNEFILE} ${dial} ${REWEIGHTFILE2} ${OPT_RATIO} ${OPTPLOT}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/B2Plot -i %s -o %s -m 8 -t %s -w %d -r %s %s %s" ${ANALYSIS_DIR} ${ANAFILE3} ${PLOTFILE3} ${TUNEFILE} ${dial} ${REWEIGHTFILE3} ${OPT_RATIO} ${OPTPLOT}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/B2Plot -i %s -o %s -m 8 -t %s -w %d -r %s %s %s" ${ANALYSIS_DIR} ${ANAFILE4} ${PLOTFILE4} ${TUNEFILE} ${dial} ${REWEIGHTFILE4} ${OPT_RATIO} ${OPTPLOT}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    
  done

fi


