#!/bin/sh

MCSIM=0
TWODRECON=0
THRDRECON=0
PLOT_ALL=0
DETRES=0
TRKEFF=1
PECHECK=0
HITEFF=0

USENO3DTRK=1

ONLYNSEL=1


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
NEUT_DIR="/hsm/nu/wagasci/neut/rhc/5.3.2.1"
OUTPUT_DIR="/home/t2k/nchikuma/b2_data/mc_neut"
OUTPUT_DIR2="/home/t2k/nchikuma/b2_data2/mc_neut"
DETRES_DIR="/home/t2k/nchikuma/b2_data3/sysE_detRes/data"

NEUTFILE1=`printf "%s/nd1/b2_numubar_nd1_o_%05d.nt"  $NEUT_DIR $SRUNID`
NEUTFILE2=`printf "%s/nd1/b2_numu_nd1_o_%05d.nt"     $NEUT_DIR $SRUNID`
NEUTFILE3=`printf "%s/nd5/b2_numubar_nd5_o_%05d.nt"  $NEUT_DIR $SRUNID`
NEUTFILE4=`printf "%s/nd5/b2_numu_nd5_o_%05d.nt"     $NEUT_DIR $SRUNID`
NEUTFILE5=`printf "%s/nd6/b2_numubar_nd6_o_%05d.nt"  $NEUT_DIR $SRUNID`
NEUTFILE6=`printf "%s/nd6/b2_numu_nd6_o_%05d.nt"     $NEUT_DIR $SRUNID`


ALLFILE1=`printf "%s/wallbg/numubar/ingrid_%08d_%04d_*.root"      $OUTPUT_DIR2 $RUNID $SRUNID`
ALLFILE2=`printf "%s/wallbg/numu/ingrid_%08d_%04d_*.root"         $OUTPUT_DIR2 $RUNID $SRUNID`
ALLFILE3_0=`printf "%s/ceiling/numubar/ingrid_%08d_%04d_*.root"   $OUTPUT_DIR2 $RUNID $SRUNID`
ALLFILE4_0=`printf "%s/ceiling/numu/ingrid_%08d_%04d_*.root"      $OUTPUT_DIR2 $RUNID $SRUNID`
ALLFILE3_1=`printf "%s/floor/numubar/ingrid_%08d_%04d_*.root"     $OUTPUT_DIR2 $RUNID $SRUNID`
ALLFILE4_1=`printf "%s/floor/numu/ingrid_%08d_%04d_*.root"        $OUTPUT_DIR2 $RUNID $SRUNID`
ALLFILE5_0=`printf "%s/pillar_r/numubar/ingrid_%08d_%04d_*.root"  $OUTPUT_DIR2 $RUNID $SRUNID`
ALLFILE6_0=`printf "%s/pillar_r/numu/ingrid_%08d_%04d_*.root"     $OUTPUT_DIR2 $RUNID $SRUNID`
ALLFILE5_1=`printf "%s/pillar_l/numubar/ingrid_%08d_%04d_*.root"  $OUTPUT_DIR2 $RUNID $SRUNID`
ALLFILE6_1=`printf "%s/pillar_l/numu/ingrid_%08d_%04d_*.root"     $OUTPUT_DIR2 $RUNID $SRUNID`

OUTPUTFILE1=`printf "%s/wallbg/numubar/ingrid_%08d_%04d_calib.root"      $OUTPUT_DIR2 $RUNID $SRUNID`
OUTPUTFILE2=`printf "%s/wallbg/numu/ingrid_%08d_%04d_calib.root"         $OUTPUT_DIR2 $RUNID $SRUNID`
OUTPUTFILE3_0=`printf "%s/ceiling/numubar/ingrid_%08d_%04d_calib.root"   $OUTPUT_DIR2 $RUNID $SRUNID`
OUTPUTFILE4_0=`printf "%s/ceiling/numu/ingrid_%08d_%04d_calib.root"      $OUTPUT_DIR2 $RUNID $SRUNID`
OUTPUTFILE3_1=`printf "%s/floor/numubar/ingrid_%08d_%04d_calib.root"     $OUTPUT_DIR2 $RUNID $SRUNID`
OUTPUTFILE4_1=`printf "%s/floor/numu/ingrid_%08d_%04d_calib.root"        $OUTPUT_DIR2 $RUNID $SRUNID`
OUTPUTFILE5_0=`printf "%s/pillar_r/numubar/ingrid_%08d_%04d_calib.root"  $OUTPUT_DIR2 $RUNID $SRUNID`
OUTPUTFILE6_0=`printf "%s/pillar_r/numu/ingrid_%08d_%04d_calib.root"     $OUTPUT_DIR2 $RUNID $SRUNID`
OUTPUTFILE5_1=`printf "%s/pillar_l/numubar/ingrid_%08d_%04d_calib.root"  $OUTPUT_DIR2 $RUNID $SRUNID`
OUTPUTFILE6_1=`printf "%s/pillar_l/numu/ingrid_%08d_%04d_calib.root"     $OUTPUT_DIR2 $RUNID $SRUNID`

RECONFILE1=`printf "%s/wallbg/numubar/ingrid_%08d_%04d_recon.root"      $OUTPUT_DIR2 $RUNID $SRUNID`
RECONFILE2=`printf "%s/wallbg/numu/ingrid_%08d_%04d_recon.root"         $OUTPUT_DIR2 $RUNID $SRUNID`
RECONFILE3_0=`printf "%s/ceiling/numubar/ingrid_%08d_%04d_recon.root"   $OUTPUT_DIR2 $RUNID $SRUNID`
RECONFILE4_0=`printf "%s/ceiling/numu/ingrid_%08d_%04d_recon.root"      $OUTPUT_DIR2 $RUNID $SRUNID`
RECONFILE3_1=`printf "%s/floor/numubar/ingrid_%08d_%04d_recon.root"     $OUTPUT_DIR2 $RUNID $SRUNID`
RECONFILE4_1=`printf "%s/floor/numu/ingrid_%08d_%04d_recon.root"        $OUTPUT_DIR2 $RUNID $SRUNID`
RECONFILE5_0=`printf "%s/pillar_r/numubar/ingrid_%08d_%04d_recon.root"  $OUTPUT_DIR2 $RUNID $SRUNID`
RECONFILE6_0=`printf "%s/pillar_r/numu/ingrid_%08d_%04d_recon.root"     $OUTPUT_DIR2 $RUNID $SRUNID`
RECONFILE5_1=`printf "%s/pillar_l/numubar/ingrid_%08d_%04d_recon.root"  $OUTPUT_DIR2 $RUNID $SRUNID`
RECONFILE6_1=`printf "%s/pillar_l/numu/ingrid_%08d_%04d_recon.root"     $OUTPUT_DIR2 $RUNID $SRUNID`

ANAFILE1=`printf "%s/wallbg/numubar/ingrid_%08d_%04d_anas1.root"      $OUTPUT_DIR2 $RUNID $SRUNID`
ANAFILE2=`printf "%s/wallbg/numu/ingrid_%08d_%04d_anas1.root"         $OUTPUT_DIR2 $RUNID $SRUNID`
ANAFILE3_0=`printf "%s/ceiling/numubar/ingrid_%08d_%04d_anas1.root"   $OUTPUT_DIR2 $RUNID $SRUNID`
ANAFILE4_0=`printf "%s/ceiling/numu/ingrid_%08d_%04d_anas1.root"      $OUTPUT_DIR2 $RUNID $SRUNID`
ANAFILE3_1=`printf "%s/floor/numubar/ingrid_%08d_%04d_anas1.root"     $OUTPUT_DIR2 $RUNID $SRUNID`
ANAFILE4_1=`printf "%s/floor/numu/ingrid_%08d_%04d_anas1.root"        $OUTPUT_DIR2 $RUNID $SRUNID`
ANAFILE5_0=`printf "%s/pillar_r/numubar/ingrid_%08d_%04d_anas1.root"  $OUTPUT_DIR2 $RUNID $SRUNID`
ANAFILE6_0=`printf "%s/pillar_r/numu/ingrid_%08d_%04d_anas1.root"     $OUTPUT_DIR2 $RUNID $SRUNID`
ANAFILE5_1=`printf "%s/pillar_l/numubar/ingrid_%08d_%04d_anas1.root"  $OUTPUT_DIR2 $RUNID $SRUNID`
ANAFILE6_1=`printf "%s/pillar_l/numu/ingrid_%08d_%04d_anas1.root"     $OUTPUT_DIR2 $RUNID $SRUNID`

PLOTFILE1=`printf "%s/wallbg/numubar/ingrid_%08d_%04d_plot.root"      $OUTPUT_DIR2 $RUNID $SRUNID`
PLOTFILE2=`printf "%s/wallbg/numu/ingrid_%08d_%04d_plot.root"         $OUTPUT_DIR2 $RUNID $SRUNID`
PLOTFILE3_0=`printf "%s/ceiling/numubar/ingrid_%08d_%04d_plot.root"   $OUTPUT_DIR2 $RUNID $SRUNID`
PLOTFILE4_0=`printf "%s/ceiling/numu/ingrid_%08d_%04d_plot.root"      $OUTPUT_DIR2 $RUNID $SRUNID`
PLOTFILE3_1=`printf "%s/floor/numubar/ingrid_%08d_%04d_plot.root"     $OUTPUT_DIR2 $RUNID $SRUNID`
PLOTFILE4_1=`printf "%s/floor/numu/ingrid_%08d_%04d_plot.root"        $OUTPUT_DIR2 $RUNID $SRUNID`
PLOTFILE5_0=`printf "%s/pillar_r/numubar/ingrid_%08d_%04d_plot.root"  $OUTPUT_DIR2 $RUNID $SRUNID`
PLOTFILE6_0=`printf "%s/pillar_r/numu/ingrid_%08d_%04d_plot.root"     $OUTPUT_DIR2 $RUNID $SRUNID`
PLOTFILE5_1=`printf "%s/pillar_l/numubar/ingrid_%08d_%04d_plot.root"  $OUTPUT_DIR2 $RUNID $SRUNID`
PLOTFILE6_1=`printf "%s/pillar_l/numu/ingrid_%08d_%04d_plot.root"     $OUTPUT_DIR2 $RUNID $SRUNID`

TRKEFFFILE1=`printf "%s/wallbg/numubar/ingrid_%08d_%04d_trkeff.root"      $OUTPUT_DIR2 $RUNID $SRUNID`
TRKEFFFILE2=`printf "%s/wallbg/numu/ingrid_%08d_%04d_trkeff.root"         $OUTPUT_DIR2 $RUNID $SRUNID`
TRKEFFFILE3_0=`printf "%s/ceiling/numubar/ingrid_%08d_%04d_trkeff.root"   $OUTPUT_DIR2 $RUNID $SRUNID`
TRKEFFFILE4_0=`printf "%s/ceiling/numu/ingrid_%08d_%04d_trkeff.root"      $OUTPUT_DIR2 $RUNID $SRUNID`
TRKEFFFILE3_1=`printf "%s/floor/numubar/ingrid_%08d_%04d_trkeff.root"     $OUTPUT_DIR2 $RUNID $SRUNID`
TRKEFFFILE4_1=`printf "%s/floor/numu/ingrid_%08d_%04d_trkeff.root"        $OUTPUT_DIR2 $RUNID $SRUNID`
TRKEFFFILE5_0=`printf "%s/pillar_r/numubar/ingrid_%08d_%04d_trkeff.root"  $OUTPUT_DIR2 $RUNID $SRUNID`
TRKEFFFILE6_0=`printf "%s/pillar_r/numu/ingrid_%08d_%04d_trkeff.root"     $OUTPUT_DIR2 $RUNID $SRUNID`
TRKEFFFILE5_1=`printf "%s/pillar_l/numubar/ingrid_%08d_%04d_trkeff.root"  $OUTPUT_DIR2 $RUNID $SRUNID`
TRKEFFFILE6_1=`printf "%s/pillar_l/numu/ingrid_%08d_%04d_trkeff.root"     $OUTPUT_DIR2 $RUNID $SRUNID`

PEFILE1=`printf "%s/wallbg/numubar/ingrid_%08d_%04d_pe.root"      $OUTPUT_DIR2 $RUNID $SRUNID`
PEFILE2=`printf "%s/wallbg/numu/ingrid_%08d_%04d_pe.root"         $OUTPUT_DIR2 $RUNID $SRUNID`
PEFILE3_0=`printf "%s/ceiling/numubar/ingrid_%08d_%04d_pe.root"   $OUTPUT_DIR2 $RUNID $SRUNID`
PEFILE4_0=`printf "%s/ceiling/numu/ingrid_%08d_%04d_pe.root"      $OUTPUT_DIR2 $RUNID $SRUNID`
PEFILE3_1=`printf "%s/floor/numubar/ingrid_%08d_%04d_pe.root"     $OUTPUT_DIR2 $RUNID $SRUNID`
PEFILE4_1=`printf "%s/floor/numu/ingrid_%08d_%04d_pe.root"        $OUTPUT_DIR2 $RUNID $SRUNID`
PEFILE5_0=`printf "%s/pillar_r/numubar/ingrid_%08d_%04d_pe.root"  $OUTPUT_DIR2 $RUNID $SRUNID`
PEFILE6_0=`printf "%s/pillar_r/numu/ingrid_%08d_%04d_pe.root"     $OUTPUT_DIR2 $RUNID $SRUNID`
PEFILE5_1=`printf "%s/pillar_l/numubar/ingrid_%08d_%04d_pe.root"  $OUTPUT_DIR2 $RUNID $SRUNID`
PEFILE6_1=`printf "%s/pillar_l/numu/ingrid_%08d_%04d_pe.root"     $OUTPUT_DIR2 $RUNID $SRUNID`

HITEFFFILE1=`printf "%s/wallbg/numubar/ingrid_%08d_%04d_hiteff.root"      $OUTPUT_DIR2 $RUNID $SRUNID`
HITEFFFILE2=`printf "%s/wallbg/numu/ingrid_%08d_%04d_hiteff.root"         $OUTPUT_DIR2 $RUNID $SRUNID`
HITEFFFILE3_0=`printf "%s/ceiling/numubar/ingrid_%08d_%04d_hiteff.root"   $OUTPUT_DIR2 $RUNID $SRUNID`
HITEFFFILE4_0=`printf "%s/ceiling/numu/ingrid_%08d_%04d_hiteff.root"      $OUTPUT_DIR2 $RUNID $SRUNID`
HITEFFFILE3_1=`printf "%s/floor/numubar/ingrid_%08d_%04d_hiteff.root"     $OUTPUT_DIR2 $RUNID $SRUNID`
HITEFFFILE4_1=`printf "%s/floor/numu/ingrid_%08d_%04d_hiteff.root"        $OUTPUT_DIR2 $RUNID $SRUNID`
HITEFFFILE5_0=`printf "%s/pillar_r/numubar/ingrid_%08d_%04d_hiteff.root"  $OUTPUT_DIR2 $RUNID $SRUNID`
HITEFFFILE6_0=`printf "%s/pillar_r/numu/ingrid_%08d_%04d_hiteff.root"     $OUTPUT_DIR2 $RUNID $SRUNID`
HITEFFFILE5_1=`printf "%s/pillar_l/numubar/ingrid_%08d_%04d_hiteff.root"  $OUTPUT_DIR2 $RUNID $SRUNID`
HITEFFFILE6_1=`printf "%s/pillar_l/numu/ingrid_%08d_%04d_hiteff.root"     $OUTPUT_DIR2 $RUNID $SRUNID`



TUNEFILE="~/b2_data/jnubeam/tunefile/tune_nd7_8_9.root"

##Set up
cmd="source ${WORK_DIR}/Run_At_Start_B2.sh"
echo -e "---------\n ${cmd}\n"; eval $cmd
cmd="source ${INGRID_SOFT_DIR}/cmt/setup.sh";
echo -e "---------\n ${cmd}\n"; eval $cmd
cmd="source ${MC_DIR}/setup.sh";
echo -e "---------\n ${cmd}\n"; eval $cmd


if test $MCSIM -eq 1
then
  echo "MC Simulation"

#  cmd="rm -f $ALLFILE1"
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd="rm -f $ALLFILE2"
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd="rm -f $ALLFILE3_0"
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd="rm -f $ALLFILE4_0"
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd="rm -f $ALLFILE5_0"
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd="rm -f $ALLFILE6_0"
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd="rm -f $ALLFILE3_1"
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd="rm -f $ALLFILE4_1"
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd="rm -f $ALLFILE5_1"
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd="rm -f $ALLFILE6_1"
#  echo -e "---------\n ${cmd}\n"; eval $cmd

#  # Upst Wall
#  cmd=`printf "%s/bin/Linux-g++/b2mc -m10 -f2 -i %s -o %s" $MC_DIR $NEUTFILE1 $OUTPUTFILE1`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd=`printf "%s/bin/Linux-g++/b2mc -m10 -f1 -i %s -o %s" $MC_DIR $NEUTFILE2 $OUTPUTFILE2`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  # Ceiling
#  cmd=`printf "%s/bin/Linux-g++/b2mc -m11 -f2 -i %s -o %s" $MC_DIR $NEUTFILE3 $OUTPUTFILE3_0`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd=`printf "%s/bin/Linux-g++/b2mc -m11 -f1 -i %s -o %s" $MC_DIR $NEUTFILE4 $OUTPUTFILE4_0`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  # Floor
#  cmd=`printf "%s/bin/Linux-g++/b2mc -m12 -f2 -i %s -o %s" $MC_DIR $NEUTFILE3 $OUTPUTFILE3_1`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd=`printf "%s/bin/Linux-g++/b2mc -m12 -f1 -i %s -o %s" $MC_DIR $NEUTFILE4 $OUTPUTFILE4_1`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  # Right Pillar
#  cmd=`printf "%s/bin/Linux-g++/b2mc -m13 -f2 -i %s -o %s" $MC_DIR $NEUTFILE5 $OUTPUTFILE5_0`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd=`printf "%s/bin/Linux-g++/b2mc -m13 -f1 -i %s -o %s" $MC_DIR $NEUTFILE6 $OUTPUTFILE6_0`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  # Left Pillar
#  cmd=`printf "%s/bin/Linux-g++/b2mc -m14 -f2 -i %s -o %s" $MC_DIR $NEUTFILE5 $OUTPUTFILE5_1`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd=`printf "%s/bin/Linux-g++/b2mc -m14 -f1 -i %s -o %s" $MC_DIR $NEUTFILE6 $OUTPUTFILE6_1`
#  echo -e "---------\n ${cmd}\n"; eval $cmd

fi

if test $TWODRECON -eq 1
then
  #Recon
  cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE1} ${RECONFILE1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE2} ${RECONFILE2}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE3_0} ${RECONFILE3_0}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE4_0} ${RECONFILE4_0}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE3_1} ${RECONFILE3_1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE4_1} ${RECONFILE4_1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE5_0} ${RECONFILE5_0}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE6_0} ${RECONFILE6_0}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE5_1} ${RECONFILE5_1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE6_1} ${RECONFILE6_1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
fi


if test $THRDRECON -eq 1
then
  #Three Dim Recon
  cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${RECONFILE1} ${ANAFILE1}     ${OPT_NOT3DRRK}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${RECONFILE2} ${ANAFILE2}     ${OPT_NOT3DRRK}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${RECONFILE3_0} ${ANAFILE3_0} ${OPT_NOT3DRRK}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${RECONFILE4_0} ${ANAFILE4_0} ${OPT_NOT3DRRK}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${RECONFILE3_1} ${ANAFILE3_1} ${OPT_NOT3DRRK}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${RECONFILE4_1} ${ANAFILE4_1} ${OPT_NOT3DRRK}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${RECONFILE5_0} ${ANAFILE5_0} ${OPT_NOT3DRRK}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${RECONFILE6_0} ${ANAFILE6_0} ${OPT_NOT3DRRK}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${RECONFILE5_1} ${ANAFILE5_1} ${OPT_NOT3DRRK}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${RECONFILE6_1} ${ANAFILE6_1} ${OPT_NOT3DRRK}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
fi


if test $PLOT_ALL -eq 1
then
  ##Plot 
  cmd=`printf "%s/B2Plot -i %s -o %s -m 10 %s" ${ANALYSIS_DIR} ${ANAFILE1}   ${PLOTFILE1}   ${OPT_ONLYNSEL}`; echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2Plot -i %s -o %s -m 10 %s" ${ANALYSIS_DIR} ${ANAFILE2}   ${PLOTFILE2}   ${OPT_ONLYNSEL}`; echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2Plot -i %s -o %s -m 11 %s" ${ANALYSIS_DIR} ${ANAFILE3_0} ${PLOTFILE3_0} ${OPT_ONLYNSEL}`; echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2Plot -i %s -o %s -m 11 %s" ${ANALYSIS_DIR} ${ANAFILE4_0} ${PLOTFILE4_0} ${OPT_ONLYNSEL}`; echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2Plot -i %s -o %s -m 12 %s" ${ANALYSIS_DIR} ${ANAFILE3_1} ${PLOTFILE3_1} ${OPT_ONLYNSEL}`; echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2Plot -i %s -o %s -m 12 %s" ${ANALYSIS_DIR} ${ANAFILE4_1} ${PLOTFILE4_1} ${OPT_ONLYNSEL}`; echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2Plot -i %s -o %s -m 13 %s" ${ANALYSIS_DIR} ${ANAFILE5_0} ${PLOTFILE5_0} ${OPT_ONLYNSEL}`; echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2Plot -i %s -o %s -m 13 %s" ${ANALYSIS_DIR} ${ANAFILE6_0} ${PLOTFILE6_0} ${OPT_ONLYNSEL}`; echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2Plot -i %s -o %s -m 14 %s" ${ANALYSIS_DIR} ${ANAFILE5_1} ${PLOTFILE5_1} ${OPT_ONLYNSEL}`; echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2Plot -i %s -o %s -m 14 %s" ${ANALYSIS_DIR} ${ANAFILE6_1} ${PLOTFILE6_1} ${OPT_ONLYNSEL}`; echo -e "---------\n ${cmd}\n"; eval $cmd
fi


if test $TRKEFF -eq 1
then
  cmd=`printf "%s/TrackEfficiency2 -i %s -o %s -m 10" ${ANALYSIS_DIR} ${ANAFILE1}   ${TRKEFFFILE1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TrackEfficiency2 -i %s -o %s -m 10" ${ANALYSIS_DIR} ${ANAFILE2}   ${TRKEFFFILE2}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TrackEfficiency2 -i %s -o %s -m 11" ${ANALYSIS_DIR} ${ANAFILE3_0} ${TRKEFFFILE3_0}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TrackEfficiency2 -i %s -o %s -m 11" ${ANALYSIS_DIR} ${ANAFILE4_0} ${TRKEFFFILE4_0}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TrackEfficiency2 -i %s -o %s -m 12" ${ANALYSIS_DIR} ${ANAFILE3_1} ${TRKEFFFILE3_1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TrackEfficiency2 -i %s -o %s -m 12" ${ANALYSIS_DIR} ${ANAFILE4_1} ${TRKEFFFILE4_1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TrackEfficiency2 -i %s -o %s -m 13" ${ANALYSIS_DIR} ${ANAFILE5_0} ${TRKEFFFILE5_0}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TrackEfficiency2 -i %s -o %s -m 13" ${ANALYSIS_DIR} ${ANAFILE6_0} ${TRKEFFFILE6_0}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TrackEfficiency2 -i %s -o %s -m 14" ${ANALYSIS_DIR} ${ANAFILE5_1} ${TRKEFFFILE5_1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TrackEfficiency2 -i %s -o %s -m 14" ${ANALYSIS_DIR} ${ANAFILE6_1} ${TRKEFFFILE6_1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd

fi


if test $PECHECK -eq 1
then
  cmd=`printf "%s/B2pe -i %s -o %s -m 10" ${ANALYSIS_DIR} ${ANAFILE1}   ${PEFILE1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2pe -i %s -o %s -m 10" ${ANALYSIS_DIR} ${ANAFILE2}   ${PEFILE2}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2pe -i %s -o %s -m 11" ${ANALYSIS_DIR} ${ANAFILE3_0} ${PEFILE3_0}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2pe -i %s -o %s -m 11" ${ANALYSIS_DIR} ${ANAFILE4_0} ${PEFILE4_0}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2pe -i %s -o %s -m 12" ${ANALYSIS_DIR} ${ANAFILE3_1} ${PEFILE3_1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2pe -i %s -o %s -m 12" ${ANALYSIS_DIR} ${ANAFILE4_1} ${PEFILE4_1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2pe -i %s -o %s -m 13" ${ANALYSIS_DIR} ${ANAFILE5_0} ${PEFILE5_0}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2pe -i %s -o %s -m 13" ${ANALYSIS_DIR} ${ANAFILE6_0} ${PEFILE6_0}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2pe -i %s -o %s -m 14" ${ANALYSIS_DIR} ${ANAFILE5_1} ${PEFILE5_1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2pe -i %s -o %s -m 14" ${ANALYSIS_DIR} ${ANAFILE6_1} ${PEFILE6_1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd

fi

if test $HITEFF -eq 1
then
  cmd=`printf "%s/HitEff2 -i %s -o %s -m 10" ${ANALYSIS_DIR} ${ANAFILE1}   ${HITEFFFILE1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/HitEff2 -i %s -o %s -m 10" ${ANALYSIS_DIR} ${ANAFILE2}   ${HITEFFFILE2}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/HitEff2 -i %s -o %s -m 11" ${ANALYSIS_DIR} ${ANAFILE3_0} ${HITEFFFILE3_0}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/HitEff2 -i %s -o %s -m 11" ${ANALYSIS_DIR} ${ANAFILE4_0} ${HITEFFFILE4_0}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/HitEff2 -i %s -o %s -m 12" ${ANALYSIS_DIR} ${ANAFILE3_1} ${HITEFFFILE3_1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/HitEff2 -i %s -o %s -m 12" ${ANALYSIS_DIR} ${ANAFILE4_1} ${HITEFFFILE4_1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/HitEff2 -i %s -o %s -m 13" ${ANALYSIS_DIR} ${ANAFILE5_0} ${HITEFFFILE5_0}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/HitEff2 -i %s -o %s -m 13" ${ANALYSIS_DIR} ${ANAFILE6_0} ${HITEFFFILE6_0}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/HitEff2 -i %s -o %s -m 14" ${ANALYSIS_DIR} ${ANAFILE5_1} ${HITEFFFILE5_1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/HitEff2 -i %s -o %s -m 14" ${ANALYSIS_DIR} ${ANAFILE6_1} ${HITEFFFILE6_1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd

fi



if test $DETRES -eq 1
then
  echo "Plot for sys. err. from detectro response"


  #for dial in `seq 0 39`
  for dial in `seq 19 26`
  do

    DET_ANAFILE1=`printf    "%s/dial_%02d/wallbg/numubar/ingrid_%08d_%04d_anas1.root"     $DETRES_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE2=`printf    "%s/dial_%02d/wallbg/numu/ingrid_%08d_%04d_anas1.root"        $DETRES_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE3_0=`printf  "%s/dial_%02d/ceiling/numubar/ingrid_%08d_%04d_anas1.root"    $DETRES_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE4_0=`printf  "%s/dial_%02d/ceiling/numu/ingrid_%08d_%04d_anas1.root"       $DETRES_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE3_1=`printf  "%s/dial_%02d/floor/numubar/ingrid_%08d_%04d_anas1.root"      $DETRES_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE4_1=`printf  "%s/dial_%02d/floor/numu/ingrid_%08d_%04d_anas1.root"         $DETRES_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE5_0=`printf  "%s/dial_%02d/pillar_r/numubar/ingrid_%08d_%04d_anas1.root"   $DETRES_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE6_0=`printf  "%s/dial_%02d/pillar_r/numu/ingrid_%08d_%04d_anas1.root"      $DETRES_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE5_1=`printf  "%s/dial_%02d/pillar_l/numubar/ingrid_%08d_%04d_anas1.root"   $DETRES_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE6_1=`printf  "%s/dial_%02d/pillar_l/numu/ingrid_%08d_%04d_anas1.root"      $DETRES_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE1=`printf   "%s/dial_%02d/wallbg/numubar/ingrid_%08d_%04d_plot.root"      $DETRES_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE2=`printf   "%s/dial_%02d/wallbg/numu/ingrid_%08d_%04d_plot.root"         $DETRES_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE3_0=`printf "%s/dial_%02d/ceiling/numubar/ingrid_%08d_%04d_plot.root"     $DETRES_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE4_0=`printf "%s/dial_%02d/ceiling/numu/ingrid_%08d_%04d_plot.root"        $DETRES_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE3_1=`printf "%s/dial_%02d/floor/numubar/ingrid_%08d_%04d_plot.root"       $DETRES_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE4_1=`printf "%s/dial_%02d/floor/numu/ingrid_%08d_%04d_plot.root"          $DETRES_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE5_0=`printf "%s/dial_%02d/pillar_r/numubar/ingrid_%08d_%04d_plot.root"    $DETRES_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE6_0=`printf "%s/dial_%02d/pillar_r/numu/ingrid_%08d_%04d_plot.root"       $DETRES_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE5_1=`printf "%s/dial_%02d/pillar_l/numubar/ingrid_%08d_%04d_plot.root"    $DETRES_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE6_1=`printf "%s/dial_%02d/pillar_l/numu/ingrid_%08d_%04d_plot.root"       $DETRES_DIR $dial $RUNID $SRUNID`

    if test $dial -lt 27
    then
      #Three Dim Recon
      cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x -n %d" ${ANALYSIS_DIR} ${RECONFILE1}   ${DET_ANAFILE1}   ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 10 -p" ${ANALYSIS_DIR} ${DET_ANAFILE1}   ${DET_PLOTFILE1}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "rm -f %s" ${DET_ANAFILE1}`
      echo -e "---------\n ${cmd}\n"; eval $cmd

      cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x -n %d" ${ANALYSIS_DIR} ${RECONFILE2}   ${DET_ANAFILE2}   ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 10 -p" ${ANALYSIS_DIR} ${DET_ANAFILE2}   ${DET_PLOTFILE2}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "rm -f %s" ${DET_ANAFILE2}`
      echo -e "---------\n ${cmd}\n"; eval $cmd


      cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x -n %d" ${ANALYSIS_DIR} ${RECONFILE3_0} ${DET_ANAFILE3_0} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 11 -p" ${ANALYSIS_DIR} ${DET_ANAFILE3_0} ${DET_PLOTFILE3_0}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "rm -f %s" ${DET_ANAFILE3_0}`
      echo -e "---------\n ${cmd}\n"; eval $cmd


      cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x -n %d" ${ANALYSIS_DIR} ${RECONFILE4_0} ${DET_ANAFILE4_0} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 11 -p" ${ANALYSIS_DIR} ${DET_ANAFILE4_0} ${DET_PLOTFILE4_0}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "rm -f %s" ${DET_ANAFILE4_0}`
      echo -e "---------\n ${cmd}\n"; eval $cmd

      cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x -n %d" ${ANALYSIS_DIR} ${RECONFILE3_1} ${DET_ANAFILE3_1} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 12 -p" ${ANALYSIS_DIR} ${DET_ANAFILE3_1} ${DET_PLOTFILE3_1}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "rm -f %s" ${DET_ANAFILE3_1}`
      echo -e "---------\n ${cmd}\n"; eval $cmd

      cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x -n %d" ${ANALYSIS_DIR} ${RECONFILE4_1} ${DET_ANAFILE4_1} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 12 -p" ${ANALYSIS_DIR} ${DET_ANAFILE4_1} ${DET_PLOTFILE4_1}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "rm -f %s" ${DET_ANAFILE4_1}`
      echo -e "---------\n ${cmd}\n"; eval $cmd

      cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x -n %d" ${ANALYSIS_DIR} ${RECONFILE5_0} ${DET_ANAFILE5_0} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 13 -p" ${ANALYSIS_DIR} ${DET_ANAFILE5_0} ${DET_PLOTFILE5_0}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "rm -f %s" ${DET_ANAFILE5_0}`
      echo -e "---------\n ${cmd}\n"; eval $cmd

      cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x -n %d" ${ANALYSIS_DIR} ${RECONFILE6_0} ${DET_ANAFILE6_0} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 13 -p" ${ANALYSIS_DIR} ${DET_ANAFILE6_0} ${DET_PLOTFILE6_0}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "rm -f %s" ${DET_ANAFILE6_0}`
      echo -e "---------\n ${cmd}\n"; eval $cmd

      cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x -n %d" ${ANALYSIS_DIR} ${RECONFILE5_1} ${DET_ANAFILE5_1} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 14 -p" ${ANALYSIS_DIR} ${DET_ANAFILE5_1} ${DET_PLOTFILE5_1}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "rm -f %s" ${DET_ANAFILE5_1}`
      echo -e "---------\n ${cmd}\n"; eval $cmd

      cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x -n %d" ${ANALYSIS_DIR} ${RECONFILE6_1} ${DET_ANAFILE6_1} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 14 -p" ${ANALYSIS_DIR} ${DET_ANAFILE6_1} ${DET_PLOTFILE6_1}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "rm -f %s" ${DET_ANAFILE6_1}`
      echo -e "---------\n ${cmd}\n"; eval $cmd

    else
      ##Plot 
      cmd=`printf "%s/B2Plot -i %s -o %s -m 10 -n %d -p" ${ANALYSIS_DIR} ${ANAFILE1}   ${DET_PLOTFILE1}   ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 10 -n %d -p" ${ANALYSIS_DIR} ${ANAFILE2}   ${DET_PLOTFILE2}   ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 11 -n %d -p" ${ANALYSIS_DIR} ${ANAFILE3_0} ${DET_PLOTFILE3_0} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 11 -n %d -p" ${ANALYSIS_DIR} ${ANAFILE4_0} ${DET_PLOTFILE4_0} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 12 -n %d -p" ${ANALYSIS_DIR} ${ANAFILE3_1} ${DET_PLOTFILE3_1} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 12 -n %d -p" ${ANALYSIS_DIR} ${ANAFILE4_1} ${DET_PLOTFILE4_1} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 13 -n %d -p" ${ANALYSIS_DIR} ${ANAFILE5_0} ${DET_PLOTFILE5_0} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 13 -n %d -p" ${ANALYSIS_DIR} ${ANAFILE6_0} ${DET_PLOTFILE6_0} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 14 -n %d -p" ${ANALYSIS_DIR} ${ANAFILE5_1} ${DET_PLOTFILE5_1} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 14 -n %d -p" ${ANALYSIS_DIR} ${ANAFILE6_1} ${DET_PLOTFILE6_1} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
    fi

 done
fi
