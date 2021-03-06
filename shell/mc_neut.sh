#!/bin/sh

MCSIM=0
TWODRECON=0
THRDRECON=0
PLOT_ALL=0
PLOT_MODE=0
PLOT_REWEI=0
DETRES=0
DETRES2=0
THRESHOLD=0
TRKEFF=0
PECHECK=0
HITEFF=0
SELECTION=1

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
NEUT_DIR="/hsm/nu/wagasci/neut/rhc/5.3.2.1"
OUTPUT_DIR="/home/t2k/nchikuma/b2_data2/mc_neut"
OUTPUT_DIR2_0="/home/t2k/nchikuma/b2_data2/reweight/rfg"
OUTPUT_DIR2_1="/home/t2k/nchikuma/b2_data2/reweight/rpa"
REWEIGHT_DIR="/hsm/t2k/t2k_JB/t2k_beam/kenichi/reweight/rfg"
REWEIGHT_DIR2="/hsm/t2k/t2k_JB/t2k_beam/kenichi/reweight/rfg_rpa"
DETRES_DIR="/home/t2k/nchikuma/b2_data3/sysE_detRes/data"
DETRES2_DIR="/home/t2k/nchikuma/b2_data3/sysE_detRes2/data"
THRES_DIR="/home/t2k/nchikuma/b2_data2/pe_thres/data"


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
REWEIGHTFILE5_0=`printf "%s/fe/numubar/ingrid_%08d_%04d_reweight.root"  $REWEIGHT_DIR $RUNID $SRUNID`
REWEIGHTFILE6_0=`printf "%s/fe/numu/ingrid_%08d_%04d_reweight.root"     $REWEIGHT_DIR $RUNID $SRUNID`
REWEIGHTFILE1_1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_reweight.root" $REWEIGHT_DIR2 $RUNID $SRUNID`
REWEIGHTFILE2_1=`printf "%s/h2o/numu/ingrid_%08d_%04d_reweight.root"    $REWEIGHT_DIR2 $RUNID $SRUNID`
REWEIGHTFILE3_1=`printf "%s/ch/numubar/ingrid_%08d_%04d_reweight.root"  $REWEIGHT_DIR2 $RUNID $SRUNID`
REWEIGHTFILE4_1=`printf "%s/ch/numu/ingrid_%08d_%04d_reweight.root"     $REWEIGHT_DIR2 $RUNID $SRUNID`
REWEIGHTFILE5_1=`printf "%s/fe/numubar/ingrid_%08d_%04d_reweight.root"  $REWEIGHT_DIR2 $RUNID $SRUNID`
REWEIGHTFILE6_1=`printf "%s/fe/numu/ingrid_%08d_%04d_reweight.root"     $REWEIGHT_DIR2 $RUNID $SRUNID`

TRKEFFFILE1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_trkeff.root" $OUTPUT_DIR $RUNID $SRUNID`
TRKEFFFILE2=`printf "%s/h2o/numu/ingrid_%08d_%04d_trkeff.root"    $OUTPUT_DIR $RUNID $SRUNID`
TRKEFFFILE3=`printf "%s/ch/numubar/ingrid_%08d_%04d_trkeff.root"  $OUTPUT_DIR $RUNID $SRUNID`
TRKEFFFILE4=`printf "%s/ch/numu/ingrid_%08d_%04d_trkeff.root"     $OUTPUT_DIR $RUNID $SRUNID`
TRKEFFFILE5=`printf "%s/fe/numubar/ingrid_%08d_%04d_trkeff.root"  $OUTPUT_DIR $RUNID $SRUNID`
TRKEFFFILE6=`printf "%s/fe/numu/ingrid_%08d_%04d_trkeff.root"     $OUTPUT_DIR $RUNID $SRUNID`

PEFILE1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_pe.root" $OUTPUT_DIR $RUNID $SRUNID`
PEFILE2=`printf "%s/h2o/numu/ingrid_%08d_%04d_pe.root"    $OUTPUT_DIR $RUNID $SRUNID`
PEFILE3=`printf "%s/ch/numubar/ingrid_%08d_%04d_pe.root"  $OUTPUT_DIR $RUNID $SRUNID`
PEFILE4=`printf "%s/ch/numu/ingrid_%08d_%04d_pe.root"     $OUTPUT_DIR $RUNID $SRUNID`
PEFILE5=`printf "%s/fe/numubar/ingrid_%08d_%04d_pe.root"  $OUTPUT_DIR $RUNID $SRUNID`
PEFILE6=`printf "%s/fe/numu/ingrid_%08d_%04d_pe.root"     $OUTPUT_DIR $RUNID $SRUNID`

HITEFFFILE1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_hiteff.root" $OUTPUT_DIR $RUNID $SRUNID`
HITEFFFILE2=`printf "%s/h2o/numu/ingrid_%08d_%04d_hiteff.root"    $OUTPUT_DIR $RUNID $SRUNID`
HITEFFFILE3=`printf "%s/ch/numubar/ingrid_%08d_%04d_hiteff.root"  $OUTPUT_DIR $RUNID $SRUNID`
HITEFFFILE4=`printf "%s/ch/numu/ingrid_%08d_%04d_hiteff.root"     $OUTPUT_DIR $RUNID $SRUNID`
HITEFFFILE5=`printf "%s/fe/numubar/ingrid_%08d_%04d_hiteff.root"  $OUTPUT_DIR $RUNID $SRUNID`
HITEFFFILE6=`printf "%s/fe/numu/ingrid_%08d_%04d_hiteff.root"     $OUTPUT_DIR $RUNID $SRUNID`



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
#  cmd="rm -f $ALLFILE5"
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd="rm -f $ALLFILE6"
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#
#  #MC Job
#  cmd=`printf "%s/bin/Linux-g++/b2mc -m7 -f2 -i %s -o %s" $MC_DIR $NEUTFILE1 $OUTPUTFILE1`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd=`printf "%s/bin/Linux-g++/b2mc -m7 -f1 -i %s -o %s" $MC_DIR $NEUTFILE2 $OUTPUTFILE2`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd=`printf "%s/bin/Linux-g++/b2mc -m8 -f2 -i %s -o %s" $MC_DIR $NEUTFILE3 $OUTPUTFILE3`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd=`printf "%s/bin/Linux-g++/b2mc -m8 -f1 -i %s -o %s" $MC_DIR $NEUTFILE4 $OUTPUTFILE4`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd=`printf "%s/bin/Linux-g++/b2mc -m9 -f2 -i %s -o %s" $MC_DIR $NEUTFILE5 $OUTPUTFILE5`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#  cmd=`printf "%s/bin/Linux-g++/b2mc -m9 -f1 -i %s -o %s" $MC_DIR $NEUTFILE6 $OUTPUTFILE6`
#  echo -e "---------\n ${cmd}\n"; eval $cmd

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
  cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE5} ${RECONFILE5}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TwoDimRecon -f %s -o %s" ${ANALYSIS_DIR} ${OUTPUTFILE6} ${RECONFILE6}`
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
  cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${RECONFILE5} ${ANAFILE5} ${OPT_NOT3DRRK}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${RECONFILE6} ${ANAFILE6} ${OPT_NOT3DRRK}`
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
  cmd=`printf "%s/B2Plot -i %s -o %s -m 9 -t %s %s %s" ${ANALYSIS_DIR} ${ANAFILE5} ${PLOTFILE5} ${TUNEFILE} ${OPT_RATIO} ${OPT_ONLYNSEL}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2Plot -i %s -o %s -m 9 -t %s %s %s" ${ANALYSIS_DIR} ${ANAFILE6} ${PLOTFILE6} ${TUNEFILE} ${OPT_RATIO} ${OPT_ONLYNSEL}`
  echo -e "---------\n ${cmd}\n"; eval $cmd

fi

if test $PLOT_MODE -eq 1
then
  echo "plot for each interaction mode"

  #Plot 
  for i in `seq 0 11`
  do
    inttype=$i
    PLOTFILE1_0=`printf "%s/h2o/numubar/ingrid_%08d_%04d_plot%d.root" $OUTPUT_DIR $RUNID $SRUNID $inttype`
    PLOTFILE3_0=`printf "%s/ch/numubar/ingrid_%08d_%04d_plot%d.root"  $OUTPUT_DIR $RUNID $SRUNID $inttype`
    PLOTFILE5_0=`printf "%s/fe/numubar/ingrid_%08d_%04d_plot%d.root"  $OUTPUT_DIR $RUNID $SRUNID $inttype`
    cmd=`printf "%s/B2Plot -i %s -o %s -m 7 -t %s -c %d %s" ${ANALYSIS_DIR} ${ANAFILE1} ${PLOTFILE1_0} ${TUNEFILE} ${inttype} ${OPT_RATIO}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/B2Plot -i %s -o %s -m 8 -t %s -c %d %s" ${ANALYSIS_DIR} ${ANAFILE3} ${PLOTFILE3_0} ${TUNEFILE} ${inttype} ${OPT_RATIO}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    #cmd=`printf "%s/B2Plot -i %s -o %s -m 9 -t %s -c %d %s" ${ANALYSIS_DIR} ${ANAFILE5} ${PLOTFILE5_0} ${TUNEFILE} ${inttype} ${OPT_RATIO}`
    #echo -e "---------\n ${cmd}\n"; eval $cmd
  done

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
      REWEIGHTFILE5=${REWEIGHTFILE5_0}
      REWEIGHTFILE6=${REWEIGHTFILE6_0}
      OUTPUT_DIR2=${OUTPUT_DIR2_0}
    else 
      REWEIGHTFILE1=${REWEIGHTFILE1_1}
      REWEIGHTFILE2=${REWEIGHTFILE2_1}
      REWEIGHTFILE3=${REWEIGHTFILE3_1}
      REWEIGHTFILE4=${REWEIGHTFILE4_1}
      REWEIGHTFILE5=${REWEIGHTFILE5_1}
      REWEIGHTFILE6=${REWEIGHTFILE6_1}
      OUTPUT_DIR2=${OUTPUT_DIR2_1}
    fi
  
    for j in `seq 1 22`
    do
      for k in 2 3 4
      do
        if [ $j -ne 1 -a $k -eq 3 ]
        then
          continue
        fi  
        dial=`expr $j '*' 7 + $k`

        if test $ONLY_NOMINAL -eq 1
        then
          if test $dial -ne 10
          then
            continue
          fi
        fi


        echo $dial

        OPTPLOT=""
        if test $dial -eq 10
        then
          OPTPLOT=$OPT_ONLYNSEL
        else
          OPTPLOT="-p"
        fi
    
        PLOTFILE1=`printf "%s/dial_%03d/h2o/numubar/ingrid_%08d_%04d_plot.root" $OUTPUT_DIR2 $dial $RUNID $SRUNID`
        PLOTFILE2=`printf "%s/dial_%03d/h2o/numu/ingrid_%08d_%04d_plot.root"    $OUTPUT_DIR2 $dial $RUNID $SRUNID`
        PLOTFILE3=`printf "%s/dial_%03d/ch/numubar/ingrid_%08d_%04d_plot.root"  $OUTPUT_DIR2 $dial $RUNID $SRUNID`
        PLOTFILE4=`printf "%s/dial_%03d/ch/numu/ingrid_%08d_%04d_plot.root"     $OUTPUT_DIR2 $dial $RUNID $SRUNID`

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
    
        if test $dial -eq 10
        then
          for i in `seq 0 11`
          do
            inttype=$i
            PLOTFILE1_0=`printf "%s/dial_%03d/h2o/numubar/ingrid_%08d_%04d_plot%d.root" $OUTPUT_DIR2 $dial $RUNID $SRUNID $inttype`
            PLOTFILE3_0=`printf "%s/dial_%03d/ch/numubar/ingrid_%08d_%04d_plot%d.root"  $OUTPUT_DIR2 $dial $RUNID $SRUNID $inttype`

            cmd=`printf "rm -f %s" ${PLOTFILE1_0}` echo -e "---------\n ${cmd}\n"; eval $cmd
            cmd=`printf "rm -f %s" ${PLOTFILE3_0}` echo -e "---------\n ${cmd}\n"; eval $cmd

            cmd=`printf "%s/B2Plot -i %s -o %s -m 7 -t %s -c %d -w %d -r %s %s %s" ${ANALYSIS_DIR} ${ANAFILE1} ${PLOTFILE1_0} ${TUNEFILE} ${inttype} ${dial} ${REWEIGHTFILE1} ${OPT_RATIO} ${OPTPLOT}`
            echo -e "---------\n ${cmd}\n"; eval $cmd
            cmd=`printf "%s/B2Plot -i %s -o %s -m 8 -t %s -c %d -w %d -r %s %s %s" ${ANALYSIS_DIR} ${ANAFILE3} ${PLOTFILE3_0} ${TUNEFILE} ${inttype} ${dial} ${REWEIGHTFILE3} ${OPT_RATIO} ${OPTPLOT}`
            echo -e "---------\n ${cmd}\n"; eval $cmd
          done
        fi
      done
    done
  done

fi


if test $TRKEFF -eq 1
then

  cmd=`printf "%s/TrackEfficiency2 -i %s -o %s -m 7" ${ANALYSIS_DIR} ${ANAFILE1} ${TRKEFFFILE1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TrackEfficiency2 -i %s -o %s -m 7" ${ANALYSIS_DIR} ${ANAFILE2} ${TRKEFFFILE2}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TrackEfficiency2 -i %s -o %s -m 8" ${ANALYSIS_DIR} ${ANAFILE3} ${TRKEFFFILE3}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TrackEfficiency2 -i %s -o %s -m 8" ${ANALYSIS_DIR} ${ANAFILE4} ${TRKEFFFILE4}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TrackEfficiency2 -i %s -o %s -m 9" ${ANALYSIS_DIR} ${ANAFILE5} ${TRKEFFFILE5}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/TrackEfficiency2 -i %s -o %s -m 9" ${ANALYSIS_DIR} ${ANAFILE6} ${TRKEFFFILE6}`
  echo -e "---------\n ${cmd}\n"; eval $cmd

fi

if test $PECHECK -eq 1
then

  cmd=`printf "%s/B2pe -i %s -o %s -m 7" ${ANALYSIS_DIR} ${ANAFILE1} ${PEFILE1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2pe -i %s -o %s -m 7" ${ANALYSIS_DIR} ${ANAFILE2} ${PEFILE2}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2pe -i %s -o %s -m 8" ${ANALYSIS_DIR} ${ANAFILE3} ${PEFILE3}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2pe -i %s -o %s -m 8" ${ANALYSIS_DIR} ${ANAFILE4} ${PEFILE4}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2pe -i %s -o %s -m 9" ${ANALYSIS_DIR} ${ANAFILE5} ${PEFILE5}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2pe -i %s -o %s -m 9" ${ANALYSIS_DIR} ${ANAFILE6} ${PEFILE6}`
  echo -e "---------\n ${cmd}\n"; eval $cmd

fi

if test $HITEFF -eq 1
then

  cmd=`printf "%s/HitEff2 -i %s -o %s -m 7" ${ANALYSIS_DIR} ${ANAFILE1} ${HITEFFFILE1}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/HitEff2 -i %s -o %s -m 7" ${ANALYSIS_DIR} ${ANAFILE2} ${HITEFFFILE2}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/HitEff2 -i %s -o %s -m 8" ${ANALYSIS_DIR} ${ANAFILE3} ${HITEFFFILE3}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/HitEff2 -i %s -o %s -m 8" ${ANALYSIS_DIR} ${ANAFILE4} ${HITEFFFILE4}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/HitEff2 -i %s -o %s -m 9" ${ANALYSIS_DIR} ${ANAFILE5} ${HITEFFFILE5}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/HitEff2 -i %s -o %s -m 9" ${ANALYSIS_DIR} ${ANAFILE6} ${HITEFFFILE6}`
  echo -e "---------\n ${cmd}\n"; eval $cmd

fi


if test $DETRES -eq 1
then
  echo "Plot for sys. err. from detectro response"


  #for dial in `seq 0 39`
  for dial in `seq 19 26`
  do
    
    DET_ANAFILE1=`printf  "%s/dial_%02d/h2o/numubar/ingrid_%08d_%04d_anas1.root" $DETRES_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE2=`printf  "%s/dial_%02d/h2o/numu/ingrid_%08d_%04d_anas1.root"    $DETRES_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE3=`printf  "%s/dial_%02d/ch/numubar/ingrid_%08d_%04d_anas1.root"  $DETRES_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE4=`printf  "%s/dial_%02d/ch/numu/ingrid_%08d_%04d_anas1.root"     $DETRES_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE5=`printf  "%s/dial_%02d/fe/numubar/ingrid_%08d_%04d_anas1.root"  $DETRES_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE6=`printf  "%s/dial_%02d/fe/numu/ingrid_%08d_%04d_anas1.root"     $DETRES_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE1=`printf "%s/dial_%02d/h2o/numubar/ingrid_%08d_%04d_plot.root"  $DETRES_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE2=`printf "%s/dial_%02d/h2o/numu/ingrid_%08d_%04d_plot.root"     $DETRES_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE3=`printf "%s/dial_%02d/ch/numubar/ingrid_%08d_%04d_plot.root"   $DETRES_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE4=`printf "%s/dial_%02d/ch/numu/ingrid_%08d_%04d_plot.root"      $DETRES_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE5=`printf "%s/dial_%02d/fe/numubar/ingrid_%08d_%04d_plot.root"   $DETRES_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE6=`printf "%s/dial_%02d/fe/numu/ingrid_%08d_%04d_plot.root"      $DETRES_DIR $dial $RUNID $SRUNID`

    if test $dial -lt 27
    then
      #Three Dim Recon
      cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x -n %d %s" ${ANALYSIS_DIR} ${RECONFILE1} ${DET_ANAFILE1} ${dial} ${OPT_NOT3DRRK}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 7 -t %s -p" ${ANALYSIS_DIR} ${DET_ANAFILE1} ${DET_PLOTFILE1} ${TUNEFILE}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "rm -f %s" ${DET_ANAFILE1}`
      echo -e "---------\n ${cmd}\n"; eval $cmd

      cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x -n %d %s" ${ANALYSIS_DIR} ${RECONFILE2} ${DET_ANAFILE2} ${dial} ${OPT_NOT3DRRK}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 7 -t %s -p" ${ANALYSIS_DIR} ${DET_ANAFILE2} ${DET_PLOTFILE2} ${TUNEFILE}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "rm -f %s" ${DET_ANAFILE2}`
      echo -e "---------\n ${cmd}\n"; eval $cmd

      cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x -n %d %s" ${ANALYSIS_DIR} ${RECONFILE3} ${DET_ANAFILE3} ${dial} ${OPT_NOT3DRRK}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 8 -t %s -p" ${ANALYSIS_DIR} ${DET_ANAFILE3} ${DET_PLOTFILE3} ${TUNEFILE}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "rm -f %s" ${DET_ANAFILE3}`
      echo -e "---------\n ${cmd}\n"; eval $cmd

      cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x -n %d %s" ${ANALYSIS_DIR} ${RECONFILE4} ${DET_ANAFILE4} ${dial} ${OPT_NOT3DRRK}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 8 -t %s -p" ${ANALYSIS_DIR} ${DET_ANAFILE4} ${DET_PLOTFILE4} ${TUNEFILE}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "rm -f %s" ${DET_ANAFILE4}`
      echo -e "---------\n ${cmd}\n"; eval $cmd


    else
      ##Plot 
      cmd=`printf "%s/B2Plot -i %s -o %s -m 7 -t %s -n %d -p" ${ANALYSIS_DIR} ${ANAFILE1} ${DET_PLOTFILE1} ${TUNEFILE} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 7 -t %s -n %d -p" ${ANALYSIS_DIR} ${ANAFILE2} ${DET_PLOTFILE2} ${TUNEFILE} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 8 -t %s -n %d -p" ${ANALYSIS_DIR} ${ANAFILE3} ${DET_PLOTFILE3} ${TUNEFILE} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/B2Plot -i %s -o %s -m 8 -t %s -n %d -p" ${ANALYSIS_DIR} ${ANAFILE4} ${DET_PLOTFILE4} ${TUNEFILE} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
    fi
    
  done

fi


if test $DETRES2 -eq 1
then
  echo "Plot for sys. err. from mppc noise & scinti crosstalk"


  for dial in `seq 0 6`
  do

    DET_RECONFILE1=`printf  "%s/dial_%02d/h2o/numubar/ingrid_%08d_%04d_recon.root" $DETRES2_DIR $dial $RUNID $SRUNID`
    DET_RECONFILE2=`printf  "%s/dial_%02d/h2o/numu/ingrid_%08d_%04d_recon.root"    $DETRES2_DIR $dial $RUNID $SRUNID`
    DET_RECONFILE3=`printf  "%s/dial_%02d/ch/numubar/ingrid_%08d_%04d_recon.root"  $DETRES2_DIR $dial $RUNID $SRUNID`
    DET_RECONFILE4=`printf  "%s/dial_%02d/ch/numu/ingrid_%08d_%04d_recon.root"     $DETRES2_DIR $dial $RUNID $SRUNID`
    DET_RECONFILE5=`printf  "%s/dial_%02d/fe/numubar/ingrid_%08d_%04d_recon.root"  $DETRES2_DIR $dial $RUNID $SRUNID`
    DET_RECONFILE6=`printf  "%s/dial_%02d/fe/numu/ingrid_%08d_%04d_recon.root"     $DETRES2_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE1=`printf    "%s/dial_%02d/h2o/numubar/ingrid_%08d_%04d_anas1.root" $DETRES2_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE2=`printf    "%s/dial_%02d/h2o/numu/ingrid_%08d_%04d_anas1.root"    $DETRES2_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE3=`printf    "%s/dial_%02d/ch/numubar/ingrid_%08d_%04d_anas1.root"  $DETRES2_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE4=`printf    "%s/dial_%02d/ch/numu/ingrid_%08d_%04d_anas1.root"     $DETRES2_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE5=`printf    "%s/dial_%02d/fe/numubar/ingrid_%08d_%04d_anas1.root"  $DETRES2_DIR $dial $RUNID $SRUNID`
    DET_ANAFILE6=`printf    "%s/dial_%02d/fe/numu/ingrid_%08d_%04d_anas1.root"     $DETRES2_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE1=`printf   "%s/dial_%02d/h2o/numubar/ingrid_%08d_%04d_plot.root"  $DETRES2_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE2=`printf   "%s/dial_%02d/h2o/numu/ingrid_%08d_%04d_plot.root"     $DETRES2_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE3=`printf   "%s/dial_%02d/ch/numubar/ingrid_%08d_%04d_plot.root"   $DETRES2_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE4=`printf   "%s/dial_%02d/ch/numu/ingrid_%08d_%04d_plot.root"      $DETRES2_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE5=`printf   "%s/dial_%02d/fe/numubar/ingrid_%08d_%04d_plot.root"   $DETRES2_DIR $dial $RUNID $SRUNID`
    DET_PLOTFILE6=`printf   "%s/dial_%02d/fe/numu/ingrid_%08d_%04d_plot.root"      $DETRES2_DIR $dial $RUNID $SRUNID`

    if test $dial -lt 6
    then
      #Recon
      cmd=`printf "%s/TwoDimRecon -f %s -o %s -n %d" ${ANALYSIS_DIR} ${OUTPUTFILE1} ${DET_RECONFILE1} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/TwoDimRecon -f %s -o %s -n %d" ${ANALYSIS_DIR} ${OUTPUTFILE2} ${DET_RECONFILE2} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/TwoDimRecon -f %s -o %s -n %d" ${ANALYSIS_DIR} ${OUTPUTFILE3} ${DET_RECONFILE3} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/TwoDimRecon -f %s -o %s -n %d" ${ANALYSIS_DIR} ${OUTPUTFILE4} ${DET_RECONFILE4} ${dial}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
    else
      #Recon
      cmd=`printf "%s/TwoDimRecon -f %s -o %s -x" ${ANALYSIS_DIR} ${OUTPUTFILE1} ${DET_RECONFILE1}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/TwoDimRecon -f %s -o %s -x" ${ANALYSIS_DIR} ${OUTPUTFILE2} ${DET_RECONFILE2}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/TwoDimRecon -f %s -o %s -x" ${ANALYSIS_DIR} ${OUTPUTFILE3} ${DET_RECONFILE3}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
      cmd=`printf "%s/TwoDimRecon -f %s -o %s -x" ${ANALYSIS_DIR} ${OUTPUTFILE4} ${DET_RECONFILE4}`
      echo -e "---------\n ${cmd}\n"; eval $cmd
    fi

    #Three Dim Recon
    cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${DET_RECONFILE1} ${DET_ANAFILE1} ${OPT_NOT3DRRK}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${DET_RECONFILE2} ${DET_ANAFILE2} ${OPT_NOT3DRRK}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${DET_RECONFILE3} ${DET_ANAFILE3} ${OPT_NOT3DRRK}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${DET_RECONFILE4} ${DET_ANAFILE4} ${OPT_NOT3DRRK}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    ##Plot 
    cmd=`printf "%s/B2Plot -i %s -o %s -m 7 -t %s" ${ANALYSIS_DIR} ${DET_ANAFILE1} ${DET_PLOTFILE1} ${TUNEFILE}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/B2Plot -i %s -o %s -m 7 -t %s" ${ANALYSIS_DIR} ${DET_ANAFILE2} ${DET_PLOTFILE2} ${TUNEFILE}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/B2Plot -i %s -o %s -m 8 -t %s" ${ANALYSIS_DIR} ${DET_ANAFILE3} ${DET_PLOTFILE3} ${TUNEFILE}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/B2Plot -i %s -o %s -m 8 -t %s" ${ANALYSIS_DIR} ${DET_ANAFILE4} ${DET_PLOTFILE4} ${TUNEFILE}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    
  done

fi


if test $THRESHOLD -eq 1
then
  echo "Plot for sys. err. from threshold pe"


  for dial in `seq 0 5`
  do

    THRES_RECONFILE1=`printf  "%s/dial_%d/h2o/numubar/ingrid_%08d_%04d_recon.root" $THRES_DIR $dial $RUNID $SRUNID`
    THRES_RECONFILE2=`printf  "%s/dial_%d/h2o/numu/ingrid_%08d_%04d_recon.root"    $THRES_DIR $dial $RUNID $SRUNID`
    THRES_RECONFILE3=`printf  "%s/dial_%d/ch/numubar/ingrid_%08d_%04d_recon.root"  $THRES_DIR $dial $RUNID $SRUNID`
    THRES_RECONFILE4=`printf  "%s/dial_%d/ch/numu/ingrid_%08d_%04d_recon.root"     $THRES_DIR $dial $RUNID $SRUNID`
    THRES_ANAFILE1=`printf    "%s/dial_%d/h2o/numubar/ingrid_%08d_%04d_anas1.root" $THRES_DIR $dial $RUNID $SRUNID`
    THRES_ANAFILE2=`printf    "%s/dial_%d/h2o/numu/ingrid_%08d_%04d_anas1.root"    $THRES_DIR $dial $RUNID $SRUNID`
    THRES_ANAFILE3=`printf    "%s/dial_%d/ch/numubar/ingrid_%08d_%04d_anas1.root"  $THRES_DIR $dial $RUNID $SRUNID`
    THRES_ANAFILE4=`printf    "%s/dial_%d/ch/numu/ingrid_%08d_%04d_anas1.root"     $THRES_DIR $dial $RUNID $SRUNID`
    THRES_PLOTFILE1=`printf   "%s/dial_%d/h2o/numubar/ingrid_%08d_%04d_plot.root"  $THRES_DIR $dial $RUNID $SRUNID`
    THRES_PLOTFILE2=`printf   "%s/dial_%d/h2o/numu/ingrid_%08d_%04d_plot.root"     $THRES_DIR $dial $RUNID $SRUNID`
    THRES_PLOTFILE3=`printf   "%s/dial_%d/ch/numubar/ingrid_%08d_%04d_plot.root"   $THRES_DIR $dial $RUNID $SRUNID`
    THRES_PLOTFILE4=`printf   "%s/dial_%d/ch/numu/ingrid_%08d_%04d_plot.root"      $THRES_DIR $dial $RUNID $SRUNID`

    #Recon
    cmd=`printf "%s/TwoDimRecon -f %s -o %s -t %d" ${ANALYSIS_DIR} ${OUTPUTFILE1} ${THRES_RECONFILE1} ${dial}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${THRES_RECONFILE1} ${THRES_ANAFILE1} ${OPT_NOT3DRRK}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "rm -f %s" ${THRES_RECONFILE1}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/B2Plot -i %s -o %s -m 7 -t %s %s %s" ${ANALYSIS_DIR} ${THRES_ANAFILE1} ${THRES_PLOTFILE1} ${TUNEFILE} ${OPT_RATIO} ${OPT_ONLYNSEL}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "rm -f %s" ${THRES_ANAFILE1}`
    echo -e "---------\n ${cmd}\n"; eval $cmd

    cmd=`printf "%s/TwoDimRecon -f %s -o %s -t %d" ${ANALYSIS_DIR} ${OUTPUTFILE2} ${THRES_RECONFILE2} ${dial}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${THRES_RECONFILE2} ${THRES_ANAFILE2} ${OPT_NOT3DRRK}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "rm -f %s" ${THRES_RECONFILE2}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/B2Plot -i %s -o %s -m 7 -t %s %s %s" ${ANALYSIS_DIR} ${THRES_ANAFILE2} ${THRES_PLOTFILE2} ${TUNEFILE} ${OPT_RATIO} ${OPT_ONLYNSEL}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "rm -f %s" ${THRES_ANAFILE2}`
    echo -e "---------\n ${cmd}\n"; eval $cmd

    cmd=`printf "%s/TwoDimRecon -f %s -o %s -t %d" ${ANALYSIS_DIR} ${OUTPUTFILE3} ${THRES_RECONFILE3} ${dial}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${THRES_RECONFILE3} ${THRES_ANAFILE3} ${OPT_NOT3DRRK}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "rm -f %s" ${THRES_RECONFILE3}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/B2Plot -i %s -o %s -m 8 -t %s %s %s" ${ANALYSIS_DIR} ${THRES_ANAFILE3} ${THRES_PLOTFILE3} ${TUNEFILE} ${OPT_RATIO} ${OPT_ONLYNSEL}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "rm -f %s" ${THRES_ANAFILE3}`
    echo -e "---------\n ${cmd}\n"; eval $cmd

    cmd=`printf "%s/TwoDimRecon -f %s -o %s -t %d" ${ANALYSIS_DIR} ${OUTPUTFILE4} ${THRES_RECONFILE4} ${dial}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/ThreeDimRecon -f %s -o %s -m1 -v -x %s" ${ANALYSIS_DIR} ${THRES_RECONFILE4} ${THRES_ANAFILE4} ${OPT_NOT3DRRK}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "rm -f %s" ${THRES_RECONFILE4}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/B2Plot -i %s -o %s -m 8 -t %s %s %s" ${ANALYSIS_DIR} ${THRES_ANAFILE4} ${THRES_PLOTFILE4} ${TUNEFILE} ${OPT_RATIO} ${OPT_ONLYNSEL}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
     cmd=`printf "rm -f %s" ${THRES_ANAFILE4}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
   
  done

fi


if test $SELECTION -eq 1
then
  echo "Plot with t2kReweight"


  for l in `seq 1 1`
  do
    if test $l -eq 0
    then
      REWEIGHTFILE1=${REWEIGHTFILE1_0}
      REWEIGHTFILE2=${REWEIGHTFILE2_0}
      REWEIGHTFILE3=${REWEIGHTFILE3_0}
      REWEIGHTFILE4=${REWEIGHTFILE4_0}
      REWEIGHTFILE5=${REWEIGHTFILE5_0}
      REWEIGHTFILE6=${REWEIGHTFILE6_0}
      OUTPUT_DIR2=${OUTPUT_DIR2_0}
    else 
      REWEIGHTFILE1=${REWEIGHTFILE1_1}
      REWEIGHTFILE2=${REWEIGHTFILE2_1}
      REWEIGHTFILE3=${REWEIGHTFILE3_1}
      REWEIGHTFILE4=${REWEIGHTFILE4_1}
      REWEIGHTFILE5=${REWEIGHTFILE5_1}
      REWEIGHTFILE6=${REWEIGHTFILE6_1}
      OUTPUT_DIR2=${OUTPUT_DIR2_1}
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
    
    PLOTFILE1=`printf "%s/h2o/numubar/ingrid_%08d_%04d_selec.root" $OUTPUT_DIR $RUNID $SRUNID`
    PLOTFILE2=`printf "%s/h2o/numu/ingrid_%08d_%04d_selec.root"    $OUTPUT_DIR $RUNID $SRUNID`
    PLOTFILE3=`printf "%s/ch/numubar/ingrid_%08d_%04d_selec.root"  $OUTPUT_DIR $RUNID $SRUNID`
    PLOTFILE4=`printf "%s/ch/numu/ingrid_%08d_%04d_selec.root"     $OUTPUT_DIR $RUNID $SRUNID`


    cmd=`printf "%s/Selection -i %s -o %s -m 7 -t %s -w %d -r %s %s %s" ${ANALYSIS_DIR} ${ANAFILE1} ${PLOTFILE1} ${TUNEFILE} ${dial} ${REWEIGHTFILE1} ${OPT_RATIO} ${OPTPLOT}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/Selection -i %s -o %s -m 7 -t %s -w %d -r %s %s %s" ${ANALYSIS_DIR} ${ANAFILE2} ${PLOTFILE2} ${TUNEFILE} ${dial} ${REWEIGHTFILE2} ${OPT_RATIO} ${OPTPLOT}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/Selection -i %s -o %s -m 8 -t %s -w %d -r %s %s %s" ${ANALYSIS_DIR} ${ANAFILE3} ${PLOTFILE3} ${TUNEFILE} ${dial} ${REWEIGHTFILE3} ${OPT_RATIO} ${OPTPLOT}`
    echo -e "---------\n ${cmd}\n"; eval $cmd
    cmd=`printf "%s/Selection -i %s -o %s -m 8 -t %s -w %d -r %s %s %s" ${ANALYSIS_DIR} ${ANAFILE4} ${PLOTFILE4} ${TUNEFILE} ${dial} ${REWEIGHTFILE4} ${OPT_RATIO} ${OPTPLOT}`
    echo -e "---------\n ${cmd}\n"; eval $cmd

  done #rpa/rfg
fi


