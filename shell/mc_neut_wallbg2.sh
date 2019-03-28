#!/bin/sh

if [ $# -ne 2 ]; then
  echo "Two arguments are required.            "
  echo "==============================================="
  echo "Usage:                                         "
  echo "  ./mc_neut.sh <run_id> <srun_id>"
  echo "==============================================="
  exit
fi

RUNID=$1;
SRUNID=`expr $2 '*' 250`;

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
OUTPUT_DIR="/home/t2k/nchikuma/b2_data/mc_neut/tmp"


# File
NEUTFILE1=`printf "%s/nd1/b2_numubar_nd1_o_%05d.nt"  $NEUT_DIR $RUNID`
NEUTFILE2=`printf "%s/nd1/b2_numu_nd1_o_%05d.nt"     $NEUT_DIR $RUNID`
NEUTFILE3=`printf "%s/nd5/b2_numubar_nd5_o_%05d.nt"  $NEUT_DIR $RUNID`
NEUTFILE4=`printf "%s/nd5/b2_numu_nd5_o_%05d.nt"     $NEUT_DIR $RUNID`
NEUTFILE5=`printf "%s/nd6/b2_numubar_nd6_o_%05d.nt"  $NEUT_DIR $RUNID`
NEUTFILE6=`printf "%s/nd6/b2_numu_nd6_o_%05d.nt"     $NEUT_DIR $RUNID`

TUNEFILE="~/b2_data/jnubeam/tunefile/tune_nd7_8_9.root"


##Set up
cmd="source ${WORK_DIR}/Run_At_Start_B2.sh"
echo -e "---------\n ${cmd}\n"; eval $cmd
cmd="source ${INGRID_SOFT_DIR}/cmt/setup.sh";
echo -e "---------\n ${cmd}\n"; eval $cmd
cmd="source ${MC_DIR}/setup.sh";
echo -e "---------\n ${cmd}\n"; eval $cmd


for i in `seq 0 249`
do
  OUTPUTFILE1=`printf "%s/wallbg/numubar/ingrid_%08d_%04d_calib.root"      $OUTPUT_DIR $RUNID $SRUNID`
  OUTPUTFILE2=`printf "%s/wallbg/numu/ingrid_%08d_%04d_calib.root"         $OUTPUT_DIR $RUNID $SRUNID`
  OUTPUTFILE3_0=`printf "%s/ceiling/numubar/ingrid_%08d_%04d_calib.root"   $OUTPUT_DIR $RUNID $SRUNID`
  OUTPUTFILE4_0=`printf "%s/ceiling/numu/ingrid_%08d_%04d_calib.root"      $OUTPUT_DIR $RUNID $SRUNID`
  OUTPUTFILE3_1=`printf "%s/floor/numubar/ingrid_%08d_%04d_calib.root"     $OUTPUT_DIR $RUNID $SRUNID`
  OUTPUTFILE4_1=`printf "%s/floor/numu/ingrid_%08d_%04d_calib.root"        $OUTPUT_DIR $RUNID $SRUNID`
  OUTPUTFILE5_0=`printf "%s/pillar_r/numubar/ingrid_%08d_%04d_calib.root"  $OUTPUT_DIR $RUNID $SRUNID`
  OUTPUTFILE6_0=`printf "%s/pillar_r/numu/ingrid_%08d_%04d_calib.root"     $OUTPUT_DIR $RUNID $SRUNID`
  OUTPUTFILE5_1=`printf "%s/pillar_l/numubar/ingrid_%08d_%04d_calib.root"  $OUTPUT_DIR $RUNID $SRUNID`
  OUTPUTFILE6_1=`printf "%s/pillar_l/numu/ingrid_%08d_%04d_calib.root"     $OUTPUT_DIR $RUNID $SRUNID`

  cmd="rm -f $OUTPUTFILE1  "; echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd="rm -f $OUTPUTFILE2  "; echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd="rm -f $OUTPUTFILE3_0"; echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd="rm -f $OUTPUTFILE4_0"; echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd="rm -f $OUTPUTFILE3_1"; echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd="rm -f $OUTPUTFILE4_1"; echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd="rm -f $OUTPUTFILE5_0"; echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd="rm -f $OUTPUTFILE6_0"; echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd="rm -f $OUTPUTFILE5_1"; echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd="rm -f $OUTPUTFILE6_1"; echo -e "---------\n ${cmd}\n"; eval $cmd
  
  # Upst Wall
  cmd=`printf "timeout 2m %s/bin/Linux-g++/b2mc -m10 -f2 -i %s -o %s " $MC_DIR $NEUTFILE1 $OUTPUTFILE1`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "timeout 2m %s/bin/Linux-g++/b2mc -m10 -f1 -i %s -o %s " $MC_DIR $NEUTFILE2 $OUTPUTFILE2`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  # Ceiling
  cmd=`printf "timeout 2m %s/bin/Linux-g++/b2mc -m11 -f2 -i %s -o %s " $MC_DIR $NEUTFILE3 $OUTPUTFILE3_0`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "timeout 2m %s/bin/Linux-g++/b2mc -m11 -f1 -i %s -o %s " $MC_DIR $NEUTFILE4 $OUTPUTFILE4_0`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  # Floor
  cmd=`printf "timeout 2m %s/bin/Linux-g++/b2mc -m12 -f2 -i %s -o %s " $MC_DIR $NEUTFILE3 $OUTPUTFILE3_1`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "timeout 2m %s/bin/Linux-g++/b2mc -m12 -f1 -i %s -o %s " $MC_DIR $NEUTFILE4 $OUTPUTFILE4_1`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  # Right Pillar
  cmd=`printf "timeout 2m %s/bin/Linux-g++/b2mc -m13 -f2 -i %s -o %s " $MC_DIR $NEUTFILE5 $OUTPUTFILE5_0`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "timeout 2m %s/bin/Linux-g++/b2mc -m13 -f1 -i %s -o %s " $MC_DIR $NEUTFILE6 $OUTPUTFILE6_0`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  # Left Pillar
  cmd=`printf "timeout 2m %s/bin/Linux-g++/b2mc -m14 -f2 -i %s -o %s " $MC_DIR $NEUTFILE5 $OUTPUTFILE5_1`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "timeout 2m %s/bin/Linux-g++/b2mc -m14 -f1 -i %s -o %s " $MC_DIR $NEUTFILE6 $OUTPUTFILE6_1`
  echo -e "---------\n ${cmd}\n"; eval $cmd
  
  echo $SRUNID
  SRUNID=`expr $SRUNID + 1`
done
