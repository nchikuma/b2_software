#!/bin/sh

if [ $# -ne 2 ]; then
  echo "Two arguments are required.            "
  echo "======================================="
  echo "Usage:                                 "
  echo "  ./decode_WAGASCI.sh <run_id> <srun_id>"
  echo "======================================="
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


WORK_DIR="/gpfs/fs03/t2k/beam/work/nchikuma/B2"
INGRID_SOFT_DIR="$WORK_DIR/b2_software/INGRID/INGRID/v1r1"  
WAGASCI_SOFT="$WORK_DIR/b2_software/wagasci_software/script/AnaAll.py" 
WAGASCI_SPILL_CORR="$WORK_DIR/b2_software/wagasci_software/Analysis/bin/wgSpillCorr"

#Set up
cmd="source ${WORK_DIR}/b2_software/Run_At_Start_B2.sh"
echo -e "---------\n ${cmd}\n"; eval $cmd
cmd="source ${INGRID_SOFT_DIR}/cmt/setup.sh";
echo -e "---------\n ${cmd}\n"; eval $cmd


raw_file_head=`printf "/gpfs/fs03/t2k/beam/work/nchikuma/B2/data/daqdata/run_%05d/run_%05d_%03d" ${RUNID} ${SRUNID}`
raw_file1=`printf "%s_dif_1_1_1.raw" ${raw_file_head}`
raw_file2=`printf "%s_dif_1_1_2.raw" ${raw_file_head}`


#Decode and Merge INGRID ligrary
#cmd=`printf "%s and %d %d all" ${WAGASCI_SOFT} ${RUNID} ${SRUNID}`
cmd=`printf "%s %d %d" ${WAGASCI_SOFT} ${RUNID} ${SRUNID}`
echo -e "---------\n ${cmd}\n"; eval $cmd

#
DECODE_DATA1=`printf "${WORK_DIR}/data/decode/run_%05d_%03d_dif_1_1_1_tree.root" ${RUNID} ${SRUNID}`
DECODE_DATA2=`printf "${WORK_DIR}/data/decode/run_%05d_%03d_dif_1_1_2_tree.root" ${RUNID} ${SRUNID}`
DST_DATA=`printf "${WORK_DIR}/data/data_dst_wagasci/run_%05d_%03d_inglib.root" ${RUNID} ${SRUNID}`
SPILL_DATA=`printf "${WORK_DIR}/data/data_dst_wagasci/run_%05d_%03d_spillcorr.root" ${RUNID} ${SRUNID}`
COSMIC_DATA=`printf "${WORK_DIR}/data/data_cosmic_wagasci/run_%05d_%03d_inglib.root" ${RUNID} ${SRUNID}`


#Remove the first decode file
sleep 1
if [ -f $DST_DATA -a -f $COSMIC_DATA ]; then
  echo "The following files are removed."
  echo $DECODE_DATA1
  echo $DECODE_DATA2
  rm -f $DECODE_DATA1
  rm -f $DECODE_DATA2
fi


#Spill Correction
cmd=`printf "%s -r %d -s %d" ${WAGASCI_SPILL_CORR} ${RUNID} ${SRUNID}`
echo -e "---------\n ${cmd}\n"; eval $cmd



if [ -f $DST_DATA ]; then
  #Remove the dst file (_inglib.root)
  sleep 1
  if [ -f $SPILL_DATA ]; then
    echo "The following files are removed."
    echo $DST_DATA
    rm -f $DST_DATA
  fi

fi
