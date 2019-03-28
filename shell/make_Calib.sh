#!/bin/sh

if [ $# -ne 2 ]; then
  echo "Two arguments are required.            "
  echo "======================================="
  echo "Usage:                                 "
  echo "  ./make_Calib.sh <run_id> <srun_id>"
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
INGRID_FORM_DIR="$WORK_DIR/b2_software/INGRID/ingrid_format"

#Set up
cmd="source ${WORK_DIR}/b2_software/Run_At_Start_B2.sh"
echo -e "---------\n ${cmd}\n"; eval $cmd
cmd="source ${INGRID_SOFT_DIR}/cmt/setup.sh";
echo -e "---------\n ${cmd}\n"; eval $cmd


#Calibration
cmd=`printf "%s/amd64_linux26/Calc_MPPC_new.exe -r %d -s %d -t 1 -w" ${INGRID_SOFT_DIR} ${RUNID} ${SRUNID}`
echo -e "---------\n ${cmd}\n"; eval $cmd
