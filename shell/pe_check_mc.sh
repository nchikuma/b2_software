#!/bin/sh

if [ $# -ne 2 ]; then
  echo "Two arguments are required.            "
  echo "==============================================="
  echo "Usage:                                         "
  echo "  ./pe_check_mc.sh <run_id> <srun_id>"
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


##Set up
cmd="source ${WORK_DIR}/Run_At_Start_B2.sh"
echo -e "---------\n ${cmd}\n"; eval $cmd
cmd="source ${INGRID_SOFT_DIR}/cmt/setup.sh";
echo -e "---------\n ${cmd}\n"; eval $cmd

for i in `seq 5 5`; do
  thres_pe=$i
#  cmd=`printf "%s/B2pe -r %d -s %d -m2 -t%d" ${ANALYSIS_DIR} ${RUNID} ${SRUNID} ${thres_pe}`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
  cmd=`printf "%s/B2pe -r %d -s %d -m3 -t%d" ${ANALYSIS_DIR} ${RUNID} ${SRUNID} ${thres_pe}`
  echo -e "---------\n ${cmd}\n"; eval $cmd
done
