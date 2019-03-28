#!/bin/sh

if [ $# -ne 2 ]; then
  echo "Two arguments are required.            "
  echo "==============================================="
  echo "Usage:                                         "
  echo "  ./analysis_B2cosmic_ing.sh <run_id> <srun_id>"
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

##Reconstruction for cosmic tirgger (except for B2 WAGASCI)
#cmd=`printf "%s/TwoDimRecon -r %d -s %d -cw" ${ANALYSIS_DIR} ${RUNID} ${SRUNID}`
#echo -e "---------\n ${cmd}\n"; eval $cmd

#cmd=`printf "%s/ThreeDimRecon -r %d -s %d -cw -m 1" ${ANALYSIS_DIR} ${RUNID} ${SRUNID}`
#echo -e "---------\n ${cmd}\n"; eval $cmd

###Hit Efficiency
#for i in `seq 0 7`; do
#  thres_pe=$i
#  cmd=`printf "%s/HitEfficiency -r %d -s %d -m1 -w -t%d" ${ANALYSIS_DIR} ${RUNID} ${SRUNID} ${thres_pe}`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#done

##P.E
thres_pe=4
cmd=`printf "%s/B2pe -r %d -s %d -m1 -w -t%d" ${ANALYSIS_DIR} ${RUNID} ${SRUNID} ${thres_pe}`
echo -e "---------\n ${cmd}\n"; eval $cmd

#P.E in scinti
thres_pe=0
cmd=`printf "%s/B2pe_Scinti -r %d -s %d -m1 -w -t%d" ${ANALYSIS_DIR} ${RUNID} ${SRUNID} ${thres_pe}`
echo -e "---------\n ${cmd}\n"; eval $cmd

###Hit Efficiency 2 (Distribution in Scinti)
#for i in `seq 4 4`; do
#  thres_pe=$i
#  cmd=`printf "%s/HitEfficiency_2 -r %d -s %d -m1 -w -t%d" ${ANALYSIS_DIR} ${RUNID} ${SRUNID} ${thres_pe}`
#  echo -e "---------\n ${cmd}\n"; eval $cmd
#done
